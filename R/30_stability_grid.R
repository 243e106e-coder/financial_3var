#!/usr/bin/env Rscript

source("R/00_config.R")

PANEL <- file.path(DERIVED_DIR,"panel_domestic_fin3.csv")
OUT <- file.path(RESULTS_DIR,"stability")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

if (!file.exists(PANEL)) stopf("Run R/10_build_financial_3var_input.R first.")
d <- read.csv(PANEL,check.names=FALSE,stringsAsFactors=FALSE)
need <- c("Quarter","Country",VARS,"gpr","brent")
if (!all(need %in% names(d))) stopf("Domestic panel is missing required columns.")

quarters <- unique(d$Quarter)
if (!all(COUNTRIES %in% unique(d$Country))) stopf("Domestic panel is missing sample economies.")
if (any(table(d$Country) != length(quarters))) stopf("Domestic panel is not balanced.")

N <- length(COUNTRIES); K <- length(VARS)
X <- array(NA_real_,c(length(quarters),N,K),dimnames=list(quarters,COUNTRIES,VARS))
for (i in seq_along(COUNTRIES)) {
  z <- d[d$Country==COUNTRIES[i],,drop=FALSE]
  z <- z[match(quarters,z$Quarter),,drop=FALSE]
  X[,i,] <- as.matrix(z[,VARS,drop=FALSE])
}
base <- d[d$Country==COUNTRIES[1],,drop=FALSE]
base <- base[match(quarters,base$Quarter),,drop=FALSE]
gpr <- num(base$gpr)
brent <- num(base$brent)
if (!all(is.finite(X)) || !all(is.finite(gpr))) stopf("Non-finite model inputs.")
include_brent <- all(is.finite(brent))

lag1 <- function(x) c(NA_real_,x[-length(x)])

make_star <- function(W) {
  S <- array(NA_real_,dim(X),dimnames=dimnames(X))
  for (i in seq_len(N)) for (v in seq_len(K)) {
    S[,i,v] <- as.numeric(X[,,v] %*% W[i,])
  }
  S
}

fit_country <- function(i,S) {
  Y <- X[,i,,drop=FALSE][,1,]
  Z <- S[,i,,drop=FALSE][,1,]
  rows <- 2:nrow(Y)
  D <- data.frame(const=rep(1,length(rows)))
  for (v in seq_len(K)) D[[paste0(VARS[v],"_L1")]] <- lag1(Y[,v])[rows]
  for (v in seq_len(K)) D[[paste0(VARS[v],"_star_0")]] <- Z[rows,v]
  for (v in seq_len(K)) D[[paste0(VARS[v],"_star_L1")]] <- lag1(Z[,v])[rows]
  D$gpr_0 <- gpr[rows]
  D$gpr_L1 <- lag1(gpr)[rows]
  if (include_brent) {
    D$brent_0 <- brent[rows]
    D$brent_L1 <- lag1(brent)[rows]
  }
  YY <- Y[rows,,drop=FALSE]
  ok <- complete.cases(D) & complete.cases(YY)
  D <- D[ok,,drop=FALSE]; YY <- YY[ok,,drop=FALSE]
  Xm <- as.matrix(D)
  if (nrow(Xm) <= ncol(Xm)+5L) stopf("Too few observations for %s",COUNTRIES[i])

  B <- matrix(NA_real_,ncol(Xm),K,dimnames=list(colnames(Xm),VARS))
  for (eq in seq_len(K)) {
    fit <- lm.fit(Xm,YY[,eq])
    if (any(!is.finite(fit$coefficients))) stopf("Rank deficiency: %s/%s",COUNTRIES[i],VARS[eq])
    B[,eq] <- fit$coefficients
  }

  A <- B0 <- B1 <- matrix(0,K,K,dimnames=list(VARS,VARS))
  for (eq in VARS) for (v in VARS) {
    A[eq,v] <- B[paste0(v,"_L1"),eq]
    B0[eq,v] <- B[paste0(v,"_star_0"),eq]
    B1[eq,v] <- B[paste0(v,"_star_L1"),eq]
  }
  rho <- max(Mod(eigen(A,only.values=TRUE)$values))
  list(A=A,B0=B0,B1=B1,n=nrow(YY),rho=rho,
       design_kappa=kappa(crossprod(Xm)),design_rcond=rcond(crossprod(Xm)))
}

selector <- function(i) {
  S <- matrix(0,K,N*K)
  S[,((i-1)*K+1):(i*K)] <- diag(K)
  S
}
star_map <- function(i,W) {
  R <- matrix(0,K,N*K)
  for (j in seq_len(N)) for (v in seq_len(K)) R[v,(j-1)*K+v] <- W[i,j]
  R
}

fit_global <- function(W,label) {
  S <- make_star(W)
  fits <- lapply(seq_len(N),fit_country,S=S)
  G0 <- matrix(0,N*K,N*K); G1 <- matrix(0,N*K,N*K)
  for (i in seq_len(N)) {
    rr <- ((i-1)*K+1):(i*K)
    Si <- selector(i); Ri <- star_map(i,W)
    G0[rr,] <- Si - fits[[i]]$B0 %*% Ri
    G1[rr,] <- fits[[i]]$A %*% Si + fits[[i]]$B1 %*% Ri
  }
  rc <- rcond(G0); kap <- kappa(G0)
  F <- tryCatch(solve(G0,G1),error=function(e) NULL)
  ev <- if (is.null(F)) rep(NA_complex_,N*K) else eigen(F,only.values=TRUE)$values
  rho <- if (all(is.na(ev))) NA_real_ else max(Mod(ev),na.rm=TRUE)
  local <- data.frame(
    Network=label,Country=COUNTRIES,
    N=vapply(fits,`[[`,numeric(1),"n"),
    DomesticSpectralRadius=vapply(fits,`[[`,numeric(1),"rho"),
    DomesticStable=vapply(fits,function(z) z$rho<1,logical(1)),
    DesignKappa=vapply(fits,`[[`,numeric(1),"design_kappa"),
    DesignRcond=vapply(fits,`[[`,numeric(1),"design_rcond")
  )
  list(label=label,rho=rho,rcond=rc,kappa=kap,eigen=ev,local=local)
}

networks <- lapply(WEIGHT_FILES, read_weight_matrix)
if (file.exists(TRADE_WEIGHT_PATH)) networks$legacy_trade <- read_weight_matrix(TRADE_WEIGHT_PATH)
fits <- lapply(names(networks), function(nm) fit_global(networks[[nm]],nm))
names(fits) <- names(networks)

summary <- do.call(rbind,lapply(fits,function(z) data.frame(
  Network=z$label,
  GlobalSpectralRadius=z$rho,
  Stable=is.finite(z$rho) && z$rho<1,
  G0_rcond=z$rcond,G0_kappa=z$kappa,
  MaxLocalDomesticRadius=max(z$local$DomesticSpectralRadius),
  UnstableLocalCount=sum(!z$local$DomesticStable),
  stringsAsFactors=FALSE
)))
write.csv(summary,file.path(OUT,"01_global_stability_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(fits,`[[`,"local")),
          file.path(OUT,"02_local_stability.csv"),row.names=FALSE)

eigen_rows <- do.call(rbind,lapply(fits,function(z) {
  q <- data.frame(Network=z$label,Real=Re(z$eigen),Imag=Im(z$eigen),Modulus=Mod(z$eigen))
  q[order(q$Modulus,decreasing=TRUE),,drop=FALSE]
}))
write.csv(eigen_rows,file.path(OUT,"03_global_eigenvalues.csv"),row.names=FALSE)

spec <- data.frame(
  Item=c("Countries","Variables","p","q","ContemporaneousForeignVariables",
         "GPR","Brent","Sample","MainFinancialW"),
  Value=c(paste(COUNTRIES,collapse=", "),paste(VARS,collapse=", "),"1","1","YES",
          "current + lag 1",ifelse(include_brent,"current + lag 1","not included"),
          paste0(quarters[1]," - ",tail(quarters,1)),MAIN_NETWORK)
)
write.csv(spec,file.path(OUT,"04_model_specification.csv"),row.names=FALSE)

main_row <- summary[summary$Network==MAIN_NETWORK,,drop=FALSE]
txt <- c(
  "CLEAN FINANCIAL 3-VARIABLE GVAR STABILITY GRID",
  "================================================",
  sprintf("Sample: %s - %s",quarters[1],tail(quarters,1)),
  "Variables: r, de, deq; p=1; q=1",
  sprintf("Networks: %s",paste(summary$Network,collapse=", ")),
  sprintf("Main spectral radius: %.10g",main_row$GlobalSpectralRadius),
  sprintf("Main stable (<1): %s",main_row$Stable),
  sprintf("Main G0 rcond: %.10g",main_row$G0_rcond),
  "This is an OLS diagnostic, not the final Bayesian TVP estimation."
)
writeLines(txt,file.path(OUT,"README_stability.txt"))
cat(paste(txt,collapse="\n"),"\n")
