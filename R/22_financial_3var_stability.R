#!/usr/bin/env Rscript

# p=1, q=1 OLS GVAR stability diagnostic.
# This is NOT the final Bayesian TVP estimation.

options(stringsAsFactors=FALSE,warn=1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("r","de","deq")
N <- length(COUNTRIES); K <- length(VARS)

PANEL <- "data/derived/panel_fin3_long.csv"
FIN_W <- "data/derived/financial_weights.csv"
TRADE_W <- Sys.getenv("FIN3_TRADE_WEIGHT_CSV",
                      "source_repo/8.12/Trade_Weights_14_Economies_2000_2014.csv")
OUT <- "results/stability"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...),call.=FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
lag1 <- function(x) c(NA_real_,x[-length(x)])

readW <- function(path) {
  d <- read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  names(d) <- trimws(names(d))
  rr <- toupper(trimws(as.character(d[[1]])))
  cc <- toupper(trimws(names(d)[-1]))
  if (!all(COUNTRIES %in% rr) || !all(COUNTRIES %in% cc))
    stopf("Incomplete W: %s",path)
  W <- matrix(NA_real_,N,N,dimnames=list(COUNTRIES,COUNTRIES))
  for (i in COUNTRIES) {
    r <- match(i,rr)
    for (j in COUNTRIES) W[i,j] <- num(d[[match(j,cc)+1L]][r])
  }
  if (any(!is.finite(W)) || any(W < -1e-12)) stopf("Invalid W: %s",path)
  diag(W) <- 0
  if (any(rowSums(W) <= 0)) stopf("Non-positive row sum: %s",path)
  W/rowSums(W)
}

if (!file.exists(PANEL)) stopf("Run R/20_build_financial_3var_input.R first")
d <- read.csv(PANEL,check.names=FALSE,stringsAsFactors=FALSE)

need <- c("Quarter","Country",VARS,"gpr","brent")
if (!all(need %in% names(d))) stopf("Panel missing required columns")
quarters <- unique(d$Quarter)
if (any(table(d$Country) != length(quarters))) stopf("Panel is not balanced")

X <- array(NA_real_,c(length(quarters),N,K),
           dimnames=list(quarters,COUNTRIES,VARS))
for (i in seq_along(COUNTRIES)) {
  z <- d[d$Country==COUNTRIES[i],,drop=FALSE]
  z <- z[match(quarters,z$Quarter),,drop=FALSE]
  X[,i,] <- as.matrix(z[,VARS,drop=FALSE])
}
base <- d[d$Country==COUNTRIES[1],,drop=FALSE]
base <- base[match(quarters,base$Quarter),,drop=FALSE]
gpr <- num(base$gpr)
brent <- num(base$brent)
include_brent <- all(is.finite(brent))
if (!all(is.finite(gpr))) stopf("GPR has missing values")

make_star <- function(W) {
  S <- array(NA_real_,dim(X),dimnames=dimnames(X))
  for (i in seq_len(N))
    for (v in seq_len(K))
      S[,i,v] <- as.numeric(X[,,v] %*% W[i,])
  S
}

fit_country <- function(i,S) {
  Y <- X[,i,,drop=FALSE][,1,]
  Z <- S[,i,,drop=FALSE][,1,]
  rows <- 2:nrow(Y)

  D <- data.frame(const=1)
  D <- D[rep(1,length(rows)),,drop=FALSE]
  rownames(D) <- NULL

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
    co <- fit$coefficients
    if (any(!is.finite(co))) stopf("Rank deficiency: %s / %s",COUNTRIES[i],VARS[eq])
    B[,eq] <- co
  }

  A <- B0 <- B1 <- matrix(0,K,K,dimnames=list(VARS,VARS))
  for (eq in VARS) for (v in VARS) {
    A[eq,v] <- B[paste0(v,"_L1"),eq]
    B0[eq,v] <- B[paste0(v,"_star_0"),eq]
    B1[eq,v] <- B[paste0(v,"_star_L1"),eq]
  }
  rho <- max(Mod(eigen(A,only.values=TRUE)$values))
  list(A=A,B0=B0,B1=B1,n=nrow(YY),rho=rho,
       design_kappa=kappa(crossprod(Xm)),
       design_rcond=rcond(crossprod(Xm)))
}

selector <- function(i) {
  S <- matrix(0,K,N*K)
  S[,((i-1)*K+1):(i*K)] <- diag(K)
  S
}
star_map <- function(i,W) {
  R <- matrix(0,K,N*K)
  for (j in seq_len(N)) for (v in seq_len(K))
    R[v,(j-1)*K+v] <- W[i,j]
  R
}

fit_global <- function(W,label) {
  S <- make_star(W)
  fits <- lapply(seq_len(N),fit_country,S=S)

  G0 <- matrix(0,N*K,N*K)
  G1 <- matrix(0,N*K,N*K)
  for (i in seq_len(N)) {
    rr <- ((i-1)*K+1):(i*K)
    Si <- selector(i)
    Ri <- star_map(i,W)
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

FIN <- readW(FIN_W)
fin <- fit_global(FIN,"GCAP_financial_2017")

trade <- NULL
if (file.exists(TRADE_W)) trade <- fit_global(readW(TRADE_W),"legacy_trade")

summary <- data.frame(
  Network=fin$label,
  GlobalSpectralRadius=fin$rho,
  Stable=is.finite(fin$rho) && fin$rho<1,
  G0_rcond=fin$rcond,
  G0_kappa=fin$kappa,
  MaxLocalDomesticRadius=max(fin$local$DomesticSpectralRadius),
  UnstableLocalCount=sum(!fin$local$DomesticStable)
)
if (!is.null(trade)) {
  summary <- rbind(summary,data.frame(
    Network=trade$label,
    GlobalSpectralRadius=trade$rho,
    Stable=is.finite(trade$rho) && trade$rho<1,
    G0_rcond=trade$rcond,
    G0_kappa=trade$kappa,
    MaxLocalDomesticRadius=max(trade$local$DomesticSpectralRadius),
    UnstableLocalCount=sum(!trade$local$DomesticStable)
  ))
}

write.csv(summary,file.path(OUT,"01_global_stability_comparison.csv"),row.names=FALSE)
write.csv(fin$local,file.path(OUT,"02_local_financial_stability.csv"),row.names=FALSE)
if (!is.null(trade)) write.csv(trade$local,file.path(OUT,"03_local_trade_stability.csv"),row.names=FALSE)

eigout <- function(x) {
  z <- data.frame(Real=Re(x$eigen),Imag=Im(x$eigen),Modulus=Mod(x$eigen))
  z[order(z$Modulus,decreasing=TRUE),,drop=FALSE]
}
write.csv(eigout(fin),file.path(OUT,"04_financial_global_eigenvalues.csv"),row.names=FALSE)
if (!is.null(trade)) write.csv(eigout(trade),file.path(OUT,"05_trade_global_eigenvalues.csv"),row.names=FALSE)

spec <- data.frame(
  Item=c("Countries","Variables","p","q","Contemporaneous foreign variables",
         "GPR","Brent","Sample","Financial W"),
  Value=c(paste(COUNTRIES,collapse=", "),paste(VARS,collapse=", "),
          "1","1","YES","current + lag 1",
          ifelse(include_brent,"current + lag 1","not included"),
          paste0(quarters[1]," - ",tail(quarters,1)),
          "GCAP 2017 Position Residency")
)
write.csv(spec,file.path(OUT,"06_model_specification.csv"),row.names=FALSE)

txt <- c(
  "3-variable financial GVAR stability diagnostic",
  paste0("Sample: ",quarters[1]," - ",tail(quarters,1)),
  "Variables: r, de, deq",
  "p=1, q=1; contemporaneous foreign variables included",
  paste0("Financial global spectral radius: ",signif(fin$rho,8)),
  paste0("Financial stable (<1): ",is.finite(fin$rho) && fin$rho<1),
  paste0("Financial G0 rcond: ",signif(fin$rcond,8)),
  paste0("Financial G0 kappa: ",signif(fin$kappa,8)),
  if (!is.null(trade))
    paste0("Trade spectral radius under same 3-variable specification: ",signif(trade$rho,8))
  else "Trade comparison unavailable",
  "Diagnostic OLS GVAR only; not the final Bayesian TVP estimation."
)
writeLines(txt,file.path(OUT,"README_stability.txt"))
cat(paste(txt,collapse="\n"),"\n")
