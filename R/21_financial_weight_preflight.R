#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE,warn=1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
FIN_W <- "data/derived/financial_weights.csv"
TRADE_W <- Sys.getenv("FIN3_TRADE_WEIGHT_CSV",
                      "source_repo/8.12/Trade_Weights_14_Economies_2000_2014.csv")
OUT <- "results/preflight"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...),call.=FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

readW <- function(path) {
  if (!file.exists(path)) stopf("Weight file not found: %s",path)
  d <- read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  names(d) <- trimws(names(d))
  rr <- toupper(trimws(as.character(d[[1]])))
  cc <- toupper(trimws(names(d)[-1]))
  if (!all(COUNTRIES %in% rr) || !all(COUNTRIES %in% cc))
    stopf("Incomplete 14-country matrix: %s",path)
  W <- matrix(NA_real_,14,14,dimnames=list(COUNTRIES,COUNTRIES))
  for (i in COUNTRIES) {
    r <- match(i,rr)
    for (j in COUNTRIES) W[i,j] <- num(d[[match(j,cc)+1L]][r])
  }
  W
}

W <- readW(FIN_W)
audit <- data.frame(
  Metric=c("Rows","Columns","Nonfinite","Negative","MaxAbsDiagonal",
           "MaxAbsRowSumMinus1","MinWeight","MaxWeight"),
  Value=c(nrow(W),ncol(W),sum(!is.finite(W)),sum(W < -1e-12,na.rm=TRUE),
          max(abs(diag(W))),max(abs(rowSums(W)-1)),min(W),max(W))
)
write.csv(audit,file.path(OUT,"01_weight_validation.csv"),row.names=FALSE)

if (any(!is.finite(W))) stopf("Financial W contains non-finite values")
if (any(W < -1e-12)) stopf("Financial W contains negative values")
if (max(abs(diag(W))) > 1e-8) stopf("Financial W diagonal is not zero")
if (max(abs(rowSums(W)-1)) > 1e-6) stopf("Financial W rows do not sum to 1")

top5 <- do.call(rbind,lapply(COUNTRIES,function(i) {
  z <- sort(W[i,setdiff(COUNTRIES,i)],decreasing=TRUE)
  k <- seq_len(min(5,length(z)))
  data.frame(Country=i,Rank=k,Counterpart=names(z)[k],Weight=as.numeric(z[k]))
}))
write.csv(top5,file.path(OUT,"02_top_counterparts.csv"),row.names=FALSE)

conc <- data.frame(
  Country=COUNTRIES,
  HHI=vapply(COUNTRIES,function(i) sum(W[i,]^2),numeric(1)),
  MaxWeight=vapply(COUNTRIES,function(i) max(W[i,]),numeric(1)),
  EffectivePartners=vapply(COUNTRIES,function(i) 1/sum(W[i,]^2),numeric(1))
)
write.csv(conc,file.path(OUT,"03_weight_concentration.csv"),row.names=FALSE)

comparison <- NULL
if (file.exists(TRADE_W)) {
  WT <- readW(TRADE_W)
  diag(WT) <- 0
  WT <- WT/rowSums(WT)
  comparison <- do.call(rbind,lapply(COUNTRIES,function(i) {
    keep <- setdiff(COUNTRIES,i)
    a <- W[i,keep]; b <- WT[i,keep]
    data.frame(
      Country=i,
      Correlation=suppressWarnings(cor(a,b)),
      L1_Distance=sum(abs(a-b)),
      MaxAbsDifference=max(abs(a-b)),
      FinancialTop=names(which.max(W[i,])),
      TradeTop=names(which.max(WT[i,]))
    )
  }))
  write.csv(comparison,file.path(OUT,"04_financial_vs_trade_weight_shift.csv"),row.names=FALSE)
}

txt <- c(
  "Financial weight preflight: PASS",
  paste0("Max |diag(W)| = ",format(max(abs(diag(W))),scientific=TRUE)),
  paste0("Max |rowsum(W)-1| = ",format(max(abs(rowSums(W)-1)),scientific=TRUE)),
  paste0("Max weight = ",signif(max(W),8)),
  paste0("Trade comparison available: ",!is.null(comparison))
)
writeLines(txt,file.path(OUT,"README_preflight.txt"))
cat(paste(txt,collapse="\n"),"\n")
