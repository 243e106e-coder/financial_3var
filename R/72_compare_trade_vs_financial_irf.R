#!/usr/bin/env Rscript

# =============================================================================
# 72_compare_trade_vs_financial_irf.R
#
# Ex-post comparison of a newly re-estimated trade-W Delta-r TVP-GVAR against
# the already accepted 2017 GCAP financial-W IRF artifact.
#
# IMPORTANT:
# - This script NEVER chooses the preferred network because it has more
#   significant results.
# - It reports stability, posterior direction/credibility and matched IRF
#   differences under two independently estimated networks.
# =============================================================================

get_env_chr <- function(name, default="") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}
stopf <- function(...) stop(sprintf(...), call.=FALSE)

TRADE_ROOT <- get_env_chr("FIN3_TRADE_IRF_ROOT", "trade_irf/formal_irf_rate_diff")
FIN_ROOT <- get_env_chr("FIN3_FINANCIAL_IRF_ROOT", "financial_ref/formal_irf_rate_diff")
OUT <- get_env_chr("FIN3_NETWORK_COMPARE_OUT", "results/network_comparison")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

need_files <- c(
  "00_formal_rate_diff_irf_gate.csv",
  "irf_posterior_summary.csv",
  "03_cumulative_reer_posterior_summary.csv"
)
for (f in need_files) {
  if (!file.exists(file.path(TRADE_ROOT,f))) stopf("Missing trade result: %s", f)
  if (!file.exists(file.path(FIN_ROOT,f))) stopf("Missing financial reference: %s", f)
}

tg <- read.csv(file.path(TRADE_ROOT,"00_formal_rate_diff_irf_gate.csv"), stringsAsFactors=FALSE)
fg <- read.csv(file.path(FIN_ROOT,"00_formal_rate_diff_irf_gate.csv"), stringsAsFactors=FALSE)
if (nrow(tg)!=1L || nrow(fg)!=1L) stopf("Malformed IRF gate.")
if (tg$Status[1]!="READY_FOR_RATE_DIFF_IRF_AUDIT") stopf("Trade IRF is not accepted.")
if (fg$Status[1]!="READY_FOR_RATE_DIFF_IRF_AUDIT") stopf("Financial reference IRF is not accepted.")

gate_compare <- data.frame(
  Metric=c(
    "MainNetwork",
    "WorstAnchorStableShare",
    "WorstAnchorG0OKShare",
    "MinValidRawIRFDraws",
    "PosteriorDrawsTotal"
  ),
  Trade=c(
    as.character(tg$MainNetwork[1]),
    as.character(tg$WorstAnchorStableShare[1]),
    as.character(tg$WorstAnchorG0OKShare[1]),
    as.character(tg$MinValidRawIRFDraws[1]),
    as.character(tg$PosteriorDrawsTotal[1])
  ),
  Financial=c(
    as.character(fg$MainNetwork[1]),
    as.character(fg$WorstAnchorStableShare[1]),
    as.character(fg$WorstAnchorG0OKShare[1]),
    as.character(fg$MinValidRawIRFDraws[1]),
    as.character(fg$PosteriorDrawsTotal[1])
  ),
  stringsAsFactors=FALSE
)
write.csv(gate_compare, file.path(OUT,"00_network_gate_comparison.csv"), row.names=FALSE)

trade_raw <- read.csv(file.path(TRADE_ROOT,"irf_posterior_summary.csv"), stringsAsFactors=FALSE)
fin_raw <- read.csv(file.path(FIN_ROOT,"irf_posterior_summary.csv"), stringsAsFactors=FALSE)
trade_reer <- read.csv(file.path(TRADE_ROOT,"03_cumulative_reer_posterior_summary.csv"), stringsAsFactors=FALSE)
fin_reer <- read.csv(file.path(FIN_ROOT,"03_cumulative_reer_posterior_summary.csv"), stringsAsFactors=FALSE)

key_cols <- c("EventID","EventSet","AnchorType","AnchorQuarter","Country","ResponseVariable","Horizon")
need_cols <- c(key_cols,"p05","median","p95","prob_positive")
for (nm in c("trade_raw","fin_raw","trade_reer","fin_reer")) {
  z <- get(nm)
  miss <- setdiff(need_cols,names(z))
  if (length(miss)) stopf("%s missing columns: %s", nm, paste(miss,collapse=","))
}

# Core cumulative REER comparison at substantively interpretable horizons.
countries_focus <- c("JP","CH","US","TR")
h_focus <- c(1,4,8,12)

summarize_reer <- function(d, label) {
  z <- d[
    d$EventSet=="CORE" &
    d$AnchorType=="EVENT_QUARTER_t0" &
    d$Country %in% countries_focus &
    d$Horizon %in% h_focus,
    ,
    drop=FALSE
  ]
  if (!nrow(z)) stopf("No core cumulative REER rows for %s.", label)
  do.call(rbind, lapply(countries_focus, function(cc) {
    x <- z[z$Country==cc,,drop=FALSE]
    data.frame(
      Network=label,
      Country=cc,
      Cells=nrow(x),
      MedianProbPositive=median(x$prob_positive,na.rm=TRUE),
      PositiveCredible90=sum(x$p05>0,na.rm=TRUE),
      NegativeCredible90=sum(x$p95<0,na.rm=TRUE),
      MedianCumulativeREER=median(x$median,na.rm=TRUE),
      stringsAsFactors=FALSE
    )
  }))
}
safe <- rbind(
  summarize_reer(trade_reer,"TRADE"),
  summarize_reer(fin_reer,"FINANCIAL_GCAP")
)
write.csv(safe, file.path(OUT,"01_core_safe_haven_direction_comparison.csv"), row.names=FALSE)

# TR raw de h=1 comparison, the anomaly identified in 08o/08p.
summarize_tr <- function(d,label) {
  z <- d[
    d$EventSet=="CORE" &
    d$AnchorType=="EVENT_QUARTER_t0" &
    d$Country=="TR" &
    d$ResponseVariable=="de" &
    d$Horizon==1,
    ,
    drop=FALSE
  ]
  if (!nrow(z)) stopf("No TR/de/h1 core rows for %s.",label)
  data.frame(
    Network=label,
    Events=nrow(z),
    MeanMedian=mean(z$median),
    MedianMedian=median(z$median),
    MinMedian=min(z$median),
    MaxMedian=max(z$median),
    MedianProbPositive=median(z$prob_positive),
    PositiveCredible90=sum(z$p05>0),
    NegativeCredible90=sum(z$p95<0),
    stringsAsFactors=FALSE
  )
}
tr <- rbind(summarize_tr(trade_raw,"TRADE"), summarize_tr(fin_raw,"FINANCIAL_GCAP"))
write.csv(tr, file.path(OUT,"02_tr_raw_reer_h1_comparison.csv"), row.names=FALSE)

# Match every raw IRF cell and measure the size of network-induced change.
merge_keys <- key_cols
a <- trade_raw[,c(merge_keys,"median","p05","p95","prob_positive")]
b <- fin_raw[,c(merge_keys,"median","p05","p95","prob_positive")]
names(a)[-(seq_along(merge_keys))] <- paste0(names(a)[-(seq_along(merge_keys))],"_Trade")
names(b)[-(seq_along(merge_keys))] <- paste0(names(b)[-(seq_along(merge_keys))],"_Financial")
m <- merge(a,b,by=merge_keys,all=FALSE)
if (!nrow(m)) stopf("No matched raw IRF cells between networks.")
m$MedianDifference_TradeMinusFinancial <- m$median_Trade - m$median_Financial
m$AbsMedianDifference <- abs(m$MedianDifference_TradeMinusFinancial)
write.csv(m, file.path(OUT,"03_matched_raw_irf_cells.csv"), row.names=FALSE)

diff_summary <- do.call(rbind,lapply(c("r","de","deq"),function(vv){
  z <- m[m$ResponseVariable==vv,,drop=FALSE]
  data.frame(
    ResponseVariable=vv,
    MatchedCells=nrow(z),
    MeanAbsMedianDifference=mean(z$AbsMedianDifference),
    MedianAbsMedianDifference=median(z$AbsMedianDifference),
    MaxAbsMedianDifference=max(z$AbsMedianDifference),
    MedianCorrelation=cor(z$median_Trade,z$median_Financial,method="spearman"),
    stringsAsFactors=FALSE
  )
}))
write.csv(diff_summary,file.path(OUT,"04_raw_irf_network_difference_summary.csv"),row.names=FALSE)

gate <- data.frame(
  Status="NETWORK_COMPARISON_COMPLETE",
  TradeNetwork=as.character(tg$MainNetwork[1]),
  FinancialNetwork=as.character(fg$MainNetwork[1]),
  TradeWorstStableShare=as.numeric(tg$WorstAnchorStableShare[1]),
  FinancialWorstStableShare=as.numeric(fg$WorstAnchorStableShare[1]),
  TradeTR_h1_Median=tr$MedianMedian[tr$Network=="TRADE"],
  FinancialTR_h1_Median=tr$MedianMedian[tr$Network=="FINANCIAL_GCAP"],
  NetworkSelectionBySignificancePerformed=FALSE,
  stringsAsFactors=FALSE
)
write.csv(gate,file.path(OUT,"00_network_comparison_gate.csv"),row.names=FALSE)

txt <- c(
  "TRADE VS FINANCIAL NETWORK COMPARISON: COMPLETE",
  "===============================================",
  sprintf("Trade network: %s", gate$TradeNetwork),
  sprintf("Financial reference: %s", gate$FinancialNetwork),
  sprintf("Worst posterior stable share: trade=%.6f; financial=%.6f",
          gate$TradeWorstStableShare, gate$FinancialWorstStableShare),
  sprintf("TR raw REER h=1 median across core events: trade=%.6f; financial=%.6f",
          gate$TradeTR_h1_Median, gate$FinancialTR_h1_Median),
  "",
  "No preferred network is selected on the basis of statistical significance.",
  "Use data provenance, coverage, structural stability and robustness jointly."
)
writeLines(txt,file.path(OUT,"README_network_comparison.txt"))
cat(paste(txt,collapse="\n"),"\n")
