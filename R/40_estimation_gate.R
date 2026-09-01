#!/usr/bin/env Rscript

source("R/00_config.R")

input <- file.path(RESULTS_DIR,"stability","01_global_stability_comparison.csv")
OUT <- file.path(RESULTS_DIR,"gate")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
if (!file.exists(input)) stopf("Run R/30_stability_grid.R first.")

d <- read.csv(input,stringsAsFactors=FALSE)
z <- d[d$Network==MAIN_NETWORK,,drop=FALSE]
if (nrow(z)!=1L) stopf("Main network missing or duplicated in stability output.")

rho_ok <- is.finite(z$GlobalSpectralRadius) && z$GlobalSpectralRadius < 1
g0_ok <- is.finite(z$G0_rcond) && z$G0_rcond >= 1e-8
local_ok <- is.finite(z$UnstableLocalCount) && z$UnstableLocalCount == 0
ready <- rho_ok && g0_ok && local_ok

status <- if (ready) "READY" else "BLOCKED"
reason <- c(
  if (!rho_ok) "global spectral radius is not below 1" else NULL,
  if (!g0_ok) "G0 reciprocal condition number is below 1e-8" else NULL,
  if (!local_ok) "at least one local domestic system is unstable" else NULL
)
if (!length(reason)) reason <- "all pre-estimation stability gates passed"

out <- data.frame(
  MainNetwork=MAIN_NETWORK,
  Status=status,
  GlobalSpectralRadius=z$GlobalSpectralRadius,
  G0_rcond=z$G0_rcond,
  UnstableLocalCount=z$UnstableLocalCount,
  Reason=paste(reason,collapse="; "),
  stringsAsFactors=FALSE
)
write.csv(out,file.path(OUT,"01_estimation_gate.csv"),row.names=FALSE)
lines <- c(
  sprintf("ESTIMATION GATE: %s",status),
  "================================",
  sprintf("Main network: %s",MAIN_NETWORK),
  sprintf("Global spectral radius: %.10g",z$GlobalSpectralRadius),
  sprintf("G0 rcond: %.10g",z$G0_rcond),
  sprintf("Unstable local systems: %d",z$UnstableLocalCount),
  sprintf("Reason: %s",paste(reason,collapse="; ")),
  if (ready) "The Bayesian TVP stage may proceed." else "Do not start the Bayesian TVP stage."
)
writeLines(lines,file.path(OUT,"ESTIMATION_GATE.txt"))
cat(paste(lines,collapse="\n"),"\n")
