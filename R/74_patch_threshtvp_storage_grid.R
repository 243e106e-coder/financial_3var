#!/usr/bin/env Rscript

# =============================================================================
# 74_patch_threshtvp_storage_grid.R
#
# Implementation-only patch for the pinned threshtvp MCMC_tvp storage grid.
#
# 10q proved that, under the 10n unrestricted contract,
#
#   thin * nsave = 1001.0000000000001
#
# leads the original source to allocate 1001 posterior rows while seq() creates
# 1002 storage times.  The final write therefore attempts row 1002 and throws
# "subscript out of bounds".
#
# This patch changes ONLY storage bookkeeping:
#   1. define one integer nstore = round(thin * nsave);
#   2. use nstore for length.out and every MCMC_tvp storage dimension;
#   3. require the generated storage grid to have exactly nstore unique points.
#
# It does NOT change likelihood, priors, state equations, TVS, SV, thresholds,
# kappa0, random-number draws before storage, or the MCMC target.
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/74_patch_threshtvp_storage_grid.R <MCMC_function.R>")
}

path <- normalizePath(args[[1]], mustWork=TRUE)
txt <- paste(readLines(path, warn=FALSE), collapse="\n")

count_fixed <- function(haystack, needle) {
  m <- gregexpr(needle, haystack, fixed=TRUE)[[1]]
  if (length(m)==1L && m[[1]] == -1L) return(0L)
  length(m)
}

replace_once <- function(x, old, new, label) {
  n <- count_fixed(x, old)
  if (n != 1L) {
    stop(sprintf(
      "Patch contract failed for %s: expected exactly 1 occurrence, found %d.",
      label, n
    ))
  }
  sub(old, new, x, fixed=TRUE)
}

grid_old <- "thin_out <- round(seq(nburn,ntot,length.out = thin*nsave))"
grid_new <- paste(
  "nstore <- as.integer(round(thin*nsave))",
  "if (!is.finite(nstore) || nstore < 1L) stop(\"Invalid integer posterior storage count\")",
  "thin_out <- round(seq(nburn,ntot,length.out = nstore))",
  paste0(
    "if (length(thin_out) != nstore || length(unique(thin_out)) != nstore) ",
    "stop(\"Posterior storage grid length/uniqueness contract failed\")"
  ),
  sep="\n    "
)

txt <- replace_once(
  txt, grid_old, grid_new,
  "original MCMC_tvp thin_out definition"
)

storage_pairs <- list(
  c("H_store <- matrix(NA,thin*nsave,T)",
    "H_store <- matrix(NA,nstore,T)"),
  c("ALPHA_store <- array(NA,c(thin*nsave,T,K_))",
    "ALPHA_store <- array(NA,c(nstore,T,K_))"),
  c("V0_store <- array(NA,c(thin*nsave,K_))",
    "V0_store <- array(NA,c(nstore,K_))"),
  c("svparms_store <- matrix(NA,thin*nsave,3)",
    "svparms_store <- matrix(NA,nstore,3)"),
  c("D_store <- array(NA,c(thin*nsave,K_,T))",
    "D_store <- array(NA,c(nstore,K_,T))"),
  c("thresholds_store <- array(NA,c(thin*nsave,K_))",
    "thresholds_store <- array(NA,c(nstore,K_))"),
  c("Omega_store <- array(NA,c(thin*nsave,K_,T))",
    "Omega_store <- array(NA,c(nstore,K_,T))"),
  c("omega_store <- matrix(NA,thin*nsave,K_)",
    "omega_store <- matrix(NA,nstore,K_)"),
  c("sigma2_store <- matrix(NA,thin*nsave,1)",
    "sigma2_store <- matrix(NA,nstore,1)"),
  c("thrshprior_store <- matrix(NA,thin*nsave,K_)",
    "thrshprior_store <- matrix(NA,nstore,K_)"),
  c("kappa_store <- matrix(NA,thin*nsave, 1)",
    "kappa_store <- matrix(NA,nstore, 1)")
)

for (i in seq_along(storage_pairs)) {
  txt <- replace_once(
    txt,
    storage_pairs[[i]][[1]],
    storage_pairs[[i]][[2]],
    sprintf("storage allocation %02d", i)
  )
}

# Hard post-patch checks: the dangerous original expressions must be gone from
# the original MCMC_tvp storage definitions, and the integer contract present.
if (count_fixed(txt, grid_old) != 0L) {
  stop("Original floating thin_out expression remains after patch.")
}
for (p in storage_pairs) {
  if (count_fixed(txt, p[[1]]) != 0L) {
    stop(sprintf("Original floating storage allocation remains: %s", p[[1]]))
  }
}
if (count_fixed(txt, "nstore <- as.integer(round(thin*nsave))") != 1L) {
  stop("Integer nstore definition missing or duplicated.")
}

writeLines(strsplit(txt, "\n", fixed=TRUE)[[1]], path, useBytes=TRUE)

manifest <- data.frame(
  Item=c(
    "TechnicalStatus",
    "PatchType",
    "GridDefinitionPatched",
    "StorageAllocationsPatched",
    "StatisticalModelChanged",
    "FormalModelChanged",
    "IRFAuthorized"
  ),
  Value=c(
    "PASS",
    "INTEGERIZE_ORIGINAL_MCMC_TVP_STORAGE_COUNT",
    "TRUE",
    as.character(length(storage_pairs)),
    "FALSE",
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors=FALSE
)

out <- Sys.getenv(
  "FIN3_10R_PATCH_MANIFEST",
  "results/10r_storage_patch_manifest.csv"
)
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
write.csv(manifest, out, row.names=FALSE)

cat("10r storage-grid patch: PASS\n")
cat("patched file:", path, "\n")
cat("storage allocations patched:", length(storage_pairs), "\n")
