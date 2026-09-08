#!/usr/bin/env Rscript

# =============================================================================
# 73_ttvp_storage_grid_offbyone_audit.R
#
# Diagnostic-only audit of the original pinned threshtvp storage arithmetic.
#
# No MCMC is run. No model, prior, data, or formal result is changed.
#
# The pinned MCMC_tvp source uses the same floating-point quantity `thin*nsave`
# in two different contexts:
#
#   thin_out <- round(seq(nburn, ntot, length.out = thin*nsave))
#   H_store  <- matrix(NA, thin*nsave, T)
#   ... analogous posterior storage objects ...
#
# R's dimension values are integer-valued, whereas seq() rounds a fractional
# length.out upward.  This audit reproduces those exact arithmetic operations
# for the successful 10p short contract and failed 10o long contract.
# =============================================================================

dir.create("results/10q_storage_grid_audit", recursive=TRUE, showWarnings=FALSE)
OUT <- "results/10q_storage_grid_audit"

contracts <- list(
  `10p_short` = list(
    burn = 1000L,
    nsave = 3000L,
    stored = 300L
  ),
  `10o_long` = list(
    burn = 5000L,
    nsave = 15000L,
    stored = 1000L
  )
)

audit_one <- function(name, burn, nsave, stored) {
  thin <- (stored + 1L) / nsave
  raw_target <- thin * nsave
  ntot <- burn + nsave

  # Reproduce the original source's dimension arithmetic with a tiny T=2
  # placeholder, so no large objects are allocated.
  H_store <- matrix(NA_real_, raw_target, 2L)
  A_store <- array(NA_real_, c(raw_target, 2L, 2L))

  # Reproduce the exact original storage schedule expression.
  thin_out <- round(seq(
    burn,
    ntot,
    length.out = raw_target
  ))

  # Simulate only the storage counter, not MCMC.
  ithin <- 0L
  first_overflow_iteration <- NA_integer_
  first_overflow_storage_index <- NA_integer_

  for (irep in seq_len(ntot)) {
    if (irep %in% thin_out) {
      ithin <- ithin + 1L
      if (is.na(first_overflow_iteration) && ithin > nrow(H_store)) {
        first_overflow_iteration <- irep
        first_overflow_storage_index <- ithin
      }
    }
  }

  # Integer-safe counterpart: define storage count once as an integer and use
  # it consistently for both sequence and storage dimensions.
  nstored_safe <- stored + 1L
  thin_out_safe <- round(seq(
    burn,
    ntot,
    length.out = nstored_safe
  ))
  H_safe <- matrix(NA_real_, nstored_safe, 2L)

  data.frame(
    Contract = name,
    Burn = burn,
    NSave = nsave,
    RequestedStrictPostBurn = stored,
    IntendedRawStored = stored + 1L,
    Thin = sprintf("%.17g", thin),
    ThinTimesNSave = sprintf("%.17g", raw_target),
    OriginalAllocatedRows_H = nrow(H_store),
    OriginalAllocatedRows_A = dim(A_store)[1],
    OriginalThinOutLength = length(thin_out),
    OriginalThinOutUniqueLength = length(unique(thin_out)),
    OriginalStorageWrites = ithin,
    OriginalOverflow = !is.na(first_overflow_iteration),
    FirstOverflowIteration = first_overflow_iteration,
    FirstOverflowStorageIndex = first_overflow_storage_index,
    SafeAllocatedRows = nrow(H_safe),
    SafeThinOutLength = length(thin_out_safe),
    SafeOverflow = length(thin_out_safe) > nrow(H_safe),
    stringsAsFactors=FALSE
  )
}

rows <- do.call(
  rbind,
  lapply(names(contracts), function(nm) {
    z <- contracts[[nm]]
    audit_one(nm, z$burn, z$nsave, z$stored)
  })
)

write.csv(
  rows,
  file.path(OUT, "01_storage_contract_comparison.csv"),
  row.names=FALSE
)

short <- rows[rows$Contract == "10p_short", , drop=FALSE]
long  <- rows[rows$Contract == "10o_long", , drop=FALSE]

if (nrow(short) != 1L || nrow(long) != 1L) {
  stop("Internal 10q contract lookup failed.")
}

if (!isTRUE(short$OriginalOverflow) && isTRUE(long$OriginalOverflow)) {
  decision <- "FLOATING_LENGTH_OUT_STORAGE_OVERFLOW_PROVEN_FOR_10O_CONTRACT"
} else if (isTRUE(short$OriginalOverflow) && isTRUE(long$OriginalOverflow)) {
  decision <- "STORAGE_OVERFLOW_REPRODUCES_IN_BOTH_CONTRACTS__REASSESS"
} else if (!isTRUE(short$OriginalOverflow) && !isTRUE(long$OriginalOverflow)) {
  decision <- "STORAGE_OVERFLOW_NOT_REPRODUCED__PROCEED_TO_SEED_BY_LENGTH_DIAGNOSTIC"
} else {
  decision <- "SHORT_ONLY_STORAGE_OVERFLOW__UNEXPECTED__INSPECT_R_RUNTIME"
}

decision_rows <- data.frame(
  Item = c(
    "TechnicalStatus",
    "RVersion",
    "ShortContractOriginalOverflow",
    "LongContractOriginalOverflow",
    "LongContractAllocatedRows",
    "LongContractStorageWrites",
    "LongContractFirstOverflowIteration",
    "LongContractFirstOverflowStorageIndex",
    "IntegerSafeLongOverflow",
    "Decision",
    "FormalModelChanged",
    "IRFAuthorized"
  ),
  Value = c(
    "PASS",
    R.version.string,
    as.character(short$OriginalOverflow),
    as.character(long$OriginalOverflow),
    as.character(long$OriginalAllocatedRows_H),
    as.character(long$OriginalStorageWrites),
    as.character(long$FirstOverflowIteration),
    as.character(long$FirstOverflowStorageIndex),
    as.character(long$SafeOverflow),
    decision,
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors=FALSE
)

write.csv(
  decision_rows,
  file.path(OUT, "02_decision.csv"),
  row.names=FALSE
)

cat("\n===== 10q STORAGE GRID AUDIT =====\n")
print(rows, row.names=FALSE)
cat("\n===== DECISION =====\n")
print(decision_rows, row.names=FALSE)

# Technical gates only. A scientifically unfavorable result is still a valid
# diagnostic and must not be converted into a workflow failure.
stopifnot(
  identical(decision_rows$Value[decision_rows$Item=="TechnicalStatus"], "PASS"),
  identical(decision_rows$Value[decision_rows$Item=="FormalModelChanged"], "FALSE"),
  identical(decision_rows$Value[decision_rows$Item=="IRFAuthorized"], "FALSE")
)
