#!/usr/bin/env Rscript

# =============================================================================
# 66_no_tr_independent_seed_mode_audit.R
#
# Independent-seed mode audit for the Delta-r formal TVP-GVAR.
#
# PURPOSE
# -------
# 08k showed that most of the Delta-r posterior converged after targeted
# long-chain rescue, but NO and TR retained hard convergence failures.
#
# 08l re-estimates ALL FOUR chains for NO and TR from a NEW seed base while
# keeping the statistical specification unchanged.  This script then asks:
#
#   1) Do the new NO/TR chains satisfy the hard R61 thresholds?
#   2) Do chain medians still separate into visibly different posterior modes?
#   3) Did the country-level convergence diagnostics improve relative to 08k?
#
# This is a convergence / mode diagnostic only.  It does not select or discard
# individual chains.  A country passes only as a complete four-chain set.
# =============================================================================

source("R/00_config.R")

get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}

get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("Environment variable %s is not numeric: %s", name, z)
  out
}

PARTS_ROOT <- get_env_chr("FIN3_PARTS_ROOT", "posterior_parts")
NEW_DIAG_DIR <- get_env_chr("FIN3_CONVERGENCE_DIR", file.path(RESULTS_DIR, "formal_tvp"))
PRIOR_DIAG_DIR <- get_env_chr("FIN3_PRIOR_CONVERGENCE_DIR", "prior_08k_convergence")
TARGETS <- c("NO", "TR")

NCHAINS <- 4L
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 2000))
RHAT_HARD <- get_env_num("FIN3_RHAT_HARD", 1.05)
ESS_HARD <- get_env_num("FIN3_ESS_HARD", 100)
SEPARATION_RATIO_THRESHOLD <- get_env_num("FIN3_MODE_SEPARATION_RATIO", 4)

OUT <- file.path(RESULTS_DIR, "independent_seed_mode_audit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

new_param_file <- file.path(NEW_DIAG_DIR, "01_mcmc_parameter_diagnostics.csv")
new_country_file <- file.path(NEW_DIAG_DIR, "02_mcmc_country_summary.csv")

if (!file.exists(new_param_file)) stopf("Missing new R61 parameter diagnostics: %s", new_param_file)
if (!file.exists(new_country_file)) stopf("Missing new R61 country summary: %s", new_country_file)

new_param <- read.csv(new_param_file, stringsAsFactors = FALSE, check.names = FALSE)
new_country <- read.csv(new_country_file, stringsAsFactors = FALSE, check.names = FALSE)

need_param <- c("Country","Parameter","Finite","ConstantAcrossAllChains","Rhat","ESS_Bulk","ESS_Tail")
if (!all(need_param %in% names(new_param))) {
  stopf("Malformed new R61 parameter diagnostics.")
}

need_country <- c("Country","MaxRhat","RhatTargetShare","MinESSBulk","MinESSTail")
if (!all(need_country %in% names(new_country))) {
  stopf("Malformed new R61 country summary.")
}

if (!all(TARGETS %in% new_country$Country)) stopf("New R61 summary is missing NO/TR.")

# -----------------------------------------------------------------------------
# 1. Read the complete four-chain sets for NO and TR
# -----------------------------------------------------------------------------

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+[.]rds$",
  recursive = TRUE,
  full.names = TRUE
)

target_files <- files[vapply(files, function(f) {
  z <- readRDS(f)
  as.character(z$meta$country) %in% TARGETS
}, logical(1))]

if (length(target_files) != length(TARGETS) * NCHAINS) {
  stopf("Expected exactly 8 NO/TR posterior RDS files; found %d.", length(target_files))
}

parts <- setNames(vector("list", length(TARGETS)), TARGETS)
index_rows <- list()
ii <- 0L

for (cc in TARGETS) {
  parts[[cc]] <- vector("list", NCHAINS)
  for (ch in seq_len(NCHAINS)) {
    cand <- target_files[vapply(target_files, function(f) {
      z <- readRDS(f)
      identical(as.character(z$meta$country), cc) &&
        as.integer(z$meta$chain) == ch
    }, logical(1))]
    if (length(cand) != 1L) stopf("Expected one %s chain %d posterior; found %d.", cc, ch, length(cand))
    z <- readRDS(cand)
    if (as.integer(z$meta$stored_draws) != EXPECTED_STORED) {
      stopf("%s chain %d stored draw mismatch.", cc, ch)
    }
    if (basename(as.character(z$meta$source_panel)) != "panel_domestic_fin3_rate_diff.csv") {
      stopf("%s chain %d is not from the Delta-r panel.", cc, ch)
    }
    parts[[cc]][[ch]] <- z

    ii <- ii + 1L
    index_rows[[ii]] <- data.frame(
      Country = cc,
      Chain = ch,
      Seed = as.integer(z$meta$seed),
      Burn = as.integer(z$meta$burn),
      KeepInternal = as.integer(z$meta$keep_internal),
      StoredDraws = as.integer(z$meta$stored_draws),
      SourcePanel = basename(as.character(z$meta$source_panel)),
      File = cand,
      stringsAsFactors = FALSE
    )
  }
}

index <- do.call(rbind, index_rows)
write.csv(index, file.path(OUT, "00_no_tr_new_seed_posterior_index.csv"), row.names = FALSE)

# New independent-seed chains must not accidentally collapse onto one seed.
if (anyDuplicated(index$Seed)) {
  stopf("Duplicate actual MCMC seeds detected within the NO/TR independent-seed audit.")
}

# -----------------------------------------------------------------------------
# 2. Chain-location separation audit
# -----------------------------------------------------------------------------

safe_mad <- function(x) {
  out <- stats::mad(x, center = stats::median(x), constant = 1, na.rm = TRUE)
  if (!is.finite(out)) NA_real_ else out
}

detail_rows <- list()
rr <- 0L

for (cc in TARGETS) {
  mons <- lapply(parts[[cc]], function(z) as.matrix(z$monitor))
  ref_names <- colnames(mons[[1]])

  if (is.null(ref_names) || !length(ref_names)) stopf("No monitor parameters for %s.", cc)
  if (!all(vapply(mons, function(x) identical(colnames(x), ref_names), logical(1)))) {
    stopf("Monitor parameter ordering differs across new-seed %s chains.", cc)
  }

  for (pp in seq_along(ref_names)) {
    pname <- ref_names[pp]
    vals <- lapply(mons, function(x) as.numeric(x[, pp]))

    med <- vapply(vals, stats::median, numeric(1), na.rm = TRUE)
    mn <- vapply(vals, mean, numeric(1), na.rm = TRUE)
    sdv <- vapply(vals, stats::sd, numeric(1), na.rm = TRUE)
    md <- vapply(vals, safe_mad, numeric(1))

    median_range <- max(med) - min(med)
    within_mad <- stats::median(md[is.finite(md) & md > 1e-12], na.rm = TRUE)
    if (!is.finite(within_mad)) within_mad <- NA_real_

    sep_ratio <- if (is.finite(within_mad) && within_mad > 0) {
      abs(median_range) / within_mad
    } else {
      NA_real_
    }

    diag_row <- new_param[
      new_param$Country == cc & new_param$Parameter == pname,
      ,
      drop = FALSE
    ]
    if (nrow(diag_row) != 1L) {
      stopf("Could not match R61 diagnostic for %s / %s.", cc, pname)
    }

    hard_issue <- (
      is.finite(diag_row$Rhat[1]) && diag_row$Rhat[1] > RHAT_HARD
    ) || (
      is.finite(diag_row$ESS_Bulk[1]) && diag_row$ESS_Bulk[1] < ESS_HARD
    ) || (
      is.finite(diag_row$ESS_Tail[1]) && diag_row$ESS_Tail[1] < ESS_HARD
    )

    strong_sep <- is.finite(sep_ratio) && sep_ratio >= SEPARATION_RATIO_THRESHOLD

    rr <- rr + 1L
    detail_rows[[rr]] <- data.frame(
      Country = cc,
      Parameter = pname,
      Chain1Median = med[1],
      Chain2Median = med[2],
      Chain3Median = med[3],
      Chain4Median = med[4],
      Chain1Mean = mn[1],
      Chain2Mean = mn[2],
      Chain3Mean = mn[3],
      Chain4Mean = mn[4],
      Chain1SD = sdv[1],
      Chain2SD = sdv[2],
      Chain3SD = sdv[3],
      Chain4SD = sdv[4],
      ChainMedianRange = median_range,
      MedianWithinChainMAD = within_mad,
      ChainMedianSeparationRatio = sep_ratio,
      Rhat = diag_row$Rhat[1],
      ESS_Bulk = diag_row$ESS_Bulk[1],
      ESS_Tail = diag_row$ESS_Tail[1],
      HardMCMCIssue = hard_issue,
      StrongChainLocationSeparation = strong_sep,
      ModeSplitSuspect = hard_issue && strong_sep,
      stringsAsFactors = FALSE
    )
  }
}

detail <- do.call(rbind, detail_rows)
write.csv(detail, file.path(OUT, "01_no_tr_parameter_mode_audit.csv"), row.names = FALSE)

# Rank parameters for human inspection.
ranked <- detail[
  order(
    detail$Country,
    -as.integer(detail$ModeSplitSuspect),
    -ifelse(is.finite(detail$Rhat), detail$Rhat, -Inf),
    -ifelse(is.finite(detail$ChainMedianSeparationRatio),
            detail$ChainMedianSeparationRatio, -Inf)
  ),
  ,
  drop = FALSE
]
write.csv(ranked, file.path(OUT, "02_no_tr_ranked_mode_suspects.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Country-level before/after comparison with 08k
# -----------------------------------------------------------------------------

prior_country_file <- file.path(PRIOR_DIAG_DIR, "02_mcmc_country_summary.csv")
prior_available <- file.exists(prior_country_file)

if (prior_available) {
  prior_country <- read.csv(prior_country_file, stringsAsFactors = FALSE, check.names = FALSE)
  prior_country <- prior_country[prior_country$Country %in% TARGETS, , drop = FALSE]
} else {
  prior_country <- data.frame()
}

country_rows <- lapply(TARGETS, function(cc) {
  n <- new_country[new_country$Country == cc, , drop = FALSE]
  d <- detail[detail$Country == cc, , drop = FALSE]
  p <- if (nrow(prior_country)) prior_country[prior_country$Country == cc, , drop = FALSE] else data.frame()

  new_hard_issue <- any(
    (is.finite(d$Rhat) & d$Rhat > RHAT_HARD) |
    (is.finite(d$ESS_Bulk) & d$ESS_Bulk < ESS_HARD) |
    (is.finite(d$ESS_Tail) & d$ESS_Tail < ESS_HARD)
  )
  mode_suspects <- sum(d$ModeSplitSuspect, na.rm = TRUE)

  data.frame(
    Country = cc,
    Prior08kAvailable = nrow(p) == 1L,
    Prior08kMaxRhat = if (nrow(p) == 1L) p$MaxRhat[1] else NA_real_,
    New08lMaxRhat = n$MaxRhat[1],
    DeltaMaxRhat_NewMinusPrior = if (nrow(p) == 1L) n$MaxRhat[1] - p$MaxRhat[1] else NA_real_,
    Prior08kMinESSBulk = if (nrow(p) == 1L) p$MinESSBulk[1] else NA_real_,
    New08lMinESSBulk = n$MinESSBulk[1],
    Prior08kMinESSTail = if (nrow(p) == 1L) p$MinESSTail[1] else NA_real_,
    New08lMinESSTail = n$MinESSTail[1],
    NewHardMCMCIssue = new_hard_issue,
    NewModeSplitSuspectParameters = mode_suspects,
    IndependentSeedModeAuditPass = !new_hard_issue && mode_suspects == 0L,
    stringsAsFactors = FALSE
  )
})

country_cmp <- do.call(rbind, country_rows)
write.csv(country_cmp, file.path(OUT, "03_no_tr_country_before_after.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Gate
# -----------------------------------------------------------------------------

all_pass <- all(country_cmp$IndependentSeedModeAuditPass)

status <- if (all_pass) {
  "INDEPENDENT_SEED_MODE_AUDIT_PASS"
} else {
  "PERSISTENT_MODE_OR_MIXING_PROBLEM"
}

gate <- data.frame(
  Status = status,
  Countries = paste(TARGETS, collapse = ","),
  ChainsPerCountry = NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  RhatHardLimit = RHAT_HARD,
  ESSHardLimit = ESS_HARD,
  SeparationRatioThreshold = SEPARATION_RATIO_THRESHOLD,
  NO_MaxRhat = country_cmp$New08lMaxRhat[country_cmp$Country == "NO"],
  NO_MinESSBulk = country_cmp$New08lMinESSBulk[country_cmp$Country == "NO"],
  NO_MinESSTail = country_cmp$New08lMinESSTail[country_cmp$Country == "NO"],
  NO_ModeSplitSuspects = country_cmp$NewModeSplitSuspectParameters[country_cmp$Country == "NO"],
  TR_MaxRhat = country_cmp$New08lMaxRhat[country_cmp$Country == "TR"],
  TR_MinESSBulk = country_cmp$New08lMinESSBulk[country_cmp$Country == "TR"],
  TR_MinESSTail = country_cmp$New08lMinESSTail[country_cmp$Country == "TR"],
  TR_ModeSplitSuspects = country_cmp$NewModeSplitSuspectParameters[country_cmp$Country == "TR"],
  IndividualChainsDiscarded = FALSE,
  SpecificationChanged = FALSE,
  stringsAsFactors = FALSE
)

write.csv(gate, file.path(OUT, "00_independent_seed_mode_audit_gate.csv"), row.names = FALSE)

readme <- c(
  sprintf("08l NO/TR INDEPENDENT-SEED MODE AUDIT: %s", status),
  "============================================================",
  "",
  "This audit never discards an individual chain.",
  "NO and TR are each evaluated as complete four-chain independent-seed sets.",
  "",
  sprintf("Rhat hard limit: %.3f", RHAT_HARD),
  sprintf("ESS hard limit: %.0f", ESS_HARD),
  sprintf("Strong chain-location separation threshold: %.2f within-chain MAD units", SEPARATION_RATIO_THRESHOLD),
  "",
  sprintf(
    "NO: MaxRhat=%.6f; MinBulkESS=%.2f; MinTailESS=%.2f; mode suspects=%d",
    gate$NO_MaxRhat, gate$NO_MinESSBulk, gate$NO_MinESSTail, gate$NO_ModeSplitSuspects
  ),
  sprintf(
    "TR: MaxRhat=%.6f; MinBulkESS=%.2f; MinTailESS=%.2f; mode suspects=%d",
    gate$TR_MaxRhat, gate$TR_MinESSBulk, gate$TR_MinESSTail, gate$TR_ModeSplitSuspects
  ),
  "",
  if (all_pass) {
    paste0(
      "Interpretation: under a completely new seed base, neither NO nor TR retains ",
      "a hard R61 failure or a hard-failure-plus-chain-location split. This is ",
      "evidence against the previous bad chain being a persistent posterior mode."
    )
  } else {
    paste0(
      "Interpretation: at least one of NO/TR still shows hard convergence failure ",
      "under the independent seed base. Do not cherry-pick chains and do not proceed ",
      "to formal IRFs; the country/equation prior or latent-threshold structure needs ",
      "a dedicated diagnostic."
    )
  }
)

writeLines(readme, file.path(OUT, "README_08l_independent_seed_mode_audit.txt"))
cat(paste(readme, collapse = "\n"), "\n")

if (!all_pass) {
  print(country_cmp)
  bad <- ranked[ranked$HardMCMCIssue | ranked$ModeSplitSuspect, , drop = FALSE]
  if (nrow(bad)) print(utils::head(bad, 30))
  stop("08l independent-seed mode audit failed.")
}
