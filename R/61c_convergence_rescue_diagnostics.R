#!/usr/bin/env Rscript

# =============================================================================
# 61c_convergence_rescue_diagnostics.R
#
# Long-chain convergence rescue diagnostics for a selected subset of countries.
# This is intentionally separate from R/61_formal_mcmc_diagnostics.R:
#   * R/61 expects all 14 countries.
#   * This script expects only FIN3_RESCUE_COUNTRIES (default US,CN,JP,CH,TR).
#
# It computes:
#   - rank-normalized split-Rhat
#   - bulk ESS
#   - tail ESS
#   - country summaries
#   - parameter-family summaries
#   - chain-level mean / sd / quantiles
#   - between-chain mean separation diagnostics
#   - worst-parameter ranking
#
# IMPORTANT:
# This script does NOT relax the formal convergence thresholds and does NOT
# generate IRFs. It is a diagnostic experiment asking whether longer chains
# alone are sufficient to resolve the formal convergence failure.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("posterior", quietly = TRUE)) {
  stopf("Package 'posterior' is required for convergence rescue diagnostics.")
}

get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("Environment variable %s is not numeric: %s", name, z)
  out
}

get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}

split_country_list <- function(x) {
  z <- toupper(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
  z <- unique(z[nzchar(z)])
  if (!length(z)) stopf("FIN3_RESCUE_COUNTRIES is empty.")
  z
}

parameter_family <- function(x) {
  if (grepl("__sv_", x, fixed = TRUE)) return("SV")
  if (grepl("__omega__", x, fixed = TRUE)) return("OMEGA")
  if (grepl("__threshold__", x, fixed = TRUE)) return("THRESHOLD")
  if (grepl("__V0__", x, fixed = TRUE)) return("V0")
  if (grepl("__coef__", x, fixed = TRUE)) return("EVENT_COEF")
  "OTHER"
}

PARTS_ROOT <- get_env_chr("FIN3_PARTS_ROOT", "posterior_parts")
RESCUE_COUNTRIES <- split_country_list(
  get_env_chr("FIN3_RESCUE_COUNTRIES", "US,CN,JP,CH,TR")
)
NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 3000))

RHAT_TARGET <- get_env_num("FIN3_RHAT_TARGET", 1.01)
RHAT_HARD <- get_env_num("FIN3_RHAT_HARD", 1.05)
ESS_TARGET <- get_env_num("FIN3_ESS_TARGET", 400)
ESS_HARD <- get_env_num("FIN3_ESS_HARD", 100)
MIN_RHAT_TARGET_SHARE <- get_env_num("FIN3_MIN_RHAT_TARGET_SHARE", 0.95)
MIN_ESS_TARGET_SHARE <- get_env_num("FIN3_MIN_ESS_TARGET_SHARE", 0.90)

if (!all(RESCUE_COUNTRIES %in% COUNTRIES)) {
  stopf(
    "Unknown rescue countries: %s",
    paste(setdiff(RESCUE_COUNTRIES, COUNTRIES), collapse = ", ")
  )
}
if (NCHAINS < 2L) stopf("At least two chains are required.")
if (EXPECTED_STORED < 200L) stopf("Expected stored draws is implausibly small.")

OUT <- file.path(RESULTS_DIR, "convergence_rescue")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(PARTS_ROOT)) {
  stopf("Posterior-parts root not found: %s", PARTS_ROOT)
}

files_all <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(files_all)) stopf("No formal TVP posterior RDS files found.")

parts_index_all <- do.call(rbind, lapply(files_all, function(f) {
  z <- readRDS(f)
  data.frame(
    File = f,
    Country = toupper(as.character(z$meta$country)),
    Chain = as.integer(z$meta$chain),
    StoredDraws = as.integer(z$meta$stored_draws),
    Burn = as.integer(z$meta$burn),
    KeepInternal = as.integer(z$meta$keep_internal),
    Seed = as.integer(z$meta$seed),
    Network = as.character(z$meta$main_network),
    stringsAsFactors = FALSE
  )
}))

parts_index <- parts_index_all[
  parts_index_all$Country %in% RESCUE_COUNTRIES,
  ,
  drop = FALSE
]

if (anyDuplicated(paste(parts_index$Country, parts_index$Chain, sep = "||"))) {
  stopf("Duplicate country-chain posterior parts detected.")
}

expected_grid <- expand.grid(
  Country = RESCUE_COUNTRIES,
  Chain = seq_len(NCHAINS),
  stringsAsFactors = FALSE
)
expected_keys <- paste(expected_grid$Country, expected_grid$Chain, sep = "||")
present_keys <- paste(parts_index$Country, parts_index$Chain, sep = "||")
missing_keys <- setdiff(expected_keys, present_keys)

if (length(missing_keys)) {
  write.csv(
    parts_index,
    file.path(OUT, "00_parts_index_partial.csv"),
    row.names = FALSE
  )
  stopf(
    "Missing %d rescue posterior part(s): %s",
    length(missing_keys),
    paste(missing_keys, collapse = ", ")
  )
}

if (nrow(parts_index) != length(RESCUE_COUNTRIES) * NCHAINS) {
  stopf(
    "Unexpected rescue part count: got %d, expected %d.",
    nrow(parts_index),
    length(RESCUE_COUNTRIES) * NCHAINS
  )
}

if (any(parts_index$StoredDraws != EXPECTED_STORED)) {
  stopf("At least one chain has an unexpected stored-draw count.")
}
if (!all(parts_index$Network == MAIN_NETWORK)) {
  stopf("At least one rescue chain uses a different financial network.")
}

parts_index <- parts_index[
  order(match(parts_index$Country, RESCUE_COUNTRIES), parts_index$Chain),
  ,
  drop = FALSE
]
write.csv(parts_index, file.path(OUT, "00_parts_index.csv"), row.names = FALSE)

diag_rows <- list()
chain_rows <- list()
rr <- 0L
cr <- 0L

for (cc in RESCUE_COUNTRIES) {
  cc_chains <- lapply(seq_len(NCHAINS), function(ch) {
    f <- parts_index$File[
      parts_index$Country == cc & parts_index$Chain == ch
    ]
    if (length(f) != 1L) stopf("Missing/duplicate %s chain %d.", cc, ch)
    z <- readRDS(f)
    m <- as.matrix(z$monitor)
    if (nrow(m) != EXPECTED_STORED) {
      stopf(
        "%s chain %d monitor rows=%d, expected=%d.",
        cc, ch, nrow(m), EXPECTED_STORED
      )
    }
    m
  })

  ref_names <- colnames(cc_chains[[1]])
  if (is.null(ref_names) || !length(ref_names)) {
    stopf("No monitored parameters for %s.", cc)
  }
  if (!all(vapply(
    cc_chains,
    function(x) identical(colnames(x), ref_names),
    logical(1)
  ))) {
    stopf("Monitor parameter names differ across chains for %s.", cc)
  }

  for (pp in seq_along(ref_names)) {
    par_name <- ref_names[pp]
    mat <- do.call(cbind, lapply(cc_chains, function(x) x[, pp]))

    finite <- all(is.finite(mat))
    chain_sd <- if (finite) apply(mat, 2, stats::sd) else rep(NA_real_, NCHAINS)
    constant <- finite && all(chain_sd < 1e-12)

    if (!finite || constant) {
      rh <- eb <- et <- NA_real_
      mean_range <- median_within_sd <- separation <- NA_real_
    } else {
      rh <- tryCatch(posterior::rhat(mat), error = function(e) NA_real_)
      eb <- tryCatch(posterior::ess_bulk(mat), error = function(e) NA_real_)
      et <- tryCatch(posterior::ess_tail(mat), error = function(e) NA_real_)

      chain_means <- colMeans(mat)
      mean_range <- diff(range(chain_means))
      median_within_sd <- stats::median(chain_sd)
      separation <- mean_range / max(median_within_sd, 1e-12)
    }

    rr <- rr + 1L
    diag_rows[[rr]] <- data.frame(
      Country = cc,
      Parameter = par_name,
      Family = parameter_family(par_name),
      Chains = NCHAINS,
      DrawsPerChain = EXPECTED_STORED,
      Finite = finite,
      ConstantAcrossAllChains = constant,
      Rhat = rh,
      ESS_Bulk = eb,
      ESS_Tail = et,
      BetweenChainMeanRange = mean_range,
      MedianWithinChainSD = median_within_sd,
      MeanRangeOverMedianWithinSD = separation,
      PassRhatTarget = if (is.finite(rh)) rh <= RHAT_TARGET else NA,
      PassRhatHard = if (is.finite(rh)) rh <= RHAT_HARD else NA,
      PassESSBulkTarget = if (is.finite(eb)) eb >= ESS_TARGET else NA,
      PassESSTailTarget = if (is.finite(et)) et >= ESS_TARGET else NA,
      PassESSBulkHard = if (is.finite(eb)) eb >= ESS_HARD else NA,
      PassESSTailHard = if (is.finite(et)) et >= ESS_HARD else NA,
      stringsAsFactors = FALSE
    )

    for (ch in seq_len(NCHAINS)) {
      vals <- mat[, ch]
      cr <- cr + 1L
      chain_rows[[cr]] <- data.frame(
        Country = cc,
        Parameter = par_name,
        Family = parameter_family(par_name),
        Chain = ch,
        Mean = if (all(is.finite(vals))) mean(vals) else NA_real_,
        SD = if (all(is.finite(vals))) stats::sd(vals) else NA_real_,
        Q05 = if (all(is.finite(vals))) unname(stats::quantile(vals, 0.05)) else NA_real_,
        Q50 = if (all(is.finite(vals))) unname(stats::quantile(vals, 0.50)) else NA_real_,
        Q95 = if (all(is.finite(vals))) unname(stats::quantile(vals, 0.95)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  rm(cc_chains)
  invisible(gc())
}

diagnostics <- do.call(rbind, diag_rows)
chain_summary <- do.call(rbind, chain_rows)

write.csv(
  diagnostics,
  file.path(OUT, "01_rescue_parameter_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  chain_summary,
  file.path(OUT, "05_chain_parameter_summary.csv"),
  row.names = FALSE
)

active <- diagnostics[
  diagnostics$Finite &
    !diagnostics$ConstantAcrossAllChains &
    is.finite(diagnostics$Rhat) &
    is.finite(diagnostics$ESS_Bulk) &
    is.finite(diagnostics$ESS_Tail),
  ,
  drop = FALSE
]
if (!nrow(active)) stopf("No active convergence diagnostics available.")

country_rows <- lapply(RESCUE_COUNTRIES, function(cc) {
  x <- active[active$Country == cc, , drop = FALSE]
  data.frame(
    Country = cc,
    ActiveDiagnostics = nrow(x),
    MaxRhat = max(x$Rhat),
    RhatTargetShare = mean(x$Rhat <= RHAT_TARGET),
    RhatHardShare = mean(x$Rhat <= RHAT_HARD),
    MinESSBulk = min(x$ESS_Bulk),
    MinESSTail = min(x$ESS_Tail),
    ESSBulkTargetShare = mean(x$ESS_Bulk >= ESS_TARGET),
    ESSTailTargetShare = mean(x$ESS_Tail >= ESS_TARGET),
    MaxMeanSeparationRatio = max(x$MeanRangeOverMedianWithinSD, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
country_summary <- do.call(rbind, country_rows)
write.csv(
  country_summary,
  file.path(OUT, "02_rescue_country_summary.csv"),
  row.names = FALSE
)

families <- unique(active$Family)
family_rows <- lapply(families, function(ff) {
  x <- active[active$Family == ff, , drop = FALSE]
  data.frame(
    Family = ff,
    ActiveDiagnostics = nrow(x),
    MaxRhat = max(x$Rhat),
    MedianRhat = stats::median(x$Rhat),
    RhatTargetShare = mean(x$Rhat <= RHAT_TARGET),
    RhatHardShare = mean(x$Rhat <= RHAT_HARD),
    MinESSBulk = min(x$ESS_Bulk),
    MinESSTail = min(x$ESS_Tail),
    ESSBulkTargetShare = mean(x$ESS_Bulk >= ESS_TARGET),
    ESSTailTargetShare = mean(x$ESS_Tail >= ESS_TARGET),
    MedianMeanSeparationRatio = stats::median(
      x$MeanRangeOverMedianWithinSD,
      na.rm = TRUE
    ),
    MaxMeanSeparationRatio = max(
      x$MeanRangeOverMedianWithinSD,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
})
family_summary <- do.call(rbind, family_rows)
family_summary <- family_summary[
  order(family_summary$MaxRhat, decreasing = TRUE),
  ,
  drop = FALSE
]
write.csv(
  family_summary,
  file.path(OUT, "03_rescue_family_summary.csv"),
  row.names = FALSE
)

worst <- active[
  order(
    -active$Rhat,
    active$ESS_Bulk,
    active$ESS_Tail,
    -active$MeanRangeOverMedianWithinSD
  ),
  ,
  drop = FALSE
]
worst <- head(worst, 100L)
write.csv(
  worst,
  file.path(OUT, "06_worst_parameters.csv"),
  row.names = FALSE
)

rhat_hard_ok <- all(active$Rhat <= RHAT_HARD)
rhat_target_share <- mean(active$Rhat <= RHAT_TARGET)
rhat_target_ok <- rhat_target_share >= MIN_RHAT_TARGET_SHARE

ess_bulk_hard_ok <- all(active$ESS_Bulk >= ESS_HARD)
ess_tail_hard_ok <- all(active$ESS_Tail >= ESS_HARD)
ess_bulk_share <- mean(active$ESS_Bulk >= ESS_TARGET)
ess_tail_share <- mean(active$ESS_Tail >= ESS_TARGET)
ess_target_ok <- (
  ess_bulk_share >= MIN_ESS_TARGET_SHARE &&
  ess_tail_share >= MIN_ESS_TARGET_SHARE
)

ready <- (
  rhat_hard_ok &&
  rhat_target_ok &&
  ess_bulk_hard_ok &&
  ess_tail_hard_ok &&
  ess_target_ok
)

hard_fail <- (
  !rhat_hard_ok ||
  !ess_bulk_hard_ok ||
  !ess_tail_hard_ok
)

status <- if (ready) {
  "LONG_CHAIN_RESCUE_PASSED"
} else if (hard_fail) {
  "HARD_CONVERGENCE_FAILURE_PERSISTS"
} else {
  "HARD_LIMITS_PASS_SOFT_TARGETS_FAIL"
}

decision <- if (ready) {
  "LONGER_CHAINS_ARE_SUFFICIENT_FOR_THIS_RESCUE_SUBSET"
} else if (hard_fail) {
  paste0(
    "CHAIN_LENGTH_ALONE_NOT_YET_SUFFICIENT; ",
    "REVIEW_THRESHOLD_TVP_SHRINKAGE_OR_PARAMETER-SPECIFIC_MIXING_BEFORE_FULL_RERUN"
  )
} else {
  paste0(
    "LONG_CHAINS_REMOVE_HARD_FAILURES_BUT_PUBLICATION_TARGETS_ARE_NOT_YET_MET; ",
    "REVIEW_WORST_PARAMETERS_BEFORE_FULL_RERUN"
  )
}

reasons <- c(
  if (!rhat_hard_ok) sprintf("at least one Rhat > %.3f", RHAT_HARD) else NULL,
  if (!rhat_target_ok) sprintf(
    "Rhat<=%.3f share %.4f < required %.4f",
    RHAT_TARGET, rhat_target_share, MIN_RHAT_TARGET_SHARE
  ) else NULL,
  if (!ess_bulk_hard_ok) sprintf("at least one bulk ESS < %.0f", ESS_HARD) else NULL,
  if (!ess_tail_hard_ok) sprintf("at least one tail ESS < %.0f", ESS_HARD) else NULL,
  if (!ess_target_ok) sprintf(
    "ESS target shares bulk/tail %.4f/%.4f below required %.4f",
    ess_bulk_share, ess_tail_share, MIN_ESS_TARGET_SHARE
  ) else NULL
)
if (!length(reasons)) reasons <- "all rescue convergence gates passed"

gate <- data.frame(
  Status = status,
  Decision = decision,
  Countries = paste(RESCUE_COUNTRIES, collapse = ","),
  CountryCount = length(RESCUE_COUNTRIES),
  ChainsPerCountry = NCHAINS,
  BurnPerChain = unique(parts_index$Burn)[1],
  InternalPostBurnIterationsPerChain = unique(parts_index$KeepInternal)[1],
  StoredDrawsPerChain = EXPECTED_STORED,
  ActiveDiagnostics = nrow(active),
  MaxRhat = max(active$Rhat),
  RhatTarget = RHAT_TARGET,
  RhatHardLimit = RHAT_HARD,
  RhatTargetShare = rhat_target_share,
  MinRhatTargetShareRequired = MIN_RHAT_TARGET_SHARE,
  MinESSBulk = min(active$ESS_Bulk),
  MinESSTail = min(active$ESS_Tail),
  ESSTarget = ESS_TARGET,
  ESSHardMinimum = ESS_HARD,
  ESSBulkTargetShare = ess_bulk_share,
  ESSTailTargetShare = ess_tail_share,
  MinESSTargetShareRequired = MIN_ESS_TARGET_SHARE,
  MaxMeanSeparationRatio = max(
    active$MeanRangeOverMedianWithinSD,
    na.rm = TRUE
  ),
  Reason = paste(reasons, collapse = "; "),
  stringsAsFactors = FALSE
)

write.csv(
  gate,
  file.path(OUT, "00_rescue_convergence_gate.csv"),
  row.names = FALSE
)

overall <- data.frame(
  Metric = c(
    "Countries",
    "ChainsPerCountry",
    "StoredDrawsPerChain",
    "ActiveDiagnostics",
    "MaxRhat",
    "ShareRhatAtOrBelowTarget",
    "ShareRhatAtOrBelowHardLimit",
    "MinESSBulk",
    "MinESSTail",
    "ShareESSBulkAtOrAboveTarget",
    "ShareESSTailAtOrAboveTarget",
    "ShareESSBulkAtOrAboveHardMinimum",
    "ShareESSTailAtOrAboveHardMinimum",
    "MaxMeanSeparationRatio"
  ),
  Value = c(
    length(RESCUE_COUNTRIES),
    NCHAINS,
    EXPECTED_STORED,
    nrow(active),
    max(active$Rhat),
    mean(active$Rhat <= RHAT_TARGET),
    mean(active$Rhat <= RHAT_HARD),
    min(active$ESS_Bulk),
    min(active$ESS_Tail),
    mean(active$ESS_Bulk >= ESS_TARGET),
    mean(active$ESS_Tail >= ESS_TARGET),
    mean(active$ESS_Bulk >= ESS_HARD),
    mean(active$ESS_Tail >= ESS_HARD),
    max(active$MeanRangeOverMedianWithinSD, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  overall,
  file.path(OUT, "04_rescue_overall_summary.csv"),
  row.names = FALSE
)

readme <- c(
  sprintf("CONVERGENCE RESCUE: %s", status),
  "==============================================",
  sprintf("Countries: %s", paste(RESCUE_COUNTRIES, collapse = ", ")),
  sprintf("Chains per country: %d", NCHAINS),
  sprintf("Burn per chain: %d", gate$BurnPerChain),
  sprintf(
    "Internal post-burn iterations per chain: %d",
    gate$InternalPostBurnIterationsPerChain
  ),
  sprintf("Stored draws per chain: %d", EXPECTED_STORED),
  sprintf("Max rank-normalized split-Rhat: %.6f", gate$MaxRhat),
  sprintf(
    "Rhat <= %.3f share: %.6f",
    RHAT_TARGET,
    gate$RhatTargetShare
  ),
  sprintf("Min bulk ESS: %.2f", gate$MinESSBulk),
  sprintf("Min tail ESS: %.2f", gate$MinESSTail),
  sprintf(
    "Bulk/Tail ESS >= %.0f shares: %.6f / %.6f",
    ESS_TARGET,
    gate$ESSBulkTargetShare,
    gate$ESSTailTargetShare
  ),
  sprintf(
    "Maximum between-chain mean-range / median within-chain SD: %.4f",
    gate$MaxMeanSeparationRatio
  ),
  sprintf("Decision: %s", decision),
  sprintf("Reason: %s", gate$Reason),
  "",
  "This rescue experiment preserves the formal model, priors, kappa0,",
  "financial network, variable transformations, event design, and sampler.",
  "Only chain length and stored draw count are increased.",
  "",
  "Interpretation:",
  "- LONG_CHAIN_RESCUE_PASSED: longer chains are sufficient for this subset.",
  "- HARD_LIMITS_PASS_SOFT_TARGETS_FAIL: hard failures disappear, but publication targets remain unmet.",
  "- HARD_CONVERGENCE_FAILURE_PERSISTS: do not run all 14 countries with longer chains yet;",
  "  inspect 03_rescue_family_summary.csv and 06_worst_parameters.csv first."
)

writeLines(readme, file.path(OUT, "README_convergence_rescue.txt"))
cat(paste(readme, collapse = "\n"), "\n")
