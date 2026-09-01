#!/usr/bin/env Rscript

# =============================================================================
# 61_formal_mcmc_diagnostics.R
#
# Aggregate all 14 countries x multiple formal TVP MCMC chains.
# Compute modern rank-normalized split-Rhat, bulk ESS and tail ESS using
# the posterior package.
#
# This script does NOT compute IRFs. It gates MCMC convergence first.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("posterior", quietly = TRUE)) {
  stopf("Package 'posterior' is required for formal MCMC diagnostics.")
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

PARTS_ROOT <- get_env_chr("FIN3_PARTS_ROOT", "posterior_parts")
NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 1500))

RHAT_TARGET <- get_env_num("FIN3_RHAT_TARGET", 1.01)
RHAT_HARD <- get_env_num("FIN3_RHAT_HARD", 1.05)
ESS_TARGET <- get_env_num("FIN3_ESS_TARGET", 400)
ESS_HARD <- get_env_num("FIN3_ESS_HARD", 100)

MIN_RHAT_TARGET_SHARE <- get_env_num("FIN3_MIN_RHAT_TARGET_SHARE", 0.95)
MIN_ESS_TARGET_SHARE <- get_env_num("FIN3_MIN_ESS_TARGET_SHARE", 0.90)

OUT <- file.path(RESULTS_DIR, "formal_tvp")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(PARTS_ROOT)) stopf("Posterior-parts root not found: %s", PARTS_ROOT)

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(files)) stopf("No formal TVP country-chain RDS files found.")

parts_index <- do.call(rbind, lapply(files, function(f) {
  z <- readRDS(f)
  data.frame(
    File = f,
    Country = as.character(z$meta$country),
    Chain = as.integer(z$meta$chain),
    StoredDraws = as.integer(z$meta$stored_draws),
    Burn = as.integer(z$meta$burn),
    KeepInternal = as.integer(z$meta$keep_internal),
    Seed = as.integer(z$meta$seed),
    Network = as.character(z$meta$main_network),
    stringsAsFactors = FALSE
  )
}))

if (anyDuplicated(paste(parts_index$Country, parts_index$Chain, sep = "||"))) {
  stopf("Duplicate country-chain posterior parts detected.")
}

expected_grid <- expand.grid(
  Country = COUNTRIES,
  Chain = seq_len(NCHAINS),
  stringsAsFactors = FALSE
)
check_grid <- merge(
  expected_grid,
  parts_index[, c("Country", "Chain")],
  by = c("Country", "Chain"),
  all.x = TRUE
)
# merge does not directly reveal unmatched because the keys remain; use explicit.
present_keys <- paste(parts_index$Country, parts_index$Chain, sep = "||")
expected_keys <- paste(expected_grid$Country, expected_grid$Chain, sep = "||")
missing_keys <- setdiff(expected_keys, present_keys)

if (length(missing_keys)) {
  write.csv(parts_index, file.path(OUT, "00_parts_index_partial.csv"), row.names = FALSE)
  stopf(
    "Missing %d formal posterior part(s): %s",
    length(missing_keys),
    paste(missing_keys, collapse = ", ")
  )
}

if (any(parts_index$StoredDraws != EXPECTED_STORED)) {
  stopf(
    "Stored draw count differs from expected %d.",
    EXPECTED_STORED
  )
}

if (!all(parts_index$Network == MAIN_NETWORK)) {
  stopf("At least one formal part uses a different financial network.")
}

write.csv(
  parts_index[order(match(parts_index$Country, COUNTRIES), parts_index$Chain), ],
  file.path(OUT, "00_parts_index.csv"),
  row.names = FALSE
)

# =============================================================================
# 1. Parameter-level R-hat / ESS
# =============================================================================

diag_rows <- list()
rr <- 0L

country_summary_rows <- list()
cc_pos <- 0L

for (cc in COUNTRIES) {
  cc_files <- parts_index$File[parts_index$Country == cc]
  cc_chains <- lapply(seq_len(NCHAINS), function(ch) {
    f <- parts_index$File[
      parts_index$Country == cc & parts_index$Chain == ch
    ]
    if (length(f) != 1L) stopf("Missing/duplicate %s chain %d", cc, ch)
    z <- readRDS(f)
    as.matrix(z$monitor)
  })

  ref_names <- colnames(cc_chains[[1]])
  if (is.null(ref_names) || !length(ref_names)) {
    stopf("No monitored parameters for %s.", cc)
  }

  if (!all(vapply(cc_chains, function(x) identical(colnames(x), ref_names), logical(1)))) {
    stopf("Monitor parameter names differ across chains for %s.", cc)
  }
  if (!all(vapply(cc_chains, nrow, integer(1)) == EXPECTED_STORED)) {
    stopf("Monitor draw counts differ across chains for %s.", cc)
  }

  for (pp in seq_along(ref_names)) {
    mat <- do.call(cbind, lapply(cc_chains, function(x) x[, pp]))

    finite <- all(is.finite(mat))
    chain_sd <- apply(mat, 2, stats::sd)
    constant <- finite && all(chain_sd < 1e-12)

    if (!finite || constant) {
      rh <- eb <- et <- NA_real_
    } else {
      rh <- tryCatch(posterior::rhat(mat), error = function(e) NA_real_)
      eb <- tryCatch(posterior::ess_bulk(mat), error = function(e) NA_real_)
      et <- tryCatch(posterior::ess_tail(mat), error = function(e) NA_real_)
    }

    rr <- rr + 1L
    diag_rows[[rr]] <- data.frame(
      Country = cc,
      Parameter = ref_names[pp],
      Chains = NCHAINS,
      DrawsPerChain = EXPECTED_STORED,
      Finite = finite,
      ConstantAcrossAllChains = constant,
      Rhat = rh,
      ESS_Bulk = eb,
      ESS_Tail = et,
      PassRhatTarget = if (is.finite(rh)) rh <= RHAT_TARGET else NA,
      PassRhatHard = if (is.finite(rh)) rh <= RHAT_HARD else NA,
      PassESSBulkTarget = if (is.finite(eb)) eb >= ESS_TARGET else NA,
      PassESSTailTarget = if (is.finite(et)) et >= ESS_TARGET else NA,
      PassESSBulkHard = if (is.finite(eb)) eb >= ESS_HARD else NA,
      PassESSTailHard = if (is.finite(et)) et >= ESS_HARD else NA,
      stringsAsFactors = FALSE
    )
  }

  cc_diag <- do.call(rbind, diag_rows)[
    vapply(diag_rows, function(x) identical(x$Country[1], cc), logical(1)),
    ,
    drop = FALSE
  ]
  # The construction above is awkward after accumulating all countries.
  # Recompute the country's slice safely from the cumulative object.
  tmp_all <- do.call(rbind, diag_rows)
  cc_diag <- tmp_all[tmp_all$Country == cc, , drop = FALSE]
  active <- cc_diag[cc_diag$Finite & !cc_diag$ConstantAcrossAllChains, , drop = FALSE]

  cc_pos <- cc_pos + 1L
  country_summary_rows[[cc_pos]] <- data.frame(
    Country = cc,
    MonitoredParameters = nrow(cc_diag),
    ActiveDiagnostics = nrow(active),
    ConstantParameters = sum(cc_diag$ConstantAcrossAllChains),
    NonFiniteParameters = sum(!cc_diag$Finite),
    MaxRhat = if (nrow(active)) max(active$Rhat, na.rm = TRUE) else NA_real_,
    RhatTargetShare = if (nrow(active)) mean(active$Rhat <= RHAT_TARGET, na.rm = TRUE) else NA_real_,
    RhatHardShare = if (nrow(active)) mean(active$Rhat <= RHAT_HARD, na.rm = TRUE) else NA_real_,
    MinESSBulk = if (nrow(active)) min(active$ESS_Bulk, na.rm = TRUE) else NA_real_,
    MinESSTail = if (nrow(active)) min(active$ESS_Tail, na.rm = TRUE) else NA_real_,
    ESSBulkTargetShare = if (nrow(active)) mean(active$ESS_Bulk >= ESS_TARGET, na.rm = TRUE) else NA_real_,
    ESSTailTargetShare = if (nrow(active)) mean(active$ESS_Tail >= ESS_TARGET, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )

  rm(cc_chains)
  invisible(gc())
}

diagnostics <- do.call(rbind, diag_rows)
country_summary <- do.call(rbind, country_summary_rows)

write.csv(
  diagnostics,
  file.path(OUT, "01_mcmc_parameter_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  country_summary,
  file.path(OUT, "02_mcmc_country_summary.csv"),
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

overall <- data.frame(
  Metric = c(
    "Countries",
    "ChainsPerCountry",
    "InternalPostBurnIterationsTotalPerCountry",
    "StoredDrawsTotalPerCountry",
    "MonitoredParameters",
    "ActiveNonconstantDiagnostics",
    "ConstantParameters",
    "NonFiniteParameters",
    "MaxRhat",
    "ShareRhatAtOrBelowTarget",
    "ShareRhatAtOrBelowHardLimit",
    "MinESSBulk",
    "MinESSTail",
    "ShareESSBulkAtOrAboveTarget",
    "ShareESSTailAtOrAboveTarget",
    "ShareESSBulkAtOrAboveHardMinimum",
    "ShareESSTailAtOrAboveHardMinimum"
  ),
  Value = c(
    length(COUNTRIES),
    NCHAINS,
    unique(parts_index$KeepInternal)[1] * NCHAINS,
    EXPECTED_STORED * NCHAINS,
    nrow(diagnostics),
    nrow(active),
    sum(diagnostics$ConstantAcrossAllChains),
    sum(!diagnostics$Finite),
    max(active$Rhat),
    mean(active$Rhat <= RHAT_TARGET),
    mean(active$Rhat <= RHAT_HARD),
    min(active$ESS_Bulk),
    min(active$ESS_Tail),
    mean(active$ESS_Bulk >= ESS_TARGET),
    mean(active$ESS_Tail >= ESS_TARGET),
    mean(active$ESS_Bulk >= ESS_HARD),
    mean(active$ESS_Tail >= ESS_HARD)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  overall,
  file.path(OUT, "03_mcmc_overall_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# 2. Aggregate posterior time-variation probabilities across chains
# =============================================================================

tv_rows <- list()
tv_pos <- 0L

for (cc in COUNTRIES) {
  cc_parts <- lapply(seq_len(NCHAINS), function(ch) {
    f <- parts_index$File[
      parts_index$Country == cc & parts_index$Chain == ch
    ]
    readRDS(f)
  })

  ref <- cc_parts[[1]]$tv_probability
  if (is.null(dim(ref)) || length(dim(ref)) != 3L) {
    stopf("Malformed time-variation probability array for %s.", cc)
  }

  avg <- ref * 0
  for (ch in seq_len(NCHAINS)) {
    if (!identical(dim(cc_parts[[ch]]$tv_probability), dim(ref))) {
      stopf("Time-variation array dimensions differ across chains for %s.", cc)
    }
    avg <- avg + cc_parts[[ch]]$tv_probability / NCHAINS
  }

  grid <- expand.grid(
    Quarter = dimnames(avg)[[1]],
    Term = dimnames(avg)[[2]],
    Equation = dimnames(avg)[[3]],
    stringsAsFactors = FALSE
  )
  grid$Country <- cc
  grid$PosteriorTimeVariationProbability <- as.vector(avg)

  tv_pos <- tv_pos + 1L
  tv_rows[[tv_pos]] <- grid[, c(
    "Country", "Equation", "Quarter", "Term",
    "PosteriorTimeVariationProbability"
  )]

  rm(cc_parts, avg)
  invisible(gc())
}

tv_df <- do.call(rbind, tv_rows)
write.csv(
  tv_df,
  file.path(OUT, "04_time_variation_probability.csv"),
  row.names = FALSE
)

# =============================================================================
# 3. Formal convergence gate
# =============================================================================

parts_ok <- (
  nrow(parts_index) == length(COUNTRIES) * NCHAINS &&
  all(parts_index$StoredDraws == EXPECTED_STORED)
)

finite_ok <- all(diagnostics$Finite)
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
  parts_ok &&
  finite_ok &&
  rhat_hard_ok &&
  rhat_target_ok &&
  ess_bulk_hard_ok &&
  ess_tail_hard_ok &&
  ess_target_ok
)

status <- if (ready) "READY_FOR_FORMAL_IRF" else "FAIL"

reasons <- c(
  if (!parts_ok) "country-chain part grid incomplete" else NULL,
  if (!finite_ok) "non-finite monitored MCMC draws detected" else NULL,
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
if (!length(reasons)) reasons <- "all formal MCMC convergence gates passed"

gate <- data.frame(
  Status = status,
  Countries = length(COUNTRIES),
  ChainsPerCountry = NCHAINS,
  BurnPerChain = unique(parts_index$Burn)[1],
  InternalPostBurnIterationsPerChain = unique(parts_index$KeepInternal)[1],
  TotalInternalPostBurnIterationsPerCountry = unique(parts_index$KeepInternal)[1] * NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  TotalStoredDrawsPerCountry = EXPECTED_STORED * NCHAINS,
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
  Reason = paste(reasons, collapse = "; "),
  stringsAsFactors = FALSE
)

write.csv(
  gate,
  file.path(OUT, "00_formal_mcmc_gate.csv"),
  row.names = FALSE
)

readme <- c(
  sprintf("FORMAL BAYESIAN TVP-GVAR MCMC: %s", status),
  "==============================================",
  sprintf("Countries: %d", length(COUNTRIES)),
  sprintf("Chains per country: %d", NCHAINS),
  sprintf("Burn per chain: %d", gate$BurnPerChain),
  sprintf("Internal post-burn iterations per chain: %d", gate$InternalPostBurnIterationsPerChain),
  sprintf(
    "Total internal post-burn iterations per country: %d",
    gate$TotalInternalPostBurnIterationsPerCountry
  ),
  sprintf("Stored draws per chain: %d", EXPECTED_STORED),
  sprintf("Total stored draws per country: %d", gate$TotalStoredDrawsPerCountry),
  sprintf("Max rank-normalized split-Rhat: %.6f", gate$MaxRhat),
  sprintf("Rhat <= %.3f share: %.6f", RHAT_TARGET, rhat_target_share),
  sprintf("Min bulk ESS: %.2f", gate$MinESSBulk),
  sprintf("Min tail ESS: %.2f", gate$MinESSTail),
  sprintf("Bulk ESS >= %.0f share: %.6f", ESS_TARGET, ess_bulk_share),
  sprintf("Tail ESS >= %.0f share: %.6f", ESS_TARGET, ess_tail_share),
  sprintf("Reason: %s", gate$Reason),
  "",
  "Diagnostics use rank-normalized split-Rhat, bulk ESS and tail ESS.",
  "Constant monitored quantities are reported but excluded from Rhat/ESS gating.",
  "Passing this gate is required before structural GPR IRFs are generated."
)

writeLines(readme, file.path(OUT, "README_formal_mcmc.txt"))
cat(paste(readme, collapse = "\n"), "\n")
