#!/usr/bin/env Rscript

# Audit the fixed 500-draw TVP-GVAR validation run.
# This script does NOT estimate the model. It verifies that the smoke-test
# output corresponds exactly to the selected hyperparameter specification and
# that the original stability gates are still satisfied.

SMOKE_DIR <- file.path("results", "tvp_smoke")
OUT_DIR <- file.path("results", "tvp_selected_validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)

get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("Environment variable %s is not numeric: %s", name, z)
  out
}

EXPECTED_DRAWS <- as.integer(get_env_num("TVP_VALIDATION_EXPECTED_DRAWS", 500))
EXPECTED_STATE_SCALE <- get_env_num("TVP_VALIDATION_EXPECTED_STATE_SCALE", 1e-6)
EXPECTED_PRIOR_SCALE <- get_env_num("TVP_VALIDATION_EXPECTED_PRIOR_SCALE", 0.02)
MIN_STABLE_SHARE <- get_env_num("TVP_MIN_STABLE_SHARE", 0.90)
MIN_G0_OK_SHARE <- 0.99
MAX_SPLIT_DIFF <- get_env_num("TVP_MAX_SPLIT_DIFF", 0.10)
MAX_ABS_STD_BETA <- get_env_num("TVP_MAX_ABS_STD_BETA", 50)

GATE_FILE <- file.path(SMOKE_DIR, "00_tvp_smoke_gate.csv")
QUARTER_FILE <- file.path(SMOKE_DIR, "01_posterior_stability_by_quarter.csv")
DRAW_FILE <- file.path(SMOKE_DIR, "03_posterior_stability_by_draw.csv")
LOCAL_FILE <- file.path(SMOKE_DIR, "05_local_posterior_stability.csv")

if (!file.exists(GATE_FILE)) {
  stopf("Missing %s. The TVP smoke test did not produce its gate output.", GATE_FILE)
}

g <- read.csv(GATE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(g) != 1L) stopf("Expected exactly one row in %s.", GATE_FILE)

required <- c(
  "MainNetwork", "Status", "Draws", "Seed", "StateScale", "PriorScale", "Sample",
  "FiniteDrawQuarterShare", "PosteriorStableShare", "G0OKShare",
  "FirstHalfStableShare", "SecondHalfStableShare", "SplitDifference",
  "MaxAbsStandardizedBeta", "MinStableShareRequired", "Reason"
)
missing_cols <- setdiff(required, names(g))
if (length(missing_cols)) stopf("Smoke gate is missing columns: %s", paste(missing_cols, collapse = ", "))

near <- function(x, y, tol = 1e-12) is.finite(x) && is.finite(y) && abs(x - y) <= tol * max(1, abs(x), abs(y))

checks <- c(
  GateStatusREADY = identical(trimws(as.character(g$Status[1])), "READY"),
  DrawsExactlyExpected = as.integer(g$Draws[1]) == EXPECTED_DRAWS,
  StateScaleExactlySelected = near(as.numeric(g$StateScale[1]), EXPECTED_STATE_SCALE),
  PriorScaleExactlySelected = near(as.numeric(g$PriorScale[1]), EXPECTED_PRIOR_SCALE),
  FiniteDrawQuarterShareIsOne = near(as.numeric(g$FiniteDrawQuarterShare[1]), 1),
  PosteriorStableSharePass = as.numeric(g$PosteriorStableShare[1]) >= MIN_STABLE_SHARE,
  G0OKSharePass = as.numeric(g$G0OKShare[1]) >= MIN_G0_OK_SHARE,
  SplitDifferencePass = as.numeric(g$SplitDifference[1]) <= MAX_SPLIT_DIFF,
  MaxAbsStandardizedBetaPass = as.numeric(g$MaxAbsStandardizedBeta[1]) <= MAX_ABS_STD_BETA
)
checks[is.na(checks)] <- FALSE

# Diagnostics only: do not introduce new rejection thresholds here.
worst_country <- NA_character_
worst_country_share <- NA_real_
worst_country_rho95 <- NA_real_
if (file.exists(LOCAL_FILE)) {
  z <- read.csv(LOCAL_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(z) && all(c("Country", "StableShare") %in% names(z))) {
    ii <- which.min(z$StableShare)
    worst_country <- as.character(z$Country[ii])
    worst_country_share <- as.numeric(z$StableShare[ii])
    if ("RhoP95" %in% names(z)) worst_country_rho95 <- as.numeric(z$RhoP95[ii])
  }
}

worst_quarter <- NA_character_
worst_quarter_share <- NA_real_
worst_quarter_rho95 <- NA_real_
if (file.exists(QUARTER_FILE)) {
  z <- read.csv(QUARTER_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(z) && all(c("Quarter", "StableShare") %in% names(z))) {
    ii <- which.min(z$StableShare)
    worst_quarter <- as.character(z$Quarter[ii])
    worst_quarter_share <- as.numeric(z$StableShare[ii])
    if ("RhoP95" %in% names(z)) worst_quarter_rho95 <- as.numeric(z$RhoP95[ii])
  }
}

n_draws_observed <- NA_integer_
all_quarters_stable_draws <- NA_integer_
zero_stable_draws <- NA_integer_
if (file.exists(DRAW_FILE)) {
  z <- read.csv(DRAW_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(z) && all(c("Draw", "StableShare") %in% names(z))) {
    n_draws_observed <- length(unique(z$Draw))
    all_quarters_stable_draws <- sum(z$StableShare >= 1 - 1e-15, na.rm = TRUE)
    zero_stable_draws <- sum(z$StableShare <= 1e-15, na.rm = TRUE)
  }
}

ready <- all(checks)
final_status <- if (ready) "READY_FOR_FORMAL_TVP" else "FAIL_VALIDATION"
failed_names <- names(checks)[!checks]
reason <- if (ready) {
  "Fixed selected specification passed all original TVP smoke-test gates at 500 draws."
} else {
  paste0("Failed checks: ", paste(failed_names, collapse = "; "))
}

out <- data.frame(
  FinalStatus = final_status,
  MainNetwork = as.character(g$MainNetwork[1]),
  DrawsExpected = EXPECTED_DRAWS,
  DrawsReported = as.integer(g$Draws[1]),
  Seed = as.integer(g$Seed[1]),
  StateScaleExpected = EXPECTED_STATE_SCALE,
  StateScaleReported = as.numeric(g$StateScale[1]),
  PriorScaleExpected = EXPECTED_PRIOR_SCALE,
  PriorScaleReported = as.numeric(g$PriorScale[1]),
  Sample = as.character(g$Sample[1]),
  FiniteDrawQuarterShare = as.numeric(g$FiniteDrawQuarterShare[1]),
  PosteriorStableShare = as.numeric(g$PosteriorStableShare[1]),
  G0OKShare = as.numeric(g$G0OKShare[1]),
  FirstHalfStableShare = as.numeric(g$FirstHalfStableShare[1]),
  SecondHalfStableShare = as.numeric(g$SecondHalfStableShare[1]),
  SplitDifference = as.numeric(g$SplitDifference[1]),
  MaxAbsStandardizedBeta = as.numeric(g$MaxAbsStandardizedBeta[1]),
  WorstCountry = worst_country,
  WorstCountryStableShare = worst_country_share,
  WorstCountryRhoP95 = worst_country_rho95,
  WorstQuarter = worst_quarter,
  WorstQuarterStableShare = worst_quarter_share,
  WorstQuarterRhoP95 = worst_quarter_rho95,
  DrawsObservedInDrawDiagnostic = n_draws_observed,
  AllQuartersStableDraws = all_quarters_stable_draws,
  ZeroStableDraws = zero_stable_draws,
  Reason = reason,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(OUT_DIR, "00_selected_spec_validation.csv"), row.names = FALSE)

check_table <- data.frame(Check = names(checks), Pass = unname(checks), stringsAsFactors = FALSE)
write.csv(check_table, file.path(OUT_DIR, "01_selected_spec_checks.csv"), row.names = FALSE)

# Preserve a self-contained copy of the underlying smoke-test outputs.
for (f in list.files(SMOKE_DIR, full.names = TRUE)) {
  if (file.info(f)$isdir) next
  file.copy(f, file.path(OUT_DIR, basename(f)), overwrite = TRUE)
}

fmt <- function(x, digits = 6) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

lines <- c(
  sprintf("SELECTED TVP SPECIFICATION VALIDATION: %s", final_status),
  "============================================================",
  sprintf("Network: %s", g$MainNetwork[1]),
  sprintf("Sample: %s", g$Sample[1]),
  sprintf("Draws: %d (expected %d)", as.integer(g$Draws[1]), EXPECTED_DRAWS),
  sprintf("Seed: %d", as.integer(g$Seed[1])),
  sprintf("State scale: %.8g (selected %.8g)", as.numeric(g$StateScale[1]), EXPECTED_STATE_SCALE),
  sprintf("Prior scale: %.8g (selected %.8g)", as.numeric(g$PriorScale[1]), EXPECTED_PRIOR_SCALE),
  sprintf("Posterior stable share: %s [required >= %.2f]", fmt(as.numeric(g$PosteriorStableShare[1])), MIN_STABLE_SHARE),
  sprintf("G0 OK share: %s [required >= %.2f]", fmt(as.numeric(g$G0OKShare[1])), MIN_G0_OK_SHARE),
  sprintf("First-half / second-half stable share: %s / %s", fmt(as.numeric(g$FirstHalfStableShare[1])), fmt(as.numeric(g$SecondHalfStableShare[1]))),
  sprintf("Split difference: %s [required <= %.2f]", fmt(as.numeric(g$SplitDifference[1])), MAX_SPLIT_DIFF),
  sprintf("Max standardized |beta|: %s [required <= %.2f]", fmt(as.numeric(g$MaxAbsStandardizedBeta[1])), MAX_ABS_STD_BETA),
  sprintf("Worst country: %s; stable share: %s; rho P95: %s", worst_country, fmt(worst_country_share), fmt(worst_country_rho95)),
  sprintf("Worst quarter: %s; stable share: %s; rho P95: %s", worst_quarter, fmt(worst_quarter_share), fmt(worst_quarter_rho95)),
  sprintf("Draw-level diagnostic: observed=%s; stable in all quarters=%s; stable in zero quarters=%s",
          n_draws_observed, all_quarters_stable_draws, zero_stable_draws),
  "",
  sprintf("Decision: %s", reason),
  "",
  "Interpretation:",
  "- This is a validation of the selected smoke-test hyperparameters, not the final publication posterior.",
  "- No stability threshold is relaxed relative to the original smoke-test gate.",
  "- READY_FOR_FORMAL_TVP means the selected architecture is robust enough to proceed to the formal Bayesian TVP-GVAR stage."
)
writeLines(lines, file.path(OUT_DIR, "README_selected_tvp_validation.txt"))
cat(paste(lines, collapse = "\n"), "\n")

if (!ready) quit(save = "no", status = 2L)
