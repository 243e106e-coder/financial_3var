#!/usr/bin/env Rscript

# =============================================================================
# 65_rate_diff_formal_decision.R
#
# Final decision gate for the formal Delta-r TVP-GVAR re-estimation.
#
# This script does NOT estimate the model and does NOT generate IRFs.
# It combines:
#   1) same-sample static rate-level vs Delta-r evidence from R/15,
#   2) formal four-chain MCMC convergence from R/61,
#   3) posterior GVAR dynamic stability from R/63,
#   4) exact posterior lineage: every country-chain must come from the
#      Delta-r panel.
#
# A READY result means:
#   - the Delta-r specification is statistically and dynamically admissible
#     for the next IRF stage;
#   - R/62 must still be adapted so the interest-rate response is cumulated
#     before it is interpreted as a rate-level response.
#
# It is NOT permission to use the old rate-level R/62 output unchanged.
# =============================================================================

source("R/00_config.R")

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
CONV_DIR <- get_env_chr("FIN3_CONVERGENCE_DIR", file.path(RESULTS_DIR, "formal_tvp"))
DYN_DIR <- get_env_chr("FIN3_DYNAMIC_DIR", file.path(RESULTS_DIR, "dynamic_stability_source"))
TRANS_DIR <- get_env_chr("FIN3_TRANSFORM_DIR", file.path(RESULTS_DIR, "variable_transform"))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 2000))
MIN_ANCHOR_STABLE <- get_env_num("FIN3_RATE_DIFF_MIN_ANCHOR_STABLE", 0.90)
EXPECTED_PANEL_BASENAME <- get_env_chr(
  "FIN3_RATE_DIFF_PANEL_BASENAME",
  "panel_domestic_fin3_rate_diff.csv"
)

OUT <- file.path(RESULTS_DIR, "rate_diff_decision")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Static same-sample comparison
# -----------------------------------------------------------------------------

cmp_file <- file.path(TRANS_DIR, "09_same_sample_rate_spec_comparison.csv")
if (!file.exists(cmp_file)) stopf("Missing transformation comparison: %s", cmp_file)

cmp <- read.csv(cmp_file, stringsAsFactors = FALSE, check.names = FALSE)
need_cmp <- c(
  "Network",
  "SampleStart", "SampleEnd", "T",
  "GlobalSpectralRadius_RateLevel", "Stable_RateLevel",
  "GlobalSpectralRadius_RateDiff", "Stable_RateDiff",
  "MaxLocalDomesticRadius_RateLevel",
  "MaxLocalDomesticRadius_RateDiff"
)
if (!all(need_cmp %in% names(cmp))) {
  stopf("Malformed same-sample transformation comparison.")
}

main_cmp <- cmp[cmp$Network == MAIN_NETWORK, , drop = FALSE]
if (nrow(main_cmp) != 1L) {
  stopf("Expected one MAIN_NETWORK row in same-sample comparison.")
}

static_level_rho <- as.numeric(main_cmp$GlobalSpectralRadius_RateLevel[1])
static_diff_rho <- as.numeric(main_cmp$GlobalSpectralRadius_RateDiff[1])
static_diff_stable <- isTRUE(as.logical(main_cmp$Stable_RateDiff[1]))
static_improves <- is.finite(static_level_rho) &&
  is.finite(static_diff_rho) &&
  static_diff_rho < static_level_rho

# -----------------------------------------------------------------------------
# 2. Formal four-chain MCMC convergence
# -----------------------------------------------------------------------------

conv_file <- file.path(CONV_DIR, "00_formal_mcmc_gate.csv")
if (!file.exists(conv_file)) stopf("Missing formal convergence gate: %s", conv_file)

conv <- read.csv(conv_file, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(conv) != 1L || !"Status" %in% names(conv)) {
  stopf("Malformed formal convergence gate.")
}
mcmc_ready <- identical(trimws(as.character(conv$Status[1])), "READY_FOR_FORMAL_IRF")

# -----------------------------------------------------------------------------
# 3. Posterior dynamic stability
# -----------------------------------------------------------------------------

dyn_file <- file.path(DYN_DIR, "00_dynamic_stability_source_gate.csv")
if (!file.exists(dyn_file)) stopf("Missing posterior dynamic stability gate: %s", dyn_file)

dyn <- read.csv(dyn_file, stringsAsFactors = FALSE, check.names = FALSE)
need_dyn <- c(
  "Status", "MainNetwork", "Countries", "Chains",
  "StoredDrawsPerChain",
  "BaselineMeanStableShare", "BaselineMinStableShare"
)
if (nrow(dyn) != 1L || !all(need_dyn %in% names(dyn))) {
  stopf("Malformed posterior dynamic stability gate.")
}

dyn_integrity <- identical(
  trimws(as.character(dyn$Status[1])),
  "DIAGNOSTIC_COMPLETE"
)
dyn_network_ok <- identical(
  trimws(as.character(dyn$MainNetwork[1])),
  MAIN_NETWORK
)
dyn_grid_ok <- as.integer(dyn$Countries[1]) == length(COUNTRIES) &&
  as.integer(dyn$Chains[1]) == 4L &&
  as.integer(dyn$StoredDrawsPerChain[1]) == EXPECTED_STORED

mean_stable <- as.numeric(dyn$BaselineMeanStableShare[1])
min_stable <- as.numeric(dyn$BaselineMinStableShare[1])
dynamic_ready <- is.finite(min_stable) && min_stable >= MIN_ANCHOR_STABLE

# -----------------------------------------------------------------------------
# 4. Exact posterior lineage: all 14 x 4 RDS files must use the Delta-r panel
# -----------------------------------------------------------------------------

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+[.]rds$",
  recursive = TRUE,
  full.names = TRUE
)

expected_n <- length(COUNTRIES) * 4L
if (length(files) != expected_n) {
  stopf("Expected exactly %d posterior RDS files; found %d.", expected_n, length(files))
}

lineage_rows <- vector("list", length(files))
for (ii in seq_along(files)) {
  z <- readRDS(files[ii])
  if (is.null(z$meta)) stopf("Missing meta in posterior: %s", files[ii])

  cc <- as.character(z$meta$country)
  ch <- as.integer(z$meta$chain)
  stored <- as.integer(z$meta$stored_draws)
  src <- as.character(z$meta$source_panel)

  lineage_rows[[ii]] <- data.frame(
    Country = cc,
    Chain = ch,
    StoredDraws = stored,
    SourcePanel = src,
    SourcePanelBasename = basename(src),
    IsRateDiffPanel = identical(basename(src), EXPECTED_PANEL_BASENAME),
    File = files[ii],
    stringsAsFactors = FALSE
  )
}
lineage <- do.call(rbind, lineage_rows)

key <- paste(lineage$Country, lineage$Chain, sep = "||")
expected_key <- paste(
  rep(COUNTRIES, each = 4L),
  rep(1:4, times = length(COUNTRIES)),
  sep = "||"
)

lineage_grid_ok <- !anyDuplicated(key) &&
  setequal(key, expected_key) &&
  all(lineage$StoredDraws == EXPECTED_STORED)

lineage_panel_ok <- all(lineage$IsRateDiffPanel)

write.csv(
  lineage,
  file.path(OUT, "02_rate_diff_posterior_lineage.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Decision table
# -----------------------------------------------------------------------------

checks <- data.frame(
  Check = c(
    "StaticRateDiffStable",
    "StaticRateDiffImprovesOnSameSampleLevel",
    "FormalFourChainMCMCReady",
    "DynamicDiagnosticIntegrity",
    "DynamicNetworkMatchesMainNetwork",
    "DynamicPosteriorGrid14x4",
    "EveryAnchorStableShareAtLeastThreshold",
    "PosteriorGridComplete14x4",
    "EveryPosteriorUsesRateDiffPanel"
  ),
  Pass = c(
    static_diff_stable && is.finite(static_diff_rho) && static_diff_rho < 1,
    static_improves,
    mcmc_ready,
    dyn_integrity,
    dyn_network_ok,
    dyn_grid_ok,
    dynamic_ready,
    lineage_grid_ok,
    lineage_panel_ok
  ),
  Detail = c(
    sprintf("Delta-r same-sample static rho=%.8f", static_diff_rho),
    sprintf(
      "Same-sample level rho=%.8f; Delta-r rho=%.8f; delta=%.8f",
      static_level_rho, static_diff_rho, static_diff_rho - static_level_rho
    ),
    sprintf("R61 Status=%s", as.character(conv$Status[1])),
    sprintf("R63 Status=%s", as.character(dyn$Status[1])),
    sprintf("R63 MainNetwork=%s; configured=%s", as.character(dyn$MainNetwork[1]), MAIN_NETWORK),
    sprintf(
      "Countries=%s; Chains=%s; Stored=%s",
      as.character(dyn$Countries[1]),
      as.character(dyn$Chains[1]),
      as.character(dyn$StoredDrawsPerChain[1])
    ),
    sprintf(
      "Baseline mean stable share=%.6f; minimum anchor stable share=%.6f; required minimum=%.6f",
      mean_stable, min_stable, MIN_ANCHOR_STABLE
    ),
    sprintf("Posterior files=%d; expected=%d", nrow(lineage), expected_n),
    sprintf(
      "Rate-diff panel files=%d/%d; expected basename=%s",
      sum(lineage$IsRateDiffPanel), nrow(lineage), EXPECTED_PANEL_BASENAME
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  checks,
  file.path(OUT, "01_rate_diff_formal_checks.csv"),
  row.names = FALSE
)

all_ready <- all(checks$Pass)

status <- if (all_ready) {
  "RATE_DIFF_READY_FOR_IRF_RECODE"
} else {
  "RATE_DIFF_NOT_READY"
}

decision <- data.frame(
  Status = status,
  MainNetwork = MAIN_NETWORK,
  SampleStart = as.character(main_cmp$SampleStart[1]),
  SampleEnd = as.character(main_cmp$SampleEnd[1]),
  SameSampleT = as.integer(main_cmp$T[1]),
  StaticRateLevelRho = static_level_rho,
  StaticRateDiffRho = static_diff_rho,
  StaticRhoImprovement = static_level_rho - static_diff_rho,
  MCMCStatus = as.character(conv$Status[1]),
  DynamicDiagnosticStatus = as.character(dyn$Status[1]),
  BaselineMeanStableShare = mean_stable,
  BaselineMinStableShare = min_stable,
  RequiredMinAnchorStableShare = MIN_ANCHOR_STABLE,
  PosteriorFiles = nrow(lineage),
  RateDiffPosteriorFiles = sum(lineage$IsRateDiffPanel),
  FormalIRFGenerated = FALSE,
  NextStep = if (all_ready) {
    paste0(
      "Adapt R62 for Delta-r: keep de/deq treatment unchanged; ",
      "cumulate the r-response across horizons before interpreting it as a ",
      "rate-level response; then run the formal GPR IRF gate."
    )
  } else {
    "Do not generate formal IRFs. Inspect failed checks and posterior stability diagnostics first."
  },
  stringsAsFactors = FALSE
)

write.csv(
  decision,
  file.path(OUT, "00_rate_diff_formal_gate.csv"),
  row.names = FALSE
)

readme <- c(
  sprintf("08j RATE-DIFFERENCE FORMAL DECISION: %s", status),
  "============================================================",
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf(
    "Same-sample static rho: rate level %.8f -> Delta-r %.8f",
    static_level_rho, static_diff_rho
  ),
  sprintf("Formal MCMC gate: %s", as.character(conv$Status[1])),
  sprintf("Posterior mean stable share: %.6f", mean_stable),
  sprintf("Posterior minimum anchor stable share: %.6f", min_stable),
  sprintf("Required minimum anchor stable share: %.6f", MIN_ANCHOR_STABLE),
  sprintf(
    "Posterior lineage: %d/%d files use %s",
    sum(lineage$IsRateDiffPanel), nrow(lineage), EXPECTED_PANEL_BASENAME
  ),
  "",
  "Interpretation:",
  "- A READY result does not generate or approve the old rate-level IRF code.",
  "- Delta-r is the modeled interest-rate variable in this specification.",
  "- Its horizon-by-horizon response must be cumulatively summed before it is interpreted as a rate-level response.",
  "- de remains REER_DLOG and may be cumulated for a cumulative REER log-level effect.",
  "- deq remains EQ_RETURN and must not be relabeled as a log-price response without upstream verification.",
  "",
  "Outputs:",
  "- 00_rate_diff_formal_gate.csv",
  "- 01_rate_diff_formal_checks.csv",
  "- 02_rate_diff_posterior_lineage.csv"
)

writeLines(readme, file.path(OUT, "README_08j_rate_diff_formal_decision.txt"))
cat(paste(readme, collapse = "\n"), "\n")

if (!all_ready) {
  print(checks[!checks$Pass, , drop = FALSE])
  stop("08j Delta-r formal decision gate failed.")
}
