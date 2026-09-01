#!/usr/bin/env Rscript

# Summarize the 3x3 TVP-GVAR hyperparameter grid produced by
# .github/workflows/04_tvp_hyperparameter_grid.yml
#
# This script does NOT estimate the TVP-GVAR itself. It collects the gate,
# quarter-level, and country-level diagnostics produced by
# R/50_tvp_gvar_smoke_test.R for each state-scale / prior-scale combination.

GRID_DIR <- Sys.getenv("TVP_GRID_DIR", "results/tvp_grid")
OUT_DIR <- GRID_DIR

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
msg <- function(fmt, ...) cat(sprintf(fmt, ...), "\n")

num_or_na <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  out <- suppressWarnings(as.numeric(x[1]))
  if (is.finite(out)) out else NA_real_
}

chr_or_na <- function(x) {
  if (length(x) == 0L) return(NA_character_)
  out <- as.character(x[1])
  if (nzchar(out)) out else NA_character_
}

if (!dir.exists(GRID_DIR)) stopf("Grid directory does not exist: %s", GRID_DIR)

spec_dirs <- list.dirs(GRID_DIR, recursive = FALSE, full.names = TRUE)
# Exclude accidental non-spec directories if the script is rerun after outputs exist.
spec_dirs <- spec_dirs[file.exists(file.path(spec_dirs, "RUN_METADATA.csv"))]
if (!length(spec_dirs)) stopf("No TVP grid specifications found under %s", GRID_DIR)

rows <- vector("list", length(spec_dirs))
local_parts <- list()
quarter_parts <- list()
local_pos <- 0L
quarter_pos <- 0L

for (ii in seq_along(spec_dirs)) {
  d <- spec_dirs[ii]
  spec_id <- basename(d)

  meta_path <- file.path(d, "RUN_METADATA.csv")
  gate_path <- file.path(d, "00_tvp_smoke_gate.csv")
  local_path <- file.path(d, "05_local_posterior_stability.csv")
  quarter_path <- file.path(d, "01_posterior_stability_by_quarter.csv")

  meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(meta) != 1L) stopf("Malformed RUN_METADATA.csv in %s", d)

  state_scale <- num_or_na(meta$StateScale)
  prior_scale <- num_or_na(meta$PriorScale)
  exit_code <- suppressWarnings(as.integer(meta$ExitCode[1]))
  if (!is.finite(exit_code)) exit_code <- NA_integer_

  gate_exists <- file.exists(gate_path)
  if (gate_exists) {
    gate <- read.csv(gate_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(gate) != 1L) stopf("Malformed gate file in %s", d)
  } else {
    gate <- data.frame(stringsAsFactors = FALSE)
  }

  gate_value <- function(name, numeric = FALSE) {
    if (!gate_exists || !name %in% names(gate)) return(if (numeric) NA_real_ else NA_character_)
    if (numeric) num_or_na(gate[[name]]) else chr_or_na(gate[[name]])
  }

  status <- gate_value("Status")
  if (is.na(status)) status <- "ERROR"

  operational_status <- if (identical(exit_code, 0L)) {
    "COMPLETED_READY"
  } else if (identical(exit_code, 2L) && gate_exists) {
    "COMPLETED_MODEL_FAIL"
  } else {
    "OPERATIONAL_ERROR"
  }

  # Summaries that make it easier to see whether JP/SG or particular quarters
  # remain the binding source of instability.
  worst_country <- NA_character_
  worst_country_stable_share <- NA_real_
  if (file.exists(local_path)) {
    z <- read.csv(local_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(z)) {
      z$SpecID <- spec_id
      z$StateScale <- state_scale
      z$PriorScale <- prior_scale
      local_pos <- local_pos + 1L
      local_parts[[local_pos]] <- z
      ok <- is.finite(suppressWarnings(as.numeric(z$StableShare)))
      if (any(ok)) {
        zz <- z[ok, , drop = FALSE]
        jj <- which.min(as.numeric(zz$StableShare))
        worst_country <- as.character(zz$Country[jj])
        worst_country_stable_share <- as.numeric(zz$StableShare[jj])
      }
    }
  }

  worst_quarter <- NA_character_
  worst_quarter_stable_share <- NA_real_
  worst_quarter_rho_p95 <- NA_real_
  if (file.exists(quarter_path)) {
    q <- read.csv(quarter_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(q)) {
      q$SpecID <- spec_id
      q$StateScale <- state_scale
      q$PriorScale <- prior_scale
      quarter_pos <- quarter_pos + 1L
      quarter_parts[[quarter_pos]] <- q
      ok <- is.finite(suppressWarnings(as.numeric(q$StableShare)))
      if (any(ok)) {
        qq <- q[ok, , drop = FALSE]
        jj <- which.min(as.numeric(qq$StableShare))
        worst_quarter <- as.character(qq$Quarter[jj])
        worst_quarter_stable_share <- as.numeric(qq$StableShare[jj])
        if ("RhoP95" %in% names(qq)) worst_quarter_rho_p95 <- num_or_na(qq$RhoP95[jj])
      }
    }
  }

  rows[[ii]] <- data.frame(
    SpecID = spec_id,
    StateScale = state_scale,
    PriorScale = prior_scale,
    ExitCode = exit_code,
    OperationalStatus = operational_status,
    GateStatus = status,
    Draws = gate_value("Draws", TRUE),
    Sample = gate_value("Sample"),
    FiniteDrawQuarterShare = gate_value("FiniteDrawQuarterShare", TRUE),
    PosteriorStableShare = gate_value("PosteriorStableShare", TRUE),
    G0OKShare = gate_value("G0OKShare", TRUE),
    FirstHalfStableShare = gate_value("FirstHalfStableShare", TRUE),
    SecondHalfStableShare = gate_value("SecondHalfStableShare", TRUE),
    SplitDifference = gate_value("SplitDifference", TRUE),
    MaxAbsStandardizedBeta = gate_value("MaxAbsStandardizedBeta", TRUE),
    MinStableShareRequired = gate_value("MinStableShareRequired", TRUE),
    MinG0RcondRequired = gate_value("MinG0RcondRequired", TRUE),
    WorstCountry = worst_country,
    WorstCountryStableShare = worst_country_stable_share,
    WorstQuarter = worst_quarter,
    WorstQuarterStableShare = worst_quarter_stable_share,
    WorstQuarterRhoP95 = worst_quarter_rho_p95,
    Reason = gate_value("Reason"),
    stringsAsFactors = FALSE
  )
}

grid <- do.call(rbind, rows)

grid$ReadyCandidate <- grid$GateStatus == "READY" & grid$ExitCode == 0L

# Selection principle:
# Among specifications that genuinely pass all gates, prefer the LEAST
# restrictive one: first the largest state innovation scale (more genuine TVP),
# then the largest prior scale (more diffuse prior). Posterior stability and
# split consistency break remaining ties. This avoids selecting a nearly
# constant-parameter model merely because it has the highest stability share.
rank_order <- order(
  !grid$ReadyCandidate,
  -grid$StateScale,
  -grid$PriorScale,
  -grid$PosteriorStableShare,
  grid$SplitDifference,
  na.last = TRUE
)
grid <- grid[rank_order, , drop = FALSE]
row.names(grid) <- NULL

grid$Recommended <- FALSE
ready_idx <- which(grid$ReadyCandidate)
if (length(ready_idx)) grid$Recommended[ready_idx[1]] <- TRUE

write.csv(grid, file.path(OUT_DIR, "00_tvp_hyperparameter_grid_summary.csv"), row.names = FALSE)
write.csv(grid[grid$ReadyCandidate, , drop = FALSE],
          file.path(OUT_DIR, "01_ready_candidates.csv"), row.names = FALSE)

if (length(local_parts)) {
  local_all <- do.call(rbind, local_parts)
  row.names(local_all) <- NULL
  write.csv(local_all, file.path(OUT_DIR, "02_local_stability_by_spec.csv"), row.names = FALSE)
}

if (length(quarter_parts)) {
  quarter_all <- do.call(rbind, quarter_parts)
  row.names(quarter_all) <- NULL
  write.csv(quarter_all, file.path(OUT_DIR, "03_quarter_stability_by_spec.csv"), row.names = FALSE)
}

n_ready <- sum(grid$ReadyCandidate, na.rm = TRUE)
n_operational_error <- sum(grid$OperationalStatus == "OPERATIONAL_ERROR", na.rm = TRUE)

if (n_ready > 0L) {
  rec <- grid[which(grid$Recommended)[1], , drop = FALSE]
  recommendation <- c(
    "RECOMMENDED TVP HYPERPARAMETER SPECIFICATION",
    "========================================",
    sprintf("SpecID: %s", rec$SpecID),
    sprintf("StateScale: %.8g", rec$StateScale),
    sprintf("PriorScale: %.8g", rec$PriorScale),
    sprintf("PosteriorStableShare: %.6f", rec$PosteriorStableShare),
    sprintf("G0OKShare: %.6f", rec$G0OKShare),
    sprintf("SplitDifference: %.6f", rec$SplitDifference),
    sprintf("Worst country: %s (stable share %.6f)", rec$WorstCountry, rec$WorstCountryStableShare),
    sprintf("Worst quarter: %s (stable share %.6f)", rec$WorstQuarter, rec$WorstQuarterStableShare),
    "",
    "Selection rule:",
    "- Must pass the original TVP smoke-test gate without lowering thresholds.",
    "- Among passing specifications, prefer the largest StateScale, then largest PriorScale.",
    "- Stability share and split consistency are tie-breakers.",
    "- This favors genuine time variation while avoiding unnecessary shrinkage."
  )
} else {
  best <- grid[order(-grid$PosteriorStableShare, grid$SplitDifference, na.last = TRUE), , drop = FALSE]
  best <- if (nrow(best)) best[1, , drop = FALSE] else NULL
  recommendation <- c(
    "NO TVP GRID SPECIFICATION PASSED ALL GATES",
    "==========================================",
    sprintf("Specifications evaluated: %d", nrow(grid)),
    sprintf("Operational errors: %d", n_operational_error),
    if (!is.null(best)) sprintf("Highest posterior stable share observed: %.6f (%s)",
                                best$PosteriorStableShare, best$SpecID) else NULL,
    "",
    "Do not lower the 0.90 stability threshold merely to obtain a green Action.",
    "Inspect 00_tvp_hyperparameter_grid_summary.csv, especially WorstCountry and WorstQuarter,",
    "before changing model architecture or adding stronger shrinkage."
  )
}

writeLines(recommendation, file.path(OUT_DIR, "RECOMMENDED_SPEC.txt"))

readme <- c(
  "TVP-GVAR 3x3 HYPERPARAMETER GRID",
  "================================",
  sprintf("Specifications evaluated: %d", nrow(grid)),
  sprintf("READY candidates: %d", n_ready),
  sprintf("Operational errors: %d", n_operational_error),
  "",
  "Grid:",
  "- StateScale: 1e-5, 1e-6, 1e-7",
  "- PriorScale: 4, 1, 0.25",
  "",
  "Important:",
  "- Exit code 2 from R/50_tvp_gvar_smoke_test.R is treated as a MODEL FAIL, not a workflow crash.",
  "- Exit codes other than 0 or the expected gate-fail code 2 are marked OPERATIONAL_ERROR.",
  "- The original posterior stability threshold is preserved.",
  "- The recommended specification is only a smoke-test candidate, not the final publication model.",
  "",
  "Files:",
  "- 00_tvp_hyperparameter_grid_summary.csv: one row per hyperparameter combination",
  "- 01_ready_candidates.csv: only specifications passing all gates",
  "- 02_local_stability_by_spec.csv: country-level diagnostics across specifications",
  "- 03_quarter_stability_by_spec.csv: quarter-level diagnostics across specifications",
  "- RECOMMENDED_SPEC.txt: mechanically selected least-restrictive passing candidate"
)
writeLines(readme, file.path(OUT_DIR, "README_tvp_hyperparameter_grid.txt"))

msg("TVP hyperparameter grid summary complete.")
msg("READY candidates: %d / %d", n_ready, nrow(grid))
cat(paste(recommendation, collapse = "\n"), "\n")
