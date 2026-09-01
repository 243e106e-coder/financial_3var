#!/usr/bin/env Rscript

GRID_DIR <- Sys.getenv("TVP_REFINEMENT_DIR", "results/tvp_prior_refinement")
OUT_DIR <- GRID_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
msg <- function(fmt, ...) cat(sprintf(fmt, ...), "\n")

num_or_na <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  z <- suppressWarnings(as.numeric(x[1]))
  if (is.finite(z)) z else NA_real_
}
chr_or_na <- function(x) {
  if (length(x) == 0L) return(NA_character_)
  z <- as.character(x[1])
  if (!is.na(z) && nzchar(z)) z else NA_character_
}

if (!dir.exists(GRID_DIR)) stopf("Missing grid directory: %s", GRID_DIR)
spec_dirs <- list.dirs(GRID_DIR, recursive = FALSE, full.names = TRUE)
spec_dirs <- spec_dirs[file.exists(file.path(spec_dirs, "RUN_METADATA.csv"))]
if (!length(spec_dirs)) stopf("No refinement specifications found.")

rows <- vector("list", length(spec_dirs))
local_parts <- list()
quarter_parts <- list()
lp <- 0L
qp <- 0L

for (ii in seq_along(spec_dirs)) {
  d <- spec_dirs[ii]
  spec_id <- basename(d)

  meta <- read.csv(file.path(d, "RUN_METADATA.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(meta) != 1L) stopf("Malformed metadata in %s", d)

  state_scale <- num_or_na(meta$StateScale)
  prior_scale <- num_or_na(meta$PriorScale)
  exit_code <- suppressWarnings(as.integer(meta$ExitCode[1]))
  if (!is.finite(exit_code)) exit_code <- NA_integer_

  gate_path <- file.path(d, "00_tvp_smoke_gate.csv")
  gate_exists <- file.exists(gate_path)
  gate <- if (gate_exists) {
    read.csv(gate_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else data.frame(stringsAsFactors = FALSE)

  gv <- function(name, numeric = FALSE) {
    if (!gate_exists || !name %in% names(gate)) {
      return(if (numeric) NA_real_ else NA_character_)
    }
    if (numeric) num_or_na(gate[[name]]) else chr_or_na(gate[[name]])
  }

  gate_status <- gv("Status")
  if (is.na(gate_status)) gate_status <- "ERROR"

  op_status <- if (identical(exit_code, 0L)) {
    "COMPLETED_READY"
  } else if (identical(exit_code, 2L) && gate_exists) {
    "COMPLETED_MODEL_FAIL"
  } else {
    "OPERATIONAL_ERROR"
  }

  worst_country <- NA_character_
  worst_country_share <- NA_real_
  local_path <- file.path(d, "05_local_posterior_stability.csv")
  if (file.exists(local_path)) {
    z <- read.csv(local_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(z)) {
      z$SpecID <- spec_id
      z$StateScale <- state_scale
      z$PriorScale <- prior_scale
      lp <- lp + 1L
      local_parts[[lp]] <- z
      ss <- suppressWarnings(as.numeric(z$StableShare))
      ok <- is.finite(ss)
      if (any(ok)) {
        jj <- which.min(ss[ok])
        zz <- z[ok, , drop = FALSE]
        worst_country <- as.character(zz$Country[jj])
        worst_country_share <- as.numeric(zz$StableShare[jj])
      }
    }
  }

  worst_quarter <- NA_character_
  worst_quarter_share <- NA_real_
  quarter_path <- file.path(d, "01_posterior_stability_by_quarter.csv")
  if (file.exists(quarter_path)) {
    q <- read.csv(quarter_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(q)) {
      q$SpecID <- spec_id
      q$StateScale <- state_scale
      q$PriorScale <- prior_scale
      qp <- qp + 1L
      quarter_parts[[qp]] <- q
      ss <- suppressWarnings(as.numeric(q$StableShare))
      ok <- is.finite(ss)
      if (any(ok)) {
        jj <- which.min(ss[ok])
        qq <- q[ok, , drop = FALSE]
        worst_quarter <- as.character(qq$Quarter[jj])
        worst_quarter_share <- as.numeric(qq$StableShare[jj])
      }
    }
  }

  rows[[ii]] <- data.frame(
    SpecID = spec_id,
    StateScale = state_scale,
    PriorScale = prior_scale,
    ExitCode = exit_code,
    OperationalStatus = op_status,
    GateStatus = gate_status,
    Draws = gv("Draws", TRUE),
    Sample = gv("Sample"),
    PosteriorStableShare = gv("PosteriorStableShare", TRUE),
    G0OKShare = gv("G0OKShare", TRUE),
    FirstHalfStableShare = gv("FirstHalfStableShare", TRUE),
    SecondHalfStableShare = gv("SecondHalfStableShare", TRUE),
    SplitDifference = gv("SplitDifference", TRUE),
    MaxAbsStandardizedBeta = gv("MaxAbsStandardizedBeta", TRUE),
    WorstCountry = worst_country,
    WorstCountryStableShare = worst_country_share,
    WorstQuarter = worst_quarter,
    WorstQuarterStableShare = worst_quarter_share,
    Reason = gv("Reason"),
    stringsAsFactors = FALSE
  )
}

grid <- do.call(rbind, rows)
grid$ReadyCandidate <- (
  grid$GateStatus == "READY" &
  grid$ExitCode == 0L &
  is.finite(grid$PosteriorStableShare) &
  grid$PosteriorStableShare >= 0.90
)

# Prefer the least restrictive passing specification:
# larger StateScale first, then larger PriorScale.
ord <- order(
  !grid$ReadyCandidate,
  -grid$StateScale,
  -grid$PriorScale,
  -grid$PosteriorStableShare,
  grid$SplitDifference,
  na.last = TRUE
)
grid <- grid[ord, , drop = FALSE]
row.names(grid) <- NULL
grid$Recommended <- FALSE
ri <- which(grid$ReadyCandidate)
if (length(ri)) grid$Recommended[ri[1]] <- TRUE

write.csv(grid,
          file.path(OUT_DIR, "00_tvp_prior_refinement_summary.csv"),
          row.names = FALSE)
write.csv(grid[grid$ReadyCandidate, , drop = FALSE],
          file.path(OUT_DIR, "01_ready_candidates.csv"),
          row.names = FALSE)

if (length(local_parts)) {
  x <- do.call(rbind, local_parts)
  row.names(x) <- NULL
  write.csv(x, file.path(OUT_DIR, "02_local_stability_by_spec.csv"),
            row.names = FALSE)
}
if (length(quarter_parts)) {
  x <- do.call(rbind, quarter_parts)
  row.names(x) <- NULL
  write.csv(x, file.path(OUT_DIR, "03_quarter_stability_by_spec.csv"),
            row.names = FALSE)
}

n_ready <- sum(grid$ReadyCandidate, na.rm = TRUE)
n_error <- sum(grid$OperationalStatus == "OPERATIONAL_ERROR", na.rm = TRUE)

if (n_ready > 0L) {
  rec <- grid[which(grid$Recommended)[1], , drop = FALSE]
  txt <- c(
    "RECOMMENDED SECOND-ROUND TVP SPECIFICATION",
    "==========================================",
    sprintf("SpecID: %s", rec$SpecID),
    sprintf("StateScale: %.8g", rec$StateScale),
    sprintf("PriorScale: %.8g", rec$PriorScale),
    sprintf("PosteriorStableShare: %.6f", rec$PosteriorStableShare),
    sprintf("G0OKShare: %.6f", rec$G0OKShare),
    sprintf("SplitDifference: %.6f", rec$SplitDifference),
    sprintf("Worst country: %s (%.6f)",
            rec$WorstCountry, rec$WorstCountryStableShare),
    sprintf("Worst quarter: %s (%.6f)",
            rec$WorstQuarter, rec$WorstQuarterStableShare),
    "",
    "This candidate passed the unchanged 0.90 stability gate.",
    "It is the least restrictive passing specification in this refinement grid.",
    "It remains a smoke-test candidate, not the final publication model."
  )
} else {
  finite_grid <- grid[is.finite(grid$PosteriorStableShare), , drop = FALSE]
  best <- if (nrow(finite_grid)) {
    finite_grid[order(-finite_grid$PosteriorStableShare,
                      finite_grid$SplitDifference,
                      na.last = TRUE), , drop = FALSE][1, , drop = FALSE]
  } else NULL

  txt <- c(
    "NO SECOND-ROUND SPECIFICATION PASSED ALL GATES",
    "==============================================",
    sprintf("Specifications evaluated: %d", nrow(grid)),
    sprintf("Operational errors: %d", n_error),
    if (!is.null(best))
      sprintf("Highest stable share: %.6f (%s)",
              best$PosteriorStableShare, best$SpecID)
    else NULL,
    "",
    "Do NOT lower the 0.90 stability threshold.",
    "Do NOT automatically shrink PriorScale below 0.02.",
    "If 0.02 still fails materially, move to a formal stability-restricted",
    "Bayesian TVP-GVAR design rather than continuing mechanical tuning."
  )
}

writeLines(txt, file.path(OUT_DIR, "RECOMMENDED_SPEC.txt"))

readme <- c(
  "TVP-GVAR SECOND-ROUND PRIOR REFINEMENT",
  "======================================",
  sprintf("Specifications evaluated: %d", nrow(grid)),
  sprintf("READY candidates: %d", n_ready),
  sprintf("Operational errors: %d", n_error),
  "",
  "Focused grid:",
  "- StateScale: 1e-6, 1e-7",
  "- PriorScale: 0.10, 0.05, 0.02",
  "",
  "The original 0.90 posterior-stability gate is preserved.",
  "Exit code 2 is treated as model-gate failure, not workflow failure.",
  "Other nonzero exit codes are operational errors.",
  "This is a smoke test, not final MCMC estimation."
)
writeLines(readme, file.path(OUT_DIR, "README_tvp_prior_refinement.txt"))

msg("Prior-refinement summary complete.")
msg("READY candidates: %d / %d", n_ready, nrow(grid))
cat(paste(txt, collapse = "\n"), "\n")
