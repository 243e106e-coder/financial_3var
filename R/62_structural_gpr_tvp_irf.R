#!/usr/bin/env Rscript

# =============================================================================
# 62_structural_gpr_tvp_irf.R
#
# Time-specific posterior GPR IRFs for the formal financial TVP-GVAR.
#
# IDENTIFICATION
# --------------
# GPR is an observed global/exogenous driver in every local country equation,
# not an endogenous domestic financial variable. Therefore a positive GPR shock
# is identified directly through the estimated gpr_0 / gpr_L1 coefficient
# vectors. No arbitrary Cholesky ordering of r/de/deq is needed for this shock.
#
# At each event anchor t and posterior draw d:
#
#   G0_{t,d} x_h = G1_{t,d} x_{h-1}
#                  + c0_{t,d} * shock at h=0
#                  + c1_{t,d} * shock at h=1
#
# with time-t coefficients frozen over the response horizon (standard
# time-specific TVP-IRF convention).
#
# GPR normalization:
#   panel uses LN_GPR_QMEAN, so a +10% GPR shock equals log(1.10).
#
# Only stable posterior draws with:
#   spectral radius(F_{t,d}) < 1
#   G0 rcond >= threshold
# are used for reported IRF quantiles.
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
EVENT_FILE <- get_env_chr(
  "FIN3_EVENT_CALENDAR",
  file.path(RESULTS_DIR, "events", "00_event_calendar.csv")
)

NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 1500))
HORIZON <- as.integer(get_env_num("FIN3_IRF_HORIZON", 12))
GPR_SHOCK_PCT <- get_env_num("FIN3_GPR_SHOCK_PCT", 10)
MIN_STABLE_SHARE <- get_env_num("FIN3_IRF_MIN_STABLE_SHARE", 0.90)
MIN_G0_RCOND <- get_env_num("FIN3_IRF_MIN_G0_RCOND", 1e-10)
IRF_CORES <- as.integer(get_env_num("FIN3_IRF_CORES", 2))

if (NCHAINS < 2L) stopf("Formal IRF requires at least two chains.")
if (EXPECTED_STORED < 200L) stopf("Expected stored draws too small.")
if (HORIZON < 1L) stopf("FIN3_IRF_HORIZON must be >=1.")
if (!(GPR_SHOCK_PCT > 0)) stopf("FIN3_GPR_SHOCK_PCT must be positive.")
if (!(MIN_STABLE_SHARE > 0 && MIN_STABLE_SHARE <= 1)) {
  stopf("FIN3_IRF_MIN_STABLE_SHARE must be in (0,1].")
}
if (IRF_CORES < 1L) IRF_CORES <- 1L

CONV_GATE <- file.path(CONV_DIR, "00_formal_mcmc_gate.csv")
if (!file.exists(CONV_GATE)) stopf("Formal MCMC convergence gate missing: %s", CONV_GATE)

cg <- read.csv(CONV_GATE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(cg) != 1L || !"Status" %in% names(cg)) stopf("Malformed formal MCMC gate.")
if (!identical(trimws(cg$Status[1]), "READY_FOR_FORMAL_IRF")) {
  stopf("Formal MCMC gate is not READY_FOR_FORMAL_IRF.")
}

if (!file.exists(EVENT_FILE)) stopf("Missing event calendar: %s", EVENT_FILE)

OUT <- file.path(RESULTS_DIR, "formal_irf")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Discover and load compact country-chain posterior parts
# =============================================================================

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(files)) stopf("No formal TVP posterior parts found.")

parts <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)

for (cc in COUNTRIES) {
  parts[[cc]] <- vector("list", NCHAINS)
  for (ch in seq_len(NCHAINS)) {
    candidates <- files[
      grepl(
        sprintf("formal_tvp_%s_chain%d\\.rds$", cc, ch),
        files
      )
    ]
    if (length(candidates) != 1L) {
      stopf("Expected one posterior part for %s chain %d; found %d.", cc, ch, length(candidates))
    }
    z <- readRDS(candidates)
    if (z$meta$stored_draws != EXPECTED_STORED) {
      stopf("Stored draw mismatch: %s chain %d", cc, ch)
    }
    parts[[cc]][[ch]] <- z
  }
}

ref <- parts[[COUNTRIES[1]]][[1]]
anchors <- ref$anchors
KEY_TERMS <- ref$meta$key_terms

required_key_terms <- c(
  paste0(VARS, "_L1"),
  paste0(VARS, "_star_0"),
  paste0(VARS, "_star_L1"),
  "gpr_0", "gpr_L1"
)
if (!all(required_key_terms %in% KEY_TERMS)) {
  stopf("Formal posterior parts do not contain all GPR-IRF coefficient blocks.")
}

# Ensure all country/chain parts use identical anchor ordering and key terms.
for (cc in COUNTRIES) {
  for (ch in seq_len(NCHAINS)) {
    z <- parts[[cc]][[ch]]
    if (!identical(z$meta$key_terms, KEY_TERMS)) {
      stopf("Key-term ordering mismatch: %s chain %d", cc, ch)
    }
    if (!identical(
      paste(z$anchors$EventID, z$anchors$AnchorType, z$anchors$AnchorQuarter),
      paste(anchors$EventID, anchors$AnchorType, anchors$AnchorQuarter)
    )) {
      stopf("Event-anchor ordering mismatch: %s chain %d", cc, ch)
    }
  }
}

N <- length(COUNTRIES)
K <- length(VARS)
NK <- N * K
DRAWS_TOTAL <- NCHAINS * EXPECTED_STORED

W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])

selector <- function(i) {
  S <- matrix(0, K, NK)
  S[, ((i - 1L) * K + 1L):(i * K)] <- diag(K)
  S
}

star_map <- function(i) {
  R <- matrix(0, K, NK)
  for (j in seq_len(N)) {
    for (v in seq_len(K)) {
      R[v, (j - 1L) * K + v] <- W[i, j]
    }
  }
  R
}

SELECTORS <- lapply(seq_len(N), selector)
STAR_MAPS <- lapply(seq_len(N), star_map)

term_dom <- paste0(VARS, "_L1")
term_star0 <- paste0(VARS, "_star_0")
term_star1 <- paste0(VARS, "_star_L1")

idx_dom <- match(term_dom, KEY_TERMS)
idx_star0 <- match(term_star0, KEY_TERMS)
idx_star1 <- match(term_star1, KEY_TERMS)
idx_gpr0 <- match("gpr_0", KEY_TERMS)
idx_gpr1 <- match("gpr_L1", KEY_TERMS)

if (any(is.na(c(idx_dom, idx_star0, idx_star1, idx_gpr0, idx_gpr1)))) {
  stopf("Failed to map formal event coefficient terms.")
}

# GPR is a logged index in the clean panel.
if (!grepl("^LN_", toupper(GPR_COLUMN))) {
  stopf(
    "Expected logged GPR column for percent-shock normalization; configured GPR_COLUMN=%s",
    GPR_COLUMN
  )
}
shock_size <- log1p(GPR_SHOCK_PCT / 100)

# =============================================================================
# 2. Process one event anchor
# =============================================================================

qsum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(
      p05 = NA_real_, p16 = NA_real_, median = NA_real_,
      p84 = NA_real_, p95 = NA_real_,
      mean = NA_real_, prob_positive = NA_real_, n = 0
    ))
  }
  qq <- stats::quantile(x, c(.05, .16, .50, .84, .95), names = FALSE)
  c(
    p05 = qq[1],
    p16 = qq[2],
    median = qq[3],
    p84 = qq[4],
    p95 = qq[5],
    mean = mean(x),
    prob_positive = mean(x > 0),
    n = length(x)
  )
}

process_anchor <- function(a) {
  meta <- anchors[a, , drop = FALSE]

  irf <- array(
    NA_real_,
    c(DRAWS_TOTAL, HORIZON + 1L, NK),
    dimnames = list(
      Draw = seq_len(DRAWS_TOTAL),
      Horizon = 0:HORIZON,
      Variable = as.vector(
        outer(COUNTRIES, VARS, paste, sep = "__")
      )
    )
  )

  rho_vec <- rep(NA_real_, DRAWS_TOTAL)
  rcond_vec <- rep(NA_real_, DRAWS_TOTAL)
  local_rho <- matrix(
    NA_real_,
    DRAWS_TOTAL,
    N,
    dimnames = list(NULL, COUNTRIES)
  )

  global_draw <- 0L

  for (ch in seq_len(NCHAINS)) {
    for (dd in seq_len(EXPECTED_STORED)) {
      global_draw <- global_draw + 1L

      G0 <- matrix(0, NK, NK)
      G1 <- matrix(0, NK, NK)
      c0 <- numeric(NK)
      c1 <- numeric(NK)

      for (i in seq_len(N)) {
        cc <- COUNTRIES[i]
        coef <- parts[[cc]][[ch]]$event_coef[dd, a, , , drop = FALSE]
        coef <- matrix(
          coef,
          nrow = K,
          ncol = length(KEY_TERMS),
          dimnames = list(VARS, KEY_TERMS)
        )

        Ai <- coef[, idx_dom, drop = FALSE]
        B0i <- coef[, idx_star0, drop = FALSE]
        B1i <- coef[, idx_star1, drop = FALSE]

        local_rho[global_draw, i] <- max(
          Mod(eigen(Ai, only.values = TRUE)$values)
        )

        rr <- ((i - 1L) * K + 1L):(i * K)
        G0[rr, ] <- SELECTORS[[i]] - B0i %*% STAR_MAPS[[i]]
        G1[rr, ] <- Ai %*% SELECTORS[[i]] + B1i %*% STAR_MAPS[[i]]

        c0[rr] <- coef[, idx_gpr0]
        c1[rr] <- coef[, idx_gpr1]
      }

      rc <- tryCatch(rcond(G0), error = function(e) NA_real_)
      rcond_vec[global_draw] <- rc

      if (!is.finite(rc) || rc < MIN_G0_RCOND) next

      impact_map <- tryCatch(solve(G0), error = function(e) NULL)
      if (is.null(impact_map) || any(!is.finite(impact_map))) next

      Fmat <- impact_map %*% G1
      rho <- tryCatch(
        max(Mod(eigen(Fmat, only.values = TRUE)$values)),
        error = function(e) NA_real_
      )
      rho_vec[global_draw] <- rho

      if (!is.finite(rho) || rho >= 1) next

      # h=0: contemporaneous GPR loading.
      xh <- as.numeric(impact_map %*% (c0 * shock_size))
      irf[global_draw, 1L, ] <- xh

      # h=1: dynamic propagation plus lagged GPR loading.
      if (HORIZON >= 1L) {
        xh <- as.numeric(
          Fmat %*% xh +
          impact_map %*% (c1 * shock_size)
        )
        irf[global_draw, 2L, ] <- xh
      }

      # h>=2: no further GPR innovation; propagate the system.
      if (HORIZON >= 2L) {
        for (hh in 2:HORIZON) {
          xh <- as.numeric(Fmat %*% xh)
          irf[global_draw, hh + 1L, ] <- xh
        }
      }
    }
  }

  stable <- is.finite(rho_vec) & rho_vec < 1
  g0_ok <- is.finite(rcond_vec) & rcond_vec >= MIN_G0_RCOND
  valid <- stable & g0_ok

  # Posterior IRF summary.
  rows <- vector("list", NK * (HORIZON + 1L))
  pos <- 0L

  for (i in seq_len(N)) {
    for (v in seq_len(K)) {
      gv <- (i - 1L) * K + v
      for (hh in 0:HORIZON) {
        pos <- pos + 1L
        qs <- qsum(irf[valid, hh + 1L, gv])
        rows[[pos]] <- data.frame(
          EventID = meta$EventID,
          EventSet = meta$EventSet,
          EventLabel = meta$EventLabel,
          ShockFamily = meta$ShockFamily,
          EventQuarter = meta$EventQuarter,
          AnchorType = meta$AnchorType,
          AnchorQuarter = meta$AnchorQuarter,
          Country = COUNTRIES[i],
          ResponseVariable = VARS[v],
          Horizon = hh,
          p05 = qs["p05"],
          p16 = qs["p16"],
          median = qs["median"],
          p84 = qs["p84"],
          p95 = qs["p95"],
          mean = qs["mean"],
          prob_positive = qs["prob_positive"],
          StablePosteriorDraws = qs["n"],
          TotalPosteriorDraws = DRAWS_TOTAL,
          GPRShockPct = GPR_SHOCK_PCT,
          GPRShockLogPoints = shock_size,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  irf_summary <- do.call(rbind, rows)

  rho_finite <- rho_vec[is.finite(rho_vec)]
  rc_finite <- rcond_vec[is.finite(rcond_vec)]

  anchor_stability <- data.frame(
    EventID = meta$EventID,
    EventSet = meta$EventSet,
    AnchorType = meta$AnchorType,
    AnchorQuarter = meta$AnchorQuarter,
    TotalDraws = DRAWS_TOTAL,
    FiniteRhoShare = mean(is.finite(rho_vec)),
    StableShare = mean(stable),
    G0OKShare = mean(g0_ok),
    ValidIRFShare = mean(valid),
    RhoMedian = if (length(rho_finite)) median(rho_finite) else NA_real_,
    RhoP95 = if (length(rho_finite)) unname(quantile(rho_finite, .95)) else NA_real_,
    RhoMax = if (length(rho_finite)) max(rho_finite) else NA_real_,
    MinG0Rcond = if (length(rc_finite)) min(rc_finite) else NA_real_,
    stringsAsFactors = FALSE
  )

  local_rows <- do.call(rbind, lapply(seq_len(N), function(i) {
    rr <- local_rho[, i]
    rr <- rr[is.finite(rr)]
    data.frame(
      EventID = meta$EventID,
      EventSet = meta$EventSet,
      AnchorType = meta$AnchorType,
      AnchorQuarter = meta$AnchorQuarter,
      Country = COUNTRIES[i],
      LocalStableShare = if (length(rr)) mean(rr < 1) else NA_real_,
      LocalRhoMedian = if (length(rr)) median(rr) else NA_real_,
      LocalRhoP95 = if (length(rr)) unname(quantile(rr, .95)) else NA_real_,
      LocalRhoMax = if (length(rr)) max(rr) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  # Draw-by-draw cumulative REER log-change effect, computed BEFORE quantiles.
  cumulative_reer_rows <- vector("list", N * (HORIZON + 1L))
  cp <- 0L
  de_idx <- match("de", VARS)

  for (i in seq_len(N)) {
    gv <- (i - 1L) * K + de_idx
    z <- irf[valid, , gv, drop = FALSE]
    z <- matrix(z, nrow = sum(valid), ncol = HORIZON + 1L)
    if (nrow(z)) {
      zcum <- t(apply(z, 1, cumsum))
    } else {
      zcum <- matrix(NA_real_, 0, HORIZON + 1L)
    }

    for (hh in 0:HORIZON) {
      cp <- cp + 1L
      qs <- qsum(if (nrow(zcum)) zcum[, hh + 1L] else numeric())
      cumulative_reer_rows[[cp]] <- data.frame(
        EventID = meta$EventID,
        EventSet = meta$EventSet,
        AnchorType = meta$AnchorType,
        AnchorQuarter = meta$AnchorQuarter,
        Country = COUNTRIES[i],
        ResponseVariable = "de_cumulative_log_reer",
        Horizon = hh,
        p05 = qs["p05"],
        p16 = qs["p16"],
        median = qs["median"],
        p84 = qs["p84"],
        p95 = qs["p95"],
        mean = qs["mean"],
        prob_positive = qs["prob_positive"],
        StablePosteriorDraws = qs["n"],
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    irf_summary = irf_summary,
    stability = anchor_stability,
    local_stability = local_rows,
    cumulative_reer = do.call(rbind, cumulative_reer_rows)
  )
}

# On GitHub Linux runners mclapply allows two anchors to be processed in
# parallel while sharing loaded posterior parts copy-on-write.
anchor_ids <- seq_len(nrow(anchors))
if (.Platform$OS.type == "unix" && IRF_CORES > 1L) {
  ans <- parallel::mclapply(
    anchor_ids,
    process_anchor,
    mc.cores = min(IRF_CORES, length(anchor_ids)),
    mc.preschedule = TRUE
  )
} else {
  ans <- lapply(anchor_ids, process_anchor)
}

irf_summary <- do.call(rbind, lapply(ans, `[[`, "irf_summary"))
anchor_stability <- do.call(rbind, lapply(ans, `[[`, "stability"))
local_stability <- do.call(rbind, lapply(ans, `[[`, "local_stability"))
cum_reer <- do.call(rbind, lapply(ans, `[[`, "cumulative_reer"))

write.csv(
  irf_summary,
  file.path(OUT, "irf_posterior_summary.csv"),
  row.names = FALSE
)
write.csv(
  anchor_stability,
  file.path(OUT, "01_irf_stability_by_anchor.csv"),
  row.names = FALSE
)
write.csv(
  local_stability,
  file.path(OUT, "02_local_stability_by_anchor_country.csv"),
  row.names = FALSE
)
write.csv(
  cum_reer,
  file.path(OUT, "03_cumulative_reer_posterior_summary.csv"),
  row.names = FALSE
)

# Core-event compact impact table.
impact <- irf_summary[
  irf_summary$EventSet == "CORE" &
  irf_summary$AnchorType == "EVENT_QUARTER_t0" &
  irf_summary$Horizon == 0,
  ,
  drop = FALSE
]
write.csv(
  impact,
  file.path(OUT, "04_core_event_impact_irf.csv"),
  row.names = FALSE
)

# =============================================================================
# 3. Formal IRF stability gate
# =============================================================================

stable_ok <- all(anchor_stability$StableShare >= MIN_STABLE_SHARE)
g0_ok <- all(anchor_stability$G0OKShare >= 0.99)
finite_ok <- all(anchor_stability$FiniteRhoShare >= 0.99)
draws_ok <- all(irf_summary$StablePosteriorDraws >= 200)

ready <- stable_ok && g0_ok && finite_ok && draws_ok
status <- if (ready) "READY_FOR_IRF_AUDIT" else "FAIL"

reasons <- c(
  if (!stable_ok) sprintf(
    "at least one anchor stable share < %.3f",
    MIN_STABLE_SHARE
  ) else NULL,
  if (!g0_ok) "at least one anchor G0-OK share < 0.99" else NULL,
  if (!finite_ok) "at least one anchor finite-rho share < 0.99" else NULL,
  if (!draws_ok) "at least one IRF cell has fewer than 200 valid stable draws" else NULL
)
if (!length(reasons)) reasons <- "all formal GPR-IRF gates passed"

gate <- data.frame(
  Status = status,
  MainNetwork = MAIN_NETWORK,
  PosteriorDrawsTotal = DRAWS_TOTAL,
  Chains = NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  EventAnchors = nrow(anchors),
  Horizon = HORIZON,
  GPRShockPct = GPR_SHOCK_PCT,
  GPRShockLogPoints = shock_size,
  MinStableShareRequired = MIN_STABLE_SHARE,
  WorstAnchorStableShare = min(anchor_stability$StableShare),
  WorstAnchorG0OKShare = min(anchor_stability$G0OKShare),
  WorstAnchorFiniteRhoShare = min(anchor_stability$FiniteRhoShare),
  MinValidIRFDraws = min(irf_summary$StablePosteriorDraws),
  Reason = paste(reasons, collapse = "; "),
  stringsAsFactors = FALSE
)

write.csv(
  gate,
  file.path(OUT, "00_formal_irf_gate.csv"),
  row.names = FALSE
)

readme <- c(
  sprintf("FORMAL TIME-SPECIFIC GPR TVP-GVAR IRF: %s", status),
  "======================================================",
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf("Posterior draws pooled across chains: %d", DRAWS_TOTAL),
  sprintf("Event anchors: %d", nrow(anchors)),
  sprintf("Horizon: 0-%d quarters", HORIZON),
  sprintf("Positive GPR shock: %.2f%% = %.8f log points", GPR_SHOCK_PCT, shock_size),
  sprintf("Worst event-anchor stable share: %.6f", gate$WorstAnchorStableShare),
  sprintf("Worst event-anchor G0-OK share: %.6f", gate$WorstAnchorG0OKShare),
  sprintf("Minimum valid stable draws in any IRF cell: %d", gate$MinValidIRFDraws),
  sprintf("Reason: %s", gate$Reason),
  "",
  "Identification:",
  "- GPR is external/global in the local country equations.",
  "- A positive GPR innovation is injected through gpr_0 and gpr_L1.",
  "- No domestic Cholesky ordering is imposed for this GPR shock.",
  "- Brent receives no simultaneous innovation in this experiment.",
  "- TVP coefficients are frozen at the selected event anchor over the response horizon.",
  "",
  "Outputs:",
  "- irf_posterior_summary.csv  [input to R/55 and R/56]",
  "- 01_irf_stability_by_anchor.csv",
  "- 02_local_stability_by_anchor_country.csv",
  "- 03_cumulative_reer_posterior_summary.csv",
  "- 04_core_event_impact_irf.csv",
  "- 00_formal_irf_gate.csv"
)

writeLines(readme, file.path(OUT, "README_formal_irf.txt"))
cat(paste(readme, collapse = "\n"), "\n")
