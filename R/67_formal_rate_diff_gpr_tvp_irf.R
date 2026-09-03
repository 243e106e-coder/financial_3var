#!/usr/bin/env Rscript

# =============================================================================
# 67_formal_rate_diff_gpr_tvp_irf.R
#
# Formal time-specific GPR IRFs for the accepted Delta-r financial TVP-GVAR.
#
# This is the Delta-r successor to R/62_structural_gpr_tvp_irf.R.
# R/62 is deliberately left unchanged for auditability.
#
# IDENTIFICATION / DYNAMICS
# -------------------------
# Exactly as in R/62:
# - GPR is an observed global/exogenous driver in every local country equation.
# - The positive GPR innovation enters through gpr_0 and gpr_L1.
# - No Cholesky ordering is imposed among r / de / deq for the GPR shock.
# - At each event anchor and posterior draw, time-t coefficients are frozen over
#   the response horizon.
# - G0, G1 and F = solve(G0) %*% G1 are reconstructed draw by draw.
# - Reported IRFs use only draws with G0 rcond above the threshold and rho(F)<1.
#
# DELTA-r INTERPRETATION
# ----------------------
# In the accepted formal specification:
#   r   = Delta RATE_LEVEL
#   de  = REER_DLOG
#   deq = EQ_RETURN
#
# Therefore this script reports BOTH:
# 1) raw r IRF = response of Delta r;
# 2) exact draw-level cumulative sum of the Delta-r IRF, interpreted as the
#    cumulative change in the interest-rate level relative to baseline.
#
# IMPORTANT: cumulative summaries are computed DRAW BY DRAW BEFORE posterior
# quantiles are taken.  We do NOT sum posterior medians or quantiles.
#
# de is handled analogously: its cumulative draw-level sum is a cumulative
# REER log-level effect. deq remains an equity-return response and is not
# relabelled as an equity-price level effect.
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
RATE_DIFF_DECISION_DIR <- get_env_chr(
  "FIN3_RATE_DIFF_DECISION_DIR",
  file.path(RESULTS_DIR, "rate_diff_decision")
)
EVENT_FILE <- get_env_chr(
  "FIN3_EVENT_CALENDAR",
  file.path(RESULTS_DIR, "events", "00_event_calendar.csv")
)
EXPECTED_PANEL_BASENAME <- get_env_chr(
  "FIN3_RATE_DIFF_PANEL_BASENAME",
  "panel_domestic_fin3_rate_diff.csv"
)

NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 2000))
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

# =============================================================================
# 0. Formal prerequisite gates
# =============================================================================

CONV_GATE <- file.path(CONV_DIR, "00_formal_mcmc_gate.csv")
if (!file.exists(CONV_GATE)) stopf("Formal MCMC convergence gate missing: %s", CONV_GATE)

cg <- read.csv(CONV_GATE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(cg) != 1L || !"Status" %in% names(cg)) stopf("Malformed formal MCMC gate.")
if (!identical(trimws(as.character(cg$Status[1])), "READY_FOR_FORMAL_IRF")) {
  stopf("Formal MCMC gate is not READY_FOR_FORMAL_IRF.")
}

RATE_DIFF_GATE <- file.path(RATE_DIFF_DECISION_DIR, "00_rate_diff_formal_gate.csv")
if (!file.exists(RATE_DIFF_GATE)) {
  stopf("Accepted Delta-r decision gate missing: %s", RATE_DIFF_GATE)
}
rdg <- read.csv(RATE_DIFF_GATE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(rdg) != 1L || !"Status" %in% names(rdg)) stopf("Malformed Delta-r decision gate.")
if (!identical(trimws(as.character(rdg$Status[1])), "RATE_DIFF_READY_FOR_IRF_RECODE")) {
  stopf("Delta-r decision is not RATE_DIFF_READY_FOR_IRF_RECODE.")
}
if ("BaselineMinStableShare" %in% names(rdg)) {
  if (!is.finite(rdg$BaselineMinStableShare[1]) || rdg$BaselineMinStableShare[1] < MIN_STABLE_SHARE) {
    stopf("Accepted Delta-r dynamic-stability share is below the formal IRF threshold.")
  }
}

if (!file.exists(EVENT_FILE)) stopf("Missing event calendar: %s", EVENT_FILE)

OUT <- file.path(RESULTS_DIR, "formal_irf_rate_diff")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Discover and validate the accepted 14 x 4 posterior lineage
# =============================================================================

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

expected_n <- length(COUNTRIES) * NCHAINS
if (length(files) != expected_n) {
  stopf("Expected exactly %d formal posterior RDS files; found %d.", expected_n, length(files))
}

parts <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)
lineage_rows <- vector("list", expected_n)
li <- 0L

for (cc in COUNTRIES) {
  parts[[cc]] <- vector("list", NCHAINS)
  for (ch in seq_len(NCHAINS)) {
    candidates <- files[grepl(sprintf("formal_tvp_%s_chain%d\\.rds$", cc, ch), files)]
    if (length(candidates) != 1L) {
      stopf("Expected one posterior part for %s chain %d; found %d.", cc, ch, length(candidates))
    }

    z <- readRDS(candidates)

    if (!identical(as.character(z$meta$country), cc)) {
      stopf("Country metadata mismatch: expected %s chain %d.", cc, ch)
    }
    if (as.integer(z$meta$chain) != ch) {
      stopf("Chain metadata mismatch: %s chain %d.", cc, ch)
    }
    if (as.integer(z$meta$stored_draws) != EXPECTED_STORED) {
      stopf("Stored draw mismatch: %s chain %d.", cc, ch)
    }
    if (!identical(basename(as.character(z$meta$source_panel)), EXPECTED_PANEL_BASENAME)) {
      stopf("Non-Delta-r posterior detected: %s chain %d.", cc, ch)
    }

    parts[[cc]][[ch]] <- z

    li <- li + 1L
    lineage_rows[[li]] <- data.frame(
      Country = cc,
      Chain = ch,
      Seed = as.integer(z$meta$seed),
      Burn = as.integer(z$meta$burn),
      KeepInternal = as.integer(z$meta$keep_internal),
      StoredDraws = as.integer(z$meta$stored_draws),
      SourcePanel = basename(as.character(z$meta$source_panel)),
      File = candidates,
      stringsAsFactors = FALSE
    )
  }
}

lineage <- do.call(rbind, lineage_rows)
write.csv(lineage, file.path(OUT, "00_posterior_lineage.csv"), row.names = FALSE)

key <- paste(lineage$Country, lineage$Chain, sep = "||")
expected_key <- paste(
  rep(COUNTRIES, each = NCHAINS),
  rep(seq_len(NCHAINS), times = length(COUNTRIES)),
  sep = "||"
)
if (anyDuplicated(key) || !setequal(key, expected_key)) {
  stopf("Formal posterior lineage is not an exact 14 x %d grid.", NCHAINS)
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

r_idx <- match("r", VARS)
de_idx <- match("de", VARS)
deq_idx <- match("deq", VARS)
if (any(is.na(c(r_idx, de_idx, deq_idx)))) {
  stopf("Expected VARS to contain r, de, deq.")
}

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

if (!grepl("^LN_", toupper(GPR_COLUMN))) {
  stopf(
    "Expected logged GPR column for percent-shock normalization; configured GPR_COLUMN=%s",
    GPR_COLUMN
  )
}
shock_size <- log1p(GPR_SHOCK_PCT / 100)

# =============================================================================
# 2. Posterior helper functions
# =============================================================================

qsum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(
      p05 = NA_real_, p16 = NA_real_, median = NA_real_,
      p84 = NA_real_, p95 = NA_real_, mean = NA_real_,
      prob_positive = NA_real_, n = 0
    ))
  }
  qq <- stats::quantile(x, c(.05, .16, .50, .84, .95), names = FALSE)
  c(
    p05 = qq[1], p16 = qq[2], median = qq[3],
    p84 = qq[4], p95 = qq[5], mean = mean(x),
    prob_positive = mean(x > 0), n = length(x)
  )
}

summary_row <- function(meta, country, variable, horizon, qs, scale_label) {
  data.frame(
    EventID = meta$EventID,
    EventSet = meta$EventSet,
    EventLabel = meta$EventLabel,
    ShockFamily = meta$ShockFamily,
    EventQuarter = meta$EventQuarter,
    AnchorType = meta$AnchorType,
    AnchorQuarter = meta$AnchorQuarter,
    Country = country,
    ResponseVariable = variable,
    ResponseScale = scale_label,
    Horizon = horizon,
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

cumulate_draw_paths <- function(z) {
  if (!nrow(z)) return(matrix(NA_real_, 0, ncol(z)))
  t(apply(z, 1, cumsum))
}

# =============================================================================
# 3. Process one event anchor
# =============================================================================

process_anchor <- function(a) {
  meta <- anchors[a, , drop = FALSE]

  irf <- array(
    NA_real_,
    c(DRAWS_TOTAL, HORIZON + 1L, NK),
    dimnames = list(
      Draw = seq_len(DRAWS_TOTAL),
      Horizon = 0:HORIZON,
      Variable = NULL
    )
  )

  rho_vec <- rep(NA_real_, DRAWS_TOTAL)
  rcond_vec <- rep(NA_real_, DRAWS_TOTAL)
  local_rho <- matrix(NA_real_, DRAWS_TOTAL, N, dimnames = list(NULL, COUNTRIES))

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

        local_rho[global_draw, i] <- max(Mod(eigen(Ai, only.values = TRUE)$values))

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

      # h = 0: contemporaneous GPR loading.
      xh <- as.numeric(impact_map %*% (c0 * shock_size))
      irf[global_draw, 1L, ] <- xh

      # h = 1: propagation plus lagged GPR loading.
      if (HORIZON >= 1L) {
        xh <- as.numeric(Fmat %*% xh + impact_map %*% (c1 * shock_size))
        irf[global_draw, 2L, ] <- xh
      }

      # h >= 2: no additional GPR innovation.
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

  # ---------------------------------------------------------------------------
  # Raw IRFs: r = Delta r; de = REER_DLOG; deq = EQ_RETURN
  # ---------------------------------------------------------------------------

  raw_rows <- vector("list", NK * (HORIZON + 1L))
  pos <- 0L

  scale_map <- c(
    r = "DELTA_RATE_LEVEL",
    de = "REER_DLOG",
    deq = "EQ_RETURN"
  )

  for (i in seq_len(N)) {
    for (v in seq_len(K)) {
      gv <- (i - 1L) * K + v
      vv <- VARS[v]
      for (hh in 0:HORIZON) {
        pos <- pos + 1L
        qs <- qsum(irf[valid, hh + 1L, gv])
        raw_rows[[pos]] <- summary_row(
          meta, COUNTRIES[i], vv, hh, qs, unname(scale_map[vv])
        )
      }
    }
  }
  irf_summary <- do.call(rbind, raw_rows)

  # ---------------------------------------------------------------------------
  # Exact draw-level cumulative interest-rate LEVEL change from Delta r.
  # ---------------------------------------------------------------------------

  cumulative_rate_rows <- vector("list", N * (HORIZON + 1L))
  cp <- 0L

  for (i in seq_len(N)) {
    gv <- (i - 1L) * K + r_idx
    z <- irf[valid, , gv, drop = FALSE]
    z <- matrix(z, nrow = sum(valid), ncol = HORIZON + 1L)
    zcum <- cumulate_draw_paths(z)

    for (hh in 0:HORIZON) {
      cp <- cp + 1L
      qs <- qsum(if (nrow(zcum)) zcum[, hh + 1L] else numeric())
      cumulative_rate_rows[[cp]] <- summary_row(
        meta,
        COUNTRIES[i],
        "r_cumulative_level_change",
        hh,
        qs,
        "CUMULATIVE_RATE_LEVEL_CHANGE_FROM_DELTA_R"
      )
    }
  }
  cumulative_rate <- do.call(rbind, cumulative_rate_rows)

  # ---------------------------------------------------------------------------
  # Exact draw-level cumulative REER log-level effect from REER_DLOG.
  # ---------------------------------------------------------------------------

  cumulative_reer_rows <- vector("list", N * (HORIZON + 1L))
  cp <- 0L

  for (i in seq_len(N)) {
    gv <- (i - 1L) * K + de_idx
    z <- irf[valid, , gv, drop = FALSE]
    z <- matrix(z, nrow = sum(valid), ncol = HORIZON + 1L)
    zcum <- cumulate_draw_paths(z)

    for (hh in 0:HORIZON) {
      cp <- cp + 1L
      qs <- qsum(if (nrow(zcum)) zcum[, hh + 1L] else numeric())
      cumulative_reer_rows[[cp]] <- summary_row(
        meta,
        COUNTRIES[i],
        "de_cumulative_log_reer",
        hh,
        qs,
        "CUMULATIVE_REER_LOG_LEVEL_EFFECT_FROM_DLOG"
      )
    }
  }
  cumulative_reer <- do.call(rbind, cumulative_reer_rows)

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

  list(
    irf_summary = irf_summary,
    cumulative_rate = cumulative_rate,
    cumulative_reer = cumulative_reer,
    stability = anchor_stability,
    local_stability = local_rows
  )
}

# =============================================================================
# 4. Run all event anchors
# =============================================================================

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
cum_rate <- do.call(rbind, lapply(ans, `[[`, "cumulative_rate"))
cum_reer <- do.call(rbind, lapply(ans, `[[`, "cumulative_reer"))
anchor_stability <- do.call(rbind, lapply(ans, `[[`, "stability"))
local_stability <- do.call(rbind, lapply(ans, `[[`, "local_stability"))

write.csv(irf_summary, file.path(OUT, "irf_posterior_summary.csv"), row.names = FALSE)
write.csv(anchor_stability, file.path(OUT, "01_irf_stability_by_anchor.csv"), row.names = FALSE)
write.csv(local_stability, file.path(OUT, "02_local_stability_by_anchor_country.csv"), row.names = FALSE)
write.csv(cum_reer, file.path(OUT, "03_cumulative_reer_posterior_summary.csv"), row.names = FALSE)

impact <- irf_summary[
  irf_summary$EventSet == "CORE" &
    irf_summary$AnchorType == "EVENT_QUARTER_t0" &
    irf_summary$Horizon == 0,
  ,
  drop = FALSE
]
write.csv(impact, file.path(OUT, "04_core_event_impact_raw_irf.csv"), row.names = FALSE)

write.csv(cum_rate, file.path(OUT, "05_cumulative_rate_level_posterior_summary.csv"), row.names = FALSE)

# Compact cumulative tables at economically interpretable horizons.
keep_h <- intersect(c(0L, 1L, 4L, 8L, 12L), 0:HORIZON)
core_cum_rate <- cum_rate[
  cum_rate$EventSet == "CORE" &
    cum_rate$AnchorType == "EVENT_QUARTER_t0" &
    cum_rate$Horizon %in% keep_h,
  ,
  drop = FALSE
]
write.csv(core_cum_rate, file.path(OUT, "06_core_event_cumulative_rate_level.csv"), row.names = FALSE)

core_cum_reer <- cum_reer[
  cum_reer$EventSet == "CORE" &
    cum_reer$AnchorType == "EVENT_QUARTER_t0" &
    cum_reer$Horizon %in% keep_h,
  ,
  drop = FALSE
]
write.csv(core_cum_reer, file.path(OUT, "07_core_event_cumulative_reer.csv"), row.names = FALSE)

# Variable interpretation manifest is explicit and machine-readable.
scale_manifest <- data.frame(
  OutputVariable = c("r", "r_cumulative_level_change", "de", "de_cumulative_log_reer", "deq"),
  Meaning = c(
    "Quarterly change in the interest-rate level (Delta r)",
    "Cumulative interest-rate level change implied by Delta-r IRFs",
    "REER log change (REER_DLOG)",
    "Cumulative REER log-level effect implied by REER_DLOG IRFs",
    "Equity return"
  ),
  CumulationMethod = c(
    "NONE",
    "DRAW_LEVEL_CUMSUM_BEFORE_POSTERIOR_QUANTILES",
    "NONE",
    "DRAW_LEVEL_CUMSUM_BEFORE_POSTERIOR_QUANTILES",
    "NONE"
  ),
  PublicationWarning = c(
    "Do not label raw r as a rate-level deviation in the Delta-r specification.",
    "This is the appropriate rate-level cumulative response for the Delta-r model.",
    "Positive de follows the configured REER convention.",
    "This is a cumulative log-level REER effect, not a percentage unless explicitly converted.",
    "Do not call this a cumulative equity-price effect without verifying the upstream return definition."
  ),
  stringsAsFactors = FALSE
)
write.csv(scale_manifest, file.path(OUT, "08_irf_scale_manifest.csv"), row.names = FALSE)

# =============================================================================
# 5. Formal full-posterior IRF gate
# =============================================================================

stable_ok <- all(anchor_stability$StableShare >= MIN_STABLE_SHARE)
g0_ok <- all(anchor_stability$G0OKShare >= 0.99)
finite_ok <- all(anchor_stability$FiniteRhoShare >= 0.99)
raw_draws_ok <- all(irf_summary$StablePosteriorDraws >= 200)
cum_rate_draws_ok <- all(cum_rate$StablePosteriorDraws >= 200)
cum_reer_draws_ok <- all(cum_reer$StablePosteriorDraws >= 200)
lineage_ok <- nrow(lineage) == expected_n && all(lineage$SourcePanel == EXPECTED_PANEL_BASENAME)

ready <- stable_ok && g0_ok && finite_ok && raw_draws_ok &&
  cum_rate_draws_ok && cum_reer_draws_ok && lineage_ok

status <- if (ready) "READY_FOR_RATE_DIFF_IRF_AUDIT" else "FAIL"

reasons <- c(
  if (!stable_ok) sprintf("at least one anchor stable share < %.3f", MIN_STABLE_SHARE) else NULL,
  if (!g0_ok) "at least one anchor G0-OK share < 0.99" else NULL,
  if (!finite_ok) "at least one anchor finite-rho share < 0.99" else NULL,
  if (!raw_draws_ok) "at least one raw IRF cell has fewer than 200 valid stable draws" else NULL,
  if (!cum_rate_draws_ok) "at least one cumulative-rate cell has fewer than 200 valid stable draws" else NULL,
  if (!cum_reer_draws_ok) "at least one cumulative-REER cell has fewer than 200 valid stable draws" else NULL,
  if (!lineage_ok) "posterior lineage is not the accepted 56-file Delta-r grid" else NULL
)
if (!length(reasons)) reasons <- "all formal Delta-r GPR-IRF gates passed"

gate <- data.frame(
  Status = status,
  MainNetwork = MAIN_NETWORK,
  PosteriorPanel = EXPECTED_PANEL_BASENAME,
  PosteriorFiles = nrow(lineage),
  PosteriorDrawsTotal = DRAWS_TOTAL,
  Chains = NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  EventAnchors = nrow(anchors),
  Horizon = HORIZON,
  GPRShockPct = GPR_SHOCK_PCT,
  GPRShockLogPoints = shock_size,
  RateMode = "difference",
  RateLevelCumulation = "DRAW_LEVEL_CUMSUM_BEFORE_POSTERIOR_QUANTILES",
  REERCumulation = "DRAW_LEVEL_CUMSUM_BEFORE_POSTERIOR_QUANTILES",
  EquityPriceCumulationPerformed = FALSE,
  MinStableShareRequired = MIN_STABLE_SHARE,
  WorstAnchorStableShare = min(anchor_stability$StableShare),
  WorstAnchorG0OKShare = min(anchor_stability$G0OKShare),
  WorstAnchorFiniteRhoShare = min(anchor_stability$FiniteRhoShare),
  MinValidRawIRFDraws = min(irf_summary$StablePosteriorDraws),
  MinValidCumulativeRateDraws = min(cum_rate$StablePosteriorDraws),
  MinValidCumulativeREERDraws = min(cum_reer$StablePosteriorDraws),
  Reason = paste(reasons, collapse = "; "),
  stringsAsFactors = FALSE
)
write.csv(gate, file.path(OUT, "00_formal_rate_diff_irf_gate.csv"), row.names = FALSE)

readme <- c(
  sprintf("FORMAL DELTA-r TIME-SPECIFIC GPR TVP-GVAR IRF: %s", status),
  "============================================================",
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf("Posterior panel: %s", EXPECTED_PANEL_BASENAME),
  sprintf("Accepted posterior files: %d", nrow(lineage)),
  sprintf("Posterior draws pooled across chains: %d", DRAWS_TOTAL),
  sprintf("Event anchors: %d", nrow(anchors)),
  sprintf("Horizon: 0-%d quarters", HORIZON),
  sprintf("Positive GPR shock: %.2f%% = %.8f log points", GPR_SHOCK_PCT, shock_size),
  sprintf("Worst full-posterior event-anchor stable share: %.6f", gate$WorstAnchorStableShare),
  sprintf("Worst event-anchor G0-OK share: %.6f", gate$WorstAnchorG0OKShare),
  "",
  "Variable interpretation:",
  "- r   = Delta interest-rate level. Raw r IRF is NOT a rate-level deviation.",
  "- r_cumulative_level_change is computed draw by draw and is the rate-level cumulative response.",
  "- de  = REER_DLOG.",
  "- de_cumulative_log_reer is computed draw by draw and is the cumulative REER log-level effect.",
  "- deq = equity return; no equity-price cumulation is imposed.",
  "",
  "Identification:",
  "- GPR is external/global in the local country equations.",
  "- Positive GPR innovation is injected through gpr_0 and gpr_L1.",
  "- No domestic Cholesky ordering is imposed for the GPR shock.",
  "- Brent receives no simultaneous innovation in this experiment.",
  "- TVP coefficients are frozen at the selected event anchor over the response horizon.",
  "",
  "Outputs:",
  "- irf_posterior_summary.csv [raw Delta-r / REER_DLOG / EQ_RETURN IRFs; input to R/55 and R/56]",
  "- 01_irf_stability_by_anchor.csv",
  "- 02_local_stability_by_anchor_country.csv",
  "- 03_cumulative_reer_posterior_summary.csv",
  "- 04_core_event_impact_raw_irf.csv",
  "- 05_cumulative_rate_level_posterior_summary.csv",
  "- 06_core_event_cumulative_rate_level.csv",
  "- 07_core_event_cumulative_reer.csv",
  "- 08_irf_scale_manifest.csv",
  "- 00_posterior_lineage.csv",
  "- 00_formal_rate_diff_irf_gate.csv"
)
writeLines(readme, file.path(OUT, "README_formal_rate_diff_irf.txt"))
cat(paste(readme, collapse = "\n"), "\n")

if (!ready) stop("Formal Delta-r GPR TVP-GVAR IRF gate failed.")
