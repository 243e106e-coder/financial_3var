#!/usr/bin/env Rscript

# =============================================================================
# 68_formal_10l_staticblock_gpr_tvp_irf.R
#
# Formal GPR-shock IRFs for the accepted 10l 14-country financial TVP-GVAR.
#
# ACCEPTED 10L MEAN-SPECIFICATION
# -------------------------------
# Countries: AU BR CA CH CN EA UK JP KR NO SG TR US ZA
# Variables:
#   r   = Delta short-term interest-rate level
#   de  = REER log change
#   deq = equity return
# p = 1, q = 1
# TVS = TRUE
# SV = TRUE
# Static contemporaneous foreign-star terms:
#   r_star_0, de_star_0, deq_star_0
# Dynamic/TVP terms:
#   r_L1, de_L1, deq_L1,
#   r_star_L1, de_star_L1, deq_star_L1,
#   gpr_0, gpr_L1, brent_0, brent_L1
# Network: trade_2000_2012
#
# IDENTIFICATION
# --------------
# GPR is an observed global/exogenous driver in every local country equation.
# A one-time positive innovation to log GPR at t enters through gpr_0 at h=0
# and, because the local model contains gpr_L1, through gpr_L1 at h=1.
# No Cholesky ordering among r/de/deq is imposed for this GPR experiment.
#
# GLOBAL DYNAMICS
# ---------------
# For every posterior draw and event anchor:
#   G0(t) x_t = G1(t) x_{t-1} + c0(t) gpr_t + c1(t) gpr_{t-1} + ...
#   F(t)      = solve(G0(t)) %*% G1(t)
# Coefficients are frozen at the selected event anchor over the IRF horizon.
# Only posterior draws satisfying both
#   rcond(G0) >= MIN_G0_RCOND
#   spectral_radius(F) < 1
# enter reported IRF summaries.
#
# TRANSFORMATIONS
# ---------------
# r is Delta rate level, so the rate-level response is obtained by exact
# draw-level cumulative summation before posterior quantiles.
# de is REER log change, so cumulative REER log-level effects are also summed
# draw by draw. A nonlinear REER percent-level effect is then computed
# draw by draw as 100 * (exp(cumulative_log_effect) - 1).
# deq is NOT cumulated here because the upstream equity-return unit is not
# redefined by this script.
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
FORMAL_GATE_DIR <- get_env_chr("FIN3_CONVERGENCE_DIR", file.path(RESULTS_DIR, "formal_tvp"))
TENL_AUDIT_DIR <- get_env_chr("FIN3_10L_AUDIT_DIR", file.path(RESULTS_DIR, "10l_staticblock_sv_audit"))
OUT <- get_env_chr("FIN3_10M_OUT", file.path(RESULTS_DIR, "formal_irf_10m"))

EXPECTED_PANEL <- get_env_chr(
  "FIN3_FORMAL_PANEL_BASENAME",
  "panel_domestic_fin3_rate_diff.csv"
)
EXPECTED_NETWORK <- get_env_chr("FIN3_EXPECTED_NETWORK", "trade_2000_2012")
EXPECTED_RESTRICTION <- "STATIC_BLOCK_ALL_CONTEMPORANEOUS_FOREIGN_STAR"
EXPECTED_STATIC_TERMS <- c("r_star_0", "de_star_0", "deq_star_0")

NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_BURN <- as.integer(get_env_num("FIN3_EXPECTED_BURN", 12500))
EXPECTED_KEEP <- as.integer(get_env_num("FIN3_EXPECTED_KEEP", 50000))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 2000))
EXPECTED_SV_INNER <- as.integer(get_env_num("FIN3_EXPECTED_SV_INNER_BURNIN", 4))

HORIZON <- as.integer(get_env_num("FIN3_IRF_HORIZON", 12))
GPR_SHOCK_PCT <- get_env_num("FIN3_GPR_SHOCK_PCT", 10)
MIN_G0_RCOND <- get_env_num("FIN3_IRF_MIN_G0_RCOND", 1e-10)
MIN_STABLE_SHARE <- get_env_num("FIN3_IRF_MIN_STABLE_SHARE", 0.90)
MIN_VALID_SHARE <- get_env_num("FIN3_IRF_MIN_VALID_SHARE", 0.90)
MIN_G0_OK_SHARE <- get_env_num("FIN3_IRF_MIN_G0_OK_SHARE", 0.99)
MIN_FINITE_RHO_SHARE <- get_env_num("FIN3_IRF_MIN_FINITE_RHO_SHARE", 0.99)
IRF_CORES <- as.integer(get_env_num("FIN3_IRF_CORES", 2))

if (NCHAINS != 4L) stopf("10m requires exactly four chains per country.")
if (EXPECTED_STORED < 1000L) stopf("10m expected stored draws must be >= 1000 per chain.")
if (HORIZON < 1L) stopf("FIN3_IRF_HORIZON must be >= 1.")
if (!is.finite(GPR_SHOCK_PCT) || GPR_SHOCK_PCT <= 0) stopf("GPR shock percent must be positive.")
if (!(MIN_STABLE_SHARE > 0 && MIN_STABLE_SHARE <= 1)) stopf("Invalid MIN_STABLE_SHARE.")
if (!(MIN_VALID_SHARE > 0 && MIN_VALID_SHARE <= 1)) stopf("Invalid MIN_VALID_SHARE.")
if (!(MIN_G0_OK_SHARE > 0 && MIN_G0_OK_SHARE <= 1)) stopf("Invalid MIN_G0_OK_SHARE.")
if (!(MIN_FINITE_RHO_SHARE > 0 && MIN_FINITE_RHO_SHARE <= 1)) stopf("Invalid MIN_FINITE_RHO_SHARE.")
if (IRF_CORES < 1L) IRF_CORES <- 1L

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. Hard prerequisite: exact accepted 10l gate
# =============================================================================

formal_gate_file <- file.path(FORMAL_GATE_DIR, "00_formal_mcmc_gate.csv")
tenl_gate_file <- file.path(TENL_AUDIT_DIR, "07_final_10l_gate.csv")

if (!file.exists(formal_gate_file)) {
  stopf("Missing accepted 10l formal MCMC gate: %s", formal_gate_file)
}
if (!file.exists(tenl_gate_file)) {
  stopf("Missing accepted final 10l gate: %s", tenl_gate_file)
}

fg <- read.csv(formal_gate_file, stringsAsFactors = FALSE, check.names = FALSE)
lg <- read.csv(tenl_gate_file, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(fg) != 1L || !"Status" %in% names(fg)) stopf("Malformed formal MCMC gate.")
if (!identical(trimws(as.character(fg$Status[1])), "READY_FOR_FORMAL_IRF")) {
  stopf("Formal MCMC gate is not READY_FOR_FORMAL_IRF.")
}
if (nrow(lg) != 1L || !"FinalStatus" %in% names(lg)) stopf("Malformed final 10l gate.")
if (!identical(trimws(as.character(lg$FinalStatus[1])), "READY_FOR_FORMAL_IRF")) {
  stopf("Final 10l gate is not READY_FOR_FORMAL_IRF.")
}

required_final_cols <- c(
  "All56ChainDegeneracyAuditsPass",
  "All56StaticBlockContractsPass",
  "StaticCoefficientHardFailCount",
  "All126StaticCoefficientDiagnosticsPassHardGate"
)
if (!all(required_final_cols %in% names(lg))) {
  stopf("Final 10l gate is missing required integrity fields.")
}
if (!isTRUE(as.logical(lg$All56ChainDegeneracyAuditsPass[1]))) {
  stopf("10l chain-degeneracy gate is not PASS.")
}
if (!isTRUE(as.logical(lg$All56StaticBlockContractsPass[1]))) {
  stopf("10l static-block contract gate is not PASS.")
}
if (as.integer(lg$StaticCoefficientHardFailCount[1]) != 0L) {
  stopf("10l static coefficient hard-fail count is nonzero.")
}
if (!isTRUE(as.logical(lg$All126StaticCoefficientDiagnosticsPassHardGate[1]))) {
  stopf("10l static coefficient formal diagnostics are not all PASS.")
}

# =============================================================================
# 1. Discover, validate, and compact the exact 14 x 4 posterior lineage
# =============================================================================

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[1-4][.]rds$",
  recursive = TRUE,
  full.names = TRUE
)
expected_n <- length(COUNTRIES) * NCHAINS
if (length(files) != expected_n) {
  stopf("Expected exactly %d 10l posterior RDS files; found %d.", expected_n, length(files))
}

parts <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)
lineage_rows <- vector("list", expected_n)
static_contract_rows <- list()
li <- 0L
si <- 0L

reference_anchors <- NULL
reference_key_terms <- NULL

for (cc in COUNTRIES) {
  parts[[cc]] <- vector("list", NCHAINS)

  for (ch in seq_len(NCHAINS)) {
    candidates <- files[
      grepl(sprintf("formal_tvp_%s_chain%d[.]rds$", cc, ch), files)
    ]
    if (length(candidates) != 1L) {
      stopf("Expected exactly one 10l posterior for %s chain %d; found %d.", cc, ch, length(candidates))
    }

    z <- readRDS(candidates)

    if (!identical(as.character(z$meta$country), cc)) stopf("Country metadata mismatch: %s chain %d", cc, ch)
    if (as.integer(z$meta$chain) != ch) stopf("Chain metadata mismatch: %s chain %d", cc, ch)
    if (as.integer(z$meta$burn) != EXPECTED_BURN) stopf("Burn mismatch: %s chain %d", cc, ch)
    if (as.integer(z$meta$keep_internal) != EXPECTED_KEEP) stopf("Keep mismatch: %s chain %d", cc, ch)
    if (as.integer(z$meta$stored_draws) != EXPECTED_STORED) stopf("Stored draws mismatch: %s chain %d", cc, ch)
    if (!isTRUE(z$meta$sv_on)) stopf("SV must be TRUE: %s chain %d", cc, ch)
    if (!isTRUE(z$meta$tvs_on)) stopf("TVS must be TRUE: %s chain %d", cc, ch)
    if (abs(as.numeric(z$meta$prior_B1) - 2) > 1e-12) stopf("B1 mismatch: %s chain %d", cc, ch)
    if (abs(as.numeric(z$meta$prior_B2) - 1) > 1e-12) stopf("B2 mismatch: %s chain %d", cc, ch)
    if (abs(as.numeric(z$meta$kappa0) - 1e-4) > 1e-12) stopf("kappa0 mismatch: %s chain %d", cc, ch)
    if (!identical(basename(as.character(z$meta$source_panel)), EXPECTED_PANEL)) {
      stopf("Panel mismatch: %s chain %d", cc, ch)
    }
    if (!identical(as.character(z$meta$main_network), EXPECTED_NETWORK)) {
      stopf("Network mismatch: %s chain %d", cc, ch)
    }
    if (!identical(as.character(z$meta$restriction_id), EXPECTED_RESTRICTION)) {
      stopf("Static-block restriction mismatch: %s chain %d", cc, ch)
    }
    if (!identical(as.character(z$meta$static_terms), EXPECTED_STATIC_TERMS)) {
      stopf("Static-term metadata mismatch: %s chain %d", cc, ch)
    }
    if (as.integer(z$meta$sv_inner_burnin) != EXPECTED_SV_INNER) {
      stopf("SV inner-burnin mismatch: %s chain %d", cc, ch)
    }

    if (is.null(z$event_coef) || is.null(z$tv_probability) || is.null(z$anchors)) {
      stopf("Posterior lacks event_coef/tv_probability/anchors: %s chain %d", cc, ch)
    }
    if (dim(z$event_coef)[1] != EXPECTED_STORED) {
      stopf("event_coef draw dimension mismatch: %s chain %d", cc, ch)
    }

    if (is.null(reference_anchors)) {
      reference_anchors <- z$anchors
      reference_key_terms <- as.character(z$meta$key_terms)
    } else {
      this_anchor_key <- paste(z$anchors$EventID, z$anchors$AnchorType, z$anchors$AnchorQuarter, sep = "||")
      ref_anchor_key <- paste(reference_anchors$EventID, reference_anchors$AnchorType, reference_anchors$AnchorQuarter, sep = "||")
      if (!identical(this_anchor_key, ref_anchor_key)) {
        stopf("Event-anchor ordering mismatch: %s chain %d", cc, ch)
      }
      if (!identical(as.character(z$meta$key_terms), reference_key_terms)) {
        stopf("Key-term ordering mismatch: %s chain %d", cc, ch)
      }
    }

    # Re-audit the actual posterior object: static contemporaneous foreign-star
    # terms must have zero TV probability and identical draw values at every anchor.
    for (eq in VARS) {
      for (nm in EXPECTED_STATIC_TERMS) {
        p <- as.numeric(z$tv_probability[, nm, eq])
        max_tv <- max(abs(p))
        base <- as.numeric(z$event_coef[, 1, eq, nm])
        max_anchor_diff <- 0
        for (a in seq_len(dim(z$event_coef)[2])) {
          max_anchor_diff <- max(
            max_anchor_diff,
            max(abs(as.numeric(z$event_coef[, a, eq, nm]) - base))
          )
        }
        posterior_sd <- stats::sd(base)
        ok <- (
          is.finite(max_tv) && max_tv <= 1e-14 &&
          is.finite(max_anchor_diff) && max_anchor_diff <= 1e-8 &&
          all(is.finite(base)) && is.finite(posterior_sd) && posterior_sd > 1e-12
        )
        si <- si + 1L
        static_contract_rows[[si]] <- data.frame(
          Country = cc,
          Chain = ch,
          Equation = eq,
          StaticTerm = nm,
          MaxAbsTVProbability = max_tv,
          MaxAnchorDifference = max_anchor_diff,
          PosteriorSD = posterior_sd,
          Status = ifelse(ok, "PASS", "FAIL"),
          stringsAsFactors = FALSE
        )
        if (!ok) stopf("10m posterior static-block re-audit failed: %s chain %d %s %s", cc, ch, eq, nm)
      }
    }

    li <- li + 1L
    lineage_rows[[li]] <- data.frame(
      Country = cc,
      Chain = ch,
      Seed = as.integer(z$meta$seed),
      Burn = as.integer(z$meta$burn),
      KeepInternal = as.integer(z$meta$keep_internal),
      StoredDraws = as.integer(z$meta$stored_draws),
      SourcePanel = basename(as.character(z$meta$source_panel)),
      Network = as.character(z$meta$main_network),
      Restriction = as.character(z$meta$restriction_id),
      StaticTerms = paste(as.character(z$meta$static_terms), collapse = ";"),
      SV = isTRUE(z$meta$sv_on),
      TVS = isTRUE(z$meta$tvs_on),
      SVInnerBurnin = as.integer(z$meta$sv_inner_burnin),
      File = candidates,
      stringsAsFactors = FALSE
    )

    # Keep only the objects required for mean-dynamic IRFs.
    parts[[cc]][[ch]] <- list(
      event_coef = z$event_coef,
      key_terms = as.character(z$meta$key_terms)
    )
    rm(z)
  }
}

lineage <- do.call(rbind, lineage_rows)
static_contract <- do.call(rbind, static_contract_rows)
write.csv(lineage, file.path(OUT, "00_posterior_lineage.csv"), row.names = FALSE)
write.csv(static_contract, file.path(OUT, "00a_static_block_posterior_reaudit.csv"), row.names = FALSE)

key <- paste(lineage$Country, lineage$Chain, sep = "||")
expected_key <- paste(
  rep(COUNTRIES, each = NCHAINS),
  rep(seq_len(NCHAINS), times = length(COUNTRIES)),
  sep = "||"
)
if (anyDuplicated(key) || !setequal(key, expected_key)) {
  stopf("10m posterior lineage is not an exact 14 x 4 grid.")
}
if (nrow(static_contract) != 56L * 9L || any(static_contract$Status != "PASS")) {
  stopf("10m posterior static-block re-audit is not complete PASS.")
}

anchors <- reference_anchors
KEY_TERMS <- reference_key_terms

required_key_terms <- c(
  paste0(VARS, "_L1"),
  paste0(VARS, "_star_0"),
  paste0(VARS, "_star_L1"),
  "gpr_0", "gpr_L1",
  "brent_0", "brent_L1"
)
if (!all(required_key_terms %in% KEY_TERMS)) {
  stopf(
    "10l posterior is missing required mean-dynamic terms: %s",
    paste(setdiff(required_key_terms, KEY_TERMS), collapse = ", ")
  )
}

N <- length(COUNTRIES)
K <- length(VARS)
NK <- N * K
DRAWS_TOTAL <- NCHAINS * EXPECTED_STORED

r_idx <- match("r", VARS)
de_idx <- match("de", VARS)
deq_idx <- match("deq", VARS)
if (any(is.na(c(r_idx, de_idx, deq_idx)))) stopf("Expected r, de and deq in VARS.")

if (!identical(MAIN_NETWORK, EXPECTED_NETWORK)) {
  stopf("Runtime MAIN_NETWORK=%s does not equal accepted 10l network=%s.", MAIN_NETWORK, EXPECTED_NETWORK)
}
W <- read_weight_matrix(network_path(MAIN_NETWORK))

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

idx_dom <- match(paste0(VARS, "_L1"), KEY_TERMS)
idx_star0 <- match(paste0(VARS, "_star_0"), KEY_TERMS)
idx_star1 <- match(paste0(VARS, "_star_L1"), KEY_TERMS)
idx_gpr0 <- match("gpr_0", KEY_TERMS)
idx_gpr1 <- match("gpr_L1", KEY_TERMS)

if (any(is.na(c(idx_dom, idx_star0, idx_star1, idx_gpr0, idx_gpr1)))) {
  stopf("Failed to map accepted 10l coefficient blocks.")
}

if (!grepl("^LN_", toupper(GPR_COLUMN))) {
  stopf("10m percent normalization requires logged GPR; GPR_COLUMN=%s", GPR_COLUMN)
}
shock_size <- log1p(GPR_SHOCK_PCT / 100)

# =============================================================================
# 2. Posterior summary helpers
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
    EventID = as.character(meta$EventID),
    EventSet = as.character(meta$EventSet),
    EventLabel = as.character(meta$EventLabel),
    ShockFamily = as.character(meta$ShockFamily),
    EventQuarter = as.character(meta$EventQuarter),
    AnchorType = as.character(meta$AnchorType),
    AnchorQuarter = as.character(meta$AnchorQuarter),
    Country = country,
    ResponseVariable = variable,
    ResponseScale = scale_label,
    Horizon = horizon,
    p05 = unname(qs["p05"]),
    p16 = unname(qs["p16"]),
    median = unname(qs["median"]),
    p84 = unname(qs["p84"]),
    p95 = unname(qs["p95"]),
    mean = unname(qs["mean"]),
    prob_positive = unname(qs["prob_positive"]),
    ValidPosteriorDraws = as.integer(unname(qs["n"])),
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
# 3. One event-anchor IRF
# =============================================================================

process_anchor <- function(a) {
  meta <- anchors[a, , drop = FALSE]

  irf <- array(
    NA_real_,
    dim = c(DRAWS_TOTAL, HORIZON + 1L, NK),
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
        arr <- parts[[cc]][[ch]]$event_coef
        coef <- arr[dd, a, , , drop = FALSE]
        coef <- matrix(
          coef,
          nrow = K,
          ncol = length(KEY_TERMS),
          dimnames = list(VARS, KEY_TERMS)
        )

        Ai <- coef[, idx_dom, drop = FALSE]
        B0i <- coef[, idx_star0, drop = FALSE]
        B1i <- coef[, idx_star1, drop = FALSE]

        local_rho[global_draw, i] <- tryCatch(
          max(Mod(eigen(Ai, only.values = TRUE)$values)),
          error = function(e) NA_real_
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

      G0inv <- tryCatch(solve(G0), error = function(e) NULL)
      if (is.null(G0inv) || any(!is.finite(G0inv))) next

      Fmat <- G0inv %*% G1
      rho <- tryCatch(
        max(Mod(eigen(Fmat, only.values = TRUE)$values)),
        error = function(e) NA_real_
      )
      rho_vec[global_draw] <- rho
      if (!is.finite(rho) || rho >= 1) next

      # h=0: one-time log-GPR innovation enters contemporaneously.
      xh <- as.numeric(G0inv %*% (c0 * shock_size))
      irf[global_draw, 1L, ] <- xh

      # h=1: endogenous propagation plus the mechanical gpr_L1 loading.
      if (HORIZON >= 1L) {
        xh <- as.numeric(Fmat %*% xh + G0inv %*% (c1 * shock_size))
        irf[global_draw, 2L, ] <- xh
      }

      # h>=2: no additional GPR innovation.
      if (HORIZON >= 2L) {
        for (hh in 2:HORIZON) {
          xh <- as.numeric(Fmat %*% xh)
          irf[global_draw, hh + 1L, ] <- xh
        }
      }
    }
  }

  finite_rho <- is.finite(rho_vec)
  stable <- finite_rho & rho_vec < 1
  g0_ok <- is.finite(rcond_vec) & rcond_vec >= MIN_G0_RCOND
  valid <- stable & g0_ok

  scale_map <- c(
    r = "DELTA_RATE_LEVEL",
    de = "REER_DLOG",
    deq = "EQ_RETURN"
  )

  raw_rows <- vector("list", NK * (HORIZON + 1L))
  pos <- 0L
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
  raw <- do.call(rbind, raw_rows)

  # Exact draw-level cumulative interest-rate level response.
  rate_rows <- vector("list", N * (HORIZON + 1L))
  rp <- 0L
  for (i in seq_len(N)) {
    gv <- (i - 1L) * K + r_idx
    z <- matrix(irf[valid, , gv, drop = FALSE], nrow = sum(valid), ncol = HORIZON + 1L)
    zcum <- cumulate_draw_paths(z)
    for (hh in 0:HORIZON) {
      rp <- rp + 1L
      qs <- qsum(if (nrow(zcum)) zcum[, hh + 1L] else numeric())
      rate_rows[[rp]] <- summary_row(
        meta, COUNTRIES[i], "r_cumulative_level_change", hh, qs,
        "CUMULATIVE_RATE_LEVEL_CHANGE_FROM_DELTA_R"
      )
    }
  }
  cum_rate <- do.call(rbind, rate_rows)

  # Exact draw-level cumulative REER log-level response and nonlinear percent response.
  reer_log_rows <- vector("list", N * (HORIZON + 1L))
  reer_pct_rows <- vector("list", N * (HORIZON + 1L))
  ep <- 0L
  for (i in seq_len(N)) {
    gv <- (i - 1L) * K + de_idx
    z <- matrix(irf[valid, , gv, drop = FALSE], nrow = sum(valid), ncol = HORIZON + 1L)
    zcum <- cumulate_draw_paths(z)
    zpct <- if (nrow(zcum)) 100 * (exp(zcum) - 1) else matrix(NA_real_, 0, HORIZON + 1L)

    for (hh in 0:HORIZON) {
      ep <- ep + 1L
      qlog <- qsum(if (nrow(zcum)) zcum[, hh + 1L] else numeric())
      qpct <- qsum(if (nrow(zpct)) zpct[, hh + 1L] else numeric())
      reer_log_rows[[ep]] <- summary_row(
        meta, COUNTRIES[i], "de_cumulative_log_reer", hh, qlog,
        "CUMULATIVE_REER_LOG_LEVEL_EFFECT"
      )
      reer_pct_rows[[ep]] <- summary_row(
        meta, COUNTRIES[i], "de_cumulative_reer_percent", hh, qpct,
        "CUMULATIVE_REER_PERCENT_LEVEL_EFFECT"
      )
    }
  }
  cum_reer_log <- do.call(rbind, reer_log_rows)
  cum_reer_pct <- do.call(rbind, reer_pct_rows)

  rho_finite_values <- rho_vec[is.finite(rho_vec)]
  rc_finite_values <- rcond_vec[is.finite(rcond_vec)]

  stability <- data.frame(
    EventID = as.character(meta$EventID),
    EventSet = as.character(meta$EventSet),
    EventLabel = as.character(meta$EventLabel),
    ShockFamily = as.character(meta$ShockFamily),
    AnchorType = as.character(meta$AnchorType),
    AnchorQuarter = as.character(meta$AnchorQuarter),
    TotalDraws = DRAWS_TOTAL,
    FiniteRhoShare = mean(finite_rho),
    StableShare = mean(stable),
    G0OKShare = mean(g0_ok),
    ValidIRFShare = mean(valid),
    ValidIRFDraws = sum(valid),
    RhoMedian = if (length(rho_finite_values)) median(rho_finite_values) else NA_real_,
    RhoP95 = if (length(rho_finite_values)) unname(quantile(rho_finite_values, .95)) else NA_real_,
    RhoMax = if (length(rho_finite_values)) max(rho_finite_values) else NA_real_,
    MinG0Rcond = if (length(rc_finite_values)) min(rc_finite_values) else NA_real_,
    stringsAsFactors = FALSE
  )

  local_rows <- do.call(rbind, lapply(seq_len(N), function(i) {
    rr <- local_rho[, i]
    finite <- is.finite(rr)
    vals <- rr[finite]
    data.frame(
      EventID = as.character(meta$EventID),
      EventSet = as.character(meta$EventSet),
      AnchorType = as.character(meta$AnchorType),
      AnchorQuarter = as.character(meta$AnchorQuarter),
      Country = COUNTRIES[i],
      LocalFiniteShare = mean(finite),
      LocalStableShare = if (length(vals)) mean(vals < 1) else NA_real_,
      LocalRhoMedian = if (length(vals)) median(vals) else NA_real_,
      LocalRhoP95 = if (length(vals)) unname(quantile(vals, .95)) else NA_real_,
      LocalRhoMax = if (length(vals)) max(vals) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(
    raw = raw,
    cum_rate = cum_rate,
    cum_reer_log = cum_reer_log,
    cum_reer_pct = cum_reer_pct,
    stability = stability,
    local_stability = local_rows
  )
}

# =============================================================================
# 4. Run all accepted event anchors
# =============================================================================

anchor_ids <- seq_len(nrow(anchors))
if (.Platform$OS.type == "unix" && IRF_CORES > 1L && length(anchor_ids) > 1L) {
  ans <- parallel::mclapply(
    anchor_ids,
    process_anchor,
    mc.cores = min(IRF_CORES, length(anchor_ids)),
    mc.preschedule = TRUE
  )
} else {
  ans <- lapply(anchor_ids, process_anchor)
}

irf_summary <- do.call(rbind, lapply(ans, `[[`, "raw"))
cum_rate <- do.call(rbind, lapply(ans, `[[`, "cum_rate"))
cum_reer_log <- do.call(rbind, lapply(ans, `[[`, "cum_reer_log"))
cum_reer_pct <- do.call(rbind, lapply(ans, `[[`, "cum_reer_pct"))
anchor_stability <- do.call(rbind, lapply(ans, `[[`, "stability"))
local_stability <- do.call(rbind, lapply(ans, `[[`, "local_stability"))

write.csv(irf_summary, file.path(OUT, "01_irf_posterior_summary.csv"), row.names = FALSE)
write.csv(anchor_stability, file.path(OUT, "02_irf_stability_by_anchor.csv"), row.names = FALSE)
write.csv(local_stability, file.path(OUT, "03_local_stability_by_anchor_country.csv"), row.names = FALSE)
write.csv(cum_rate, file.path(OUT, "04_cumulative_rate_level_posterior_summary.csv"), row.names = FALSE)
write.csv(cum_reer_log, file.path(OUT, "05_cumulative_reer_log_posterior_summary.csv"), row.names = FALSE)
write.csv(cum_reer_pct, file.path(OUT, "06_cumulative_reer_percent_posterior_summary.csv"), row.names = FALSE)

# Publication-oriented compact tables at standard horizons.
keep_h <- intersect(c(0L, 1L, 4L, 8L, 12L), 0:HORIZON)
core_raw <- irf_summary[
  irf_summary$EventSet == "CORE" &
    irf_summary$AnchorType == "EVENT_QUARTER_t0" &
    irf_summary$Horizon %in% keep_h,
  , drop = FALSE
]
core_rate <- cum_rate[
  cum_rate$EventSet == "CORE" &
    cum_rate$AnchorType == "EVENT_QUARTER_t0" &
    cum_rate$Horizon %in% keep_h,
  , drop = FALSE
]
core_reer_pct <- cum_reer_pct[
  cum_reer_pct$EventSet == "CORE" &
    cum_reer_pct$AnchorType == "EVENT_QUARTER_t0" &
    cum_reer_pct$Horizon %in% keep_h,
  , drop = FALSE
]

write.csv(core_raw, file.path(OUT, "07_core_event_raw_irf.csv"), row.names = FALSE)
write.csv(core_rate, file.path(OUT, "08_core_event_cumulative_rate_level.csv"), row.names = FALSE)
write.csv(core_reer_pct, file.path(OUT, "09_core_event_cumulative_reer_percent.csv"), row.names = FALSE)

# Safe-haven input table: metrics only, no arbitrary classification threshold.
safe_haven_inputs <- core_reer_pct[
  core_reer_pct$Horizon %in% intersect(c(1L, 4L, 8L, 12L), keep_h),
  c(
    "EventID", "EventLabel", "ShockFamily", "AnchorQuarter", "Country", "Horizon",
    "p05", "median", "p95", "prob_positive", "ValidPosteriorDraws"
  ),
  drop = FALSE
]
names(safe_haven_inputs)[names(safe_haven_inputs) == "p05"] <- "REER_AppreciationPct_p05"
names(safe_haven_inputs)[names(safe_haven_inputs) == "median"] <- "REER_AppreciationPct_Median"
names(safe_haven_inputs)[names(safe_haven_inputs) == "p95"] <- "REER_AppreciationPct_p95"
names(safe_haven_inputs)[names(safe_haven_inputs) == "prob_positive"] <- "PosteriorProb_REER_Appreciation"
write.csv(safe_haven_inputs, file.path(OUT, "10_safe_haven_reer_inputs.csv"), row.names = FALSE)

scale_manifest <- data.frame(
  OutputVariable = c(
    "r",
    "r_cumulative_level_change",
    "de",
    "de_cumulative_log_reer",
    "de_cumulative_reer_percent",
    "deq"
  ),
  Meaning = c(
    "Quarterly change in the short-term interest-rate level",
    "Cumulative interest-rate level change implied by Delta-r IRFs",
    "REER log change",
    "Cumulative REER log-level effect implied by REER_DLOG IRFs",
    "Cumulative REER percent-level effect: 100*(exp(cumulative log effect)-1)",
    "Equity return"
  ),
  CumulationMethod = c(
    "NONE",
    "DRAW_LEVEL_CUMSUM_BEFORE_POSTERIOR_QUANTILES",
    "NONE",
    "DRAW_LEVEL_CUMSUM_BEFORE_POSTERIOR_QUANTILES",
    "DRAW_LEVEL_CUMSUM_THEN_NONLINEAR_PERCENT_CONVERSION_BEFORE_POSTERIOR_QUANTILES",
    "NONE"
  ),
  PublicationWarning = c(
    "Do not label raw r as a rate-level deviation.",
    "Use this series when discussing the level of the policy/short rate relative to baseline.",
    "Positive de follows the configured REER appreciation convention.",
    "This is a cumulative log-level REER effect.",
    "Positive values mean REER appreciation under the configured convention.",
    "Do not cumulate as an equity-price level effect unless the upstream return unit is separately verified."
  ),
  stringsAsFactors = FALSE
)
write.csv(scale_manifest, file.path(OUT, "11_irf_scale_manifest.csv"), row.names = FALSE)

identification_manifest <- data.frame(
  Item = c(
    "Shock",
    "ShockNormalization",
    "ContemporaneousEndogenousOrdering",
    "GPRImpactLoading",
    "GPRLagLoading",
    "BrentSimultaneousShock",
    "TVPCoefficientTreatmentOverHorizon",
    "ContemporaneousForeignStarTreatment",
    "PosteriorStabilityFilter"
  ),
  Value = c(
    "Positive one-time innovation to log GPR",
    sprintf("%.6f percent = %.12f log points", GPR_SHOCK_PCT, shock_size),
    "NONE_FOR_GPR_EXPERIMENT",
    "gpr_0",
    "gpr_L1 at h=1 only",
    "ZERO",
    "Frozen at each accepted event anchor",
    "Static in 10l: r_star_0;de_star_0;deq_star_0",
    sprintf("rcond(G0)>=%g and spectral_radius(F)<1", MIN_G0_RCOND)
  ),
  stringsAsFactors = FALSE
)
write.csv(identification_manifest, file.path(OUT, "12_identification_manifest.csv"), row.names = FALSE)

# =============================================================================
# 5. Formal 10m IRF integrity gate
# =============================================================================

finite_summary <- function(x) {
  all(is.finite(x$p05)) &&
    all(is.finite(x$median)) &&
    all(is.finite(x$p95)) &&
    all(x$p05 <= x$median) &&
    all(x$median <= x$p95)
}

stable_ok <- all(anchor_stability$StableShare >= MIN_STABLE_SHARE)
valid_ok <- all(anchor_stability$ValidIRFShare >= MIN_VALID_SHARE)
g0_ok <- all(anchor_stability$G0OKShare >= MIN_G0_OK_SHARE)
finite_rho_ok <- all(anchor_stability$FiniteRhoShare >= MIN_FINITE_RHO_SHARE)
raw_numeric_ok <- finite_summary(irf_summary)
rate_numeric_ok <- finite_summary(cum_rate)
reer_log_numeric_ok <- finite_summary(cum_reer_log)
reer_pct_numeric_ok <- finite_summary(cum_reer_pct)
lineage_ok <- (
  nrow(lineage) == 56L &&
  all(lineage$SourcePanel == EXPECTED_PANEL) &&
  all(lineage$Network == EXPECTED_NETWORK) &&
  all(lineage$Restriction == EXPECTED_RESTRICTION) &&
  all(lineage$StaticTerms == paste(EXPECTED_STATIC_TERMS, collapse = ";")) &&
  all(lineage$SV) &&
  all(lineage$TVS) &&
  all(lineage$SVInnerBurnin == EXPECTED_SV_INNER)
)
static_reaudit_ok <- (
  nrow(static_contract) == 504L &&
  all(static_contract$Status == "PASS")
)

ready <- all(c(
  stable_ok,
  valid_ok,
  g0_ok,
  finite_rho_ok,
  raw_numeric_ok,
  rate_numeric_ok,
  reer_log_numeric_ok,
  reer_pct_numeric_ok,
  lineage_ok,
  static_reaudit_ok
))

reasons <- c(
  if (!lineage_ok) "10l posterior lineage contract failed" else NULL,
  if (!static_reaudit_ok) "10l static-block posterior re-audit failed" else NULL,
  if (!stable_ok) sprintf("at least one anchor StableShare < %.3f", MIN_STABLE_SHARE) else NULL,
  if (!valid_ok) sprintf("at least one anchor ValidIRFShare < %.3f", MIN_VALID_SHARE) else NULL,
  if (!g0_ok) sprintf("at least one anchor G0OKShare < %.3f", MIN_G0_OK_SHARE) else NULL,
  if (!finite_rho_ok) sprintf("at least one anchor FiniteRhoShare < %.3f", MIN_FINITE_RHO_SHARE) else NULL,
  if (!raw_numeric_ok) "raw IRF posterior summaries failed numeric integrity" else NULL,
  if (!rate_numeric_ok) "cumulative rate posterior summaries failed numeric integrity" else NULL,
  if (!reer_log_numeric_ok) "cumulative REER log summaries failed numeric integrity" else NULL,
  if (!reer_pct_numeric_ok) "cumulative REER percent summaries failed numeric integrity" else NULL
)
if (!length(reasons)) reasons <- "all formal 10m GPR-IRF integrity gates passed"

status <- if (ready) "READY_FOR_FORMAL_10M_IRF_AUDIT" else "FAIL"

gate <- data.frame(
  Status = status,
  SourcePosteriorGate = as.character(lg$FinalStatus[1]),
  Countries = length(COUNTRIES),
  PosteriorFiles = nrow(lineage),
  ChainsPerCountry = NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  GlobalPosteriorDrawsPerAnchor = DRAWS_TOTAL,
  EventAnchors = nrow(anchors),
  Horizon = HORIZON,
  GPRShockPct = GPR_SHOCK_PCT,
  GPRShockLogPoints = shock_size,
  Network = EXPECTED_NETWORK,
  Panel = EXPECTED_PANEL,
  StaticRestriction = EXPECTED_RESTRICTION,
  StaticTerms = paste(EXPECTED_STATIC_TERMS, collapse = ";"),
  SV = TRUE,
  TVS = TRUE,
  SVInnerBurnin = EXPECTED_SV_INNER,
  MinStableShareRequired = MIN_STABLE_SHARE,
  MinValidShareRequired = MIN_VALID_SHARE,
  MinG0OKShareRequired = MIN_G0_OK_SHARE,
  MinFiniteRhoShareRequired = MIN_FINITE_RHO_SHARE,
  WorstAnchorStableShare = min(anchor_stability$StableShare),
  WorstAnchorValidIRFShare = min(anchor_stability$ValidIRFShare),
  WorstAnchorG0OKShare = min(anchor_stability$G0OKShare),
  WorstAnchorFiniteRhoShare = min(anchor_stability$FiniteRhoShare),
  MaxAnchorRhoP95 = max(anchor_stability$RhoP95, na.rm = TRUE),
  MinAnchorG0Rcond = min(anchor_stability$MinG0Rcond, na.rm = TRUE),
  StaticPosteriorReauditChecks = nrow(static_contract),
  StaticPosteriorReauditFailures = sum(static_contract$Status != "PASS"),
  RawNumericIntegrity = raw_numeric_ok,
  CumulativeRateNumericIntegrity = rate_numeric_ok,
  CumulativeREERLogNumericIntegrity = reer_log_numeric_ok,
  CumulativeREERPercentNumericIntegrity = reer_pct_numeric_ok,
  Reason = paste(reasons, collapse = "; "),
  stringsAsFactors = FALSE
)
write.csv(gate, file.path(OUT, "00_formal_10m_irf_gate.csv"), row.names = FALSE)

readme <- c(
  sprintf("FORMAL 10M GPR-SHOCK TVP-GVAR IRF: %s", status),
  "============================================================",
  sprintf("Accepted 10l posterior gate: %s", as.character(lg$FinalStatus[1])),
  sprintf("Countries: %d", length(COUNTRIES)),
  sprintf("Posterior files: %d", nrow(lineage)),
  sprintf("Chains per country: %d", NCHAINS),
  sprintf("Stored draws per chain: %d", EXPECTED_STORED),
  sprintf("Global posterior draws per anchor: %d", DRAWS_TOTAL),
  sprintf("Network: %s", EXPECTED_NETWORK),
  sprintf("Panel: %s", EXPECTED_PANEL),
  sprintf("Event anchors: %d", nrow(anchors)),
  sprintf("IRF horizon: 0-%d quarters", HORIZON),
  sprintf("Positive GPR shock: %.6f%% = %.12f log points", GPR_SHOCK_PCT, shock_size),
  sprintf("Worst anchor stable share: %.6f", gate$WorstAnchorStableShare),
  sprintf("Worst anchor valid-IRF share: %.6f", gate$WorstAnchorValidIRFShare),
  sprintf("Worst anchor G0-OK share: %.6f", gate$WorstAnchorG0OKShare),
  "",
  "Identification:",
  "- GPR is an observed global/exogenous driver.",
  "- No Cholesky ordering among r/de/deq is imposed for the GPR experiment.",
  "- gpr_0 loads at h=0 and gpr_L1 loads at h=1.",
  "- Brent receives no simultaneous innovation.",
  "- TVP coefficients are frozen at each event anchor over the response horizon.",
  "- r_star_0/de_star_0/deq_star_0 remain static exactly as accepted in 10l.",
  "",
  "Transformation:",
  "- raw r is Delta rate; use cumulative r output for rate-level discussion.",
  "- raw de is REER log change.",
  "- cumulative REER percent is computed draw by draw after log-level cumulation.",
  "- deq is not mechanically cumulated into an equity-price level.",
  "",
  "Primary outputs:",
  "- 01_irf_posterior_summary.csv",
  "- 02_irf_stability_by_anchor.csv",
  "- 04_cumulative_rate_level_posterior_summary.csv",
  "- 05_cumulative_reer_log_posterior_summary.csv",
  "- 06_cumulative_reer_percent_posterior_summary.csv",
  "- 07_core_event_raw_irf.csv",
  "- 08_core_event_cumulative_rate_level.csv",
  "- 09_core_event_cumulative_reer_percent.csv",
  "- 10_safe_haven_reer_inputs.csv",
  "- 11_irf_scale_manifest.csv",
  "- 12_identification_manifest.csv",
  "- 00_formal_10m_irf_gate.csv"
)
writeLines(readme, file.path(OUT, "README_formal_10m_irf.txt"))
cat(paste(readme, collapse = "\n"), "\n")

if (!ready) {
  stop("Formal 10m IRF integrity gate failed. Do not interpret or publish IRFs.")
}
