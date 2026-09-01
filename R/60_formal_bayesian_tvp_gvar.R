#!/usr/bin/env Rscript

# =============================================================================
# 60_formal_bayesian_tvp_gvar.R
#
# One-country / one-chain formal Bayesian TVP-GVAR estimator.
#
# MODEL
# -----
# For country i and each domestic response variable k in {r, de, deq}:
#
#   y_{i,k,t} = x_{i,t}' beta_{i,k,t} + epsilon_{i,k,t}
#
# where x contains:
#   constant
#   domestic lag-1:              r_L1, de_L1, deq_L1
#   contemporaneous foreign*:    r_star_0, de_star_0, deq_star_0
#   foreign* lag-1:              r_star_L1, de_star_L1, deq_star_L1
#   GPR current + lag-1:         gpr_0, gpr_L1
#   Brent current + lag-1:       brent_0, brent_L1
#
# p = 1, q = 1, fixed MAIN_NETWORK financial weights.
#
# FORMAL TVP PRIOR / SAMPLER
# --------------------------
# Uses the latent-threshold TVP sampler in threshtvp with:
#   * stochastic volatility,
#   * latent threshold time variation (TVS=TRUE),
#   * Normal-Gamma shrinkage on the time-invariant component,
#   * MCMC burn-in + post-burn iterations,
#   * multiple independent chains (handled by GitHub Actions matrix).
#
# The threshtvp package is pinned by the workflow to:
#   gregorkastner/threshtvp@c0133f262fcf0f935805168f8b7def02b43dc6d5
#
# IMPORTANT
# ---------
# 1) This is the formal Bayesian replacement for R/50 smoke-test draws.
# 2) It does NOT use the smoke-test PriorScale=.02 as the publication prior.
# 3) Predictors and each response are standardized before MCMC; event-date
#    coefficient draws are transformed back to ORIGINAL DATA UNITS.
# 4) The baseline variable transformations remain:
#       r   = RATE_LEVEL
#       de  = REER_DLOG
#       deq = EQ_RETURN
# 5) GPR is treated as an observed global/exogenous driver. The later GPR IRF
#    therefore does not require a Cholesky ordering of domestic residuals.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("threshtvp", quietly = TRUE)) {
  stopf("Package 'threshtvp' is required.")
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

COUNTRY <- toupper(get_env_chr("FIN3_COUNTRY", "AU"))
CHAIN_ID <- as.integer(get_env_num("FIN3_CHAIN_ID", 1))
BURN <- as.integer(get_env_num("FIN3_BURN", 5000))
KEEP <- as.integer(get_env_num("FIN3_KEEP", 7500))
STORED_PER_CHAIN <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 1500))
SEED_BASE <- as.integer(get_env_num("FIN3_SEED_BASE", 20260910))

PRIOR_B1 <- get_env_num("FIN3_TTVP_B1", 2)
PRIOR_B2 <- get_env_num("FIN3_TTVP_B2", 1)
KAPPA0 <- get_env_num("FIN3_TTVP_KAPPA0", 1e-7)

SV_ON <- tolower(get_env_chr("FIN3_SV_ON", "true")) %in% c("1", "true", "yes", "y")
TVS_ON <- tolower(get_env_chr("FIN3_TVS_ON", "true")) %in% c("1", "true", "yes", "y")

PANEL <- get_env_chr(
  "FIN3_FORMAL_PANEL",
  file.path(DERIVED_DIR, "panel_domestic_fin3.csv")
)
EVENT_FILE <- get_env_chr(
  "FIN3_EVENT_CALENDAR",
  file.path(RESULTS_DIR, "events", "00_event_calendar.csv")
)

if (!COUNTRY %in% COUNTRIES) {
  stopf("FIN3_COUNTRY must be one of: %s", paste(COUNTRIES, collapse = ", "))
}
if (CHAIN_ID < 1L) stopf("FIN3_CHAIN_ID must be >= 1.")
if (BURN < 100L) stopf("FIN3_BURN must be at least 100.")
if (KEEP < 500L) stopf("FIN3_KEEP must be at least 500.")
if (STORED_PER_CHAIN < 200L) stopf("FIN3_STORED_PER_CHAIN must be at least 200.")
if (STORED_PER_CHAIN >= KEEP) {
  stopf("FIN3_STORED_PER_CHAIN must be smaller than FIN3_KEEP.")
}
if (!(PRIOR_B1 > 0) || !(PRIOR_B2 > 0)) stopf("TTVP B1/B2 must be positive.")
if (!(KAPPA0 > 0)) stopf("TTVP kappa0 must be positive.")

# threshtvp stores thin * KEEP draws and its first stored iteration lies exactly
# at the burn boundary. Request one additional stored draw, then drop that
# boundary observation so the final artifact contains exactly STORED_PER_CHAIN
# strictly post-burn draws.
THIN_FRACTION <- (STORED_PER_CHAIN + 1) / KEEP
if (!(THIN_FRACTION > 0 && THIN_FRACTION <= 1)) {
  stopf("Invalid threshtvp thin fraction.")
}

SEED <- SEED_BASE + 10000L * CHAIN_ID + 101L * match(COUNTRY, COUNTRIES)
set.seed(SEED)

if (!file.exists(PANEL)) stopf("Missing formal panel: %s", PANEL)
if (!file.exists(EVENT_FILE)) stopf("Missing event calendar: %s", EVENT_FILE)

OUT <- file.path(
  RESULTS_DIR,
  "formal_tvp_parts",
  COUNTRY,
  paste0("chain_", CHAIN_ID)
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load panel and build the exact accepted p=1, q=1 design
# =============================================================================

d <- read.csv(PANEL, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("Quarter", "Country", VARS, "gpr", "brent")
if (!all(need %in% names(d))) {
  stopf("Formal panel is missing: %s", paste(setdiff(need, names(d)), collapse = ", "))
}

d$Country <- toupper(trimws(as.character(d$Country)))
d$Quarter <- toupper(trimws(as.character(d$Quarter)))

quarters <- unique(d$Quarter)
if (!all(COUNTRIES %in% unique(d$Country))) {
  stopf("Formal panel is missing sample economies.")
}
if (any(table(d$Country) != length(quarters))) {
  stopf("Formal panel is not balanced.")
}

N <- length(COUNTRIES)
K <- length(VARS)
Tn <- length(quarters)

Xglobal <- array(
  NA_real_,
  c(Tn, N, K),
  dimnames = list(quarters, COUNTRIES, VARS)
)

for (i in seq_along(COUNTRIES)) {
  z <- d[d$Country == COUNTRIES[i], , drop = FALSE]
  z <- z[match(quarters, z$Quarter), , drop = FALSE]
  Xglobal[, i, ] <- as.matrix(z[, VARS, drop = FALSE])
}

base <- d[d$Country == COUNTRIES[1], , drop = FALSE]
base <- base[match(quarters, base$Quarter), , drop = FALSE]
gpr <- num(base$gpr)
brent <- num(base$brent)

if (!all(is.finite(Xglobal)) || !all(is.finite(gpr)) || !all(is.finite(brent))) {
  stopf("Non-finite values in formal model inputs.")
}

W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])
country_i <- match(COUNTRY, COUNTRIES)

STAR <- array(
  NA_real_,
  c(Tn, N, K),
  dimnames = list(quarters, COUNTRIES, VARS)
)

for (i in seq_len(N)) {
  for (v in seq_len(K)) {
    STAR[, i, v] <- as.numeric(Xglobal[, , v] %*% W[i, ])
  }
}

lag1 <- function(x) c(NA_real_, x[-length(x)])

Yraw <- Xglobal[, country_i, , drop = FALSE][, 1, ]
Zraw <- STAR[, country_i, , drop = FALSE][, 1, ]
rows <- 2:nrow(Yraw)

D <- data.frame(const = rep(1, length(rows)), check.names = FALSE)
for (v in seq_len(K)) {
  D[[paste0(VARS[v], "_L1")]] <- lag1(Yraw[, v])[rows]
}
for (v in seq_len(K)) {
  D[[paste0(VARS[v], "_star_0")]] <- Zraw[rows, v]
}
for (v in seq_len(K)) {
  D[[paste0(VARS[v], "_star_L1")]] <- lag1(Zraw[, v])[rows]
}
D$gpr_0 <- gpr[rows]
D$gpr_L1 <- lag1(gpr)[rows]
D$brent_0 <- brent[rows]
D$brent_L1 <- lag1(brent)[rows]

Y <- Yraw[rows, , drop = FALSE]
model_quarters <- quarters[rows]

ok <- complete.cases(D) & complete.cases(Y)
D <- D[ok, , drop = FALSE]
Y <- Y[ok, , drop = FALSE]
model_quarters <- model_quarters[ok]

if (nrow(D) <= ncol(D) + 10L) {
  stopf("Too few observations in formal design for %s.", COUNTRY)
}

TERMS <- colnames(D)
EXPECTED_TERMS <- c(
  "const",
  paste0(VARS, "_L1"),
  paste0(VARS, "_star_0"),
  paste0(VARS, "_star_L1"),
  "gpr_0", "gpr_L1",
  "brent_0", "brent_L1"
)
if (!identical(TERMS, EXPECTED_TERMS)) {
  stopf(
    "Unexpected formal design ordering. Found: %s",
    paste(TERMS, collapse = ", ")
  )
}

# =============================================================================
# 2. Event-anchor table
# =============================================================================

ev <- read.csv(EVENT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
ev_need <- c(
  "EventID", "EventSet", "EventLabel", "ShockFamily",
  "EventQuarter", "AnchorQuarter_t0", "AnchorQuarter_t1", "EventWindowReady"
)
if (!all(ev_need %in% names(ev))) stopf("Malformed event calendar.")
if (any(!ev$EventWindowReady)) stopf("Event calendar contains non-ready event windows.")

anchors <- do.call(rbind, lapply(seq_len(nrow(ev)), function(i) {
  data.frame(
    EventID = ev$EventID[i],
    EventSet = ev$EventSet[i],
    EventLabel = ev$EventLabel[i],
    ShockFamily = ev$ShockFamily[i],
    EventQuarter = ev$EventQuarter[i],
    AnchorType = c("EVENT_QUARTER_t0", "NEXT_QUARTER_t1"),
    AnchorQuarter = c(ev$AnchorQuarter_t0[i], ev$AnchorQuarter_t1[i]),
    stringsAsFactors = FALSE
  )
}))

anchors$AnchorIndex <- match(anchors$AnchorQuarter, model_quarters)
if (any(is.na(anchors$AnchorIndex))) {
  bad <- paste(
    anchors$EventID[is.na(anchors$AnchorIndex)],
    anchors$AnchorQuarter[is.na(anchors$AnchorIndex)],
    sep = "@",
    collapse = ", "
  )
  stopf("Event anchors outside formal model sample: %s", bad)
}

KEY_TERMS <- c(
  paste0(VARS, "_L1"),
  paste0(VARS, "_star_0"),
  paste0(VARS, "_star_L1"),
  "gpr_0", "gpr_L1",
  "brent_0", "brent_L1"
)
key_index <- match(KEY_TERMS, TERMS)
if (any(is.na(key_index))) stopf("Missing key IRF term(s) in formal design.")

# =============================================================================
# 3. Standardize predictors once; standardize each response equation separately
# =============================================================================

Xm_raw <- as.matrix(D)
x_center <- rep(0, ncol(Xm_raw)); names(x_center) <- TERMS
x_scale <- rep(1, ncol(Xm_raw)); names(x_scale) <- TERMS
Xm <- Xm_raw

for (j in seq_len(ncol(Xm))) {
  if (TERMS[j] == "const") next
  x_center[j] <- mean(Xm[, j])
  x_scale[j] <- stats::sd(Xm[, j])
  if (!is.finite(x_scale[j]) || x_scale[j] < 1e-10) {
    stopf("Near-constant predictor in formal design: %s", TERMS[j])
  }
  Xm[, j] <- (Xm[, j] - x_center[j]) / x_scale[j]
}

# Constant remains exactly 1.
Xm[, "const"] <- 1

raw_kappa <- tryCatch(kappa(crossprod(Xm_raw)), error = function(e) Inf)
scaled_kappa <- tryCatch(kappa(crossprod(Xm)), error = function(e) Inf)

# =============================================================================
# 4. Run three formal equation samplers, one at a time
# =============================================================================

n_anchor <- nrow(anchors)
n_key <- length(KEY_TERMS)

event_coef <- array(
  NA_real_,
  c(STORED_PER_CHAIN, n_anchor, K, n_key),
  dimnames = list(
    Draw = seq_len(STORED_PER_CHAIN),
    Anchor = paste(anchors$EventID, anchors$AnchorType, sep = "__"),
    Equation = VARS,
    Term = KEY_TERMS
  )
)

event_logvar <- array(
  NA_real_,
  c(STORED_PER_CHAIN, n_anchor, K),
  dimnames = list(
    Draw = seq_len(STORED_PER_CHAIN),
    Anchor = paste(anchors$EventID, anchors$AnchorType, sep = "__"),
    Equation = VARS
  )
)

tv_probability <- array(
  NA_real_,
  c(length(model_quarters), length(TERMS), K),
  dimnames = list(
    Quarter = model_quarters,
    Term = TERMS,
    Equation = VARS
  )
)

monitor_blocks <- vector("list", K)
equation_summary <- vector("list", K)

for (eq in seq_len(K)) {
  eq_name <- VARS[eq]
  y_raw <- as.numeric(Y[, eq])

  y_center <- mean(y_raw)
  y_scale <- stats::sd(y_raw)
  if (!is.finite(y_scale) || y_scale < 1e-10) {
    stopf("Near-constant response: %s/%s", COUNTRY, eq_name)
  }

  y_std <- (y_raw - y_center) / y_scale

  msg(
    "Formal TVP MCMC: country=%s equation=%s chain=%d burn=%d keep=%d stored=%d",
    COUNTRY, eq_name, CHAIN_ID, BURN, KEEP, STORED_PER_CHAIN
  )

  fit <- threshtvp::estimate_tvp(
    Y = matrix(y_std, ncol = 1),
    X = Xm,
    save = KEEP,
    burn = BURN,
    priorbtheta = list(
      B_1 = PRIOR_B1,
      B_2 = PRIOR_B2,
      kappa0 = KAPPA0
    ),
    priorb0 = list(
      a_tau = 0.1,
      c_tau = 0.01,
      d_tau = 0.01
    ),
    priorsig = c(0.01, 0.01),
    priorphi = c(2, 2),
    priormu = c(0, 10^2),
    h0prior = "stationary",
    grid.length = 150,
    thrsh.pct = 0.1,
    thrsh.pct.high = 1,
    sv_on = SV_ON,
    TVS = TVS_ON,
    cons.mod = FALSE,
    p = 1,
    thin = THIN_FRACTION,
    CPU = 1,
    approx = FALSE,
    sim.kappa0 = FALSE
  )

  post <- fit$posterior

  if (is.null(post$A) || length(dim(post$A)) != 3L) {
    stopf("Unexpected threshtvp posterior A structure: %s/%s", COUNTRY, eq_name)
  }

  nstored_raw <- dim(post$A)[1]
  if (nstored_raw != STORED_PER_CHAIN + 1L) {
    stopf(
      "Expected %d raw stored draws but threshtvp returned %d for %s/%s.",
      STORED_PER_CHAIN + 1L, nstored_raw, COUNTRY, eq_name
    )
  }

  # Drop the exact burn-boundary stored draw.
  keep_draw <- 2:nstored_raw

  Astd <- post$A[keep_draw, , , drop = FALSE]
  if (dim(Astd)[1] != STORED_PER_CHAIN ||
      dim(Astd)[2] != length(model_quarters) ||
      dim(Astd)[3] != length(TERMS)) {
    stopf("Unexpected coefficient-array dimensions after boundary-drop.")
  }

  term_names <- dimnames(post$A)[[3]]
  if (is.null(term_names)) term_names <- TERMS
  if (!identical(term_names, TERMS)) {
    stopf(
      "threshtvp coefficient term order differs from formal design for %s/%s.",
      COUNTRY, eq_name
    )
  }

  # Transform selected slopes back to original units:
  # beta_orig_j = y_sd * beta_std_j / x_sd_j.
  for (a in seq_len(n_anchor)) {
    tt <- anchors$AnchorIndex[a]
    for (kterm in seq_along(KEY_TERMS)) {
      jj <- key_index[kterm]
      event_coef[, a, eq, kterm] <- (
        y_scale * Astd[, tt, jj] / x_scale[jj]
      )
    }
  }

  # Stochastic volatility is estimated on standardized y. If h is log variance,
  # original-unit log variance adds 2*log(y_sd).
  Hdraw <- post$H[keep_draw, , drop = FALSE]
  if (nrow(Hdraw) != STORED_PER_CHAIN ||
      ncol(Hdraw) != length(model_quarters)) {
    stopf("Unexpected stochastic-volatility draw dimensions.")
  }
  for (a in seq_len(n_anchor)) {
    tt <- anchors$AnchorIndex[a]
    event_logvar[, a, eq] <- Hdraw[, tt] + 2 * log(y_scale)
  }

  # Posterior probability that the latent threshold classifies a coefficient
  # as being in its high-variation regime at each date.
  Ddyn <- post$D_dyn[keep_draw, , , drop = FALSE]
  if (length(dim(Ddyn)) != 3L ||
      dim(Ddyn)[2] != length(TERMS) ||
      dim(Ddyn)[3] != length(model_quarters)) {
    stopf("Unexpected D_dyn dimensions.")
  }
  tv_probability[, , eq] <- t(apply(Ddyn, c(2, 3), mean, na.rm = TRUE))

  # ---------------------------------------------------------------------------
  # Monitor parameters for multi-chain R-hat / ESS.
  # Keep hyperparameters plus representative event-date coefficients.
  # ---------------------------------------------------------------------------

  mon <- data.frame(stringsAsFactors = FALSE)

  svp <- post$svparms[keep_draw, , drop = FALSE]
  if (ncol(svp) >= 3L) {
    mon[[paste0(eq_name, "__sv_mu")]] <- svp[, 1]
    mon[[paste0(eq_name, "__sv_phi")]] <- svp[, 2]
    mon[[paste0(eq_name, "__sv_sigma")]] <- svp[, 3]
  }

  omega <- post$omega[keep_draw, , drop = FALSE]
  thresholds <- post$thresholds[keep_draw, , drop = FALSE]
  v0 <- post$V0[keep_draw, , drop = FALSE]

  if (ncol(omega) != length(TERMS) ||
      ncol(thresholds) != length(TERMS) ||
      ncol(v0) != length(TERMS)) {
    stopf("Unexpected threshtvp hyperparameter dimensions.")
  }

  # Monitor all dynamic-state scales, but only a compact subset of thresholds
  # and Normal-Gamma initial-state scales to control artifact size.
  for (jj in seq_along(TERMS)) {
    nm <- TERMS[jj]
    mon[[paste0(eq_name, "__omega__", nm)]] <- omega[, jj]
  }

  compact_terms <- unique(c(
    paste0(eq_name, "_L1"),
    paste0(eq_name, "_star_0"),
    "gpr_0",
    "gpr_L1"
  ))
  compact_index <- match(compact_terms, TERMS)
  compact_index <- compact_index[!is.na(compact_index)]

  for (jj in compact_index) {
    nm <- TERMS[jj]
    mon[[paste0(eq_name, "__threshold__", nm)]] <- thresholds[, jj]
    mon[[paste0(eq_name, "__V0__", nm)]] <- v0[, jj]
  }

  # Representative coefficient draws at three core event dates.
  monitor_quarters <- intersect(
    c("2008Q3", "2020Q1", "2024Q2"),
    model_quarters
  )
  own_lag <- paste0(eq_name, "_L1")
  monitor_terms <- c(own_lag, "gpr_0")

  for (qq in monitor_quarters) {
    tt <- match(qq, model_quarters)
    for (nm in monitor_terms) {
      jj <- match(nm, TERMS)
      vals <- y_scale * Astd[, tt, jj] / x_scale[jj]
      mon[[paste0(eq_name, "__coef__", qq, "__", nm)]] <- vals
    }
  }

  monitor_blocks[[eq]] <- mon

  equation_summary[[eq]] <- data.frame(
    Country = COUNTRY,
    Chain = CHAIN_ID,
    Equation = eq_name,
    RawN = length(y_raw),
    Predictors = ncol(Xm),
    RawDesignKappa = raw_kappa,
    StandardizedDesignKappa = scaled_kappa,
    YMean = y_center,
    YSD = y_scale,
    StoredDraws = STORED_PER_CHAIN,
    Burn = BURN,
    KeepInternal = KEEP,
    Seed = SEED,
    SV = SV_ON,
    TVS = TVS_ON,
    stringsAsFactors = FALSE
  )

  rm(fit, post, Astd, Ddyn, Hdraw, omega, thresholds, v0, mon)
  invisible(gc())
}

monitor <- do.call(cbind, monitor_blocks)
if (nrow(monitor) != STORED_PER_CHAIN) {
  stopf("Final monitor matrix has unexpected draw count.")
}
if (any(!vapply(monitor, is.numeric, logical(1)))) {
  stopf("Non-numeric column in monitor matrix.")
}

equation_summary <- do.call(rbind, equation_summary)

# =============================================================================
# 5. Save compact formal posterior part
# =============================================================================

pkg_desc <- utils::packageDescription("threshtvp")
remote_sha <- if (!is.null(pkg_desc$RemoteSha)) as.character(pkg_desc$RemoteSha) else NA_character_
pkg_version <- as.character(utils::packageVersion("threshtvp"))

part <- list(
  meta = list(
    country = COUNTRY,
    chain = CHAIN_ID,
    seed = SEED,
    burn = BURN,
    keep_internal = KEEP,
    stored_draws = STORED_PER_CHAIN,
    threshtvp_thin_fraction = THIN_FRACTION,
    threshtvp_version = pkg_version,
    threshtvp_remote_sha = remote_sha,
    sv_on = SV_ON,
    tvs_on = TVS_ON,
    prior_B1 = PRIOR_B1,
    prior_B2 = PRIOR_B2,
    kappa0 = KAPPA0,
    main_network = MAIN_NETWORK,
    variables = VARS,
    model_quarters = model_quarters,
    terms = TERMS,
    key_terms = KEY_TERMS,
    source_panel = PANEL
  ),
  anchors = anchors,
  event_coef = event_coef,
  event_logvar = event_logvar,
  monitor = as.matrix(monitor),
  tv_probability = tv_probability,
  equation_summary = equation_summary
)

rds_path <- file.path(
  OUT,
  sprintf("formal_tvp_%s_chain%d.rds", COUNTRY, CHAIN_ID)
)
saveRDS(part, rds_path, compress = "xz")

write.csv(
  equation_summary,
  file.path(OUT, "01_equation_summary.csv"),
  row.names = FALSE
)

monitor_summary <- data.frame(
  Parameter = colnames(part$monitor),
  Mean = vapply(seq_len(ncol(part$monitor)), function(j) mean(part$monitor[, j]), numeric(1)),
  SD = vapply(seq_len(ncol(part$monitor)), function(j) stats::sd(part$monitor[, j]), numeric(1)),
  Min = vapply(seq_len(ncol(part$monitor)), function(j) min(part$monitor[, j]), numeric(1)),
  Max = vapply(seq_len(ncol(part$monitor)), function(j) max(part$monitor[, j]), numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(
  monitor_summary,
  file.path(OUT, "02_chain_monitor_summary.csv"),
  row.names = FALSE
)

tv_rows <- do.call(rbind, lapply(seq_len(K), function(eq) {
  out <- expand.grid(
    Quarter = model_quarters,
    Term = TERMS,
    stringsAsFactors = FALSE
  )
  out$Country <- COUNTRY
  out$Chain <- CHAIN_ID
  out$Equation <- VARS[eq]
  out$TimeVariationProbability <- as.vector(tv_probability[, , eq])
  out[, c(
    "Country", "Chain", "Equation", "Quarter", "Term",
    "TimeVariationProbability"
  )]
}))
write.csv(
  tv_rows,
  file.path(OUT, "03_time_variation_probability.csv"),
  row.names = FALSE
)

manifest <- data.frame(
  Item = c(
    "Country", "Chain", "Seed", "BurnPerChain", "PostBurnIterationsPerChain",
    "StoredDrawsPerChain", "SV", "LatentThresholdTVP", "MainNetwork",
    "Sample", "Variables", "Predictors", "threshtvpVersion",
    "threshtvpRemoteSha", "RDS"
  ),
  Value = c(
    COUNTRY,
    CHAIN_ID,
    SEED,
    BURN,
    KEEP,
    STORED_PER_CHAIN,
    SV_ON,
    TVS_ON,
    MAIN_NETWORK,
    paste0(model_quarters[1], " - ", tail(model_quarters, 1)),
    paste(VARS, collapse = ", "),
    paste(TERMS, collapse = ", "),
    pkg_version,
    remote_sha,
    rds_path
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(OUT, "00_chain_manifest.csv"),
  row.names = FALSE
)

readme <- c(
  "FORMAL BAYESIAN TVP-GVAR COUNTRY/CHAIN PART",
  "===========================================",
  sprintf("Country: %s", COUNTRY),
  sprintf("Chain: %d", CHAIN_ID),
  sprintf("Seed: %d", SEED),
  sprintf("Sample: %s - %s", model_quarters[1], tail(model_quarters, 1)),
  sprintf("Variables: %s", paste(VARS, collapse = ", ")),
  sprintf("Main financial network: %s", MAIN_NETWORK),
  sprintf("Burn: %d", BURN),
  sprintf("Internal post-burn MCMC iterations: %d", KEEP),
  sprintf("Stored strictly-post-burn draws: %d", STORED_PER_CHAIN),
  sprintf("Stochastic volatility: %s", SV_ON),
  sprintf("Latent-threshold TVP: %s", TVS_ON),
  sprintf("threshtvp version: %s", pkg_version),
  sprintf("threshtvp RemoteSha: %s", remote_sha),
  "",
  "The RDS stores only event-date coefficient draws, event stochastic volatility,",
  "compact convergence monitors, and posterior time-variation probabilities.",
  "The full T x K coefficient paths are intentionally not retained to control",
  "artifact size. They can be regenerated exactly from the chain specification."
)

writeLines(readme, file.path(OUT, "README_formal_chain.txt"))
cat(paste(readme, collapse = "\n"), "\n")
