#!/usr/bin/env Rscript

# Bayesian TVP-GVAR smoke test for the clean 14-economy / 3-variable system.
#
# Purpose:
#   * verify that a time-varying-parameter version of the accepted p=1, q=1
#     financial GVAR can be filtered and sampled without numerical failure;
#   * draw coefficient paths with a random-walk TVP state equation using
#     Kalman filtering + FFBS (base R only);
#   * rebuild G0_t and G1_t for every posterior draw and quarter;
#   * measure posterior global/local stability and G0 invertibility;
#   * produce a mechanical READY / FAIL smoke-test gate.
#
# Important:
#   This is NOT the final publication Bayesian TVP-GVAR. Hyperparameters are
#   fixed/empirical-Bayes and country equations are sampled independently.
#   It is a deliberately lightweight architecture test before the final model.

source("R/00_config.R")

PANEL <- file.path(DERIVED_DIR, "panel_domestic_fin3.csv")
PRE_GATE <- file.path(RESULTS_DIR, "gate", "01_estimation_gate.csv")
OUT <- file.path(RESULTS_DIR, "tvp_smoke")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("Environment variable %s is not numeric: %s", name, z)
  out
}

NDRAWS <- as.integer(get_env_num("TVP_SMOKE_DRAWS", 60))
SEED <- as.integer(get_env_num("TVP_SMOKE_SEED", 20260901))
STATE_SCALE <- get_env_num("TVP_STATE_SCALE", 1e-5)
PRIOR_SCALE <- get_env_num("TVP_PRIOR_SCALE", 4)
RIDGE <- get_env_num("TVP_RIDGE", 1e-8)
MIN_STABLE_SHARE <- get_env_num("TVP_MIN_STABLE_SHARE", 0.90)
MIN_G0_RCOND <- get_env_num("TVP_MIN_G0_RCOND", 1e-10)
MAX_SPLIT_DIFF <- get_env_num("TVP_MAX_SPLIT_DIFF", 0.10)
MAX_ABS_STD_BETA <- get_env_num("TVP_MAX_ABS_STD_BETA", 50)

if (NDRAWS < 20L) stopf("TVP_SMOKE_DRAWS must be at least 20.")
if (!(STATE_SCALE > 0)) stopf("TVP_STATE_SCALE must be positive.")
if (!(PRIOR_SCALE > 0)) stopf("TVP_PRIOR_SCALE must be positive.")
if (!(RIDGE > 0)) stopf("TVP_RIDGE must be positive.")
if (!(MIN_STABLE_SHARE > 0 && MIN_STABLE_SHARE <= 1)) stopf("TVP_MIN_STABLE_SHARE must be in (0,1].")
set.seed(SEED)

if (!file.exists(PANEL)) stopf("Run R/10_build_financial_3var_input.R first.")
if (!file.exists(PRE_GATE)) stopf("Run R/40_estimation_gate.R first.")

pre <- read.csv(PRE_GATE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(pre) != 1L || !"Status" %in% names(pre)) stopf("Malformed pre-estimation gate file.")
if (!identical(trimws(pre$Status[1]), "READY")) {
  stopf("Pre-estimation gate is not READY. Do not run TVP smoke test.")
}

W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])

d <- read.csv(PANEL, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("Quarter", "Country", VARS, "gpr", "brent")
if (!all(need %in% names(d))) stopf("Domestic panel is missing required columns.")

quarters <- unique(d$Quarter)
if (!all(COUNTRIES %in% unique(d$Country))) stopf("Domestic panel is missing sample economies.")
if (any(table(d$Country) != length(quarters))) stopf("Domestic panel is not balanced.")

N <- length(COUNTRIES)
K <- length(VARS)
Tn <- length(quarters)
if (Tn < 40L) stopf("TVP smoke test requires at least 40 quarters.")

X <- array(NA_real_, c(Tn, N, K), dimnames = list(quarters, COUNTRIES, VARS))
for (i in seq_along(COUNTRIES)) {
  z <- d[d$Country == COUNTRIES[i], , drop = FALSE]
  z <- z[match(quarters, z$Quarter), , drop = FALSE]
  X[, i, ] <- as.matrix(z[, VARS, drop = FALSE])
}
base <- d[d$Country == COUNTRIES[1], , drop = FALSE]
base <- base[match(quarters, base$Quarter), , drop = FALSE]
gpr <- num(base$gpr)
brent <- num(base$brent)
if (!all(is.finite(X)) || !all(is.finite(gpr)) || !all(is.finite(brent))) {
  stopf("Non-finite model inputs in final panel.")
}

lag1 <- function(x) c(NA_real_, x[-length(x)])

make_star <- function(W) {
  S <- array(NA_real_, dim(X), dimnames = dimnames(X))
  for (i in seq_len(N)) for (v in seq_len(K)) {
    S[, i, v] <- as.numeric(X[, , v] %*% W[i, ])
  }
  S
}
STAR <- make_star(W)

selector <- function(i) {
  S <- matrix(0, K, N * K)
  S[, ((i - 1L) * K + 1L):(i * K)] <- diag(K)
  S
}
star_map <- function(i, W) {
  R <- matrix(0, K, N * K)
  for (j in seq_len(N)) for (v in seq_len(K)) {
    R[v, (j - 1L) * K + v] <- W[i, j]
  }
  R
}
SELECTORS <- lapply(seq_len(N), selector)
STAR_MAPS <- lapply(seq_len(N), star_map, W = W)

sym <- function(A) (A + t(A)) / 2

safe_solve <- function(A, ridge = 0) {
  A <- sym(A)
  if (ridge > 0) A <- A + diag(ridge, nrow(A))
  out <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) {
    out <- tryCatch(qr.solve(A, diag(nrow(A)), tol = 1e-12), error = function(e) NULL)
  }
  if (is.null(out) || any(!is.finite(out))) stopf("Matrix inversion failed in TVP filter.")
  out
}

rmvnorm_psd <- function(mu, Sigma) {
  Sigma <- sym(Sigma)
  ee <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  as.numeric(mu + ee$vectors %*% (sqrt(vals) * rnorm(length(mu))))
}

build_design <- function(i) {
  Y <- X[, i, , drop = FALSE][, 1, ]
  Z <- STAR[, i, , drop = FALSE][, 1, ]
  rows <- 2:nrow(Y)

  D <- data.frame(const = rep(1, length(rows)), check.names = FALSE)
  for (v in seq_len(K)) D[[paste0(VARS[v], "_L1")]] <- lag1(Y[, v])[rows]
  for (v in seq_len(K)) D[[paste0(VARS[v], "_star_0")]] <- Z[rows, v]
  for (v in seq_len(K)) D[[paste0(VARS[v], "_star_L1")]] <- lag1(Z[, v])[rows]
  D$gpr_0 <- gpr[rows]
  D$gpr_L1 <- lag1(gpr)[rows]
  D$brent_0 <- brent[rows]
  D$brent_L1 <- lag1(brent)[rows]

  YY <- Y[rows, , drop = FALSE]
  ok <- complete.cases(D) & complete.cases(YY)
  D <- D[ok, , drop = FALSE]
  YY <- YY[ok, , drop = FALSE]
  qq <- quarters[rows][ok]

  if (nrow(D) <= ncol(D) + 10L) stopf("Too few observations for TVP design: %s", COUNTRIES[i])
  list(D = D, Y = YY, Quarter = qq)
}

standardize_predictors <- function(D) {
  Xm <- as.matrix(D)
  centers <- rep(0, ncol(Xm)); names(centers) <- colnames(Xm)
  scales <- rep(1, ncol(Xm)); names(scales) <- colnames(Xm)

  for (j in seq_len(ncol(Xm))) {
    if (colnames(Xm)[j] == "const") next
    centers[j] <- mean(Xm[, j])
    scales[j] <- sd(Xm[, j])
    if (!is.finite(scales[j]) || scales[j] < 1e-10) {
      stopf("Near-constant predictor in TVP design: %s", colnames(Xm)[j])
    }
    Xm[, j] <- (Xm[, j] - centers[j]) / scales[j]
  }
  list(X = Xm, center = centers, scale = scales)
}

# Filter a single country-equation once. FFBS draws can then reuse the stored
# smoothing matrices, which keeps the smoke test fast enough for GitHub Actions.
fit_tvp_filter <- function(country_i, eq_i, design) {
  std <- standardize_predictors(design$D)
  Xm <- std$X
  y <- as.numeric(design$Y[, eq_i])
  P <- ncol(Xm)
  M <- nrow(Xm)

  ols <- lm.fit(Xm, y)
  if (any(!is.finite(ols$coefficients))) {
    stopf("Rank deficiency in standardized TVP design: %s/%s", COUNTRIES[country_i], VARS[eq_i])
  }
  beta0 <- as.numeric(ols$coefficients)
  resid <- as.numeric(ols$residuals)
  sse <- sum(resid^2)
  df <- M - P
  if (df <= 5L) stopf("Insufficient residual degrees of freedom: %s/%s", COUNTRIES[country_i], VARS[eq_i])

  XtX <- crossprod(Xm)
  C0 <- PRIOR_SCALE * safe_solve(XtX, ridge = RIDGE)
  Wstate <- STATE_SCALE * diag(P)

  m <- matrix(NA_real_, M, P, dimnames = list(design$Quarter, colnames(Xm)))
  a <- matrix(NA_real_, M, P, dimnames = list(design$Quarter, colnames(Xm)))
  C_store <- array(NA_real_, c(P, P, M))
  R_store <- array(NA_real_, c(P, P, M))

  m_prev <- beta0
  C_prev <- C0
  for (tt in seq_len(M)) {
    at <- m_prev
    Rt <- sym(C_prev + Wstate)
    xt <- as.numeric(Xm[tt, ])
    rx <- as.numeric(Rt %*% xt)
    Qt <- 1 + sum(xt * rx)
    if (!is.finite(Qt) || Qt <= 1e-12) stopf("Invalid Kalman forecast variance: %s/%s", COUNTRIES[country_i], VARS[eq_i])
    err <- y[tt] - sum(xt * at)
    mt <- at + rx * (err / Qt)
    Ct <- sym(Rt - tcrossprod(rx) / Qt)

    a[tt, ] <- at
    m[tt, ] <- mt
    C_store[, , tt] <- Ct
    R_store[, , tt] <- Rt
    m_prev <- mt
    C_prev <- Ct
  }

  J <- if (M > 1L) array(NA_real_, c(P, P, M - 1L)) else array(0, c(P, P, 0L))
  H <- if (M > 1L) array(NA_real_, c(P, P, M - 1L)) else array(0, c(P, P, 0L))
  if (M > 1L) {
    for (tt in seq_len(M - 1L)) {
      Jt <- C_store[, , tt] %*% safe_solve(R_store[, , tt + 1L], ridge = RIDGE)
      Ht <- sym(C_store[, , tt] - Jt %*% C_store[, , tt])
      J[, , tt] <- Jt
      H[, , tt] <- Ht
    }
  }

  raw_kappa <- tryCatch(kappa(crossprod(as.matrix(design$D))), error = function(e) Inf)
  scaled_kappa <- tryCatch(kappa(XtX), error = function(e) Inf)

  list(
    country = COUNTRIES[country_i],
    equation = VARS[eq_i],
    quarter = design$Quarter,
    terms = colnames(Xm),
    center = std$center,
    scale = std$scale,
    m = m,
    a = a,
    J = J,
    H = H,
    C_last = C_store[, , M],
    ig_shape = df / 2,
    ig_rate = max(sse / 2, 1e-12),
    sigma2_ols = max(sse / df, 1e-12),
    raw_kappa = raw_kappa,
    scaled_kappa = scaled_kappa,
    max_abs_filter_mean = max(abs(m)),
    n = M,
    p = P
  )
}

draw_ffbs <- function(fit) {
  M <- nrow(fit$m)
  P <- ncol(fit$m)
  sigma2 <- 1 / rgamma(1, shape = fit$ig_shape, rate = fit$ig_rate)
  if (!is.finite(sigma2) || sigma2 <= 0) stopf("Invalid sigma2 draw in FFBS.")

  beta <- matrix(NA_real_, M, P, dimnames = dimnames(fit$m))
  beta[M, ] <- rmvnorm_psd(fit$m[M, ], sigma2 * fit$C_last)
  if (M > 1L) {
    for (tt in (M - 1L):1L) {
      h <- fit$m[tt, ] + fit$J[, , tt] %*% (beta[tt + 1L, ] - fit$a[tt + 1L, ])
      beta[tt, ] <- rmvnorm_psd(h, sigma2 * fit$H[, , tt])
    }
  }
  attr(beta, "sigma2") <- sigma2
  beta
}

# Build and filter all 42 country-equations.
designs <- lapply(seq_len(N), build_design)
model_quarters <- designs[[1]]$Quarter
if (!all(vapply(designs, function(z) identical(z$Quarter, model_quarters), logical(1)))) {
  stopf("Country TVP designs do not share the same quarterly sample.")
}
M <- length(model_quarters)

filters <- vector("list", N * K)
idx_fk <- function(i, eq) (i - 1L) * K + eq

for (i in seq_len(N)) {
  msg("Filtering %s ...", COUNTRIES[i])
  for (eq in seq_len(K)) {
    filters[[idx_fk(i, eq)]] <- fit_tvp_filter(i, eq, designs[[i]])
  }
}

filter_diag <- do.call(rbind, lapply(filters, function(f) data.frame(
  Country = f$country,
  Equation = f$equation,
  N = f$n,
  Predictors = f$p,
  Sigma2_OLS = f$sigma2_ols,
  RawDesignKappa = f$raw_kappa,
  StandardizedDesignKappa = f$scaled_kappa,
  MaxAbsFilteredMean = f$max_abs_filter_mean,
  stringsAsFactors = FALSE
)))
write.csv(filter_diag, file.path(OUT, "04_equation_filter_diagnostics.csv"), row.names = FALSE)

# One row per draw-quarter. At 60 draws and ~100 quarters this is small enough
# to retain as an audit trail.
stab_list <- vector("list", NDRAWS)
draw_diag <- vector("list", NDRAWS)
local_rho <- matrix(NA_real_, NDRAWS * M, N, dimnames = list(NULL, COUNTRIES))
local_pos <- 0L

term_dom <- paste0(VARS, "_L1")
term_star0 <- paste0(VARS, "_star_0")
term_star1 <- paste0(VARS, "_star_L1")

for (dd in seq_len(NDRAWS)) {
  if (dd == 1L || dd %% 10L == 0L || dd == NDRAWS) msg("Posterior draw %d / %d", dd, NDRAWS)

  Aarr <- array(0, c(M, N, K, K))
  B0arr <- array(0, c(M, N, K, K))
  B1arr <- array(0, c(M, N, K, K))
  max_abs_std <- 0
  sigma2_draws <- numeric(N * K)

  for (i in seq_len(N)) for (eq in seq_len(K)) {
    ff <- filters[[idx_fk(i, eq)]]
    path_std <- draw_ffbs(ff)
    sigma2_draws[idx_fk(i, eq)] <- attr(path_std, "sigma2")
    max_abs_std <- max(max_abs_std, abs(path_std))

    path_orig <- path_std
    for (jj in seq_len(ncol(path_orig))) {
      if (ff$terms[jj] != "const") path_orig[, jj] <- path_std[, jj] / ff$scale[jj]
    }

    for (v in seq_len(K)) {
      Aarr[, i, eq, v] <- path_orig[, match(term_dom[v], ff$terms)]
      B0arr[, i, eq, v] <- path_orig[, match(term_star0[v], ff$terms)]
      B1arr[, i, eq, v] <- path_orig[, match(term_star1[v], ff$terms)]
    }
  }

  rho_vec <- rep(NA_real_, M)
  rcond_vec <- rep(NA_real_, M)
  local_draw <- matrix(NA_real_, M, N)

  for (tt in seq_len(M)) {
    G0 <- matrix(0, N * K, N * K)
    G1 <- matrix(0, N * K, N * K)

    for (i in seq_len(N)) {
      Ai <- matrix(Aarr[tt, i, , ], K, K)
      B0i <- matrix(B0arr[tt, i, , ], K, K)
      B1i <- matrix(B1arr[tt, i, , ], K, K)
      local_draw[tt, i] <- max(Mod(eigen(Ai, only.values = TRUE)$values))

      rr <- ((i - 1L) * K + 1L):(i * K)
      Si <- SELECTORS[[i]]
      Ri <- STAR_MAPS[[i]]
      G0[rr, ] <- Si - B0i %*% Ri
      G1[rr, ] <- Ai %*% Si + B1i %*% Ri
    }

    rc <- tryCatch(rcond(G0), error = function(e) NA_real_)
    Fmat <- if (is.finite(rc) && rc > 0) tryCatch(solve(G0, G1), error = function(e) NULL) else NULL
    rho <- if (is.null(Fmat)) NA_real_ else tryCatch(max(Mod(eigen(Fmat, only.values = TRUE)$values)), error = function(e) NA_real_)
    rho_vec[tt] <- rho
    rcond_vec[tt] <- rc
  }

  rows <- (local_pos + 1L):(local_pos + M)
  local_rho[rows, ] <- local_draw
  local_pos <- local_pos + M

  stable <- is.finite(rho_vec) & rho_vec < 1
  g0_ok <- is.finite(rcond_vec) & rcond_vec >= MIN_G0_RCOND

  stab_list[[dd]] <- data.frame(
    Draw = dd,
    Quarter = model_quarters,
    GlobalSpectralRadius = rho_vec,
    G0_rcond = rcond_vec,
    Stable = stable,
    G0_OK = g0_ok,
    stringsAsFactors = FALSE
  )

  draw_diag[[dd]] <- data.frame(
    Draw = dd,
    StableShare = mean(stable),
    G0OKShare = mean(g0_ok),
    MaxGlobalSpectralRadius = if (any(is.finite(rho_vec))) max(rho_vec, na.rm = TRUE) else NA_real_,
    MinG0_rcond = if (any(is.finite(rcond_vec))) min(rcond_vec, na.rm = TRUE) else NA_real_,
    MaxAbsStandardizedBeta = max_abs_std,
    MedianSigma2Draw = median(sigma2_draws),
    stringsAsFactors = FALSE
  )
}

stab <- do.call(rbind, stab_list)
draw_df <- do.call(rbind, draw_diag)
write.csv(stab, file.path(OUT, "02_posterior_stability_draw_quarter.csv"), row.names = FALSE)
write.csv(draw_df, file.path(OUT, "03_posterior_stability_by_draw.csv"), row.names = FALSE)

# Quarter-level posterior stability summary.
qsplit <- split(stab, factor(stab$Quarter, levels = model_quarters))
quarter_summary <- do.call(rbind, lapply(qsplit, function(z) {
  rr <- z$GlobalSpectralRadius[is.finite(z$GlobalSpectralRadius)]
  rc <- z$G0_rcond[is.finite(z$G0_rcond)]
  data.frame(
    Quarter = z$Quarter[1],
    StableShare = mean(z$Stable),
    G0OKShare = mean(z$G0_OK),
    RhoMedian = if (length(rr)) median(rr) else NA_real_,
    RhoP05 = if (length(rr)) unname(quantile(rr, 0.05)) else NA_real_,
    RhoP95 = if (length(rr)) unname(quantile(rr, 0.95)) else NA_real_,
    RhoMax = if (length(rr)) max(rr) else NA_real_,
    MinG0_rcond = if (length(rc)) min(rc) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(quarter_summary, file.path(OUT, "01_posterior_stability_by_quarter.csv"), row.names = FALSE)

# Country-level local stability summary.
local_summary <- do.call(rbind, lapply(seq_len(N), function(i) {
  rr <- local_rho[, i]
  rr <- rr[is.finite(rr)]
  data.frame(
    Country = COUNTRIES[i],
    StableShare = if (length(rr)) mean(rr < 1) else NA_real_,
    RhoMedian = if (length(rr)) median(rr) else NA_real_,
    RhoP95 = if (length(rr)) unname(quantile(rr, 0.95)) else NA_real_,
    RhoMax = if (length(rr)) max(rr) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(local_summary, file.path(OUT, "05_local_posterior_stability.csv"), row.names = FALSE)

# Mechanical smoke-test gate.
finite_share <- mean(is.finite(stab$GlobalSpectralRadius) & is.finite(stab$G0_rcond))
stable_share <- mean(stab$Stable)
g0_ok_share <- mean(stab$G0_OK)
max_abs_beta <- max(draw_df$MaxAbsStandardizedBeta, na.rm = TRUE)

half <- floor(NDRAWS / 2)
first_draws <- seq_len(half)
second_draws <- (half + 1L):NDRAWS
first_share <- mean(stab$Stable[stab$Draw %in% first_draws])
second_share <- mean(stab$Stable[stab$Draw %in% second_draws])
split_diff <- abs(first_share - second_share)

finite_ok <- is.finite(finite_share) && finite_share == 1
stable_ok <- is.finite(stable_share) && stable_share >= MIN_STABLE_SHARE
g0_ok <- is.finite(g0_ok_share) && g0_ok_share >= 0.99
split_ok <- is.finite(split_diff) && split_diff <= MAX_SPLIT_DIFF
beta_ok <- is.finite(max_abs_beta) && max_abs_beta <= MAX_ABS_STD_BETA
ready <- finite_ok && stable_ok && g0_ok && split_ok && beta_ok
status <- if (ready) "READY" else "FAIL"

reasons <- c(
  if (!finite_ok) sprintf("finite draw-quarter share %.4f < 1", finite_share) else NULL,
  if (!stable_ok) sprintf("posterior stable share %.4f < %.4f", stable_share, MIN_STABLE_SHARE) else NULL,
  if (!g0_ok) sprintf("G0-OK share %.4f < 0.99", g0_ok_share) else NULL,
  if (!split_ok) sprintf("first/second-half stable-share difference %.4f > %.4f", split_diff, MAX_SPLIT_DIFF) else NULL,
  if (!beta_ok) sprintf("max standardized |beta| %.4f > %.4f", max_abs_beta, MAX_ABS_STD_BETA) else NULL
)
if (!length(reasons)) reasons <- "all TVP smoke-test gates passed"

smoke_gate <- data.frame(
  MainNetwork = MAIN_NETWORK,
  Status = status,
  Draws = NDRAWS,
  Seed = SEED,
  StateScale = STATE_SCALE,
  PriorScale = PRIOR_SCALE,
  Sample = paste0(model_quarters[1], " - ", tail(model_quarters, 1)),
  FiniteDrawQuarterShare = finite_share,
  PosteriorStableShare = stable_share,
  G0OKShare = g0_ok_share,
  FirstHalfStableShare = first_share,
  SecondHalfStableShare = second_share,
  SplitDifference = split_diff,
  MaxAbsStandardizedBeta = max_abs_beta,
  MinStableShareRequired = MIN_STABLE_SHARE,
  MinG0RcondRequired = MIN_G0_RCOND,
  Reason = paste(reasons, collapse = "; "),
  stringsAsFactors = FALSE
)
write.csv(smoke_gate, file.path(OUT, "00_tvp_smoke_gate.csv"), row.names = FALSE)

lines <- c(
  sprintf("TVP-GVAR SMOKE TEST: %s", status),
  "========================================",
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf("Sample: %s - %s", model_quarters[1], tail(model_quarters, 1)),
  sprintf("Variables: %s; p=1; q=1", paste(VARS, collapse = ", ")),
  sprintf("Posterior draws: %d; seed: %d", NDRAWS, SEED),
  sprintf("TVP state scale: %.8g", STATE_SCALE),
  sprintf("Prior scale: %.8g", PRIOR_SCALE),
  sprintf("Finite draw-quarter share: %.6f", finite_share),
  sprintf("Posterior stable share (rho < 1): %.6f [required >= %.4f]", stable_share, MIN_STABLE_SHARE),
  sprintf("G0 OK share (rcond >= %.3g): %.6f", MIN_G0_RCOND, g0_ok_share),
  sprintf("Stable share first half / second half: %.6f / %.6f", first_share, second_share),
  sprintf("Split difference: %.6f [required <= %.4f]", split_diff, MAX_SPLIT_DIFF),
  sprintf("Max standardized |beta|: %.6f [required <= %.4f]", max_abs_beta, MAX_ABS_STD_BETA),
  sprintf("Reason: %s", paste(reasons, collapse = "; ")),
  "",
  "Interpretation:",
  "- This is a conditional-hyperparameter Bayesian TVP architecture smoke test.",
  "- FFBS draws are independent conditional draws, not an MCMC chain; R-hat/ESS are therefore not applicable here.",
  "- READY means the accepted 3-variable GVAR architecture can proceed to the final Bayesian TVP specification.",
  "- FAIL means inspect quarter/draw/local-stability diagnostics before increasing the model complexity.",
  "- Do not report these smoke-test posterior draws as final paper estimates."
)
writeLines(lines, file.path(OUT, "README_tvp_smoke.txt"))
cat(paste(lines, collapse = "\n"), "\n")

if (!ready) quit(save = "no", status = 2L)
