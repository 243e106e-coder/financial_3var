#!/usr/bin/env Rscript

# =============================================================================
# 70_tr_data_design_forensic_audit.R
#
# 08p: TR REER_DLOG / lagged-GPR data-design forensic audit
#
# PURPOSE
# -------
# 08o established that the TR/de/h=1 anomaly is dominated by the local
# TR/de/gpr_L1 coefficient. 08p does NOT re-estimate the Bayesian TVP-GVAR.
# It audits the already-used Delta-r model input and the exact R/60 design to
# answer:
#
# 1) Is the large gpr_L1 coefficient already present in a static OLS/FWL fit?
# 2) Is it driven by extreme TR REER_DLOG observations?
# 3) Is it driven by extreme / low-residual-variance lagged GPR observations?
# 4) Is gpr_L1 nearly collinear with the other R/60 regressors?
# 5) Which quarters contribute most to the FWL coefficient?
# 6) Does leave-one-quarter-out or dropping the top influence points collapse it?
# 7) Is the coefficient stable across broad subsamples / rolling windows?
#
# IMPORTANT
# ---------
# * No MCMC.
# * No prior changes.
# * No country deletion.
# * No weight change.
# * No automatic outlier deletion.
# * Robust / leave-out estimates are diagnostics only.
# * The exact R/60 p=1,q=1 design is reconstructed.
# =============================================================================

source("R/00_config.R")

get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}
get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  x <- suppressWarnings(as.numeric(z))
  if (!is.finite(x)) stopf("Environment variable %s is not numeric: %s", name, z)
  x
}

PANEL <- get_env_chr(
  "FIN3_08P_PANEL",
  file.path(DERIVED_DIR, "panel_domestic_fin3_rate_diff.csv")
)
SOURCE_08O_ROOT <- get_env_chr("FIN3_08O_ROOT", "source_08o")
OUT <- get_env_chr(
  "FIN3_08P_OUT",
  file.path(RESULTS_DIR, "tr_data_design_forensic")
)

SOURCE_08J_RUN <- get_env_chr("FIN3_08J_RUN_ID", "")
SOURCE_08O_RUN <- get_env_chr("FIN3_08O_RUN_ID", "")
SOURCE_08O_SHA <- get_env_chr("FIN3_08O_HEAD_SHA", "")

TARGET_COUNTRY <- "TR"
TARGET_EQUATION <- "de"
TARGET_TERM <- "gpr_L1"
GPR_SHOCK_PCT <- get_env_num("FIN3_GPR_SHOCK_PCT", 10)
ROLLING_WINDOW <- as.integer(get_env_num("FIN3_08P_ROLLING_WINDOW", 40))
ROBUST_Z_CUTOFF <- get_env_num("FIN3_08P_ROBUST_Z_CUTOFF", 3.5)
VIF_FLAG <- get_env_num("FIN3_08P_VIF_FLAG", 10)
LOO_REL_CHANGE_FLAG <- get_env_num("FIN3_08P_LOO_REL_CHANGE_FLAG", 0.25)
TOPK_REL_CHANGE_FLAG <- get_env_num("FIN3_08P_TOPK_REL_CHANGE_FLAG", 0.50)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!(TARGET_COUNTRY %in% COUNTRIES)) stopf("Target country not configured.")
if (!(TARGET_EQUATION %in% VARS)) stopf("Target equation not configured.")
if (!(GPR_SHOCK_PCT > 0)) stopf("GPR shock must be positive.")
if (!(ROLLING_WINDOW >= 30L)) stopf("Rolling window must be >=30 quarters.")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

lag1 <- function(x) c(NA_real_, x[-length(x)])

qv <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, names = FALSE, type = 7))
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  sc <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(sc) || sc < 1e-12) return(rep(NA_real_, length(x)))
  (x - med) / sc
}

ols_fit <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  fit <- stats::lm.fit(x = X, y = y)
  if (fit$rank < ncol(X)) {
    return(list(
      ok = FALSE, fit = fit, coef = rep(NA_real_, ncol(X)),
      se = rep(NA_real_, ncol(X)), vcov = matrix(NA_real_, ncol(X), ncol(X)),
      sigma = NA_real_
    ))
  }
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  rss <- sum(fit$residuals^2)
  sigma2 <- rss / df
  xtx_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(xtx_inv)) {
    return(list(
      ok = FALSE, fit = fit, coef = fit$coefficients,
      se = rep(NA_real_, ncol(X)), vcov = matrix(NA_real_, ncol(X), ncol(X)),
      sigma = sqrt(sigma2)
    ))
  }
  vc <- sigma2 * xtx_inv
  list(
    ok = TRUE,
    fit = fit,
    coef = fit$coefficients,
    se = sqrt(diag(vc)),
    vcov = vc,
    sigma = sqrt(sigma2)
  )
}

huber_irls <- function(X, y, c = 1.345, maxit = 200L, tol = 1e-9) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  init <- ols_fit(X, y)
  if (!init$ok) {
    return(list(coef = rep(NA_real_, ncol(X)), converged = FALSE, iterations = 0L))
  }
  b <- init$coef
  converged <- FALSE

  for (it in seq_len(maxit)) {
    r <- as.numeric(y - X %*% b)
    s <- stats::mad(r, center = 0, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-12) {
      converged <- TRUE
      break
    }
    u <- r / s
    w <- pmin(1, c / pmax(abs(u), 1e-12))
    sw <- sqrt(w)
    fw <- stats::lm.fit(X * sw, y * sw)
    if (fw$rank < ncol(X)) break
    b_new <- fw$coefficients
    if (max(abs(b_new - b), na.rm = TRUE) < tol) {
      b <- b_new
      converged <- TRUE
      break
    }
    b <- b_new
  }

  list(coef = b, converged = converged, iterations = it)
}

fwl_residual <- function(target, controls) {
  fit <- stats::lm.fit(x = as.matrix(controls), y = as.numeric(target))
  if (fit$rank < ncol(as.matrix(controls))) {
    stopf("FWL control design is rank deficient.")
  }
  fit$residuals
}

# -----------------------------------------------------------------------------
# 0. Validate accepted 08o source
# -----------------------------------------------------------------------------

gate08o_file <- file.path(SOURCE_08O_ROOT, "00_08o_gate.csv")
coef08o_file <- file.path(SOURCE_08O_ROOT, "09_tr_de_core_event_coefficients.csv")
guard08o_file <- file.path(SOURCE_08O_ROOT, "14_claims_guardrail.csv")

for (f in c(gate08o_file, coef08o_file, guard08o_file)) {
  if (!file.exists(f)) stopf("Missing accepted 08o artifact file: %s", f)
}

gate08o <- read.csv(gate08o_file, stringsAsFactors = FALSE, check.names = FALSE)
coef08o <- read.csv(coef08o_file, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(gate08o) != 1L || gate08o$Status[1] != "DIAGNOSTIC_COMPLETE") {
  stopf("08o gate is not DIAGNOSTIC_COMPLETE.")
}
if (gate08o$H1SourceClassification[1] != "LOCAL_TR_DE_GPR_L1_DOMINANT") {
  stopf("08p requires the accepted local TR/de/gpr_L1 diagnosis from 08o.")
}
if (isTRUE(gate08o$ReestimationPerformed[1]) ||
    isTRUE(gate08o$ModelSpecificationChanged[1]) ||
    isTRUE(gate08o$WeightMatrixChanged[1])) {
  stopf("Accepted 08o source violated read-only diagnostic contract.")
}

# -----------------------------------------------------------------------------
# 1. Reconstruct exact R/60 Delta-r design for TR
# -----------------------------------------------------------------------------

if (!file.exists(PANEL)) stopf("Missing Delta-r panel: %s", PANEL)

d <- read.csv(PANEL, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("Quarter", "Country", VARS, "gpr", "brent")
if (!all(need %in% names(d))) {
  stopf("Delta-r panel missing columns: %s", paste(setdiff(need, names(d)), collapse = ", "))
}

d$Country <- toupper(trimws(as.character(d$Country)))
d$Quarter <- toupper(trimws(as.character(d$Quarter)))

quarters <- unique(d$Quarter)
if (!all(COUNTRIES %in% unique(d$Country))) stopf("Panel missing configured countries.")
if (any(table(d$Country) != length(quarters))) stopf("Panel is not balanced.")

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
  stopf("Non-finite values in Delta-r formal inputs.")
}

W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])
country_i <- match(TARGET_COUNTRY, COUNTRIES)

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

Yraw <- Xglobal[, country_i, , drop = FALSE][, 1, ]
Zraw <- STAR[, country_i, , drop = FALSE][, 1, ]
rows <- 2:nrow(Yraw)

D <- data.frame(const = rep(1, length(rows)), check.names = FALSE)
for (v in seq_len(K)) D[[paste0(VARS[v], "_L1")]] <- lag1(Yraw[, v])[rows]
for (v in seq_len(K)) D[[paste0(VARS[v], "_star_0")]] <- Zraw[rows, v]
for (v in seq_len(K)) D[[paste0(VARS[v], "_star_L1")]] <- lag1(Zraw[, v])[rows]
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

EXPECTED_TERMS <- c(
  "const",
  paste0(VARS, "_L1"),
  paste0(VARS, "_star_0"),
  paste0(VARS, "_star_L1"),
  "gpr_0", "gpr_L1",
  "brent_0", "brent_L1"
)
if (!identical(colnames(D), EXPECTED_TERMS)) {
  stopf("08p design does not reproduce R/60 term ordering.")
}
if (!(TARGET_TERM %in% colnames(D))) stopf("Target term missing.")

y <- as.numeric(Y[, match(TARGET_EQUATION, VARS)])
X <- as.matrix(D)
term_j <- match(TARGET_TERM, colnames(X))

if (nrow(X) <= ncol(X) + 10L) stopf("Too few observations for forensic audit.")

design_manifest <- data.frame(
  Item = c(
    "Source08jRunID","Source08oRunID","Source08oHeadSHA",
    "Panel","Country","Equation","TargetTerm","Observations","Predictors",
    "SampleStart","SampleEnd","MainNetwork","GPRColumn",
    "GPRShockPct","GPRShockLogUnits"
  ),
  Value = c(
    SOURCE_08J_RUN, SOURCE_08O_RUN, SOURCE_08O_SHA,
    basename(PANEL), TARGET_COUNTRY, TARGET_EQUATION, TARGET_TERM,
    nrow(X), ncol(X), model_quarters[1], tail(model_quarters,1),
    MAIN_NETWORK, GPR_COLUMN, GPR_SHOCK_PCT, log1p(GPR_SHOCK_PCT/100)
  ),
  stringsAsFactors = FALSE
)
write.csv(design_manifest, file.path(OUT, "00_design_manifest.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Raw scale and outlier audit
# -----------------------------------------------------------------------------

raw_tbl <- data.frame(
  Quarter = model_quarters,
  de = y,
  gpr_0 = X[, "gpr_0"],
  gpr_L1 = X[, "gpr_L1"],
  brent_0 = X[, "brent_0"],
  brent_L1 = X[, "brent_L1"],
  stringsAsFactors = FALSE
)
raw_tbl$de_robust_z <- robust_z(raw_tbl$de)
raw_tbl$gpr_L1_robust_z <- robust_z(raw_tbl$gpr_L1)

raw_scale <- data.frame(
  Series = c("TR_de_REER_DLOG","gpr_0","gpr_L1"),
  N = c(length(y), nrow(X), nrow(X)),
  Mean = c(mean(y), mean(X[,"gpr_0"]), mean(X[,"gpr_L1"])),
  SD = c(sd(y), sd(X[,"gpr_0"]), sd(X[,"gpr_L1"])),
  Median = c(median(y), median(X[,"gpr_0"]), median(X[,"gpr_L1"])),
  P01 = c(qv(y,.01), qv(X[,"gpr_0"],.01), qv(X[,"gpr_L1"],.01)),
  P05 = c(qv(y,.05), qv(X[,"gpr_0"],.05), qv(X[,"gpr_L1"],.05)),
  P95 = c(qv(y,.95), qv(X[,"gpr_0"],.95), qv(X[,"gpr_L1"],.95)),
  P99 = c(qv(y,.99), qv(X[,"gpr_0"],.99), qv(X[,"gpr_L1"],.99)),
  Min = c(min(y), min(X[,"gpr_0"]), min(X[,"gpr_L1"])),
  Max = c(max(y), max(X[,"gpr_0"]), max(X[,"gpr_L1"])),
  MaxAbsRobustZ = c(
    max(abs(raw_tbl$de_robust_z), na.rm=TRUE),
    max(abs(robust_z(X[,"gpr_0"])), na.rm=TRUE),
    max(abs(raw_tbl$gpr_L1_robust_z), na.rm=TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(raw_scale, file.path(OUT, "01_raw_scale_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Full static OLS + standardized OLS
# -----------------------------------------------------------------------------

fit_full <- ols_fit(X, y)
if (!fit_full$ok) stopf("Full R/60 static design is rank deficient.")

beta_full <- unname(fit_full$coef[term_j])
se_full <- unname(fit_full$se[term_j])
t_full <- beta_full / se_full

# Exact standardized representation used conceptually by R/60 before MCMC.
Xs <- X
x_center <- rep(0, ncol(X)); names(x_center) <- colnames(X)
x_scale <- rep(1, ncol(X)); names(x_scale) <- colnames(X)
for (j in seq_len(ncol(Xs))) {
  if (colnames(Xs)[j] == "const") next
  x_center[j] <- mean(Xs[,j])
  x_scale[j] <- sd(Xs[,j])
  Xs[,j] <- (Xs[,j] - x_center[j]) / x_scale[j]
}
ys <- (y - mean(y)) / sd(y)
fit_std <- ols_fit(Xs, ys)
if (!fit_std$ok) stopf("Standardized static design is rank deficient.")
beta_std <- unname(fit_std$coef[term_j])

raw_kappa <- tryCatch(kappa(crossprod(X)), error = function(e) Inf)
scaled_kappa <- tryCatch(kappa(crossprod(Xs)), error = function(e) Inf)

# -----------------------------------------------------------------------------
# 4. FWL partialling-out for gpr_L1
# -----------------------------------------------------------------------------

Z <- X[, -term_j, drop = FALSE]
x_target <- X[, term_j]

rx <- fwl_residual(x_target, Z)
ry <- fwl_residual(y, Z)

fwl_beta <- sum(rx * ry) / sum(rx^2)
fwl_gap <- abs(fwl_beta - beta_full)

fit_x_on_others <- ols_fit(Z, x_target)
if (!fit_x_on_others$ok) stopf("gpr_L1-on-controls regression is rank deficient.")
sst_x <- sum((x_target - mean(x_target))^2)
sse_x <- sum(fit_x_on_others$fit$residuals^2)
r2_x_others <- 1 - sse_x/sst_x
vif_gpr_l1 <- 1 / pmax(1-r2_x_others, 1e-15)
partial_corr <- stats::cor(rx, ry)
partial_r2 <- partial_corr^2
resid_x_sd_ratio <- sd(rx) / sd(x_target)

# Observation-level FWL contribution:
# beta = sum(rx_i * ry_i) / sum(rx_i^2)
fwl_contrib <- (rx * ry) / sum(rx^2)
raw_tbl$gpr_L1_FWL_residual <- rx
raw_tbl$de_FWL_residual <- ry
raw_tbl$gpr_L1_FWL_residual_robust_z <- robust_z(rx)
raw_tbl$de_FWL_residual_robust_z <- robust_z(ry)
raw_tbl$FWL_beta_contribution <- fwl_contrib
raw_tbl$AbsFWLContribution <- abs(fwl_contrib)
raw_tbl$FWLContributionAbsShare <- abs(fwl_contrib) / sum(abs(fwl_contrib))

# -----------------------------------------------------------------------------
# 5. Classical influence diagnostics
# -----------------------------------------------------------------------------

lm_data <- data.frame(y = y, D, check.names = FALSE)
form <- stats::as.formula(
  paste("y ~", paste(colnames(D)[colnames(D)!="const"], collapse = " + "))
)
lm_obj <- stats::lm(form, data = lm_data)

# This formula-generated intercept must match D$const.
coef_formula <- stats::coef(lm_obj)
if (!(TARGET_TERM %in% names(coef_formula))) stopf("Target coefficient absent from lm object.")
if (abs(unname(coef_formula[TARGET_TERM]) - beta_full) > 1e-10) {
  stopf("lm formula coefficient does not match exact matrix OLS.")
}

hat <- stats::hatvalues(lm_obj)
cook <- stats::cooks.distance(lm_obj)
stud <- stats::rstudent(lm_obj)
dfb <- stats::dfbetas(lm_obj)
dfb_target <- dfb[, TARGET_TERM]

raw_tbl$Hat = hat
raw_tbl$CooksD = cook
raw_tbl$StudentizedResidual = stud
raw_tbl$DFBETA_gpr_L1 = dfb_target

# FWL outlier flags are diagnostic only; never automatic deletion.
raw_tbl$Flag_de_RobustZ <- abs(raw_tbl$de_robust_z) > ROBUST_Z_CUTOFF
raw_tbl$Flag_gprL1_RobustZ <- abs(raw_tbl$gpr_L1_robust_z) > ROBUST_Z_CUTOFF
raw_tbl$Flag_FWLx_RobustZ <- abs(raw_tbl$gpr_L1_FWL_residual_robust_z) > ROBUST_Z_CUTOFF
raw_tbl$Flag_FWLy_RobustZ <- abs(raw_tbl$de_FWL_residual_robust_z) > ROBUST_Z_CUTOFF

# Rank by multiple notions of influence.
raw_tbl$RankAbsFWLContribution <- rank(-raw_tbl$AbsFWLContribution, ties.method="min")
raw_tbl$RankCook <- rank(-raw_tbl$CooksD, ties.method="min")
raw_tbl$RankAbsDFBETA <- rank(-abs(raw_tbl$DFBETA_gpr_L1), ties.method="min")
write.csv(raw_tbl, file.path(OUT, "02_quarter_level_forensic_table.csv"), row.names = FALSE)

top_union <- raw_tbl[
  raw_tbl$RankAbsFWLContribution <= 15 |
  raw_tbl$RankCook <= 15 |
  raw_tbl$RankAbsDFBETA <= 15 |
  raw_tbl$Flag_de_RobustZ |
  raw_tbl$Flag_gprL1_RobustZ |
  raw_tbl$Flag_FWLx_RobustZ |
  raw_tbl$Flag_FWLy_RobustZ,
  ,
  drop = FALSE
]
top_union <- top_union[
  order(top_union$RankAbsFWLContribution, top_union$RankAbsDFBETA, top_union$RankCook),
  ,
  drop = FALSE
]
write.csv(top_union, file.path(OUT, "03_top_influential_quarters.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Leave-one-quarter-out coefficient audit
# -----------------------------------------------------------------------------

loo_rows <- vector("list", nrow(X))
for (i in seq_len(nrow(X))) {
  keep <- setdiff(seq_len(nrow(X)), i)
  z <- ols_fit(X[keep,,drop=FALSE], y[keep])
  b <- if (z$ok) unname(z$coef[term_j]) else NA_real_
  loo_rows[[i]] <- data.frame(
    DroppedQuarter = model_quarters[i],
    Dropped_de = y[i],
    Dropped_gpr_L1 = x_target[i],
    Dropped_CooksD = cook[i],
    Dropped_DFBETA_gpr_L1 = dfb_target[i],
    Dropped_FWLContribution = fwl_contrib[i],
    gpr_L1_Coefficient = b,
    ChangeFromFull = b - beta_full,
    AbsChangeFromFull = abs(b - beta_full),
    RelativeAbsChange = abs(b - beta_full) / pmax(abs(beta_full), 1e-15),
    stringsAsFactors = FALSE
  )
}
loo <- do.call(rbind, loo_rows)
loo <- loo[order(-loo$AbsChangeFromFull), , drop = FALSE]
write.csv(loo, file.path(OUT, "04_leave_one_quarter_out.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Pre-specified top-k influence deletion diagnostics
#    Diagnostic only. Never used as publication model selection.
# -----------------------------------------------------------------------------

topk_rows <- list()
kk <- 0L
criteria <- list(
  Cook = order(-cook),
  AbsDFBETA_gpr_L1 = order(-abs(dfb_target)),
  AbsFWLContribution = order(-abs(fwl_contrib))
)

for (criterion in names(criteria)) {
  ord <- criteria[[criterion]]
  for (kdrop in c(1L,2L,3L,5L)) {
    if (kdrop >= nrow(X) - ncol(X) - 2L) next
    drop_idx <- ord[seq_len(kdrop)]
    keep <- setdiff(seq_len(nrow(X)), drop_idx)
    z <- ols_fit(X[keep,,drop=FALSE], y[keep])
    b <- if (z$ok) unname(z$coef[term_j]) else NA_real_

    kk <- kk + 1L
    topk_rows[[kk]] <- data.frame(
      Criterion = criterion,
      K = kdrop,
      DroppedQuarters = paste(model_quarters[drop_idx], collapse = ";"),
      gpr_L1_Coefficient = b,
      ChangeFromFull = b - beta_full,
      RelativeAbsChange = abs(b-beta_full)/pmax(abs(beta_full),1e-15),
      RemainingN = length(keep),
      stringsAsFactors = FALSE
    )
  }
}
topk <- do.call(rbind, topk_rows)
write.csv(topk, file.path(OUT, "05_topk_influence_deletion_diagnostic.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Core-event-window deletion
# -----------------------------------------------------------------------------

event_coef08o <- coef08o[
  coef08o$Country == TARGET_COUNTRY &
    coef08o$Equation == TARGET_EQUATION &
    coef08o$Term == TARGET_TERM,
  ,
  drop = FALSE
]
if (!nrow(event_coef08o)) stopf("08o TR/de/gpr_L1 event coefficient rows missing.")

event_quarters <- unique(event_coef08o[, c("EventID","AnchorQuarter"), drop=FALSE])
event_quarters$qid <- quarter_id(event_quarters$AnchorQuarter)
model_qid <- quarter_id(model_quarters)

event_drop_rows <- list()
kk <- 0L
for (i in seq_len(nrow(event_quarters))) {
  q0 <- event_quarters$qid[i]
  drop_idx <- which(model_qid %in% (q0 + c(-1L,0L,1L)))
  keep <- setdiff(seq_len(nrow(X)), drop_idx)
  z <- ols_fit(X[keep,,drop=FALSE], y[keep])
  b <- if (z$ok) unname(z$coef[term_j]) else NA_real_

  kk <- kk + 1L
  event_drop_rows[[kk]] <- data.frame(
    EventID = event_quarters$EventID[i],
    AnchorQuarter = event_quarters$AnchorQuarter[i],
    Window = "t-1;t0;t+1",
    DroppedQuarters = paste(model_quarters[drop_idx], collapse = ";"),
    gpr_L1_Coefficient = b,
    ChangeFromFull = b-beta_full,
    RelativeAbsChange = abs(b-beta_full)/pmax(abs(beta_full),1e-15),
    stringsAsFactors = FALSE
  )
}
all_event_idx <- unique(unlist(lapply(event_quarters$qid, function(q0) {
  which(model_qid %in% (q0 + c(-1L,0L,1L)))
})))
keep <- setdiff(seq_len(nrow(X)), all_event_idx)
z <- ols_fit(X[keep,,drop=FALSE], y[keep])
b_all <- if (z$ok) unname(z$coef[term_j]) else NA_real_
event_drop_rows[[length(event_drop_rows)+1L]] <- data.frame(
  EventID = "ALL_SIX_CORE_WINDOWS",
  AnchorQuarter = NA_character_,
  Window = "all six t-1;t0;t+1 windows",
  DroppedQuarters = paste(model_quarters[all_event_idx], collapse = ";"),
  gpr_L1_Coefficient = b_all,
  ChangeFromFull = b_all-beta_full,
  RelativeAbsChange = abs(b_all-beta_full)/pmax(abs(beta_full),1e-15),
  stringsAsFactors = FALSE
)
event_drop <- do.call(rbind, event_drop_rows)
write.csv(event_drop, file.path(OUT, "06_core_event_window_deletion_diagnostic.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. Broad subsamples and rolling-window OLS
# -----------------------------------------------------------------------------

mid <- floor(nrow(X)/2)
subsample_specs <- list(
  Full = seq_len(nrow(X)),
  FirstHalf = seq_len(mid),
  SecondHalf = (mid+1L):nrow(X)
)

# Calendar split, if both sides have sufficient observations.
qid2013 <- quarter_id("2013Q1")
pre2013 <- which(model_qid < qid2013)
post2013 <- which(model_qid >= qid2013)
if (length(pre2013) > ncol(X)+10L && length(post2013) > ncol(X)+10L) {
  subsample_specs$Pre2013 <- pre2013
  subsample_specs$Post2013 <- post2013
}

sub_rows <- list()
kk <- 0L
for (nm in names(subsample_specs)) {
  ii <- subsample_specs[[nm]]
  z <- ols_fit(X[ii,,drop=FALSE], y[ii])
  b <- if (z$ok) unname(z$coef[term_j]) else NA_real_
  kk <- kk + 1L
  sub_rows[[kk]] <- data.frame(
    Subsample = nm,
    StartQuarter = model_quarters[min(ii)],
    EndQuarter = model_quarters[max(ii)],
    N = length(ii),
    gpr_L1_Coefficient = b,
    RelativeToFull = b / beta_full,
    stringsAsFactors = FALSE
  )
}
subsample <- do.call(rbind, sub_rows)
write.csv(subsample, file.path(OUT, "07_subsample_coefficient_diagnostic.csv"), row.names = FALSE)

roll_rows <- list()
kk <- 0L
if (ROLLING_WINDOW > ncol(X)+10L && ROLLING_WINDOW <= nrow(X)) {
  for (end_i in ROLLING_WINDOW:nrow(X)) {
    ii <- (end_i-ROLLING_WINDOW+1L):end_i
    z <- ols_fit(X[ii,,drop=FALSE], y[ii])
    b <- if (z$ok) unname(z$coef[term_j]) else NA_real_
    kk <- kk + 1L
    roll_rows[[kk]] <- data.frame(
      StartQuarter = model_quarters[min(ii)],
      EndQuarter = model_quarters[max(ii)],
      N = length(ii),
      gpr_L1_Coefficient = b,
      stringsAsFactors = FALSE
    )
  }
}
rolling <- if (length(roll_rows)) do.call(rbind, roll_rows) else data.frame()
write.csv(rolling, file.path(OUT, "08_rolling_window_coefficient_diagnostic.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. Huber IRLS robust diagnostic
# -----------------------------------------------------------------------------

hub <- huber_irls(X, y)
beta_huber <- if (hub$converged) unname(hub$coef[term_j]) else NA_real_

# -----------------------------------------------------------------------------
# 11. Compare static design evidence with accepted 08o posterior coefficient
# -----------------------------------------------------------------------------

post_medians <- event_coef08o$median
post_center <- median(post_medians, na.rm = TRUE)
post_range <- diff(range(post_medians, finite = TRUE))

shock_log <- log1p(GPR_SHOCK_PCT/100)
static_shock_effect <- beta_full * shock_log
huber_shock_effect <- beta_huber * shock_log
posterior_shock_effect <- post_center * shock_log

comparison <- data.frame(
  Quantity = c(
    "Static_full_OLS_gpr_L1",
    "Static_standardized_OLS_gpr_L1",
    "FWL_gpr_L1",
    "Huber_IRLS_gpr_L1",
    "Accepted_08o_event_posterior_median_center",
    "Accepted_08o_event_posterior_median_range",
    "Static_OLS_10pct_GPR_direct_effect",
    "Huber_10pct_GPR_direct_effect",
    "Posterior_center_10pct_GPR_direct_effect",
    "TR_de_sample_SD",
    "TR_de_sample_MaxAbs"
  ),
  Value = c(
    beta_full,
    beta_std,
    fwl_beta,
    beta_huber,
    post_center,
    post_range,
    static_shock_effect,
    huber_shock_effect,
    posterior_shock_effect,
    sd(y),
    max(abs(y))
  ),
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(OUT, "09_static_vs_posterior_scale_comparison.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 12. Summary metrics and transparent flags
# -----------------------------------------------------------------------------

max_loo_rel <- max(loo$RelativeAbsChange, na.rm = TRUE)
max_top3_rel <- max(
  topk$RelativeAbsChange[topk$K <= 3L],
  na.rm = TRUE
)

posterior_ols_rel_gap <- abs(post_center-beta_full)/pmax(abs(post_center),1e-15)
huber_ols_rel_gap <- abs(beta_huber-beta_full)/pmax(abs(beta_full),1e-15)

flag_high_vif <- is.finite(vif_gpr_l1) && vif_gpr_l1 >= VIF_FLAG
flag_single_obs <- is.finite(max_loo_rel) && max_loo_rel >= LOO_REL_CHANGE_FLAG
flag_topfew <- is.finite(max_top3_rel) && max_top3_rel >= TOPK_REL_CHANGE_FLAG
flag_robust_shift <- is.finite(huber_ols_rel_gap) && huber_ols_rel_gap >= 0.25
flag_posterior_tracks_static <- is.finite(posterior_ols_rel_gap) && posterior_ols_rel_gap <= 0.25

# A descriptive source diagnosis; thresholds are written to the gate.
diagnosis <- if (flag_high_vif && (flag_single_obs || flag_topfew)) {
  "COLLINEARITY_PLUS_INFLUENCE_RISK"
} else if (flag_high_vif) {
  "COLLINEARITY_RISK"
} else if (flag_single_obs || flag_topfew || flag_robust_shift) {
  "OBSERVATION_INFLUENCE_RISK"
} else if (flag_posterior_tracks_static) {
  "PERSISTENT_FULL_SAMPLE_STATIC_RELATION"
} else {
  "MIXED_OR_BAYESIAN_SPECIFIC_RELATION"
}

summary_metrics <- data.frame(
  Metric = c(
    "FullOLS_gpr_L1","FullOLS_SE","FullOLS_t",
    "StandardizedOLS_gpr_L1","FWL_gpr_L1","FWL_OLS_gap",
    "PartialCorrelation","PartialR2",
    "R2_gprL1_on_other_predictors","VIF_gpr_L1",
    "ResidualSD_gprL1_to_rawSD_ratio",
    "RawDesignKappa","ScaledDesignKappa",
    "HuberIRLS_gpr_L1","HuberConverged",
    "PosteriorEventMedianCenter","PosteriorEventMedianRange",
    "PosteriorVsOLSRelativeGap","HuberVsOLSRelativeGap",
    "MaxLeaveOneOutRelativeChange","MaxTop3DeletionRelativeChange",
    "MaxAbsQuarterFWLContribution",
    "Top5AbsFWLContributionShare",
    "de_RobustZ_FlagCount","gprL1_RobustZ_FlagCount",
    "FWLx_RobustZ_FlagCount","FWLy_RobustZ_FlagCount"
  ),
  Value = c(
    beta_full,se_full,t_full,
    beta_std,fwl_beta,fwl_gap,
    partial_corr,partial_r2,
    r2_x_others,vif_gpr_l1,
    resid_x_sd_ratio,
    raw_kappa,scaled_kappa,
    beta_huber,hub$converged,
    post_center,post_range,
    posterior_ols_rel_gap,huber_ols_rel_gap,
    max_loo_rel,max_top3_rel,
    max(abs(fwl_contrib)),
    sum(sort(abs(fwl_contrib), decreasing=TRUE)[1:min(5,length(fwl_contrib))]) /
      sum(abs(fwl_contrib)),
    sum(raw_tbl$Flag_de_RobustZ, na.rm=TRUE),
    sum(raw_tbl$Flag_gprL1_RobustZ, na.rm=TRUE),
    sum(raw_tbl$Flag_FWLx_RobustZ, na.rm=TRUE),
    sum(raw_tbl$Flag_FWLy_RobustZ, na.rm=TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_metrics, file.path(OUT, "10_forensic_summary_metrics.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 13. Claims guardrail + gate
# -----------------------------------------------------------------------------

guardrail <- data.frame(
  Claim = c(
    "08p re-estimates the Bayesian TVP-GVAR",
    "08p changes the financial weight matrix",
    "08p automatically deletes influential quarters",
    "FWL gpr_L1 coefficient exactly reconciles with full static OLS",
    "A large VIF proves the Bayesian posterior is wrong",
    "A large leave-one-out change proves the underlying data are erroneous",
    "The static OLS coefficient can diagnose whether the 5.11 posterior loading is rooted in the same sample/design",
    "Huber IRLS is a robustness diagnostic, not the publication model",
    "The 10 percent GPR effect uses log(1.1), consistent with logged GPR"
  ),
  Supported = c(
    FALSE,FALSE,FALSE,
    fwl_gap < 1e-10,
    FALSE,FALSE,
    TRUE,TRUE,
    grepl("^LN_", toupper(GPR_COLUMN))
  ),
  EvidenceOrLimit = c(
    "No threshtvp/MCMC estimator is called.",
    sprintf("Uses the accepted repository MAIN_NETWORK=%s only to reconstruct R/60 stars.", MAIN_NETWORK),
    "Influence deletions are reported as diagnostics only.",
    sprintf("absolute FWL-vs-OLS gap = %.12g", fwl_gap),
    "VIF diagnoses linear dependence in the static design; it does not invalidate Bayesian estimation by itself.",
    "Influence can arise from genuine crisis observations as well as data problems.",
    "Same Delta-r panel, same p=1/q=1 regressors, same financial W, same original units.",
    "Huber coefficients are never fed back into the TVP-GVAR.",
    sprintf("GPR_COLUMN=%s; +10%% shock = %.8f log units.", GPR_COLUMN, shock_log)
  ),
  stringsAsFactors = FALSE
)
write.csv(guardrail, file.path(OUT, "11_claims_guardrail.csv"), row.names = FALSE)

gate <- data.frame(
  Status = "DATA_DESIGN_AUDIT_COMPLETE",
  Source08jRunID = SOURCE_08J_RUN,
  Source08oRunID = SOURCE_08O_RUN,
  Source08oHeadSHA = SOURCE_08O_SHA,
  Panel = basename(PANEL),
  MainNetwork = MAIN_NETWORK,
  Target = "TR__de__gpr_L1",
  Observations = nrow(X),
  Predictors = ncol(X),
  FullOLS_gpr_L1 = beta_full,
  FWL_gpr_L1 = fwl_beta,
  FWL_OLS_AbsoluteGap = fwl_gap,
  HuberIRLS_gpr_L1 = beta_huber,
  PosteriorEventMedianCenter = post_center,
  PosteriorVsOLSRelativeGap = posterior_ols_rel_gap,
  VIF_gpr_L1 = vif_gpr_l1,
  R2_gprL1_on_other_predictors = r2_x_others,
  MaxLeaveOneOutRelativeChange = max_loo_rel,
  MaxTop3DeletionRelativeChange = max_top3_rel,
  HuberVsOLSRelativeGap = huber_ols_rel_gap,
  VIFThreshold = VIF_FLAG,
  LOOFlagThreshold = LOO_REL_CHANGE_FLAG,
  TopKFlagThreshold = TOPK_REL_CHANGE_FLAG,
  FlagHighVIF = flag_high_vif,
  FlagSingleQuarterInfluence = flag_single_obs,
  FlagTopFewInfluence = flag_topfew,
  FlagRobustEstimateShift = flag_robust_shift,
  FlagPosteriorTracksStaticOLS = flag_posterior_tracks_static,
  ForensicDiagnosis = diagnosis,
  ReestimationPerformed = FALSE,
  ModelSpecificationChanged = FALSE,
  WeightMatrixChanged = FALSE,
  AutomaticObservationDeletion = FALSE,
  stringsAsFactors = FALSE
)
write.csv(gate, file.path(OUT, "00_08p_gate.csv"), row.names = FALSE)

readme <- c(
  sprintf("08p TR DATA/DESIGN FORENSIC AUDIT: %s", gate$Status),
  "==========================================================",
  "",
  sprintf("Target: %s", gate$Target),
  sprintf("Panel: %s", gate$Panel),
  sprintf("Sample: %s to %s; N=%d", model_quarters[1], tail(model_quarters,1), nrow(X)),
  sprintf("Exact R/60 predictors: %d", ncol(X)),
  sprintf("Main network: %s", MAIN_NETWORK),
  "",
  sprintf("Full static OLS gpr_L1 = %.8f", beta_full),
  sprintf("FWL gpr_L1 = %.8f; reconciliation gap = %.12g", fwl_beta, fwl_gap),
  sprintf("Huber IRLS gpr_L1 = %.8f", beta_huber),
  sprintf("Accepted 08o posterior event-median center = %.8f", post_center),
  sprintf("Posterior-vs-OLS relative gap = %.6f", posterior_ols_rel_gap),
  "",
  sprintf("R2(gpr_L1 ~ other predictors) = %.6f", r2_x_others),
  sprintf("VIF(gpr_L1) = %.6f", vif_gpr_l1),
  sprintf("Residual SD(gpr_L1 | controls) / raw SD = %.6f", resid_x_sd_ratio),
  sprintf("Partial correlation = %.6f", partial_corr),
  "",
  sprintf("Max leave-one-quarter relative coefficient change = %.6f", max_loo_rel),
  sprintf("Max top-3 influence deletion relative change = %.6f", max_top3_rel),
  sprintf("Forensic diagnosis = %s", diagnosis),
  "",
  "Interpretation:",
  "- Large OLS and posterior coefficients together imply the anomaly is rooted in the same sample/design, not created only by the IRF algebra.",
  "- Large VIF or tiny residual SD ratio indicates gpr_L1 is weakly separately identified from the other regressors.",
  "- Large leave-one-out / FWL contribution concentration indicates a few quarters drive the slope.",
  "- Neither condition proves a data error. Genuine crisis observations can be influential.",
  "- No observation is automatically removed.",
  "",
  "Key outputs:",
  "- 00_08p_gate.csv",
  "- 00_design_manifest.csv",
  "- 01_raw_scale_summary.csv",
  "- 02_quarter_level_forensic_table.csv",
  "- 03_top_influential_quarters.csv",
  "- 04_leave_one_quarter_out.csv",
  "- 05_topk_influence_deletion_diagnostic.csv",
  "- 06_core_event_window_deletion_diagnostic.csv",
  "- 07_subsample_coefficient_diagnostic.csv",
  "- 08_rolling_window_coefficient_diagnostic.csv",
  "- 09_static_vs_posterior_scale_comparison.csv",
  "- 10_forensic_summary_metrics.csv",
  "- 11_claims_guardrail.csv"
)
writeLines(readme, file.path(OUT, "README_08p.txt"))

cat("08P TR DATA/DESIGN FORENSIC AUDIT: DATA_DESIGN_AUDIT_COMPLETE\n")
print(gate)
