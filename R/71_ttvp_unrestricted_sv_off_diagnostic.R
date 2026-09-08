#!/usr/bin/env Rscript

# =============================================================================
# 71_ttvp_unrestricted_sv_off_diagnostic.R
#
# Diagnostic-only unrestricted TTVP × SV numerical-isolation audit.
#
# Purpose
# -------
# Compare the accepted 10l calibration against the official replication-code
# calibration used in:
#   Crespo Cuaresma, Doppelhofer, Feldkircher & Huber (2019),
#   "Spillovers from US monetary policy: evidence from a time varying parameter
#    global vector auto-regressive model", JRSS Series A.
#
# This script DOES NOT replace the accepted 10l baseline and DOES NOT authorize
# IRFs. It is a positive-control / prior-sensitivity diagnostic.
#
# Scenarios
# ---------
# REAL_JP / REAL_SG / REAL_TR:
#   exact accepted 3-variable Trade-W design, one country, three equations.
#
# SYNTH_JP:
#   uses the exact standardized JP design matrix but generates a synthetic
#   response whose gpr_0 coefficient has known breaks:
#       pre-2008Q3:  0.0
#       2008Q3-2019Q4: -0.6
#       2020Q1 onward: +0.5
#   Only the deq equation is estimated for the synthetic positive control.
#
# Calibrations
# ------------
# current_10l:
#   B1=2, B2=1, kappa0=+1e-4, threshold high=1.0, SV mu prior sd=100.
#
# paper_code:
#   official replication-code calibration:
#   B1=3, B2=0.03, kappa0=-0.005, threshold high=1.5,
#   SV mu prior sd=10.
#   In threshtvp, negative kappa0 means coefficient-specific scaling:
#       low-regime innovation SD_j = 0.005 * OLS_SE_j.
#
# paper_code_wide:
#   same as paper_code except threshold high=3.0.
#   This is a sensitivity case, NOT claimed to be the official replication code.
#
# static vs unrestricted
# ----------------------
# static:
#   r_star_0, de_star_0, deq_star_0 are sampled in the exact static Gibbs block
#   used by 10k-fixed/10l.
#
# unrestricted:
#   no static block; all coefficients use the original pinned MCMC_tvp branch.
#   This is used only for the synthetic paper_code positive-control comparison.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("threshtvp", quietly = TRUE)) {
  stopf("Package 'threshtvp' is required.")
}

get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}
get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("%s is not numeric: %s", name, z)
  out
}

CASE <- toupper(get_env_chr("FIN3_10N_CASE", "REAL_JP"))
CALIB <- tolower(get_env_chr("FIN3_10N_CALIBRATION", "current_10l"))
MODE <- tolower(get_env_chr("FIN3_10N_MODE", "static"))
CHAIN_ID <- as.integer(get_env_num("FIN3_CHAIN_ID", 1))

BURN <- as.integer(get_env_num("FIN3_BURN", 5000))
KEEP <- as.integer(get_env_num("FIN3_KEEP", 15000))
STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 1000))
SEED_BASE <- as.integer(get_env_num("FIN3_SEED_BASE", 20267001))
SV_INNER_BURNIN <- as.integer(get_env_num("FIN3_SV_INNER_BURNIN", 4))

PANEL <- get_env_chr(
  "FIN3_FORMAL_PANEL",
  file.path(DERIVED_DIR, "panel_domestic_fin3_rate_diff.csv")
)
EVENT_FILE <- get_env_chr(
  "FIN3_EVENT_CALENDAR",
  file.path(RESULTS_DIR, "events", "00_event_calendar.csv")
)

if (!CASE %in% c("REAL_JP","REAL_SG","REAL_TR","SYNTH_JP")) {
  stopf("Unknown FIN3_10N_CASE: %s", CASE)
}
if (!CALIB %in% c("current_10l","paper_code","paper_code_wide")) {
  stopf("Unknown FIN3_10N_CALIBRATION: %s", CALIB)
}
if (!MODE %in% c("static","unrestricted")) {
  stopf("Unknown FIN3_10N_MODE: %s", MODE)
}
if (MODE == "unrestricted" && !(CASE == "SYNTH_JP" && CALIB == "paper_code")) {
  stopf("Unrestricted mode is allowed only for SYNTH_JP + paper_code.")
}
if (CHAIN_ID < 1L || CHAIN_ID > 4L) stopf("CHAIN_ID must be 1..4.")
if (BURN < 1000L || KEEP < 5000L || STORED < 500L || STORED >= KEEP) {
  stopf("Invalid diagnostic MCMC iteration contract.")
}
if (SV_INNER_BURNIN < 1L || SV_INNER_BURNIN > 20L) {
  stopf("SV inner burn-in must be 1..20.")
}

calibration <- switch(
  CALIB,
  current_10l = list(
    B1 = 2,
    B2 = 1,
    kappa0 = 1e-4,
    thr_low = 0.1,
    thr_high = 1.0,
    sv_mu_sd = 100,
    label = "ACCEPTED_10L_CALIBRATION"
  ),
  paper_code = list(
    B1 = 3,
    B2 = 0.03,
    kappa0 = -0.005,
    thr_low = 0.1,
    thr_high = 1.5,
    sv_mu_sd = 10,
    label = "JRSSA_2019_OFFICIAL_REPLICATION_CODE_CALIBRATION"
  ),
  paper_code_wide = list(
    B1 = 3,
    B2 = 0.03,
    kappa0 = -0.005,
    thr_low = 0.1,
    thr_high = 3.0,
    sv_mu_sd = 10,
    label = "JRSSA_2019_CODE_CALIBRATION_WITH_WIDER_THRESHOLD_SENSITIVITY"
  )
)

case_country <- switch(
  CASE,
  REAL_JP = "JP",
  REAL_SG = "SG",
  REAL_TR = "TR",
  SYNTH_JP = "JP"
)

scenario_id <- paste(CASE, CALIB, MODE, sep = "__")
OUT <- file.path(
  RESULTS_DIR, "10n_ttvp_recovery_parts",
  scenario_id, paste0("chain_", CHAIN_ID)
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Exact accepted p=1/q=1 design construction
# -----------------------------------------------------------------------------

if (!file.exists(PANEL)) stopf("Missing formal panel: %s", PANEL)
d <- read.csv(PANEL, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("Quarter","Country",VARS,"gpr","brent")
if (!all(need %in% names(d))) {
  stopf("Panel missing: %s", paste(setdiff(need,names(d)),collapse=", "))
}

d$Country <- toupper(trimws(as.character(d$Country)))
d$Quarter <- toupper(trimws(as.character(d$Quarter)))

quarters <- unique(d$Quarter)
if (!all(COUNTRIES %in% unique(d$Country))) stopf("Missing sample economies.")
if (any(table(d$Country) != length(quarters))) stopf("Panel is not balanced.")

N <- length(COUNTRIES)
K <- length(VARS)
Tn <- length(quarters)

Xglobal <- array(
  NA_real_,
  c(Tn,N,K),
  dimnames = list(quarters,COUNTRIES,VARS)
)
for (i in seq_along(COUNTRIES)) {
  z <- d[d$Country == COUNTRIES[i],,drop=FALSE]
  z <- z[match(quarters,z$Quarter),,drop=FALSE]
  Xglobal[,i,] <- as.matrix(z[,VARS,drop=FALSE])
}

base <- d[d$Country == COUNTRIES[1],,drop=FALSE]
base <- base[match(quarters,base$Quarter),,drop=FALSE]
gpr <- num(base$gpr)
brent <- num(base$brent)

if (!all(is.finite(Xglobal)) || !all(is.finite(gpr)) || !all(is.finite(brent))) {
  stopf("Non-finite accepted model input.")
}

if (!identical(MAIN_NETWORK, "trade_2000_2012")) {
  stopf("10n requires FIN3_NETWORK=trade_2000_2012.")
}
W <- read_weight_matrix(network_path(MAIN_NETWORK))
country_i <- match(case_country, COUNTRIES)

STAR <- array(
  NA_real_,
  c(Tn,N,K),
  dimnames = list(quarters,COUNTRIES,VARS)
)
for (i in seq_len(N)) {
  for (v in seq_len(K)) {
    STAR[,i,v] <- as.numeric(Xglobal[,,v] %*% W[i,])
  }
}

lag1 <- function(x) c(NA_real_,x[-length(x)])

Yraw <- Xglobal[,country_i,,drop=FALSE][,1,]
Zraw <- STAR[,country_i,,drop=FALSE][,1,]
rows <- 2:nrow(Yraw)

D <- data.frame(const=rep(1,length(rows)),check.names=FALSE)
for (v in seq_len(K)) D[[paste0(VARS[v],"_L1")]] <- lag1(Yraw[,v])[rows]
for (v in seq_len(K)) D[[paste0(VARS[v],"_star_0")]] <- Zraw[rows,v]
for (v in seq_len(K)) D[[paste0(VARS[v],"_star_L1")]] <- lag1(Zraw[,v])[rows]
D$gpr_0 <- gpr[rows]
D$gpr_L1 <- lag1(gpr)[rows]
D$brent_0 <- brent[rows]
D$brent_L1 <- lag1(brent)[rows]

Yreal <- Yraw[rows,,drop=FALSE]
model_quarters <- quarters[rows]

ok <- complete.cases(D) & complete.cases(Yreal)
D <- D[ok,,drop=FALSE]
Yreal <- Yreal[ok,,drop=FALSE]
model_quarters <- model_quarters[ok]

TERMS <- colnames(D)
EXPECTED_TERMS <- c(
  "const",
  paste0(VARS,"_L1"),
  paste0(VARS,"_star_0"),
  paste0(VARS,"_star_L1"),
  "gpr_0","gpr_L1","brent_0","brent_L1"
)
if (!identical(TERMS, EXPECTED_TERMS)) {
  stopf("Unexpected design term ordering.")
}
if (nrow(D) <= ncol(D) + 10L) stopf("Too few observations.")

# Standardize predictors exactly as accepted R/60.
Xm_raw <- as.matrix(D)
Xm <- Xm_raw
x_center <- rep(0,ncol(Xm)); names(x_center) <- TERMS
x_scale <- rep(1,ncol(Xm)); names(x_scale) <- TERMS
for (j in seq_len(ncol(Xm))) {
  if (TERMS[j] == "const") next
  x_center[j] <- mean(Xm[,j])
  x_scale[j] <- stats::sd(Xm[,j])
  if (!is.finite(x_scale[j]) || x_scale[j] < 1e-10) {
    stopf("Near-constant predictor: %s", TERMS[j])
  }
  Xm[,j] <- (Xm[,j] - x_center[j]) / x_scale[j]
}
Xm[,"const"] <- 1

# Static block contract used in 10k-fixed / 10l.
STATIC_TERMS <- paste0(VARS,"_star_0")
STATIC_IDX <- match(STATIC_TERMS,TERMS)
if (any(is.na(STATIC_IDX))) stopf("Static block terms not found.")
static_idx_call <- if (MODE == "static") STATIC_IDX else integer(0)

# Paper-code negative kappa0 is coefficient-specific because the sampler
# transforms it to -kappa0 * OLS_SE for each coefficient.
low_sd_from_kappa <- function(kappa0, sd_ols) {
  if (kappa0 < 0) (-kappa0) * sd_ols else rep(kappa0,length(sd_ols))
}

THIN_FRACTION <- (STORED + 1L) / KEEP
if (!(THIN_FRACTION > 0 && THIN_FRACTION <= 1)) stopf("Invalid thinning.")

# Fixed synthetic dataset seed: identical data across calibrations/chains.
make_synthetic_y <- function(Xm, quarters, terms) {
  set.seed(20261909)
  TT <- nrow(Xm)
  P <- ncol(Xm)

  beta <- matrix(0,TT,P,dimnames=list(quarters,terms))
  beta[,"const"] <- 0.05
  beta[,c("r_L1","de_L1","deq_L1")] <- matrix(
    rep(c(0.15,0.10,0.15),each=TT), nrow=TT
  )
  beta[,c("r_star_0","de_star_0","deq_star_0")] <- matrix(
    rep(c(0.04,0.03,0.05),each=TT), nrow=TT
  )
  beta[,c("r_star_L1","de_star_L1","deq_star_L1")] <- matrix(
    rep(c(0.03,0.02,0.03),each=TT), nrow=TT
  )
  beta[,"gpr_L1"] <- 0
  beta[,"brent_0"] <- 0.04
  beta[,"brent_L1"] <- -0.02

  qid <- quarter_id(quarters)
  b1 <- quarter_id("2008Q3")
  b2 <- quarter_id("2020Q1")
  beta[, "gpr_0"] <- ifelse(
    qid < b1, 0.0,
    ifelse(qid < b2, -0.60, 0.50)
  )

  # Mild stochastic volatility so the positive control still exercises SV.
  mu_h <- -2.0
  phi_h <- 0.95
  sig_h <- 0.12
  h <- numeric(TT)
  h[1] <- mu_h
  if (TT > 1L) {
    for (tt in 2:TT) {
      h[tt] <- mu_h + phi_h*(h[tt-1]-mu_h) + rnorm(1,0,sig_h)
    }
  }
  eps <- exp(h/2) * rnorm(TT)

  signal <- rowSums(Xm * beta)
  y <- signal + eps

  list(
    y_std = as.numeric(y),
    beta_true = beta,
    h_true = h,
    break_quarters = c("2008Q3","2020Q1")
  )
}

synthetic <- NULL
if (CASE == "SYNTH_JP") {
  synthetic <- make_synthetic_y(Xm,model_quarters,TERMS)
}

equations <- if (CASE == "SYNTH_JP") "deq" else VARS

# MCMC seed changes only sampling, not the synthetic dataset.
SEED <- SEED_BASE +
  100000L*match(CASE,c("REAL_JP","REAL_SG","REAL_TR","SYNTH_JP")) +
  10000L*match(CALIB,c("current_10l","paper_code","paper_code_wide")) +
  1000L*(MODE == "unrestricted") +
  CHAIN_ID
set.seed(SEED)



# =============================================================================
# 10o hard lock
# =============================================================================
# This diagnostic changes ONE switch relative to the failed unrestricted
# synthetic paper-code path: sv_on = FALSE.
#
# It is NOT a candidate formal specification, does NOT replace 10l/10n, and
# does NOT authorize IRFs.
# =============================================================================

if (!identical(CASE, "SYNTH_JP")) {
  stopf("10o hard lock requires FIN3_10N_CASE=SYNTH_JP.")
}
if (!identical(CALIB, "paper_code")) {
  stopf("10o hard lock requires FIN3_10N_CALIBRATION=paper_code.")
}
if (!identical(MODE, "unrestricted")) {
  stopf("10o hard lock requires FIN3_10N_MODE=unrestricted.")
}
if (!identical(equations, "deq")) {
  stopf("10o hard lock expected the synthetic deq equation only.")
}

OUT10O <- file.path(
  RESULTS_DIR,
  "10o_unrestricted_sv_off",
  paste0("chain_", CHAIN_ID)
)
dir.create(OUT10O, recursive=TRUE, showWarnings=FALSE)

write_status <- function(
  sampler_completed,
  error_message="",
  posterior_schema_available=FALSE,
  coefficient_path_available=FALSE,
  indicator_path_available=FALSE
) {
  z <- data.frame(
    TechnicalStatus="PASS",
    Diagnostic="10o_unrestricted_synthetic_paper_code_sv_off",
    Case=CASE,
    Calibration=CALIB,
    SamplerMode=MODE,
    SV=FALSE,
    Chain=CHAIN_ID,
    Seed=SEED,
    Burn=BURN,
    PostBurnInternal=KEEP,
    RequestedStored=STORED,
    B1=calibration$B1,
    B2=calibration$B2,
    Kappa0=calibration$kappa0,
    ThresholdLow=calibration$thr_low,
    ThresholdHigh=calibration$thr_high,
    SamplerCompleted=as.logical(sampler_completed),
    PosteriorSchemaAvailable=as.logical(posterior_schema_available),
    CoefficientPathAvailable=as.logical(coefficient_path_available),
    IndicatorPathAvailable=as.logical(indicator_path_available),
    ErrorMessage=as.character(error_message),
    FormalModelChanged=FALSE,
    IRFAuthorized=FALSE,
    stringsAsFactors=FALSE
  )
  write.csv(z, file.path(OUT10O, "00_status.csv"), row.names=FALSE)
  invisible(z)
}

# Write a pre-MCMC status file so that an estimator error still leaves a
# diagnostic artifact rather than only an Actions exit code.
write_status(
  sampler_completed=FALSE,
  error_message="MCMC_NOT_YET_RUN"
)

eq_name <- "deq"

# Synthetic positive control is already on the standardized response scale.
y_std <- synthetic$y_std
if (length(y_std) != nrow(Xm) || any(!is.finite(y_std))) {
  stopf("10o synthetic response contract failed.")
}

msg(
  paste0(
    "10o: SYNTH_JP paper_code unrestricted; chain=%d; ",
    "SV=FALSE; B1=%g B2=%g kappa0=%g threshold=[%g,%g]"
  ),
  CHAIN_ID,
  calibration$B1,
  calibration$B2,
  calibration$kappa0,
  calibration$thr_low,
  calibration$thr_high
)

# ---------------------------------------------------------------------------
# The only intended experimental change from the failed unrestricted path is:
#     sv_on = FALSE
#
# Keep all other estimation arguments identical to the current 10n R69 call.
# ---------------------------------------------------------------------------

fit <- tryCatch(
  threshtvp::estimate_tvp(
    Y = matrix(y_std,ncol=1),
    X = Xm,
    save = KEEP,
    burn = BURN,
    priorbtheta = list(
      B_1 = calibration$B1,
      B_2 = calibration$B2,
      kappa0 = calibration$kappa0
    ),
    priorb0 = list(
      a_tau=0.1,
      c_tau=0.01,
      d_tau=0.01
    ),
    priorsig = c(0.01,0.01),
    priorphi = c(2,2),
    priormu = c(0,calibration$sv_mu_sd),
    h0prior = "stationary",
    grid.length = 150,
    thrsh.pct = calibration$thr_low,
    thrsh.pct.high = calibration$thr_high,
    sv_on = FALSE,
    TVS = TRUE,
    cons.mod = FALSE,
    static_idx = integer(0),
    sv_inner_burnin = SV_INNER_BURNIN,
    p = 1,
    thin = THIN_FRACTION,
    CPU = 1,
    approx = FALSE,
    sim.kappa0 = FALSE
  ),
  error=function(e) e
)

if (inherits(fit, "error")) {
  err <- conditionMessage(fit)

  write_status(
    sampler_completed=FALSE,
    error_message=err
  )

  writeLines(
    c(
      "10o unrestricted synthetic paper_code SV=FALSE",
      paste0("chain=", CHAIN_ID),
      paste0("seed=", SEED),
      paste0("error=", err),
      "",
      "Interpretation rule:",
      "If this error still contains chol.default(SigHigh), the numerical",
      "failure is not specific to the stochastic-volatility layer."
    ),
    file.path(OUT10O, "00_error.txt")
  )

  cat("\n===== 10o SAMPLER DID NOT COMPLETE =====\n")
  cat("chain:", CHAIN_ID, "\n")
  cat("error:", err, "\n")
  cat(
    "This is a scientific/numerical diagnostic outcome; ",
    "the script exits successfully after recording it.\n",
    sep=""
  )

  quit(save="no", status=0L)
}

post <- fit$posterior
if (is.null(post) || !is.list(post)) {
  write_status(
    sampler_completed=TRUE,
    error_message="MCMC_COMPLETED_BUT_POSTERIOR_LIST_MISSING"
  )
  quit(save="no", status=0L)
}

# ---------------------------------------------------------------------------
# Posterior schema audit
# ---------------------------------------------------------------------------

schema_rows <- lapply(names(post), function(nm) {
  x <- post[[nm]]
  dx <- dim(x)
  data.frame(
    Object=nm,
    Class=paste(class(x),collapse="|"),
    HasDim=!is.null(dx),
    Dimensions=if (is.null(dx)) "" else paste(dx,collapse="x"),
    Length=length(x),
    stringsAsFactors=FALSE
  )
})
schema <- if (length(schema_rows)) do.call(rbind,schema_rows) else data.frame(
  Object=character(),Class=character(),HasDim=logical(),
  Dimensions=character(),Length=integer()
)
write.csv(
  schema,
  file.path(OUT10O, "01_posterior_schema.csv"),
  row.names=FALSE
)

# Robustly standardize a 3-D posterior array to [draw, time, term].
orient_draw_time_term <- function(x, object_name) {
  dx <- dim(x)
  if (is.null(dx) || length(dx) != 3L) {
    stopf("%s is not a 3-D posterior array.", object_name)
  }

  term_axis <- which(dx == length(TERMS))
  time_axis <- which(dx == length(model_quarters))

  if (length(term_axis) != 1L || length(time_axis) != 1L) {
    stopf(
      "%s dimensions cannot be uniquely mapped to term/time axes: %s",
      object_name,paste(dx,collapse="x")
    )
  }

  draw_axis <- setdiff(seq_len(3L),c(term_axis,time_axis))
  if (length(draw_axis) != 1L) {
    stopf("%s draw axis is ambiguous.",object_name)
  }

  aperm(x,c(draw_axis,time_axis,term_axis))
}

coef_available <- FALSE
indicator_available <- FALSE
path <- NULL

# ---------------------------------------------------------------------------
# Coefficient-path recovery from A
# ---------------------------------------------------------------------------

if (!is.null(post$A)) {
  A3 <- tryCatch(
    orient_draw_time_term(post$A,"A"),
    error=function(e) e
  )

  if (!inherits(A3,"error")) {
    jg <- match("gpr_0",TERMS)
    if (!is.na(jg)) {
      gdraw <- A3[,,jg,drop=FALSE][,,1]
      if (is.vector(gdraw)) {
        gdraw <- matrix(gdraw,nrow=dim(A3)[1])
      }

      post_mean <- colMeans(gdraw,na.rm=TRUE)
      post_sd <- apply(gdraw,2,stats::sd,na.rm=TRUE)

      truth <- synthetic$beta_true[,"gpr_0"]

      path <- data.frame(
        Chain=CHAIN_ID,
        Quarter=model_quarters,
        TrueGPR0BetaStd=as.numeric(truth),
        PosteriorMeanGPR0BetaStd=as.numeric(post_mean),
        PosteriorSDGPR0BetaStd=as.numeric(post_sd),
        stringsAsFactors=FALSE
      )
      coef_available <- TRUE
    }
  }
}

# ---------------------------------------------------------------------------
# Threshold-indicator recovery, if the unrestricted posterior exposes it.
# We do not assume a field name; try known candidates and record schema either
# way. Missing indicators do not invalidate the SV-off numerical diagnostic.
# ---------------------------------------------------------------------------

D_candidates <- c("D_dyn","D","d")
D_name <- D_candidates[D_candidates %in% names(post)]
if (length(D_name)) {
  D_name <- D_name[1]
  D3 <- tryCatch(
    orient_draw_time_term(post[[D_name]],D_name),
    error=function(e) e
  )

  if (!inherits(D3,"error")) {
    jg <- match("gpr_0",TERMS)
    dg <- D3[,,jg,drop=FALSE][,,1]
    if (is.vector(dg)) {
      dg <- matrix(dg,nrow=dim(D3)[1])
    }
    tv_prob <- colMeans(dg,na.rm=TRUE)

    if (is.null(path)) {
      path <- data.frame(
        Chain=CHAIN_ID,
        Quarter=model_quarters,
        TrueGPR0BetaStd=as.numeric(synthetic$beta_true[,"gpr_0"]),
        stringsAsFactors=FALSE
      )
    }
    path$TVProbabilityGPR0 <- as.numeric(tv_prob)
    indicator_available <- TRUE
  }
}

if (!is.null(path)) {
  write.csv(
    path,
    file.path(OUT10O, "02_gpr0_recovery_path.csv"),
    row.names=FALSE
  )
}

# ---------------------------------------------------------------------------
# Chain-level recovery summary
# ---------------------------------------------------------------------------

qid10o <- function(x) quarter_id(as.character(x))
q <- qid10o(model_quarters)

if (coef_available) {
  pre <- q < quarter_id("2008Q3")
  mid <- q >= quarter_id("2008Q3") & q < quarter_id("2020Q1")
  postseg <- q >= quarter_id("2020Q1")

  pre_beta <- mean(path$PosteriorMeanGPR0BetaStd[pre],na.rm=TRUE)
  mid_beta <- mean(path$PosteriorMeanGPR0BetaStd[mid],na.rm=TRUE)
  post_beta <- mean(path$PosteriorMeanGPR0BetaStd[postseg],na.rm=TRUE)

  rmse <- sqrt(mean(
    (path$PosteriorMeanGPR0BetaStd-path$TrueGPR0BetaStd)^2,
    na.rm=TRUE
  ))

  coefficient_order_pass <- (
    is.finite(pre_beta) &&
    is.finite(mid_beta) &&
    is.finite(post_beta) &&
    mid_beta < pre_beta - 0.15 &&
    post_beta > mid_beta + 0.25 &&
    post_beta > pre_beta + 0.10
  )
} else {
  pre_beta <- mid_beta <- post_beta <- rmse <- NA_real_
  coefficient_order_pass <- FALSE
}

if (indicator_available) {
  w1 <- abs(q-quarter_id("2008Q3")) <= 1L
  w2 <- abs(q-quarter_id("2020Q1")) <= 1L
  nonbreak <- !(w1|w2)

  break1 <- max(path$TVProbabilityGPR0[w1],na.rm=TRUE)
  break2 <- max(path$TVProbabilityGPR0[w2],na.rm=TRUE)
  nonbreak_mean <- mean(path$TVProbabilityGPR0[nonbreak],na.rm=TRUE)
} else {
  break1 <- break2 <- nonbreak_mean <- NA_real_
}

recovery <- data.frame(
  Chain=CHAIN_ID,
  SamplerCompleted=TRUE,
  SV=FALSE,
  CoefficientPathAvailable=coef_available,
  IndicatorPathAvailable=indicator_available,
  PrePosteriorMeanBeta=pre_beta,
  MidPosteriorMeanBeta=mid_beta,
  PostPosteriorMeanBeta=post_beta,
  CoefficientPathRMSE=rmse,
  CoefficientOrderPass=coefficient_order_pass,
  Break2008MaxTVProb=break1,
  Break2020MaxTVProb=break2,
  NonBreakMeanTVProb=nonbreak_mean,
  stringsAsFactors=FALSE
)
write.csv(
  recovery,
  file.path(OUT10O, "03_recovery_summary.csv"),
  row.names=FALSE
)

# Keep only the minimal posterior slices required to inspect the diagnostic.
mini <- list(
  meta=list(
    diagnostic="10o_unrestricted_synthetic_paper_code_sv_off",
    chain=CHAIN_ID,
    seed=SEED,
    sv=FALSE,
    calibration=calibration,
    terms=TERMS,
    quarters=model_quarters
  ),
  schema=schema,
  gpr0_recovery_path=path,
  recovery_summary=recovery
)
saveRDS(
  mini,
  file.path(OUT10O, "10o_minimal_result.rds"),
  compress="xz"
)

write_status(
  sampler_completed=TRUE,
  error_message="",
  posterior_schema_available=TRUE,
  coefficient_path_available=coef_available,
  indicator_path_available=indicator_available
)

cat("\n===== 10o SV-OFF DIAGNOSTIC COMPLETE =====\n")
print(read.csv(file.path(OUT10O,"00_status.csv"),check.names=FALSE))
print(recovery)
