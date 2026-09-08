#!/usr/bin/env Rscript

# =============================================================================
# 69_ttvp_recovery_chain.R
#
# Diagnostic-only TTVP recovery / calibration audit.
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
# 10r note: unrestricted mode assumes the integerized MCMC_tvp storage-grid
# implementation validated by 10q; this is an implementation-only bug fix.
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

equation_objects <- list()
summary_rows <- list()
tv_rows <- list()
truth_rows <- list()
monitor <- data.frame(.draw=seq_len(STORED))
sr <- 0L
tr <- 0L
rr <- 0L

for (eq_name in equations) {
  eq <- match(eq_name,VARS)

  if (CASE == "SYNTH_JP") {
    y_std <- synthetic$y_std
    y_center <- 0
    y_scale <- 1
  } else {
    y_raw <- as.numeric(Yreal[,eq])
    y_center <- mean(y_raw)
    y_scale <- stats::sd(y_raw)
    if (!is.finite(y_scale) || y_scale < 1e-10) {
      stopf("Near-constant response %s/%s",case_country,eq_name)
    }
    y_std <- (y_raw-y_center)/y_scale
  }

  msg(
    "10n: scenario=%s equation=%s chain=%d B1=%g B2=%g kappa0=%g thr_high=%g static=%s",
    scenario_id,eq_name,CHAIN_ID,
    calibration$B1,calibration$B2,calibration$kappa0,
    calibration$thr_high,MODE
  )

  fit <- threshtvp::estimate_tvp(
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
    sv_on = TRUE,
    TVS = TRUE,
    cons.mod = FALSE,
    static_idx = static_idx_call,
    sv_inner_burnin = SV_INNER_BURNIN,
    p = 1,
    thin = THIN_FRACTION,
    CPU = 1,
    approx = FALSE,
    sim.kappa0 = FALSE
  )

  post <- fit$posterior
  # and ending after:
  #     svp <- post$svparms[keep_draw,,drop=FALSE]
  #
  # IMPORTANT:
  # - MODE == "static": raw draw 1 is the static-block burn-boundary/storage row.
  #   Discard it, then retain exactly STORED strictly post-burn draws.
  # - MODE == "unrestricted": fit$posterior is already post-burn storage.
  #   Do NOT require or discard an extra boundary row.
  # - This changes storage bookkeeping only. It does NOT change the statistical
  #   model, priors, burn-in, internal post-burn iterations, TVS, SV, or gate.
  # =============================================================================

  post <- fit$posterior

  if (is.null(post$A) || length(dim(post$A)) != 3L) {
    stopf("Unexpected A structure.")
  }

  nstored_raw <- dim(post$A)[1]

  # Return the first dimension / row count for posterior objects that are expected
  # to share the same draw index.
  posterior_draw_count <- function(x, name) {
    if (is.null(x)) stopf("Missing posterior object: %s", name)

    dx <- dim(x)
    if (is.null(dx)) {
      stopf("Posterior object %s has no draw dimension.", name)
    }

    n <- dx[1]
    if (!is.finite(n) || n < 1L) {
      stopf("Invalid raw draw count for %s: %s", name, as.character(n))
    }
    as.integer(n)
  }

  raw_counts <- c(
    A          = posterior_draw_count(post$A,          "A"),
    D_dyn      = posterior_draw_count(post$D_dyn,      "D_dyn"),
    thresholds = posterior_draw_count(post$thresholds, "thresholds"),
    omega      = posterior_draw_count(post$omega,      "omega"),
    V0         = posterior_draw_count(post$V0,         "V0"),
    H          = posterior_draw_count(post$H,          "H"),
    svparms    = posterior_draw_count(post$svparms,    "svparms")
  )

  if (length(unique(raw_counts)) != 1L) {
    stopf(
      "Raw posterior draw-count mismatch: %s",
      paste(sprintf("%s=%d", names(raw_counts), raw_counts), collapse = ", ")
    )
  }

  if (raw_counts[["A"]] != nstored_raw) {
    stopf("Internal raw-draw bookkeeping mismatch for A.")
  }

  # -----------------------------------------------------------------------------
  # MODE-SPECIFIC storage contract
  # -----------------------------------------------------------------------------
  # static:
  #   tested static-block extension can expose one leading burn-boundary/storage
  #   row.  That row is not a strictly post-burn posterior draw and must be
  #   discarded.
  #
  # unrestricted:
  #   original pinned threshtvp branch returns posterior storage directly.  If it
  #   returns exactly STORED rows, all STORED rows are valid posterior draws.
  #   There is no additional boundary row to remove.
  # -----------------------------------------------------------------------------

  if (MODE == "static") {
    if (nstored_raw < STORED + 1L) {
      stopf(
        paste0(
          "Too few raw draws in static mode: need at least %d ",
          "(1 boundary + %d post-burn), found %d."
        ),
        STORED + 1L, STORED, nstored_raw
      )
    }

    candidates <- seq.int(2L, nstored_raw)
    boundary_dropped <- 1L

  } else if (MODE == "unrestricted") {
    # 10r storage-fixed original MCMC_tvp contract:
    #   nstore = round(thin * nsave) = STORED + 1
    # and thin_out starts exactly at irep = burn.  Therefore raw row 1 is the
    # burn-boundary storage row and must be dropped before formal diagnostics.
    expected_raw <- STORED + 1L
    if (nstored_raw != expected_raw) {
      stopf(
        paste0(
          "Storage-fixed unrestricted expected exactly %d raw draws ",
          "(1 burn-boundary + %d strictly post-burn), found %d. ",
          "The integer storage patch may be missing or the sampler contract changed."
        ),
        expected_raw, STORED, nstored_raw
      )
    }

    candidates <- seq.int(2L, nstored_raw)
    boundary_dropped <- 1L

  } else {
    stopf("Unknown sampler MODE in draw bookkeeping: %s", MODE)
  }

  n_candidates <- length(candidates)
  if (n_candidates < STORED) {
    stopf(
      "Too few eligible posterior candidates in %s mode: need %d, found %d.",
      MODE, STORED, n_candidates
    )
  }

  # Select exactly STORED draws.  If the sampler exposes more eligible storage
  # rows than requested (e.g. due to floating-point construction of the thinning
  # grid), choose a deterministic approximately-even subset.  This does not alter
  # the MCMC target.
  if (n_candidates == STORED) {
    select_pos <- seq_len(STORED)
  } else {
    select_pos <- as.integer(round(
      seq(1, n_candidates, length.out = STORED)
    ))
  }

  if (length(select_pos) != STORED ||
      length(unique(select_pos)) != STORED ||
      min(select_pos) < 1L ||
      max(select_pos) > n_candidates) {
    stopf(
      paste0(
        "Deterministic posterior selection failed in %s mode: ",
        "candidates=%d selected=%d unique=%d."
      ),
      MODE, n_candidates, length(select_pos), length(unique(select_pos))
    )
  }

  keep_draw <- candidates[select_pos]

  if (length(keep_draw) != STORED || any(diff(keep_draw) <= 0L)) {
    stopf("Posterior draw contract failed after selection in %s mode.", MODE)
  }

  if (MODE %in% c("static","unrestricted") && any(keep_draw <= 1L)) {
    stopf("%s-mode posterior selection retained the burn-boundary draw.", MODE)
  }

  msg(
    paste0(
      "10n draw bookkeeping: mode=%s; raw=%d; boundary_dropped=%d; ",
      "eligible_candidates=%d; selected=%d; first_raw_index=%d; ",
      "last_raw_index=%d"
    ),
    MODE,
    nstored_raw,
    boundary_dropped,
    n_candidates,
    length(keep_draw),
    min(keep_draw),
    max(keep_draw)
  )

  Astd       <- post$A[keep_draw,,,drop=FALSE]
  Ddyn       <- post$D_dyn[keep_draw,,,drop=FALSE]
  thresholds <- post$thresholds[keep_draw,,drop=FALSE]
  omega      <- post$omega[keep_draw,,drop=FALSE]
  V0         <- post$V0[keep_draw,,drop=FALSE]
  H          <- post$H[keep_draw,,drop=FALSE]
  svp        <- post$svparms[keep_draw,,drop=FALSE]

  if (dim(Astd)[1] != STORED ||
      dim(Ddyn)[1] != STORED ||
      nrow(thresholds) != STORED ||
      nrow(omega) != STORED ||
      nrow(V0) != STORED ||
      nrow(H) != STORED ||
      nrow(svp) != STORED) {
    stopf("Selected posterior objects do not all contain exactly STORED draws.")
  }

  if (dim(Astd)[2] != length(model_quarters) ||
      dim(Astd)[3] != length(TERMS)) {
    stopf("A dimension mismatch.")
  }
  if (dim(Ddyn)[2] != length(TERMS) ||
      dim(Ddyn)[3] != length(model_quarters)) {
    stopf("D_dyn dimension mismatch.")
  }
  if (ncol(thresholds) != length(TERMS) ||
      ncol(omega) != length(TERMS) ||
      ncol(V0) != length(TERMS)) {
    stopf("Hyperparameter dimension mismatch.")
  }

  term_names <- dimnames(post$A)[[3]]
  if (!is.null(term_names) && !identical(term_names,TERMS)) {
    stopf("Coefficient order mismatch.")
  }

  tv_prob <- t(apply(Ddyn,c(2,3),mean,na.rm=TRUE))
  colnames(tv_prob) <- TERMS
  rownames(tv_prob) <- model_quarters

  beta_mean <- apply(Astd,c(2,3),mean,na.rm=TRUE)
  beta_sd <- apply(Astd,c(2,3),stats::sd,na.rm=TRUE)
  colnames(beta_mean) <- colnames(beta_sd) <- TERMS
  rownames(beta_mean) <- rownames(beta_sd) <- model_quarters

  abs_dbeta <- array(NA_real_,dim=c(STORED,length(model_quarters),length(TERMS)))
  abs_dbeta[,1,] <- 0
  if (length(model_quarters) > 1L) {
    abs_dbeta[,2:length(model_quarters),] <- abs(
      Astd[,2:length(model_quarters),,drop=FALSE] -
      Astd[,1:(length(model_quarters)-1L),,drop=FALSE]
    )
  }

  sd_ols <- as.numeric(post$sd.OLS)
  if (length(sd_ols) != length(TERMS)) stopf("sd.OLS length mismatch.")
  low_sd <- low_sd_from_kappa(calibration$kappa0,sd_ols)

  dynamic_terms <- if (MODE == "static") setdiff(TERMS,STATIC_TERMS) else TERMS

  for (jj in seq_along(TERMS)) {
    nm <- TERMS[jj]
    threshold_med <- stats::median(thresholds[,jj],na.rm=TRUE)
    high_sd_med <- stats::median(sqrt(pmax(omega[,jj],0)),na.rm=TRUE)
    db <- as.numeric(abs_dbeta[,,jj])
    db <- db[is.finite(db)]

    sr <- sr + 1L
    summary_rows[[sr]] <- data.frame(
      Scenario=scenario_id,
      Case=CASE,
      Country=case_country,
      Calibration=CALIB,
      SamplerMode=MODE,
      Chain=CHAIN_ID,
      Equation=eq_name,
      Term=nm,
      IsStaticByDesign=(MODE=="static" && nm %in% STATIC_TERMS),
      MeanTVProbability=mean(tv_prob[,jj],na.rm=TRUE),
      MaxTVProbability=max(tv_prob[,jj],na.rm=TRUE),
      ShareQuartersTVProbGT10=mean(tv_prob[,jj]>.10,na.rm=TRUE),
      ShareQuartersTVProbGT50=mean(tv_prob[,jj]>.50,na.rm=TRUE),
      MedianThreshold=threshold_med,
      MedianHighInnovationSD=high_sd_med,
      LowInnovationSD=low_sd[jj],
      ThresholdToLowSD=ifelse(
        is.finite(low_sd[jj]) && low_sd[jj]>0,
        threshold_med/low_sd[jj],NA_real_
      ),
      MedianAbsDeltaBeta=stats::median(db,na.rm=TRUE),
      P95AbsDeltaBeta=unname(stats::quantile(db,.95,na.rm=TRUE)),
      stringsAsFactors=FALSE
    )

    for (tt in seq_along(model_quarters)) {
      tr <- tr + 1L
      tv_rows[[tr]] <- data.frame(
        Scenario=scenario_id,
        Case=CASE,
        Country=case_country,
        Calibration=CALIB,
        SamplerMode=MODE,
        Chain=CHAIN_ID,
        Equation=eq_name,
        Quarter=model_quarters[tt],
        Term=nm,
        TVProbability=tv_prob[tt,jj],
        PosteriorMeanBetaStd=beta_mean[tt,jj],
        PosteriorSDBetaStd=beta_sd[tt,jj],
        stringsAsFactors=FALSE
      )
    }
  }

  if (CASE == "SYNTH_JP") {
    jj <- match("gpr_0",TERMS)
    for (tt in seq_along(model_quarters)) {
      rr <- rr + 1L
      truth_rows[[rr]] <- data.frame(
        Scenario=scenario_id,
        Chain=CHAIN_ID,
        Quarter=model_quarters[tt],
        TrueGPR0BetaStd=synthetic$beta_true[tt,jj],
        PosteriorMeanGPR0BetaStd=beta_mean[tt,jj],
        TVProbabilityGPR0=tv_prob[tt,jj],
        stringsAsFactors=FALSE
      )
    }
  }

  # Monitor convergence of SV, scales, thresholds and key coefficients.
  if (ncol(svp) >= 3L) {
    monitor[[paste0(eq_name,"__sv_mu")]] <- svp[,1]
    monitor[[paste0(eq_name,"__sv_phi")]] <- svp[,2]
    monitor[[paste0(eq_name,"__sv_sigma")]] <- svp[,3]
  }

  for (jj in seq_along(TERMS)) {
    nm <- TERMS[jj]
    monitor[[paste0(eq_name,"__omega__",nm)]] <- omega[,jj]
  }

  key_terms <- unique(c(
    paste0(eq_name,"_L1"),
    "gpr_0","gpr_L1",
    "r_star_L1","de_star_L1","deq_star_L1"
  ))
  key_terms <- key_terms[key_terms %in% TERMS]

  for (nm in key_terms) {
    jj <- match(nm,TERMS)
    monitor[[paste0(eq_name,"__threshold__",nm)]] <- thresholds[,jj]
    monitor[[paste0(eq_name,"__V0__",nm)]] <- V0[,jj]

    for (qq in intersect(c("2008Q3","2020Q1","2024Q2"),model_quarters)) {
      tt <- match(qq,model_quarters)
      monitor[[paste0(eq_name,"__coef__",qq,"__",nm)]] <- Astd[,tt,jj]
    }
  }

  equation_objects[[eq_name]] <- list(
    tv_probability=tv_prob,
    beta_mean_std=beta_mean,
    beta_sd_std=beta_sd,
    low_innovation_sd=setNames(low_sd,TERMS),
    mean_threshold=setNames(colMeans(thresholds,na.rm=TRUE),TERMS),
    mean_omega=setNames(colMeans(omega,na.rm=TRUE),TERMS)
  )
}

monitor$.draw <- NULL

meta <- list(
  scenario=scenario_id,
  case=CASE,
  country=case_country,
  calibration=CALIB,
  calibration_label=calibration$label,
  sampler_mode=MODE,
  chain=CHAIN_ID,
  seed=SEED,
  burn=BURN,
  keep=KEEP,
  stored=STORED,
  sv_inner_burnin=SV_INNER_BURNIN,
  B1=calibration$B1,
  B2=calibration$B2,
  kappa0=calibration$kappa0,
  thr_low=calibration$thr_low,
  thr_high=calibration$thr_high,
  sv_mu_sd=calibration$sv_mu_sd,
  terms=TERMS,
  static_terms=if (MODE=="static") STATIC_TERMS else character(),
  model_quarters=model_quarters,
  network=MAIN_NETWORK,
  panel=basename(PANEL),
  synthetic_breaks=if (CASE=="SYNTH_JP") synthetic$break_quarters else character()
)

part <- list(
  meta=meta,
  equations=equation_objects,
  monitor=as.matrix(monitor),
  synthetic_truth=if (CASE=="SYNTH_JP") synthetic$beta_true else NULL
)

saveRDS(part,file.path(OUT,"10n_chain.rds"),compress="xz")

summary_df <- do.call(rbind,summary_rows)
tv_df <- do.call(rbind,tv_rows)
write.csv(summary_df,file.path(OUT,"01_chain_term_summary.csv"),row.names=FALSE)
write.csv(tv_df,file.path(OUT,"02_chain_tvp_by_quarter.csv"),row.names=FALSE)

if (length(truth_rows)) {
  write.csv(
    do.call(rbind,truth_rows),
    file.path(OUT,"03_synthetic_gpr0_recovery.csv"),
    row.names=FALSE
  )
}

manifest <- data.frame(
  Scenario=scenario_id,
  Case=CASE,
  Country=case_country,
  Calibration=CALIB,
  CalibrationLabel=calibration$label,
  SamplerMode=MODE,
  Chain=CHAIN_ID,
  B1=calibration$B1,
  B2=calibration$B2,
  Kappa0=calibration$kappa0,
  ThresholdLow=calibration$thr_low,
  ThresholdHigh=calibration$thr_high,
  SVMuPriorSD=calibration$sv_mu_sd,
  Burn=BURN,
  PostBurnInternal=KEEP,
  Stored=STORED,
  SVInnerBurnin=SV_INNER_BURNIN,
  Seed=SEED,
  Network=MAIN_NETWORK,
  Panel=basename(PANEL),
  stringsAsFactors=FALSE
)
write.csv(manifest,file.path(OUT,"00_chain_manifest.csv"),row.names=FALSE)

cat("\n10N CHAIN COMPLETE\n")
print(manifest)
