#!/usr/bin/env Rscript

# =============================================================================
# 72_ttvp_unrestricted_sv_off_traceback.R
#
# 10p diagnostic only.
# Reproduces the SYNTH_JP + paper-code + unrestricted + SV=FALSE path and
# captures the R call stack BEFORE stack unwinding.
#
# It compares three package/call paths selected by FIN3_10P_IMPLEMENTATION:
#   pristine_direct  : pinned original package, call MCMC_tvp() directly
#   extended_direct  : 10k static-block extension installed, call original
#                      MCMC_tvp() directly
#   extended_wrapper : 10k extension installed, call estimate_tvp() with
#                      static_idx=integer(0)
#
# No formal model is changed. No IRF is authorized.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("threshtvp", quietly=TRUE)) {
  stopf("Package 'threshtvp' is required.")
}

get_chr <- function(name, default="") {
  z <- trimws(Sys.getenv(name, ""))
  if (nzchar(z)) z else default
}
get_int <- function(name, default) {
  z <- get_chr(name, as.character(default))
  x <- suppressWarnings(as.integer(z))
  if (is.na(x)) stopf("%s must be integer: %s", name, z)
  x
}

IMPL <- tolower(get_chr("FIN3_10P_IMPLEMENTATION", "pristine_direct"))
if (!IMPL %in% c("pristine_direct","extended_direct","extended_wrapper")) {
  stopf("Unknown FIN3_10P_IMPLEMENTATION: %s", IMPL)
}

BURN <- get_int("FIN3_10P_BURN", 1000L)
SAVE <- get_int("FIN3_10P_SAVE", 3000L)
STORED <- get_int("FIN3_10P_STORED", 300L)
SEED <- get_int("FIN3_10P_SEED", 20268001L)
if (BURN < 10L || SAVE < 100L || STORED < 10L || STORED >= SAVE) {
  stopf("Invalid short diagnostic iteration contract.")
}
THIN <- (STORED + 1) / SAVE

PANEL <- get_chr(
  "FIN3_FORMAL_PANEL",
  file.path(DERIVED_DIR, "panel_domestic_fin3_rate_diff.csv")
)
OUT <- file.path(RESULTS_DIR, "10p_traceback", IMPL)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

writeLines(
  c(
    "FormalModelChanged=FALSE",
    "IRFAuthorized=FALSE",
    paste0("Implementation=", IMPL),
    paste0("Burn=", BURN),
    paste0("Save=", SAVE),
    paste0("Thin=", format(THIN, digits=16))
  ),
  file.path(OUT, "00_hard_lock.txt")
)

# -----------------------------------------------------------------------------
# Build the exact accepted JP p=1/q=1 Trade-W design used by R/71.
# -----------------------------------------------------------------------------
if (!file.exists(PANEL)) stopf("Missing panel: %s", PANEL)
d <- read.csv(PANEL, stringsAsFactors=FALSE, check.names=FALSE)
need <- c("Quarter","Country",VARS,"gpr","brent")
if (!all(need %in% names(d))) {
  stopf("Panel missing: %s", paste(setdiff(need,names(d)), collapse=", "))
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
  NA_real_, c(Tn,N,K),
  dimnames=list(quarters,COUNTRIES,VARS)
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
  stopf("10p requires FIN3_NETWORK=trade_2000_2012.")
}

W <- read_weight_matrix(network_path(MAIN_NETWORK))
country_i <- match("JP", COUNTRIES)
if (is.na(country_i)) stopf("JP not found in COUNTRIES.")

STAR <- array(
  NA_real_, c(Tn,N,K),
  dimnames=list(quarters,COUNTRIES,VARS)
)
for (i in seq_len(N)) {
  for (v in seq_len(K)) {
    STAR[,i,v] <- as.numeric(Xglobal[,,v] %*% W[i,])
  }
}

lag1 <- function(x) c(NA_real_, x[-length(x)])
Yraw <- Xglobal[,country_i,,drop=FALSE][,1,]
Zraw <- STAR[,country_i,,drop=FALSE][,1,]
rows <- 2:nrow(Yraw)

D <- data.frame(const=rep(1,length(rows)), check.names=FALSE)
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
model_quarters <- model_quarters[ok]

TERMS <- colnames(D)
EXPECTED_TERMS <- c(
  "const",
  paste0(VARS,"_L1"),
  paste0(VARS,"_star_0"),
  paste0(VARS,"_star_L1"),
  "gpr_0","gpr_L1","brent_0","brent_L1"
)
if (!identical(TERMS, EXPECTED_TERMS)) stopf("Unexpected design term ordering.")

Xm <- as.matrix(D)
for (j in seq_len(ncol(Xm))) {
  if (TERMS[j] == "const") next
  sc <- stats::sd(Xm[,j])
  if (!is.finite(sc) || sc < 1e-10) stopf("Near-constant predictor: %s", TERMS[j])
  Xm[,j] <- (Xm[,j] - mean(Xm[,j])) / sc
}
Xm[,"const"] <- 1

# -----------------------------------------------------------------------------
# Exact 10n/10o synthetic positive-control response.
# -----------------------------------------------------------------------------
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
  beta[,"gpr_0"] <- ifelse(qid < b1, 0, ifelse(qid < b2, -0.60, 0.50))

  mu_h <- -2
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
  as.numeric(rowSums(Xm * beta) + eps)
}
y_std <- make_synthetic_y(Xm, model_quarters, TERMS)
if (length(y_std) != nrow(Xm) || any(!is.finite(y_std))) {
  stopf("Synthetic response contract failed.")
}

# Paper replication-code calibration.
priorbtheta <- list(B_1=3, B_2=0.03, kappa0=-0.005)
priorb0 <- list(a_tau=0.1, c_tau=0.01, d_tau=0.01)
priorsig <- c(0.01,0.01)
priorphi <- c(2,2)
priormu <- c(0,10)
set.seed(SEED)

# -----------------------------------------------------------------------------
# Pre-unwind traceback capture.
# -----------------------------------------------------------------------------
summarize_obj <- function(x) {
  if (is.null(x)) return("NULL")
  cl <- paste(class(x),collapse="/")
  dm <- dim(x)
  n <- length(x)
  nm <- names(x)
  bits <- c(paste0("class=",cl), paste0("length=",n))
  if (!is.null(dm)) bits <- c(bits, paste0("dim=",paste(dm,collapse="x")))
  if (!is.null(nm)) bits <- c(bits, paste0("names=",paste(head(nm,12),collapse="|")))
  if (is.atomic(x) && length(x) <= 8L) {
    bits <- c(bits, paste0("value=",paste(as.character(x),collapse="|")))
  }
  paste(bits,collapse=";")
}

TRACE_VARS <- c(
  "irep","ithin","ntot","T","K_","K","M","nr","thin","thin_out",
  "hv","sig_eta","svdraw","ALPHA","D_t","Omega_t","d","tau2",
  "H_store","ALPHA_store","svparms_store","D_store",
  "thresholds_store","Omega_store","omega_store","V0_store",
  "sigma2_store","thrshprior_store"
)

captured <- new.env(parent=emptyenv())
captured$done <- FALSE

capture_handler <- function(e) {
  calls <- sys.calls()
  frames <- sys.frames()

  call_txt <- vapply(
    seq_along(calls),
    function(i) sprintf("%03d  %s", i, paste(deparse(calls[[i]]),collapse=" ")),
    character(1)
  )
  writeLines(
    c(
      paste0("error=", conditionMessage(e)),
      paste0("condition_class=", paste(class(e),collapse="|")),
      "",
      call_txt
    ),
    file.path(OUT,"02_call_stack_pre_unwind.txt")
  )

  rows <- list()
  rr <- 0L
  for (i in seq_along(frames)) {
    fr <- frames[[i]]
    present <- intersect(TRACE_VARS, ls(fr, all.names=TRUE))
    call_i <- if (i <= length(calls)) paste(deparse(calls[[i]]),collapse=" ") else ""
    if (!length(present)) {
      rr <- rr + 1L
      rows[[rr]] <- data.frame(
        Frame=i, Call=call_i, Variable="", Summary="",
        stringsAsFactors=FALSE
      )
    } else {
      for (nm in present) {
        rr <- rr + 1L
        val <- tryCatch(get(nm,envir=fr,inherits=FALSE), error=function(e2) e2)
        sm <- if (inherits(val,"error")) {
          paste0("READ_ERROR:",conditionMessage(val))
        } else {
          tryCatch(summarize_obj(val), error=function(e2) paste0("SUMMARY_ERROR:",conditionMessage(e2)))
        }
        rows[[rr]] <- data.frame(
          Frame=i, Call=call_i, Variable=nm, Summary=sm,
          stringsAsFactors=FALSE
        )
      }
    }
  }
  snap <- do.call(rbind, rows)
  write.csv(snap, file.path(OUT,"03_frame_snapshot_pre_unwind.csv"), row.names=FALSE)

  # A compact deepest-sampler-frame snapshot for rapid classification.
  sampler_frames <- which(vapply(
    seq_along(frames),
    function(i) "irep" %in% ls(frames[[i]], all.names=TRUE),
    logical(1)
  ))
  if (length(sampler_frames)) {
    ii <- tail(sampler_frames,1)
    fr <- frames[[ii]]
    val <- function(nm, default=NA_character_) {
      if (!exists(nm,envir=fr,inherits=FALSE)) return(default)
      z <- get(nm,envir=fr,inherits=FALSE)
      if (length(z)==1L && is.atomic(z)) as.character(z) else summarize_obj(z)
    }
    quick <- data.frame(
      Item=c("Error","SamplerFrame","irep","ithin","ntot","thin_out",
             "hv","sig_eta","svdraw","H_store","svparms_store",
             "ALPHA_store","D_store","Omega_store","sigma2_store"),
      Value=c(
        conditionMessage(e), as.character(ii), val("irep"), val("ithin"),
        val("ntot"), val("thin_out"), val("hv"), val("sig_eta"),
        val("svdraw"), val("H_store"), val("svparms_store"),
        val("ALPHA_store"), val("D_store"), val("Omega_store"),
        val("sigma2_store")
      ),
      stringsAsFactors=FALSE
    )
    write.csv(quick, file.path(OUT,"04_deepest_sampler_frame.csv"), row.names=FALSE)
  }

  captured$message <- conditionMessage(e)
  captured$done <- TRUE
}

run_sampler <- function() {
  common <- list(
    Y=matrix(y_std,ncol=1),
    X=Xm,
    priorbtheta=priorbtheta,
    priorb0=priorb0,
    priorsig=priorsig,
    priorphi=priorphi,
    priormu=priormu,
    h0prior="stationary",
    grid.length=150,
    thrsh.pct=0.1,
    thrsh.pct.high=1.5,
    sv_on=FALSE,
    TVS=TRUE,
    cons.mod=FALSE,
    thin=THIN
  )

  if (IMPL %in% c("pristine_direct","extended_direct")) {
    f <- getFromNamespace("MCMC_tvp","threshtvp")
    args <- c(
      common,
      list(
        nburn=BURN,
        nsave=SAVE,
        nr=1,
        robust=TRUE,
        a.approx=FALSE,
        sim.kappa=FALSE,
        kappa.grid=seq(1e-4,0.1,10)
      )
    )
    return(do.call(f,args))
  }

  args <- c(
    common,
    list(
      save=SAVE,
      burn=BURN,
      static_idx=integer(0),
      sv_inner_burnin=4L,
      p=1,
      CPU=1,
      approx=FALSE,
      sim.kappa0=FALSE,
      kappa0.grid=seq(1e-4,0.1,10)
    )
  )
  do.call(threshtvp::estimate_tvp,args)
}

fit <- tryCatch(
  withCallingHandlers(run_sampler(), error=capture_handler),
  error=function(e) e
)

completed <- !inherits(fit,"error")
err <- if (completed) "" else conditionMessage(fit)

if (!completed && !isTRUE(captured$done)) {
  writeLines(
    c(
      paste0("error=",err),
      "WARNING=calling handler did not capture a pre-unwind stack"
    ),
    file.path(OUT,"02_call_stack_pre_unwind.txt")
  )
}

phase <- "UNKNOWN"
if (!completed && file.exists(file.path(OUT,"04_deepest_sampler_frame.csv"))) {
  q <- read.csv(file.path(OUT,"04_deepest_sampler_frame.csv"), stringsAsFactors=FALSE)
  getv <- function(k) {
    z <- q$Value[q$Item==k]
    if (length(z)) z[1] else NA_character_
  }
  ir <- suppressWarnings(as.integer(getv("irep")))
  ith <- suppressWarnings(as.integer(getv("ithin")))
  if (is.finite(ir) && is.finite(ith)) {
    phase <- if (ith > 0L) "AT_OR_AFTER_FIRST_STORAGE" else "BEFORE_FIRST_STORAGE"
  }
}

status <- data.frame(
  TechnicalStatus="PASS",
  Diagnostic="10p_unrestricted_sv_off_traceback",
  Implementation=IMPL,
  PackageVersion=as.character(utils::packageVersion("threshtvp")),
  SamplerCompleted=completed,
  ErrorMessage=err,
  ErrorClass=if (completed) "" else paste(class(fit),collapse="|"),
  FailurePhase=phase,
  Burn=BURN,
  Save=SAVE,
  RequestedStored=STORED,
  Thin=THIN,
  Seed=SEED,
  SV=FALSE,
  TVS=TRUE,
  B1=3,
  B2=0.03,
  Kappa0=-0.005,
  ThresholdLow=0.1,
  ThresholdHigh=1.5,
  FormalModelChanged=FALSE,
  IRFAuthorized=FALSE,
  stringsAsFactors=FALSE
)
write.csv(status, file.path(OUT,"01_status.csv"), row.names=FALSE)

if (completed) {
  saveRDS(fit, file.path(OUT,"05_short_fit.rds"))
}

cat("\n===== 10p STATUS =====\n")
print(status)
cat("FormalModelChanged=FALSE\n")
cat("IRFAuthorized=FALSE\n")
