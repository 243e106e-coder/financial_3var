#!/usr/bin/env Rscript

# =============================================================================
# 70_ttvp_recovery_summarize.R
#
# Aggregates the 10n 4-chain diagnostic scenarios.
# This is a diagnostic decision layer, NOT a formal model-selection gate.
# =============================================================================

if (!requireNamespace("posterior",quietly=TRUE)) {
  stop("posterior package is required.",call.=FALSE)
}

ROOT <- Sys.getenv("FIN3_10N_PARTS_ROOT","posterior_parts")
OUT <- Sys.getenv("FIN3_10N_SUMMARY_OUT","results/10n_ttvp_recovery_audit")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

RHAT_HARD <- 1.05
ESS_HARD <- 100
RHAT_TARGET <- 1.01
ESS_TARGET <- 400

expected_scenarios <- c(
  "REAL_JP__current_10l__static",
  "REAL_JP__paper_code__static",
  "REAL_JP__paper_code_wide__static",
  "REAL_SG__current_10l__static",
  "REAL_SG__paper_code__static",
  "REAL_SG__paper_code_wide__static",
  "REAL_TR__current_10l__static",
  "REAL_TR__paper_code__static",
  "REAL_TR__paper_code_wide__static",
  "SYNTH_JP__current_10l__static",
  "SYNTH_JP__paper_code__static",
  "SYNTH_JP__paper_code_wide__static",
  "SYNTH_JP__paper_code__unrestricted"
)

files <- list.files(
  ROOT,pattern="10n_chain[.]rds$",
  recursive=TRUE,full.names=TRUE
)
if (length(files) != length(expected_scenarios)*4L) {
  stop(
    sprintf(
      "Expected %d chain files, found %d.",
      length(expected_scenarios)*4L,length(files)
    ),
    call.=FALSE
  )
}

parts <- lapply(files,readRDS)
scenario <- vapply(parts,function(z) as.character(z$meta$scenario),character(1))
chain <- vapply(parts,function(z) as.integer(z$meta$chain),integer(1))

for (sc in expected_scenarios) {
  ix <- which(scenario==sc)
  if (length(ix)!=4L || !identical(sort(chain[ix]),1:4)) {
    stop(sprintf("Incomplete scenario: %s",sc),call.=FALSE)
  }
}

# -----------------------------------------------------------------------------
# Calibration / lineage manifest
# -----------------------------------------------------------------------------

manifest <- do.call(rbind,lapply(parts,function(z){
  m <- z$meta
  data.frame(
    Scenario=m$scenario,
    Case=m$case,
    Country=m$country,
    Calibration=m$calibration,
    CalibrationLabel=m$calibration_label,
    SamplerMode=m$sampler_mode,
    Chain=m$chain,
    Seed=m$seed,
    Burn=m$burn,
    PostBurnInternal=m$keep,
    Stored=m$stored,
    SVInnerBurnin=m$sv_inner_burnin,
    B1=m$B1,
    B2=m$B2,
    Kappa0=m$kappa0,
    ThresholdLow=m$thr_low,
    ThresholdHigh=m$thr_high,
    SVMuPriorSD=m$sv_mu_sd,
    Network=m$network,
    Panel=m$panel,
    stringsAsFactors=FALSE
  )
}))
manifest <- manifest[order(manifest$Scenario,manifest$Chain),]
write.csv(manifest,file.path(OUT,"00_calibration_manifest.csv"),row.names=FALSE)

source_manifest <- data.frame(
  Source=c(
    "JRSSA_2019_official_replication_main",
    "JRSSA_2019_official_replication_BVAR",
    "JRSSA_2019_official_replication_threshtvp_tar"
  ),
  SHA256=c(
    "8d7e8501208467831794de8d05fb0b26aed26c62ca2eed27202948fe45560ce4",
    "a0226fd591fdab524bddb9ce0695ac3998e03c32920cd6fe8f9688360dee3a18",
    "d573e044c8f9218ec1ada2578b589cd6a1e95a8bf54278de6b9104627eff9661"
  ),
  Note=c(
    "01b_baseline_estimation_sp_levels.r: B1=3, B2=0.03, kappa0=-0.1/20",
    "BVAR_ttvp.r: thrsh.pct=0.1, thrsh.pct.high=1.5, priormu=c(0,10)",
    "Official replication archive bundled threshtvp 0.2 source"
  ),
  stringsAsFactors=FALSE
)
write.csv(source_manifest,file.path(OUT,"00b_replication_source_manifest.csv"),row.names=FALSE)

# -----------------------------------------------------------------------------
# Convergence diagnostics by scenario
# -----------------------------------------------------------------------------

safe_metric <- function(fun,x,bad) {
  z <- try(fun(x),silent=TRUE)
  if (inherits(z,"try-error") || length(z)!=1L || !is.finite(z)) return(bad)
  as.numeric(z)
}

diag_rows <- list()
scenario_conv <- list()
dr <- 0L
sr <- 0L

for (sc in expected_scenarios) {
  ix <- which(scenario==sc)
  zlist <- parts[ix][order(chain[ix])]
  mons <- lapply(zlist,function(z) as.matrix(z$monitor))
  pnames <- colnames(mons[[1]])

  if (!all(vapply(mons,function(m) identical(colnames(m),pnames),logical(1)))) {
    stop(sprintf("Monitor columns differ: %s",sc),call.=FALSE)
  }

  scenario_param <- list()
  pr <- 0L

  for (j in seq_along(pnames)) {
    mat <- do.call(cbind,lapply(mons,function(m)m[,j]))
    finite <- all(is.finite(mat))
    constant <- finite && (max(mat)-min(mat)<1e-12)
    active <- finite && !constant

    if (active) {
      rh <- safe_metric(posterior::rhat,mat,Inf)
      eb <- safe_metric(posterior::ess_bulk,mat,0)
      et <- safe_metric(posterior::ess_tail,mat,0)
    } else {
      rh <- eb <- et <- NA_real_
    }

    dr <- dr+1L
    diag_rows[[dr]] <- data.frame(
      Scenario=sc,
      Parameter=pnames[j],
      Active=active,
      Rhat=rh,
      ESSBulk=eb,
      ESSTail=et,
      HardPass=!active || (
        is.finite(rh) && rh<=RHAT_HARD &&
        is.finite(eb) && eb>=ESS_HARD &&
        is.finite(et) && et>=ESS_HARD
      ),
      stringsAsFactors=FALSE
    )
  }

  dsc <- do.call(rbind,diag_rows)[vapply(diag_rows,function(x)x$Scenario[1]==sc,logical(1)),]
  a <- dsc[dsc$Active,,drop=FALSE]

  sr <- sr+1L
  scenario_conv[[sr]] <- data.frame(
    Scenario=sc,
    ActiveDiagnostics=nrow(a),
    MaxRhat=if(nrow(a)) max(a$Rhat) else NA_real_,
    MinESSBulk=if(nrow(a)) min(a$ESSBulk) else NA_real_,
    MinESSTail=if(nrow(a)) min(a$ESSTail) else NA_real_,
    HardFailCount=if(nrow(a)) sum(!a$HardPass) else 0L,
    RhatTargetShare=if(nrow(a)) mean(a$Rhat<=RHAT_TARGET) else NA_real_,
    BulkTargetShare=if(nrow(a)) mean(a$ESSBulk>=ESS_TARGET) else NA_real_,
    TailTargetShare=if(nrow(a)) mean(a$ESSTail>=ESS_TARGET) else NA_real_,
    stringsAsFactors=FALSE
  )
}

diag_df <- do.call(rbind,diag_rows)
conv_df <- do.call(rbind,scenario_conv)
write.csv(diag_df,file.path(OUT,"01_convergence_parameter_diagnostics.csv"),row.names=FALSE)
write.csv(conv_df,file.path(OUT,"02_convergence_by_scenario.csv"),row.names=FALSE)

# -----------------------------------------------------------------------------
# Combine TV probabilities and posterior mean coefficient paths across chains
# -----------------------------------------------------------------------------

tv_records <- list()
lock_records <- list()
ri <- 0L
li <- 0L

for (sc in expected_scenarios) {
  ix <- which(scenario==sc)
  zlist <- parts[ix][order(chain[ix])]

  eqs <- names(zlist[[1]]$equations)
  for (eq in eqs) {
    terms <- colnames(zlist[[1]]$equations[[eq]]$tv_probability)

    for (nm in terms) {
      tv_stack <- do.call(cbind,lapply(
        zlist,
        function(z) z$equations[[eq]]$tv_probability[,nm]
      ))
      beta_stack <- do.call(cbind,lapply(
        zlist,
        function(z) z$equations[[eq]]$beta_mean_std[,nm]
      ))
      quarters <- rownames(zlist[[1]]$equations[[eq]]$tv_probability)

      for (tt in seq_along(quarters)) {
        ri <- ri+1L
        tv_records[[ri]] <- data.frame(
          Scenario=sc,
          Case=zlist[[1]]$meta$case,
          Country=zlist[[1]]$meta$country,
          Calibration=zlist[[1]]$meta$calibration,
          SamplerMode=zlist[[1]]$meta$sampler_mode,
          Equation=eq,
          Term=nm,
          Quarter=quarters[tt],
          MeanTVProbability=mean(tv_stack[tt,],na.rm=TRUE),
          MinChainTVProbability=min(tv_stack[tt,],na.rm=TRUE),
          MaxChainTVProbability=max(tv_stack[tt,],na.rm=TRUE),
          MeanPosteriorBetaStd=mean(beta_stack[tt,],na.rm=TRUE),
          SDChainPosteriorBetaMean=sd(beta_stack[tt,],na.rm=TRUE),
          stringsAsFactors=FALSE
        )
      }

      low <- vapply(
        zlist,
        function(z) unname(z$equations[[eq]]$low_innovation_sd[nm]),
        numeric(1)
      )
      thr <- vapply(
        zlist,
        function(z) unname(z$equations[[eq]]$mean_threshold[nm]),
        numeric(1)
      )
      omg <- vapply(
        zlist,
        function(z) unname(z$equations[[eq]]$mean_omega[nm]),
        numeric(1)
      )

      li <- li+1L
      lock_records[[li]] <- data.frame(
        Scenario=sc,
        Case=zlist[[1]]$meta$case,
        Country=zlist[[1]]$meta$country,
        Calibration=zlist[[1]]$meta$calibration,
        SamplerMode=zlist[[1]]$meta$sampler_mode,
        Equation=eq,
        Term=nm,
        MeanLowInnovationSD=mean(low,na.rm=TRUE),
        MeanThreshold=mean(thr,na.rm=TRUE),
        MeanHighInnovationSD=mean(sqrt(pmax(omg,0)),na.rm=TRUE),
        ThresholdToLowSD=ifelse(
          mean(low,na.rm=TRUE)>0,
          mean(thr,na.rm=TRUE)/mean(low,na.rm=TRUE),
          NA_real_
        ),
        stringsAsFactors=FALSE
      )
    }
  }
}

tv_df <- do.call(rbind,tv_records)
lock_df <- do.call(rbind,lock_records)
write.csv(tv_df,file.path(OUT,"03_tvp_probability_by_quarter.csv"),row.names=FALSE)
write.csv(lock_df,file.path(OUT,"04_threshold_lock_metrics.csv"),row.names=FALSE)

# Key-term activation summary.
is_key <- tv_df$Term %in% c(
  "gpr_0","gpr_L1",
  "r_L1","de_L1","deq_L1",
  "r_star_L1","de_star_L1","deq_star_L1"
)
key <- tv_df[is_key,,drop=FALSE]

key_summary <- do.call(rbind,lapply(
  split(key,paste(key$Scenario,key$Equation,key$Term,sep="||")),
  function(z) data.frame(
    Scenario=z$Scenario[1],
    Case=z$Case[1],
    Country=z$Country[1],
    Calibration=z$Calibration[1],
    SamplerMode=z$SamplerMode[1],
    Equation=z$Equation[1],
    Term=z$Term[1],
    MeanTVProbability=mean(z$MeanTVProbability),
    MaxTVProbability=max(z$MeanTVProbability),
    ShareQuartersTVProbGT10=mean(z$MeanTVProbability>.10),
    ShareQuartersTVProbGT50=mean(z$MeanTVProbability>.50),
    BetaPathRange=max(z$MeanPosteriorBetaStd)-min(z$MeanPosteriorBetaStd),
    stringsAsFactors=FALSE
  )
))
write.csv(key_summary,file.path(OUT,"05_key_term_activation_summary.csv"),row.names=FALSE)

# -----------------------------------------------------------------------------
# Synthetic positive-control recovery
# -----------------------------------------------------------------------------

synth <- tv_df[
  tv_df$Case=="SYNTH_JP" &
  tv_df$Equation=="deq" &
  tv_df$Term=="gpr_0",
,drop=FALSE]

qid <- function(x) {
  x <- toupper(as.character(x))
  yr <- as.integer(substr(x,1,4))
  q <- as.integer(substr(x,6,6))
  4L*yr+q
}

truth_beta <- function(q) {
  z <- qid(q)
  ifelse(
    z < qid("2008Q3"),0,
    ifelse(z < qid("2020Q1"),-0.60,0.50)
  )
}

synth$TrueGPR0BetaStd <- truth_beta(synth$Quarter)
synth$IsBreakWindow <- (
  abs(qid(synth$Quarter)-qid("2008Q3"))<=1L |
  abs(qid(synth$Quarter)-qid("2020Q1"))<=1L
)

write.csv(
  synth,
  file.path(OUT,"06_synthetic_recovery_by_quarter.csv"),
  row.names=FALSE
)

synth_summary <- do.call(rbind,lapply(split(synth,synth$Scenario),function(z){
  pre <- z[qid(z$Quarter)<qid("2008Q3"),]
  mid <- z[qid(z$Quarter)>=qid("2008Q3") & qid(z$Quarter)<qid("2020Q1"),]
  post <- z[qid(z$Quarter)>=qid("2020Q1"),]

  w1 <- z[abs(qid(z$Quarter)-qid("2008Q3"))<=1L,]
  w2 <- z[abs(qid(z$Quarter)-qid("2020Q1"))<=1L,]
  nonbreak <- z[!z$IsBreakWindow,]

  pre_beta <- mean(pre$MeanPosteriorBetaStd)
  mid_beta <- mean(mid$MeanPosteriorBetaStd)
  post_beta <- mean(post$MeanPosteriorBetaStd)

  break1 <- max(w1$MeanTVProbability)
  break2 <- max(w2$MeanTVProbability)
  false_act <- mean(nonbreak$MeanTVProbability)

  sign_order_pass <- (
    is.finite(pre_beta) && is.finite(mid_beta) && is.finite(post_beta) &&
    mid_beta < pre_beta - 0.15 &&
    post_beta > mid_beta + 0.25 &&
    post_beta > pre_beta + 0.10
  )
  tv_break_pass <- (break1>=0.50 && break2>=0.50)
  false_activation_pass <- false_act <= 0.35

  data.frame(
    Scenario=z$Scenario[1],
    Calibration=z$Calibration[1],
    SamplerMode=z$SamplerMode[1],
    Break2008MaxTVProb=break1,
    Break2020MaxTVProb=break2,
    NonBreakMeanTVProb=false_act,
    PrePosteriorMeanBeta=pre_beta,
    MidPosteriorMeanBeta=mid_beta,
    PostPosteriorMeanBeta=post_beta,
    SignOrderPass=sign_order_pass,
    TVBreakPass=tv_break_pass,
    FalseActivationPass=false_activation_pass,
    RecoveryPass=(sign_order_pass && tv_break_pass && false_activation_pass),
    stringsAsFactors=FALSE
  )
}))
write.csv(
  synth_summary,
  file.path(OUT,"07_synthetic_recovery_summary.csv"),
  row.names=FALSE
)

# -----------------------------------------------------------------------------
# Real-data calibration comparison
# -----------------------------------------------------------------------------

real_key <- key_summary[
  grepl("^REAL_",key_summary$Case) &
  key_summary$Term %in% c("gpr_0","gpr_L1"),
,drop=FALSE]

real_compare <- do.call(rbind,lapply(
  split(real_key,paste(real_key$Case,real_key$Equation,sep="||")),
  function(z) {
    getv <- function(cal,col,fun=max) {
      q <- z[z$Calibration==cal,col]
      if (!length(q)) NA_real_ else fun(q,na.rm=TRUE)
    }
    data.frame(
      Case=z$Case[1],
      Country=z$Country[1],
      Equation=z$Equation[1],
      CurrentMaxGPRTVProb=getv("current_10l","MaxTVProbability"),
      PaperCodeMaxGPRTVProb=getv("paper_code","MaxTVProbability"),
      PaperWideMaxGPRTVProb=getv("paper_code_wide","MaxTVProbability"),
      CurrentMeanGPRTVProb=getv("current_10l","MeanTVProbability",mean),
      PaperCodeMeanGPRTVProb=getv("paper_code","MeanTVProbability",mean),
      PaperWideMeanGPRTVProb=getv("paper_code_wide","MeanTVProbability",mean),
      stringsAsFactors=FALSE
    )
  }
))
write.csv(
  real_compare,
  file.path(OUT,"08_real_gpr_calibration_comparison.csv"),
  row.names=FALSE
)

# -----------------------------------------------------------------------------
# Decision logic
# -----------------------------------------------------------------------------

get_recovery <- function(sc) {
  z <- synth_summary[synth_summary$Scenario==sc,]
  if (!nrow(z)) return(NA)
  isTRUE(z$RecoveryPass[1])
}

cur_syn <- get_recovery("SYNTH_JP__current_10l__static")
paper_static_syn <- get_recovery("SYNTH_JP__paper_code__static")
paper_unres_syn <- get_recovery("SYNTH_JP__paper_code__unrestricted")
paper_wide_syn <- get_recovery("SYNTH_JP__paper_code_wide__static")

real_current_max <- max(
  real_compare$CurrentMaxGPRTVProb,
  na.rm=TRUE
)
real_paper_max <- max(
  real_compare$PaperCodeMaxGPRTVProb,
  na.rm=TRUE
)

decision_code <- if (isTRUE(paper_unres_syn) && !isTRUE(paper_static_syn)) {
  "STATICBLOCK_IMPLEMENTATION_SUSPECT__UNRESTRICTED_RECOVERS"
} else if (!isTRUE(paper_unres_syn) && !isTRUE(paper_static_syn)) {
  "SYNTHETIC_RECOVERY_FAILURE__DO_NOT_CHANGE_FORMAL_MODEL"
} else if (isTRUE(paper_static_syn) && !isTRUE(cur_syn)) {
  if (is.finite(real_paper_max) && real_paper_max > 0.10) {
    "CURRENT_10L_CALIBRATION_LOCK_IN_SUPPORTED__PAPER_CALIBRATION_OPENS_REAL_TVP"
  } else {
    "CURRENT_10L_CALIBRATION_LOCK_IN_SUPPORTED__REAL_DATA_TVP_STILL_WEAK"
  }
} else if (isTRUE(paper_static_syn) && isTRUE(cur_syn)) {
  if (is.finite(real_current_max) && real_current_max <= 0.01 &&
      is.finite(real_paper_max) && real_paper_max > 0.10) {
    "REAL_DATA_TVP_HIGHLY_PRIOR_SENSITIVE"
  } else {
    "BOTH_CALIBRATIONS_RECOVER_SYNTHETIC__REAL_DATA_RESULT_NEEDS_INTERPRETATION"
  }
} else {
  "INCONCLUSIVE_10N_DIAGNOSTIC"
}

# Technical status only checks file/parse completeness; an unfavorable recovery
# result is a scientifically valid diagnostic outcome and does not make the
# workflow a technical failure.
technical_status <- if (
  length(files)==length(expected_scenarios)*4L &&
  all(is.finite(manifest$B1)) &&
  all(is.finite(manifest$B2)) &&
  all(is.finite(manifest$Kappa0))
) "PASS" else "FAIL"

decision <- data.frame(
  Item=c(
    "TechnicalStatus",
    "CurrentSyntheticRecoveryPass",
    "PaperStaticSyntheticRecoveryPass",
    "PaperUnrestrictedSyntheticRecoveryPass",
    "PaperWideSyntheticRecoveryPass",
    "RealCurrentMaxGPRTVProbability",
    "RealPaperCodeMaxGPRTVProbability",
    "Decision",
    "FormalModelChanged",
    "IRFAuthorized"
  ),
  Value=c(
    technical_status,
    as.character(cur_syn),
    as.character(paper_static_syn),
    as.character(paper_unres_syn),
    as.character(paper_wide_syn),
    format(real_current_max,digits=8),
    format(real_paper_max,digits=8),
    decision_code,
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors=FALSE
)
write.csv(decision,file.path(OUT,"09_decision.csv"),row.names=FALSE)

cat("\n===== 10N SYNTHETIC RECOVERY =====\n")
print(synth_summary)
cat("\n===== 10N REAL GPR CALIBRATION COMPARISON =====\n")
print(real_compare)
cat("\n===== 10N DECISION =====\n")
print(decision)

if (!identical(technical_status,"PASS")) {
  stop("10n technical integrity failed.",call.=FALSE)
}
