#!/usr/bin/env Rscript

# =============================================================================
# 69_tr_h1_source_and_event_coefficient_tvp_audit.R
#
# 08o: TR h=1 structural-source decomposition + event-coefficient TVP audit
#
# NO MCMC / NO re-estimation.
# Reconstructs the accepted 56 posterior parts and reuses the exact R/67 algebra.
# =============================================================================

source("R/00_config.R")

get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, "")); if (!nzchar(z)) default else z
}
get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, "")); if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("Environment variable %s is not numeric: %s", name, z)
  out
}
qsum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(p05=NA,p16=NA,median=NA,p84=NA,p95=NA,mean=NA,prob_positive=NA,n=0))
  qq <- quantile(x, c(.05,.16,.50,.84,.95), names=FALSE)
  c(p05=qq[1],p16=qq[2],median=qq[3],p84=qq[4],p95=qq[5],mean=mean(x),prob_positive=mean(x>0),n=length(x))
}
safe_median <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else median(x) }
safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
safe_max <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }
safe_p95 <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else unname(quantile(x,.95)) }

PARTS_ROOT <- get_env_chr("FIN3_PARTS_ROOT", "posterior_parts")
SOURCE_08M_ROOT <- get_env_chr("FIN3_08M_ROOT", "source_08m")
SOURCE_08N_ROOT <- get_env_chr("FIN3_08N_ROOT", "source_08n")
OUT <- get_env_chr("FIN3_08O_OUT", file.path(RESULTS_DIR,"tr_h1_source_tvp_audit"))
NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS",4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN",2000))
GPR_SHOCK_PCT <- get_env_num("FIN3_GPR_SHOCK_PCT",10)
MIN_G0_RCOND <- get_env_num("FIN3_IRF_MIN_G0_RCOND",1e-10)
MIN_VALID_SHARE <- get_env_num("FIN3_08O_MIN_VALID_SHARE",.90)
REPRO_TOL <- get_env_num("FIN3_08O_REPRO_TOL",1e-8)
ALGEBRA_TOL <- get_env_num("FIN3_08O_ALGEBRA_TOL",1e-10)
EXPECTED_PANEL_BASENAME <- get_env_chr("FIN3_RATE_DIFF_PANEL_BASENAME","panel_domestic_fin3_rate_diff.csv")
SOURCE_08J_RUN <- get_env_chr("FIN3_08J_RUN_ID","")
SOURCE_08K_RUN <- get_env_chr("FIN3_08K_RUN_ID","")
SOURCE_08L_RUN <- get_env_chr("FIN3_08L_RUN_ID","")
SOURCE_08M_RUN <- get_env_chr("FIN3_08M_RUN_ID","")
SOURCE_08N_RUN <- get_env_chr("FIN3_08N_RUN_ID","")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

if (NCHAINS != 4L) stopf("08o expects exactly four chains.")
if (EXPECTED_STORED != 2000L) stopf("08o expects 2000 stored draws per chain.")

# ----- accepted 08m / 08n gates ------------------------------------------------
m_gate_file <- file.path(SOURCE_08M_ROOT,"formal_irf_rate_diff","00_formal_rate_diff_irf_gate.csv")
m_raw_file <- file.path(SOURCE_08M_ROOT,"formal_irf_rate_diff","irf_posterior_summary.csv")
n_gate_file <- file.path(SOURCE_08N_ROOT,"00_08n_gate.csv")
if (!file.exists(m_gate_file) || !file.exists(m_raw_file) || !file.exists(n_gate_file)) stopf("Missing accepted 08m/08n source artifacts.")
m_gate <- read.csv(m_gate_file, stringsAsFactors=FALSE, check.names=FALSE)
m_raw <- read.csv(m_raw_file, stringsAsFactors=FALSE, check.names=FALSE)
n_gate <- read.csv(n_gate_file, stringsAsFactors=FALSE, check.names=FALSE)
if (nrow(m_gate)!=1L || m_gate$Status[1] != "READY_FOR_RATE_DIFF_IRF_AUDIT") stopf("08m gate not accepted.")
if (nrow(n_gate)!=1L || n_gate$Status[1] != "DIAGNOSTIC_COMPLETE") stopf("08n gate not accepted.")
if ("TRDiagnosticStatus" %in% names(n_gate) && n_gate$TRDiagnosticStatus[1] != "TR_REER_H1_ANOMALY_CONCENTRATED") stopf("08n does not contain expected TR h=1 diagnosis.")

# ----- exact 56 posterior lineage ----------------------------------------------
files <- list.files(PARTS_ROOT, pattern="^formal_tvp_[A-Z]{2}_chain[0-9]+[.]rds$", recursive=TRUE, full.names=TRUE)
expected_n <- length(COUNTRIES)*NCHAINS
if (length(files)!=expected_n) stopf("Expected %d posterior RDS files; found %d.", expected_n, length(files))
parts <- setNames(vector("list",length(COUNTRIES)),COUNTRIES)
lineage_rows <- list()
for (cc in COUNTRIES) {
  parts[[cc]] <- vector("list",NCHAINS)
  for (ch in seq_len(NCHAINS)) {
    cand <- files[grepl(sprintf("formal_tvp_%s_chain%d[.]rds$",cc,ch),files)]
    if (length(cand)!=1L) stopf("Expected one posterior for %s chain %d.",cc,ch)
    z <- readRDS(cand)
    if (as.character(z$meta$country)!=cc || as.integer(z$meta$chain)!=ch) stopf("Posterior metadata mismatch.")
    if (as.integer(z$meta$stored_draws)!=EXPECTED_STORED) stopf("Stored draw mismatch.")
    if (basename(as.character(z$meta$source_panel))!=EXPECTED_PANEL_BASENAME) stopf("Non-Delta-r posterior detected.")
    if (!isTRUE(z$meta$tvs_on) || !isTRUE(z$meta$sv_on)) stopf("TVS/SV metadata mismatch.")
    parts[[cc]][[ch]] <- z
    lineage_rows[[length(lineage_rows)+1L]] <- data.frame(Country=cc,Chain=ch,Seed=as.integer(z$meta$seed),Burn=as.integer(z$meta$burn),KeepInternal=as.integer(z$meta$keep_internal),StoredDraws=as.integer(z$meta$stored_draws),SourcePanel=basename(as.character(z$meta$source_panel)),MainNetwork=as.character(z$meta$main_network),File=cand,stringsAsFactors=FALSE)
  }
}
lineage <- do.call(rbind,lineage_rows)
key <- paste(lineage$Country,lineage$Chain,sep="||")
expected_key <- paste(rep(COUNTRIES,each=NCHAINS),rep(seq_len(NCHAINS),times=length(COUNTRIES)),sep="||")
if (anyDuplicated(key) || !setequal(key,expected_key)) stopf("Posterior grid is not exact 14x4.")
if (any(lineage$MainNetwork!=MAIN_NETWORK)) stopf("Main network mismatch.")
write.csv(lineage,file.path(OUT,"00_posterior_lineage.csv"),row.names=FALSE)

ref <- parts[[COUNTRIES[1]]][[1]]
anchors <- ref$anchors
KEY_TERMS <- ref$meta$key_terms
ALL_TERMS <- ref$meta$terms
for (cc in COUNTRIES) for (ch in seq_len(NCHAINS)) {
  z <- parts[[cc]][[ch]]
  if (!identical(z$meta$key_terms,KEY_TERMS) || !identical(z$meta$terms,ALL_TERMS)) stopf("Term ordering mismatch: %s chain %d",cc,ch)
  if (!identical(paste(z$anchors$EventID,z$anchors$AnchorType,z$anchors$AnchorQuarter),paste(anchors$EventID,anchors$AnchorType,anchors$AnchorQuarter))) stopf("Anchor ordering mismatch: %s chain %d",cc,ch)
}
core_idx <- which(anchors$EventSet=="CORE" & anchors$AnchorType=="EVENT_QUARTER_t0")
if (length(core_idx)!=6L) stopf("Expected six core t0 anchors; found %d.",length(core_idx))
core_anchors <- anchors[core_idx,,drop=FALSE]
core_events <- as.character(core_anchors$EventID)

N <- length(COUNTRIES); K <- length(VARS); NK <- N*K; DRAWS_TOTAL <- NCHAINS*EXPECTED_STORED
tr_i <- match("TR",COUNTRIES); de_i <- match("de",VARS); target <- (tr_i-1L)*K+de_i
idx_dom <- match(paste0(VARS,"_L1"),KEY_TERMS)
idx_star0 <- match(paste0(VARS,"_star_0"),KEY_TERMS)
idx_star1 <- match(paste0(VARS,"_star_L1"),KEY_TERMS)
idx_gpr0 <- match("gpr_0",KEY_TERMS); idx_gpr1 <- match("gpr_L1",KEY_TERMS)
if (any(is.na(c(idx_dom,idx_star0,idx_star1,idx_gpr0,idx_gpr1)))) stopf("Required GVAR/GPR terms missing.")
if (!grepl("^LN_",toupper(GPR_COLUMN))) stopf("Expected logged GPR column.")
shock_size <- log1p(GPR_SHOCK_PCT/100)
W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])
selector <- function(i) { S <- matrix(0,K,NK); S[,((i-1L)*K+1L):(i*K)] <- diag(K); S }
star_map <- function(i) { R <- matrix(0,K,NK); for (j in seq_len(N)) for (v in seq_len(K)) R[v,(j-1L)*K+v] <- W[i,j]; R }
SELECTORS <- lapply(seq_len(N),selector); STAR_MAPS <- lapply(seq_len(N),star_map)
source_country <- rep(COUNTRIES,each=K); source_variable <- rep(VARS,times=N); source_label <- paste(source_country,source_variable,sep="__")

# ----- exact TR/de/h=1 decomposition ------------------------------------------
decomp_rows <- list(); prop_source_rows <- list(); lag_source_rows <- list(); cf_rows <- list(); repro_rows <- list()
for (ee in seq_along(core_idx)) {
  a <- core_idx[ee]; meta <- anchors[a,,drop=FALSE]
  total <- propagation <- lagged_gpr <- local_tr_de_lag_mapped <- local_tr_de_lag_raw <- counterfactual <- rep(NA_real_,DRAWS_TOTAL)
  prop_contrib <- matrix(NA_real_,DRAWS_TOTAL,NK,dimnames=list(NULL,source_label))
  lag_contrib <- matrix(NA_real_,DRAWS_TOTAL,NK,dimnames=list(NULL,source_label))
  gd <- 0L
  for (ch in seq_len(NCHAINS)) for (dd in seq_len(EXPECTED_STORED)) {
    gd <- gd+1L; G0 <- matrix(0,NK,NK); G1 <- matrix(0,NK,NK); c0 <- numeric(NK); c1 <- numeric(NK)
    for (i in seq_len(N)) {
      cc <- COUNTRIES[i]
      coef <- parts[[cc]][[ch]]$event_coef[dd,a,,,drop=FALSE]
      coef <- matrix(coef,nrow=K,ncol=length(KEY_TERMS),dimnames=list(VARS,KEY_TERMS))
      Ai <- coef[,idx_dom,drop=FALSE]; B0i <- coef[,idx_star0,drop=FALSE]; B1i <- coef[,idx_star1,drop=FALSE]
      rr <- ((i-1L)*K+1L):(i*K)
      G0[rr,] <- SELECTORS[[i]] - B0i %*% STAR_MAPS[[i]]
      G1[rr,] <- Ai %*% SELECTORS[[i]] + B1i %*% STAR_MAPS[[i]]
      c0[rr] <- coef[,idx_gpr0]; c1[rr] <- coef[,idx_gpr1]
    }
    rc <- tryCatch(rcond(G0),error=function(e) NA_real_)
    if (!is.finite(rc) || rc<MIN_G0_RCOND) next
    impact <- tryCatch(solve(G0),error=function(e) NULL)
    if (is.null(impact) || any(!is.finite(impact))) next
    Fmat <- impact %*% G1
    rho <- tryCatch(max(Mod(eigen(Fmat,only.values=TRUE)$values)),error=function(e) NA_real_)
    if (!is.finite(rho) || rho>=1) next
    x0 <- as.numeric(impact %*% (c0*shock_size))
    pc <- as.numeric(Fmat[target,]*x0)
    lc <- as.numeric(impact[target,]*(c1*shock_size))
    prop_contrib[gd,] <- pc; lag_contrib[gd,] <- lc
    propagation[gd] <- sum(pc); lagged_gpr[gd] <- sum(lc); total[gd] <- propagation[gd]+lagged_gpr[gd]
    local_tr_de_lag_raw[gd] <- c1[target]*shock_size
    local_tr_de_lag_mapped[gd] <- lc[target]
    counterfactual[gd] <- total[gd]-lc[target]
  }
  valid <- is.finite(total); valid_share <- mean(valid)
  if (valid_share<MIN_VALID_SHARE) stopf("%s valid share %.6f below threshold.",meta$EventID,valid_share)
  max_alg_error <- max(abs(total[valid]-propagation[valid]-lagged_gpr[valid]),na.rm=TRUE)
  lag_abs_share <- abs(lagged_gpr)/pmax(abs(propagation)+abs(lagged_gpr),1e-300)
  prop_abs_share <- abs(propagation)/pmax(abs(propagation)+abs(lagged_gpr),1e-300)
  lag_abs_total <- rowSums(abs(lag_contrib),na.rm=TRUE)
  local_lag_share <- abs(local_tr_de_lag_mapped)/pmax(lag_abs_total,1e-300)
  qt <- qsum(total[valid]); qp <- qsum(propagation[valid]); ql <- qsum(lagged_gpr[valid]); qloc <- qsum(local_tr_de_lag_mapped[valid]); qcf <- qsum(counterfactual[valid]); qraw <- qsum(local_tr_de_lag_raw[valid])
  decomp_rows[[ee]] <- data.frame(EventID=meta$EventID,EventLabel=meta$EventLabel,EventQuarter=meta$EventQuarter,AnchorQuarter=meta$AnchorQuarter,ValidDraws=sum(valid),ValidShare=valid_share,Total_p05=qt["p05"],Total_median=qt["median"],Total_p95=qt["p95"],Propagation_p05=qp["p05"],Propagation_median=qp["median"],Propagation_p95=qp["p95"],LaggedGPR_p05=ql["p05"],LaggedGPR_median=ql["median"],LaggedGPR_p95=ql["p95"],LocalTRdeGPRL1Mapped_p05=qloc["p05"],LocalTRdeGPRL1Mapped_median=qloc["median"],LocalTRdeGPRL1Mapped_p95=qloc["p95"],MedianAbsPropagationShare=safe_median(prop_abs_share[valid]),MedianAbsLaggedGPRShare=safe_median(lag_abs_share[valid]),ProbAbsLaggedGPR_GT_Propagation=safe_mean(abs(lagged_gpr[valid])>abs(propagation[valid])),MedianAbsLocalTRdeShareWithinLagBlock=safe_median(local_lag_share[valid]),CounterfactualZeroTRdeGPRL1_median=qcf["median"],MedianDifferenceActualMinusCounterfactual=safe_median(total[valid]-counterfactual[valid]),MaxAlgebraReconciliationError=max_alg_error,stringsAsFactors=FALSE)
  prop_abs_rows <- rowSums(abs(prop_contrib[valid,,drop=FALSE]),na.rm=TRUE)
  lag_abs_rows <- rowSums(abs(lag_contrib[valid,,drop=FALSE]),na.rm=TRUE)
  for (j in seq_len(NK)) {
    pz <- prop_contrib[valid,j]; lz <- lag_contrib[valid,j]; qpj <- qsum(pz); qlj <- qsum(lz)
    prop_source_rows[[length(prop_source_rows)+1L]] <- data.frame(EventID=meta$EventID,SourceCountry=source_country[j],SourceVariable=source_variable[j],SourceIndex=j,p05=qpj["p05"],median=qpj["median"],p95=qpj["p95"],MeanAbsContribution=mean(abs(pz)),MedianAbsShareWithinPropagationBlock=safe_median(abs(pz)/pmax(prop_abs_rows,1e-300)),stringsAsFactors=FALSE)
    lag_source_rows[[length(lag_source_rows)+1L]] <- data.frame(EventID=meta$EventID,SourceCountry=source_country[j],SourceVariable=source_variable[j],SourceIndex=j,p05=qlj["p05"],median=qlj["median"],p95=qlj["p95"],MeanAbsContribution=mean(abs(lz)),MedianAbsShareWithinLaggedGPRBlock=safe_median(abs(lz)/pmax(lag_abs_rows,1e-300)),IsLocalTRdeGPRL1Source=(j==target),stringsAsFactors=FALSE)
  }
  cf_rows[[ee]] <- data.frame(EventID=meta$EventID,LocalTRdeGPRL1RawShockContribution_p05=qraw["p05"],LocalTRdeGPRL1RawShockContribution_median=qraw["median"],LocalTRdeGPRL1RawShockContribution_p95=qraw["p95"],LocalTRdeGPRL1MappedContribution_median=qloc["median"],ActualTRdeH1_median=qt["median"],ZeroTRdeGPRL1CounterfactualH1_median=qcf["median"],ActualMinusCounterfactual_median=safe_median(total[valid]-counterfactual[valid]),stringsAsFactors=FALSE)
  refrow <- m_raw[m_raw$EventID==meta$EventID & m_raw$EventSet=="CORE" & m_raw$AnchorType=="EVENT_QUARTER_t0" & toupper(m_raw$Country)=="TR" & tolower(m_raw$ResponseVariable)=="de" & as.integer(m_raw$Horizon)==1L,,drop=FALSE]
  if (nrow(refrow)!=1L) stopf("Unique 08m TR/de/h1 row not found for %s.",meta$EventID)
  repro_rows[[ee]] <- data.frame(EventID=meta$EventID,RecomputedMedian=qt["median"],Accepted08mMedian=as.numeric(refrow$median[1]),AbsoluteMedianGap=abs(qt["median"]-as.numeric(refrow$median[1])),RecomputedP05=qt["p05"],Accepted08mP05=as.numeric(refrow$p05[1]),AbsoluteP05Gap=abs(qt["p05"]-as.numeric(refrow$p05[1])),RecomputedP95=qt["p95"],Accepted08mP95=as.numeric(refrow$p95[1]),AbsoluteP95Gap=abs(qt["p95"]-as.numeric(refrow$p95[1])),ValidDraws=sum(valid),stringsAsFactors=FALSE)
}

decomp <- do.call(rbind,decomp_rows); prop_sources <- do.call(rbind,prop_source_rows); lag_sources <- do.call(rbind,lag_source_rows); counterfactual <- do.call(rbind,cf_rows); reproduction <- do.call(rbind,repro_rows)
prop_sources <- prop_sources[order(prop_sources$EventID,-prop_sources$MeanAbsContribution),,drop=FALSE]
lag_sources <- lag_sources[order(lag_sources$EventID,-lag_sources$MeanAbsContribution),,drop=FALSE]
write.csv(decomp,file.path(OUT,"01_tr_h1_exact_decomposition_summary.csv"),row.names=FALSE)
write.csv(prop_sources,file.path(OUT,"02_tr_h1_propagation_source_contributions.csv"),row.names=FALSE)
write.csv(lag_sources,file.path(OUT,"03_tr_h1_lagged_gpr_source_contributions.csv"),row.names=FALSE)
write.csv(counterfactual,file.path(OUT,"04_tr_de_gpr_l1_counterfactual.csv"),row.names=FALSE)
write.csv(reproduction,file.path(OUT,"05_reconciliation_with_accepted_08m.csv"),row.names=FALSE)

# ----- event-date coefficient posterior summaries + paired differences --------
coef_summary_rows <- list(); coef_variation_rows <- list(); coef_pair_rows <- list(); pair_index <- combn(seq_along(core_idx),2,simplify=FALSE)
for (cc in COUNTRIES) for (eq in seq_len(K)) for (tt in seq_along(KEY_TERMS)) {
  drawmat <- do.call(rbind,lapply(seq_len(NCHAINS),function(ch) {
    z <- parts[[cc]][[ch]]$event_coef[,core_idx,eq,tt,drop=FALSE]
    matrix(z,nrow=EXPECTED_STORED,ncol=length(core_idx))
  }))
  colnames(drawmat) <- core_events
  event_medians <- numeric(length(core_idx)); event_width90 <- numeric(length(core_idx))
  for (ee in seq_along(core_idx)) {
    qs <- qsum(drawmat[,ee]); event_medians[ee] <- qs["median"]; event_width90[ee] <- qs["p95"]-qs["p05"]
    coef_summary_rows[[length(coef_summary_rows)+1L]] <- data.frame(Country=cc,Equation=VARS[eq],Term=KEY_TERMS[tt],EventID=core_events[ee],AnchorQuarter=core_anchors$AnchorQuarter[ee],p05=qs["p05"],p16=qs["p16"],median=qs["median"],p84=qs["p84"],p95=qs["p95"],mean=qs["mean"],prob_positive=qs["prob_positive"],Draws=qs["n"],stringsAsFactors=FALSE)
  }
  rr <- range(event_medians,finite=TRUE)
  coef_variation_rows[[length(coef_variation_rows)+1L]] <- data.frame(Country=cc,Equation=VARS[eq],Term=KEY_TERMS[tt],CoreEvents=length(core_idx),MinEventMedian=rr[1],MaxEventMedian=rr[2],AcrossEventMedianRange=rr[2]-rr[1],AcrossEventSDofMedians=sd(event_medians),MedianWithinEventWidth90=median(event_width90),RangeToMedianWidth90=(rr[2]-rr[1])/pmax(median(event_width90),1e-300),stringsAsFactors=FALSE)
  for (pp in pair_index) {
    e1 <- pp[1]; e2 <- pp[2]; dz <- drawmat[,e2]-drawmat[,e1]; qd <- qsum(dz)
    coef_pair_rows[[length(coef_pair_rows)+1L]] <- data.frame(Country=cc,Equation=VARS[eq],Term=KEY_TERMS[tt],EventA=core_events[e1],EventB=core_events[e2],DifferenceDefinition="EventB_minus_EventA",p05=qd["p05"],median=qd["median"],p95=qd["p95"],ProbDifferencePositive=qd["prob_positive"],CredibleDifference90=(is.finite(qd["p05"]) && qd["p05"]>0)||(is.finite(qd["p95"]) && qd["p95"]<0),Draws=qd["n"],stringsAsFactors=FALSE)
  }
}
coef_summary <- do.call(rbind,coef_summary_rows); coef_variation <- do.call(rbind,coef_variation_rows); coef_pairs <- do.call(rbind,coef_pair_rows)
write.csv(coef_summary,file.path(OUT,"06_event_coefficient_summary_all.csv"),row.names=FALSE)
write.csv(coef_variation,file.path(OUT,"07_event_coefficient_variation_all.csv"),row.names=FALSE)
write.csv(coef_pairs,file.path(OUT,"08_paired_event_coefficient_differences_all.csv"),row.names=FALSE)
tr_de_coef <- coef_summary[coef_summary$Country=="TR" & coef_summary$Equation=="de",,drop=FALSE]
tr_de_pairs <- coef_pairs[coef_pairs$Country=="TR" & coef_pairs$Equation=="de",,drop=FALSE]
write.csv(tr_de_coef,file.path(OUT,"09_tr_de_core_event_coefficients.csv"),row.names=FALSE)
write.csv(tr_de_pairs,file.path(OUT,"10_tr_de_paired_event_coefficient_differences.csv"),row.names=FALSE)

# ----- latent-threshold time-variation probability -----------------------------
tv_core_rows <- list(); tv_summary_rows <- list()
for (cc in COUNTRIES) {
  mq <- parts[[cc]][[1]]$meta$model_quarters; terms <- parts[[cc]][[1]]$meta$terms
  for (ch in seq_len(NCHAINS)) {
    z <- parts[[cc]][[ch]]
    if (!identical(z$meta$model_quarters,mq) || !identical(z$meta$terms,terms)) stopf("Time-variation metadata mismatch for %s.",cc)
  }
  core_q_idx <- match(core_anchors$AnchorQuarter,mq)
  if (any(is.na(core_q_idx))) stopf("Core quarter missing from tv_probability for %s.",cc)
  for (eq in seq_len(K)) for (tt in seq_along(terms)) {
    probmat <- do.call(cbind,lapply(seq_len(NCHAINS),function(ch) as.numeric(parts[[cc]][[ch]]$tv_probability[,tt,eq])))
    pooled_prob <- rowMeans(probmat,na.rm=TRUE)
    tv_summary_rows[[length(tv_summary_rows)+1L]] <- data.frame(Country=cc,Equation=VARS[eq],Term=terms[tt],SampleMeanTVProbability=safe_mean(pooled_prob),SampleMedianTVProbability=safe_median(pooled_prob),SampleP95TVProbability=safe_p95(pooled_prob),SampleMaxTVProbability=safe_max(pooled_prob),CoreEventMeanTVProbability=safe_mean(pooled_prob[core_q_idx]),CoreEventMaxTVProbability=safe_max(pooled_prob[core_q_idx]),stringsAsFactors=FALSE)
    for (ee in seq_along(core_q_idx)) tv_core_rows[[length(tv_core_rows)+1L]] <- data.frame(Country=cc,Equation=VARS[eq],Term=terms[tt],EventID=core_events[ee],AnchorQuarter=core_anchors$AnchorQuarter[ee],MeanAcrossChainsTimeVariationProbability=pooled_prob[core_q_idx[ee]],stringsAsFactors=FALSE)
  }
}
tv_core <- do.call(rbind,tv_core_rows); tv_summary <- do.call(rbind,tv_summary_rows); tv_summary <- tv_summary[order(-tv_summary$CoreEventMaxTVProbability),,drop=FALSE]
write.csv(tv_core,file.path(OUT,"11_time_variation_probability_core_events.csv"),row.names=FALSE)
write.csv(tv_summary,file.path(OUT,"12_time_variation_probability_summary.csv"),row.names=FALSE)
tr_de_tv <- tv_core[tv_core$Country=="TR" & tv_core$Equation=="de" & tv_core$Term %in% c("de_L1","de_star_0","de_star_L1","gpr_0","gpr_L1"),,drop=FALSE]
write.csv(tr_de_tv,file.path(OUT,"13_tr_de_key_term_time_variation_probability.csv"),row.names=FALSE)

# ----- gate + classification ---------------------------------------------------
repro_max_gap <- max(reproduction$AbsoluteMedianGap,reproduction$AbsoluteP05Gap,reproduction$AbsoluteP95Gap,na.rm=TRUE)
max_alg_error <- max(decomp$MaxAlgebraReconciliationError,na.rm=TRUE)
if (!is.finite(repro_max_gap) || repro_max_gap>REPRO_TOL) stopf("08m reproduction failed: max gap %.12g",repro_max_gap)
if (!is.finite(max_alg_error) || max_alg_error>ALGEBRA_TOL) stopf("Algebra reconciliation failed: max error %.12g",max_alg_error)
mean_lag_share <- mean(decomp$MedianAbsLaggedGPRShare); mean_local_share_lag <- mean(decomp$MedianAbsLocalTRdeShareWithinLagBlock); mean_prop_share <- mean(decomp$MedianAbsPropagationShare)
source_classification <- if (mean_lag_share>=.80 && mean_local_share_lag>=.80) "LOCAL_TR_DE_GPR_L1_DOMINANT" else if (mean_lag_share>=.80) "LAGGED_GPR_NETWORK_MAPPED_DOMINANT" else if (mean_prop_share>=.80) "DYNAMIC_PROPAGATION_DOMINANT" else "MIXED_H1_SOURCE"
credible_coef_pairs <- sum(coef_pairs$CredibleDifference90,na.rm=TRUE); total_coef_pairs <- nrow(coef_pairs); credible_coef_pair_share <- credible_coef_pairs/total_coef_pairs
tr_de_credible_pairs <- sum(tr_de_pairs$CredibleDifference90,na.rm=TRUE); tr_de_total_pairs <- nrow(tr_de_pairs)
gpr_l1_tr_de <- tr_de_coef[tr_de_coef$Term=="gpr_L1",,drop=FALSE]; gpr_l1_median_range <- diff(range(gpr_l1_tr_de$median,finite=TRUE))
tr_de_tv_gpr_l1 <- tr_de_tv[tr_de_tv$Term=="gpr_L1",,drop=FALSE]; tr_de_gpr_l1_tv_mean <- safe_mean(tr_de_tv_gpr_l1$MeanAcrossChainsTimeVariationProbability); tr_de_gpr_l1_tv_max <- safe_max(tr_de_tv_gpr_l1$MeanAcrossChainsTimeVariationProbability)

# Top source tables
top_prop <- do.call(rbind,lapply(split(prop_sources,prop_sources$EventID),function(x) head(x[order(-x$MeanAbsContribution),,drop=FALSE],10)))
top_lag <- do.call(rbind,lapply(split(lag_sources,lag_sources$EventID),function(x) head(x[order(-x$MeanAbsContribution),,drop=FALSE],10)))
write.csv(top_prop,file.path(OUT,"02b_top10_propagation_sources_by_event.csv"),row.names=FALSE)
write.csv(top_lag,file.path(OUT,"03b_top10_lagged_gpr_sources_by_event.csv"),row.names=FALSE)

guardrails <- data.frame(Claim=c("08o reproduces accepted 08m TR/de/h1 summaries","h1 decomposition is exact draw by draw","lagged-GPR block dominates h1 anomaly","local TR/de gpr_L1 dominates lagged-GPR block","event coefficient changes use paired posterior draws","full T-by-K coefficient paths were retained","tv_probability is directly retained from latent-threshold D_dyn classification","08o changes model, MCMC, countries or W"),Supported=c(repro_max_gap<=REPRO_TOL,max_alg_error<=ALGEBRA_TOL,mean_lag_share>=.80,mean_local_share_lag>=.80,TRUE,FALSE,TRUE,FALSE),EvidenceOrLimit=c(sprintf("max p05/median/p95 gap %.12g",repro_max_gap),sprintf("max algebra error %.12g",max_alg_error),sprintf("mean median-absolute lag share %.6f",mean_lag_share),sprintf("mean local TR/de share in lag block %.6f",mean_local_share_lag),"event_coef retains same draw at all anchors within country/chain","R/60 intentionally stores compact event draws plus tv_probability, not full paths","R/60 constructs tv_probability from D_dyn high-variation classifications","read-only diagnostic"),stringsAsFactors=FALSE)
write.csv(guardrails,file.path(OUT,"14_claims_guardrail.csv"),row.names=FALSE)

gate <- data.frame(Status="DIAGNOSTIC_COMPLETE",Source08jRunID=SOURCE_08J_RUN,Source08kRunID=SOURCE_08K_RUN,Source08lRunID=SOURCE_08L_RUN,Source08mRunID=SOURCE_08M_RUN,Source08nRunID=SOURCE_08N_RUN,MainNetwork=MAIN_NETWORK,PosteriorFiles=nrow(lineage),PosteriorDrawsPerCoreEvent=DRAWS_TOTAL,CoreEvents=length(core_events),Target="TR__de__h1",GPRShockPct=GPR_SHOCK_PCT,Max08mReproductionGap=repro_max_gap,MaxAlgebraReconciliationError=max_alg_error,MeanMedianAbsPropagationShare=mean_prop_share,MeanMedianAbsLaggedGPRShare=mean_lag_share,MeanMedianAbsLocalTRdeShareWithinLaggedGPRBlock=mean_local_share_lag,H1SourceClassification=source_classification,CrediblePairedCoefficientDifferences90=credible_coef_pairs,TotalPairedCoefficientDifferences=total_coef_pairs,CrediblePairedCoefficientDifferenceShare90=credible_coef_pair_share,TRdeCrediblePairedCoefficientDifferences90=tr_de_credible_pairs,TRdeTotalPairedCoefficientDifferences=tr_de_total_pairs,TRdeGPRL1AcrossCoreEventMedianRange=gpr_l1_median_range,TRdeGPRL1CoreEventMeanTVProbability=tr_de_gpr_l1_tv_mean,TRdeGPRL1CoreEventMaxTVProbability=tr_de_gpr_l1_tv_max,FullCoefficientPathsRetained=FALSE,ReestimationPerformed=FALSE,ModelSpecificationChanged=FALSE,WeightMatrixChanged=FALSE,stringsAsFactors=FALSE)
write.csv(gate,file.path(OUT,"00_08o_gate.csv"),row.names=FALSE)

readme <- c("08o TR h=1 SOURCE + EVENT-COEFFICIENT TVP AUDIT","=================================================",sprintf("Target: %s",gate$Target),sprintf("Accepted posterior files: %d",gate$PosteriorFiles),sprintf("Posterior draws per event: %d",gate$PosteriorDrawsPerCoreEvent),sprintf("Max gap vs accepted 08m: %.12g",repro_max_gap),sprintf("Max algebra error: %.12g",max_alg_error),sprintf("Mean propagation abs share: %.6f",mean_prop_share),sprintf("Mean lagged-GPR abs share: %.6f",mean_lag_share),sprintf("Mean local TR/de share inside lag block: %.6f",mean_local_share_lag),sprintf("Source classification: %s",source_classification),sprintf("Paired event coefficient differences credible at 90%%: %d / %d (%.6f)",credible_coef_pairs,total_coef_pairs,credible_coef_pair_share),sprintf("TR/de paired differences credible at 90%%: %d / %d",tr_de_credible_pairs,tr_de_total_pairs),sprintf("TR/de gpr_L1 event-median range: %.10g",gpr_l1_median_range),sprintf("TR/de gpr_L1 mean core-event TV probability: %.6f",tr_de_gpr_l1_tv_mean),sprintf("TR/de gpr_L1 max core-event TV probability: %.6f",tr_de_gpr_l1_tv_max),"","No MCMC or model re-estimation is performed.","Full coefficient paths are not retained; paired event coefficients and tv_probability are available.")
writeLines(readme,file.path(OUT,"README_08o.txt"))
cat("08O TR H1 SOURCE + EVENT-COEFFICIENT TVP AUDIT: DIAGNOSTIC_COMPLETE\n")
print(gate)
