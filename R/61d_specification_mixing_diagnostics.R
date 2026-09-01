#!/usr/bin/env Rscript

# =============================================================================
# 61d_specification_mixing_diagnostics.R
#
# Diagnostic comparison for four formal TVP specifications after 08c showed
# that longer chains alone do not solve convergence.
#
# Expected specs:
#   BASELINE : SV=TRUE,  B1=2, B2=1,   kappa0=1e-7
#   NOSV     : SV=FALSE, B1=2, B2=1,   kappa0=1e-7
#   SHRINK   : SV=TRUE,  B1=5, B2=0.1, kappa0=1e-7
#   KAPPA    : SV=TRUE,  B1=2, B2=1,   kappa0=1e-4
#
# This script is diagnostic only. It does not choose a publication model based
# on fit or economics. It asks which parameter block is responsible for poor
# multi-chain mixing in US/CN.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("posterior", quietly = TRUE)) {
  stopf("Package 'posterior' is required.")
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

split_csv <- function(x, upper = TRUE) {
  z <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  z <- z[nzchar(z)]
  if (upper) z <- toupper(z)
  unique(z)
}

parameter_family <- function(x) {
  if (grepl("__sv_", x, fixed = TRUE)) return("SV")
  if (grepl("__omega__", x, fixed = TRUE)) return("OMEGA")
  if (grepl("__threshold__", x, fixed = TRUE)) return("THRESHOLD")
  if (grepl("__V0__", x, fixed = TRUE)) return("V0")
  if (grepl("__coef__", x, fixed = TRUE)) return("EVENT_COEF")
  "OTHER"
}

spec_from_path <- function(path) {
  p <- gsub("\\\\", "/", path)
  m <- regexec("specdiag_parts/([^/]+)/", p)
  z <- regmatches(p, m)[[1]]
  if (length(z) < 2L) return(NA_character_)
  toupper(z[2])
}

safe_mean <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
safe_median <- function(x) if (length(x)) stats::median(x, na.rm = TRUE) else NA_real_
safe_min <- function(x) if (length(x)) min(x, na.rm = TRUE) else NA_real_
safe_max <- function(x) if (length(x)) max(x, na.rm = TRUE) else NA_real_

PARTS_ROOT <- get_env_chr("FIN3_PARTS_ROOT", "posterior_parts")
DIAG_COUNTRIES <- split_csv(get_env_chr("FIN3_DIAG_COUNTRIES", "US,CN"))
SPECS <- split_csv(get_env_chr("FIN3_DIAG_SPECS", "BASELINE,NOSV,SHRINK,KAPPA"))
NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 1000))

RHAT_TARGET <- get_env_num("FIN3_RHAT_TARGET", 1.01)
RHAT_HARD <- get_env_num("FIN3_RHAT_HARD", 1.05)
ESS_TARGET <- get_env_num("FIN3_ESS_TARGET", 400)
ESS_HARD <- get_env_num("FIN3_ESS_HARD", 100)
MIN_RHAT_TARGET_SHARE <- get_env_num("FIN3_MIN_RHAT_TARGET_SHARE", 0.95)
MIN_ESS_TARGET_SHARE <- get_env_num("FIN3_MIN_ESS_TARGET_SHARE", 0.90)

if (!all(DIAG_COUNTRIES %in% COUNTRIES)) {
  stopf("Unknown diagnostic countries: %s", paste(setdiff(DIAG_COUNTRIES, COUNTRIES), collapse=", "))
}
if (NCHAINS < 2L) stopf("Need at least two chains.")
if (EXPECTED_STORED < 200L) stopf("Stored draws per chain is implausibly small.")

OUT <- file.path(RESULTS_DIR, "specification_mixing_diagnostic")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(files)) stopf("No formal TVP posterior RDS files found under %s.", PARTS_ROOT)

index_rows <- lapply(files, function(f) {
  z <- readRDS(f)
  sp <- spec_from_path(f)
  data.frame(
    File = f,
    Spec = sp,
    Country = toupper(as.character(z$meta$country)),
    Chain = as.integer(z$meta$chain),
    StoredDraws = as.integer(z$meta$stored_draws),
    Burn = as.integer(z$meta$burn),
    KeepInternal = as.integer(z$meta$keep_internal),
    SV = isTRUE(z$meta$sv_on),
    TVS = isTRUE(z$meta$tvs_on),
    B1 = as.numeric(z$meta$prior_B1),
    B2 = as.numeric(z$meta$prior_B2),
    Kappa0 = as.numeric(z$meta$kappa0),
    Network = as.character(z$meta$main_network),
    stringsAsFactors = FALSE
  )
})
parts_index <- do.call(rbind, index_rows)
parts_index <- parts_index[
  parts_index$Country %in% DIAG_COUNTRIES & parts_index$Spec %in% SPECS,
  , drop = FALSE
]

if (!nrow(parts_index)) stopf("No diagnostic posterior parts matched requested specs/countries.")
if (any(is.na(parts_index$Spec))) stopf("Could not infer Spec from at least one posterior path.")

expected <- expand.grid(
  Spec = SPECS,
  Country = DIAG_COUNTRIES,
  Chain = seq_len(NCHAINS),
  stringsAsFactors = FALSE
)
key <- function(d) paste(d$Spec, d$Country, d$Chain, sep = "||")
missing <- setdiff(key(expected), key(parts_index))
extra <- setdiff(key(parts_index), key(expected))

write.csv(parts_index, file.path(OUT, "00_parts_index.csv"), row.names = FALSE)
if (length(missing) || length(extra)) {
  write.csv(
    data.frame(Type=c(rep("MISSING", length(missing)), rep("EXTRA", length(extra))),
               Key=c(missing, extra), stringsAsFactors=FALSE),
    file.path(OUT, "00_grid_mismatch.csv"), row.names=FALSE
  )
  stopf("Posterior grid mismatch: missing=%d extra=%d", length(missing), length(extra))
}
if (anyDuplicated(key(parts_index))) stopf("Duplicate spec-country-chain posterior parts detected.")
if (any(parts_index$StoredDraws != EXPECTED_STORED)) stopf("Unexpected stored-draw count in at least one part.")
if (!all(parts_index$Network == MAIN_NETWORK)) stopf("Network mismatch across diagnostic parts.")

# Validate that each spec actually used the intended switches.
contract <- data.frame(
  Spec = c("BASELINE","NOSV","SHRINK","KAPPA"),
  SV = c(TRUE,FALSE,TRUE,TRUE),
  TVS = c(TRUE,TRUE,TRUE,TRUE),
  B1 = c(2,2,5,2),
  B2 = c(1,1,0.1,1),
  Kappa0 = c(1e-7,1e-7,1e-7,1e-4),
  stringsAsFactors = FALSE
)
contract <- contract[contract$Spec %in% SPECS, , drop=FALSE]
for (i in seq_len(nrow(contract))) {
  x <- parts_index[parts_index$Spec == contract$Spec[i], , drop=FALSE]
  if (!nrow(x)) stopf("No rows for spec %s.", contract$Spec[i])
  ok <- x$SV == contract$SV[i] & x$TVS == contract$TVS[i] &
    abs(x$B1-contract$B1[i]) < 1e-12 & abs(x$B2-contract$B2[i]) < 1e-12 &
    abs(x$Kappa0-contract$Kappa0[i]) < max(1e-14, abs(contract$Kappa0[i])*1e-10)
  if (!all(ok)) stopf("Specification contract mismatch for %s.", contract$Spec[i])
}
write.csv(contract, file.path(OUT, "00_specification_contract.csv"), row.names=FALSE)

param_rows <- list(); chain_rows <- list(); rr <- 0L; cr <- 0L

for (sp in SPECS) {
  for (cc in DIAG_COUNTRIES) {
    chains <- lapply(seq_len(NCHAINS), function(ch) {
      f <- parts_index$File[parts_index$Spec==sp & parts_index$Country==cc & parts_index$Chain==ch]
      if (length(f) != 1L) stopf("Missing/duplicate %s/%s chain %d", sp, cc, ch)
      z <- readRDS(f)
      m <- as.matrix(z$monitor)
      if (nrow(m) != EXPECTED_STORED) stopf("Unexpected monitor rows for %s/%s/%d", sp, cc, ch)
      m
    })
    ref <- colnames(chains[[1]])
    if (is.null(ref) || !length(ref)) stopf("No monitored parameters for %s/%s", sp, cc)
    if (!all(vapply(chains, function(x) identical(colnames(x), ref), logical(1)))) {
      stopf("Monitor columns differ across chains for %s/%s", sp, cc)
    }

    for (pp in seq_along(ref)) {
      par_name <- ref[pp]
      mat <- do.call(cbind, lapply(chains, function(x) x[,pp]))
      finite <- all(is.finite(mat))
      chain_sd <- if (finite) apply(mat, 2, stats::sd) else rep(NA_real_, NCHAINS)
      constant <- finite && all(chain_sd < 1e-12)

      if (!finite || constant) {
        rh <- eb <- et <- sep_ratio <- NA_real_
      } else {
        rh <- tryCatch(posterior::rhat(mat), error=function(e) NA_real_)
        eb <- tryCatch(posterior::ess_bulk(mat), error=function(e) NA_real_)
        et <- tryCatch(posterior::ess_tail(mat), error=function(e) NA_real_)
        mean_range <- diff(range(colMeans(mat)))
        sep_ratio <- mean_range / max(stats::median(chain_sd), 1e-12)
      }

      rr <- rr + 1L
      param_rows[[rr]] <- data.frame(
        Spec=sp, Country=cc, Parameter=par_name, Family=parameter_family(par_name),
        Finite=finite, ConstantAcrossAllChains=constant,
        Rhat=rh, ESS_Bulk=eb, ESS_Tail=et, MeanSeparationRatio=sep_ratio,
        PassRhatTarget=if (is.finite(rh)) rh<=RHAT_TARGET else NA,
        PassRhatHard=if (is.finite(rh)) rh<=RHAT_HARD else NA,
        PassBulkTarget=if (is.finite(eb)) eb>=ESS_TARGET else NA,
        PassTailTarget=if (is.finite(et)) et>=ESS_TARGET else NA,
        PassBulkHard=if (is.finite(eb)) eb>=ESS_HARD else NA,
        PassTailHard=if (is.finite(et)) et>=ESS_HARD else NA,
        stringsAsFactors=FALSE
      )

      for (ch in seq_len(NCHAINS)) {
        vals <- mat[,ch]
        cr <- cr + 1L
        chain_rows[[cr]] <- data.frame(
          Spec=sp, Country=cc, Parameter=par_name, Family=parameter_family(par_name), Chain=ch,
          Mean=if(all(is.finite(vals))) mean(vals) else NA_real_,
          SD=if(all(is.finite(vals))) stats::sd(vals) else NA_real_,
          Q05=if(all(is.finite(vals))) unname(stats::quantile(vals,.05)) else NA_real_,
          Q50=if(all(is.finite(vals))) unname(stats::quantile(vals,.50)) else NA_real_,
          Q95=if(all(is.finite(vals))) unname(stats::quantile(vals,.95)) else NA_real_,
          stringsAsFactors=FALSE
        )
      }
    }
    rm(chains); invisible(gc())
  }
}

diagnostics <- do.call(rbind, param_rows)
chain_summary <- do.call(rbind, chain_rows)
write.csv(diagnostics, file.path(OUT, "01_parameter_diagnostics.csv"), row.names=FALSE)
write.csv(chain_summary, file.path(OUT, "07_chain_parameter_summary.csv"), row.names=FALSE)

active <- diagnostics[
  diagnostics$Finite & !diagnostics$ConstantAcrossAllChains &
    is.finite(diagnostics$Rhat) & is.finite(diagnostics$ESS_Bulk) & is.finite(diagnostics$ESS_Tail),
  , drop=FALSE
]
if (!nrow(active)) stopf("No active diagnostics.")

summarize_block <- function(x) {
  data.frame(
    ActiveDiagnostics=nrow(x),
    MaxRhat=max(x$Rhat), MedianRhat=stats::median(x$Rhat),
    RhatTargetShare=mean(x$Rhat<=RHAT_TARGET), RhatHardShare=mean(x$Rhat<=RHAT_HARD),
    MinESSBulk=min(x$ESS_Bulk), MinESSTail=min(x$ESS_Tail),
    ESSBulkTargetShare=mean(x$ESS_Bulk>=ESS_TARGET), ESSTailTargetShare=mean(x$ESS_Tail>=ESS_TARGET),
    ESSBulkHardShare=mean(x$ESS_Bulk>=ESS_HARD), ESSTailHardShare=mean(x$ESS_Tail>=ESS_HARD),
    HardFailCount=sum(x$Rhat>RHAT_HARD | x$ESS_Bulk<ESS_HARD | x$ESS_Tail<ESS_HARD),
    MedianMeanSeparationRatio=safe_median(x$MeanSeparationRatio),
    MaxMeanSeparationRatio=safe_max(x$MeanSeparationRatio),
    stringsAsFactors=FALSE
  )
}

country_spec <- do.call(rbind, lapply(SPECS, function(sp) {
  do.call(rbind, lapply(DIAG_COUNTRIES, function(cc) {
    x <- active[active$Spec==sp & active$Country==cc,,drop=FALSE]
    cbind(data.frame(Spec=sp, Country=cc, stringsAsFactors=FALSE), summarize_block(x))
  }))
}))
write.csv(country_spec, file.path(OUT, "02_country_spec_summary.csv"), row.names=FALSE)

family_spec <- do.call(rbind, lapply(SPECS, function(sp) {
  fams <- sort(unique(active$Family[active$Spec==sp]))
  do.call(rbind, lapply(fams, function(ff) {
    x <- active[active$Spec==sp & active$Family==ff,,drop=FALSE]
    cbind(data.frame(Spec=sp, Family=ff, stringsAsFactors=FALSE), summarize_block(x))
  }))
}))
write.csv(family_spec, file.path(OUT, "03_family_spec_summary.csv"), row.names=FALSE)

overall <- do.call(rbind, lapply(SPECS, function(sp) {
  x <- active[active$Spec==sp,,drop=FALSE]
  cbind(data.frame(Spec=sp, stringsAsFactors=FALSE), summarize_block(x))
}))

overall$HardGatePass <- with(overall,
  MaxRhat <= RHAT_HARD & MinESSBulk >= ESS_HARD & MinESSTail >= ESS_HARD)
overall$PublicationTargetPass <- with(overall,
  HardGatePass & RhatTargetShare >= MIN_RHAT_TARGET_SHARE &
    ESSBulkTargetShare >= MIN_ESS_TARGET_SHARE & ESSTailTargetShare >= MIN_ESS_TARGET_SHARE)

overall <- overall[order(
  !overall$PublicationTargetPass,
  !overall$HardGatePass,
  overall$HardFailCount,
  -overall$RhatTargetShare,
  overall$MaxRhat,
  -pmin(overall$ESSBulkTargetShare, overall$ESSTailTargetShare)
),,drop=FALSE]
overall$DiagnosticRank <- seq_len(nrow(overall))
write.csv(overall, file.path(OUT, "04_overall_spec_ranking.csv"), row.names=FALSE)

# Baseline deltas by family make the mechanism visible.
base_family <- family_spec[family_spec$Spec=="BASELINE",,drop=FALSE]
comparison <- do.call(rbind, lapply(setdiff(SPECS, "BASELINE"), function(sp) {
  x <- family_spec[family_spec$Spec==sp,,drop=FALSE]
  m <- merge(x, base_family, by="Family", suffixes=c("_Spec","_Baseline"), all=FALSE)
  if (!nrow(m)) return(NULL)
  data.frame(
    Spec=sp, Family=m$Family,
    DeltaMedianRhat=m$MedianRhat_Spec-m$MedianRhat_Baseline,
    DeltaMaxRhat=m$MaxRhat_Spec-m$MaxRhat_Baseline,
    DeltaRhatTargetShare=m$RhatTargetShare_Spec-m$RhatTargetShare_Baseline,
    DeltaMinESSBulk=m$MinESSBulk_Spec-m$MinESSBulk_Baseline,
    DeltaMinESSTail=m$MinESSTail_Spec-m$MinESSTail_Baseline,
    DeltaHardFailCount=m$HardFailCount_Spec-m$HardFailCount_Baseline,
    stringsAsFactors=FALSE
  )
}))
if (!is.null(comparison)) write.csv(comparison, file.path(OUT, "05_family_improvement_vs_baseline.csv"), row.names=FALSE)

worst <- active[order(active$Spec, -active$Rhat, active$ESS_Bulk, active$ESS_Tail),,drop=FALSE]
worst <- do.call(rbind, lapply(SPECS, function(sp) head(worst[worst$Spec==sp,,drop=FALSE], 40L)))
write.csv(worst, file.path(OUT, "06_worst_parameters_by_spec.csv"), row.names=FALSE)

best <- overall[1,,drop=FALSE]
status <- if (best$PublicationTargetPass) {
  "SPECIFICATION_CANDIDATE_PASSES_SCREEN"
} else if (best$HardGatePass) {
  "BEST_SPEC_REMOVES_HARD_FAILURES_ONLY"
} else {
  "NO_SPEC_REMOVES_HARD_FAILURES"
}

mechanism_note <- switch(best$Spec,
  "BASELINE" = "No tested simplification beats the baseline; deeper reparameterization or direct V0-prior work is indicated.",
  "NOSV" = "Turning off stochastic volatility gives the best mixing; SV-state interaction is the leading suspect.",
  "SHRINK" = "Stronger state-innovation shrinkage gives the best mixing; dynamic-state variance prior interaction is the leading suspect.",
  "KAPPA" = "Relaxing the lower-regime kappa0 gives the best mixing; near-degenerate low-variation states are the leading suspect.",
  "Best diagnostic specification identified."
)

gate <- data.frame(
  Status=status,
  BestSpec=best$Spec,
  BestSpecRank=best$DiagnosticRank,
  HardGatePass=best$HardGatePass,
  PublicationTargetPass=best$PublicationTargetPass,
  MaxRhat=best$MaxRhat,
  RhatTargetShare=best$RhatTargetShare,
  MinESSBulk=best$MinESSBulk,
  MinESSTail=best$MinESSTail,
  ESSBulkTargetShare=best$ESSBulkTargetShare,
  ESSTailTargetShare=best$ESSTailTargetShare,
  HardFailCount=best$HardFailCount,
  MechanismNote=mechanism_note,
  stringsAsFactors=FALSE
)
write.csv(gate, file.path(OUT, "00_specification_diagnostic_gate.csv"), row.names=FALSE)

readme <- c(
  sprintf("08d SPECIFICATION MIXING DIAGNOSTIC: %s", status),
  "============================================================",
  sprintf("Countries: %s", paste(DIAG_COUNTRIES, collapse=", ")),
  sprintf("Specs: %s", paste(SPECS, collapse=", ")),
  sprintf("Chains/spec/country: %d", NCHAINS),
  sprintf("Stored draws/chain: %d", EXPECTED_STORED),
  "",
  sprintf("Best diagnostic spec: %s", best$Spec),
  sprintf("Max Rhat: %.6f", best$MaxRhat),
  sprintf("Rhat <= %.3f share: %.6f", RHAT_TARGET, best$RhatTargetShare),
  sprintf("Min bulk/tail ESS: %.2f / %.2f", best$MinESSBulk, best$MinESSTail),
  sprintf("Hard failures: %d", best$HardFailCount),
  sprintf("Interpretation: %s", mechanism_note),
  "",
  "IMPORTANT:",
  "This is a sampler/specification diagnostic, not an automatic publication-model selector.",
  "Do not switch the substantive baseline only because one diagnostic spec mixes better.",
  "Use 03_family_spec_summary.csv and 05_family_improvement_vs_baseline.csv to identify",
  "whether V0, event coefficients, SV, or threshold-related blocks actually improve."
)
writeLines(readme, file.path(OUT, "README_08d_specification_diagnostic.txt"))
cat(paste(readme, collapse="\n"), "\n")
