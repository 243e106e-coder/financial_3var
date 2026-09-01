#!/usr/bin/env Rscript

# =============================================================================
# IRF SIGNIFICANCE / SIGN / PERSISTENCE AUDIT
# Financial 3-variable TVP-GVAR
#
# Response variables:
#   r   = RATE_LEVEL in baseline, or Delta RATE_LEVEL in the rate-diff robustness
#   de  = REER_DLOG
#   deq = EQ_RETURN
#
# This version is transformation-aware:
# - It can require the transformation audit gate.
# - It keeps raw IRF path sums separate from genuine cumulative level effects.
# - It does not assume EQ_RETURN is a log return unless explicitly requested.
# =============================================================================

source("R/00_config.R")

IRF_FILE <- Sys.getenv(
  "FIN3_IRF_SUMMARY",
  file.path(RESULTS_DIR, "formal_irf", "irf_posterior_summary.csv")
)

OUT <- file.path(RESULTS_DIR, "irf_audit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REER_POSITIVE_MEANING <- Sys.getenv(
  "FIN3_REER_POSITIVE_MEANING",
  "REER appreciation direction"
)

RATE_MODE <- tolower(trimws(Sys.getenv("FIN3_RATE_MODE", "level")))
if (!RATE_MODE %in% c("level", "difference")) {
  stopf("FIN3_RATE_MODE must be 'level' or 'difference'.")
}

EQ_RETURN_MODE <- tolower(trimws(Sys.getenv("FIN3_EQ_RETURN_MODE", "simple_return")))
if (!EQ_RETURN_MODE %in% c("simple_return", "log_return")) {
  stopf("FIN3_EQ_RETURN_MODE must be 'simple_return' or 'log_return'.")
}

REQUIRE_TRANSFORM_AUDIT <- tolower(trimws(
  Sys.getenv("FIN3_REQUIRE_TRANSFORM_AUDIT", "0")
)) %in% c("1", "true", "yes", "y")

TRANSFORM_GATE <- file.path(
  RESULTS_DIR,
  "variable_transform",
  "11_transformation_integrity_gate.csv"
)

if (REQUIRE_TRANSFORM_AUDIT) {
  if (!file.exists(TRANSFORM_GATE)) {
    stopf(
      "Transformation audit is required but missing: %s. Run R/15_variable_transformation_audit.R first.",
      TRANSFORM_GATE
    )
  }
  tg <- read.csv(TRANSFORM_GATE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Check", "Pass") %in% names(tg))) {
    stopf("Malformed transformation audit gate: %s", TRANSFORM_GATE)
  }
  if (!all(tg$Pass)) {
    stopf("Transformation audit integrity gate is not READY.")
  }
}

if (!file.exists(IRF_FILE)) {
  stopf(
    paste0(
      "Missing formal IRF posterior summary: %s\n",
      "Set FIN3_IRF_SUMMARY to the correct CSV path, or run the formal TVP-IRF stage first."
    ),
    IRF_FILE
  )
}

d <- read.csv(IRF_FILE, stringsAsFactors = FALSE, check.names = FALSE)

# -------------------------------------------------------------------------
# Conservative column-name resolver
# -------------------------------------------------------------------------

pick_col <- function(data, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) return(hit[1])
  if (required) {
    stopf(
      "IRF file is missing required field. Expected one of: %s",
      paste(candidates, collapse = ", ")
    )
  }
  NA_character_
}

map <- list(
  EventID = pick_col(d, c("EventID", "event_id", "event")),
  EventSet = pick_col(d, c("EventSet", "event_set"), required = FALSE),
  EventLabel = pick_col(d, c("EventLabel", "event_label"), required = FALSE),
  ShockFamily = pick_col(d, c("ShockFamily", "shock_family"), required = FALSE),
  AnchorType = pick_col(d, c("AnchorType", "anchor_type")),
  AnchorQuarter = pick_col(d, c("AnchorQuarter", "anchor_quarter", "quarter"), required = FALSE),
  Country = pick_col(d, c("Country", "country", "economy")),
  ResponseVariable = pick_col(d, c("ResponseVariable", "response_variable", "variable", "var")),
  Horizon = pick_col(d, c("Horizon", "horizon", "h")),
  Median = pick_col(d, c("median", "Median", "posterior_median", "irf_median")),
  P05 = pick_col(d, c("p05", "P05", "q05", "lower90", "lower_90")),
  P95 = pick_col(d, c("p95", "P95", "q95", "upper90", "upper_90")),
  P16 = pick_col(d, c("p16", "P16", "q16", "lower68", "lower_68"), required = FALSE),
  P84 = pick_col(d, c("p84", "P84", "q84", "upper68", "upper_68"), required = FALSE),
  ProbPositive = pick_col(
    d,
    c("prob_positive", "ProbPositive", "sign_probability", "posterior_prob_positive"),
    required = FALSE
  )
)

canon <- data.frame(
  EventID = as.character(d[[map$EventID]]),
  EventSet = if (!is.na(map$EventSet)) as.character(d[[map$EventSet]]) else NA_character_,
  EventLabel = if (!is.na(map$EventLabel)) as.character(d[[map$EventLabel]]) else NA_character_,
  ShockFamily = if (!is.na(map$ShockFamily)) as.character(d[[map$ShockFamily]]) else NA_character_,
  AnchorType = as.character(d[[map$AnchorType]]),
  AnchorQuarter = if (!is.na(map$AnchorQuarter)) as.character(d[[map$AnchorQuarter]]) else NA_character_,
  Country = toupper(trimws(as.character(d[[map$Country]]))),
  ResponseVariable = tolower(trimws(as.character(d[[map$ResponseVariable]]))),
  Horizon = num(d[[map$Horizon]]),
  median = num(d[[map$Median]]),
  p05 = num(d[[map$P05]]),
  p95 = num(d[[map$P95]]),
  p16 = if (!is.na(map$P16)) num(d[[map$P16]]) else NA_real_,
  p84 = if (!is.na(map$P84)) num(d[[map$P84]]) else NA_real_,
  prob_positive = if (!is.na(map$ProbPositive)) num(d[[map$ProbPositive]]) else NA_real_,
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------

if (any(!canon$Country %in% COUNTRIES)) {
  bad <- unique(canon$Country[!canon$Country %in% COUNTRIES])
  stopf("Unknown country codes in IRF summary: %s", paste(bad, collapse = ", "))
}

if (any(!canon$ResponseVariable %in% VARS)) {
  bad <- unique(canon$ResponseVariable[!canon$ResponseVariable %in% VARS])
  stopf(
    "Unknown response variables in IRF summary: %s. Expected exactly: %s",
    paste(bad, collapse = ", "),
    paste(VARS, collapse = ", ")
  )
}

if (any(!is.finite(canon$Horizon))) stopf("Non-finite IRF horizon.")
if (any(canon$Horizon < 0)) stopf("Negative IRF horizon.")
if (any(!is.finite(canon$median))) stopf("Non-finite posterior median.")
if (any(!is.finite(canon$p05)) || any(!is.finite(canon$p95))) {
  stopf("Non-finite p05/p95 credible-band values.")
}
if (any(canon$p05 > canon$p95)) stopf("Found p05 > p95.")
if (any(canon$median < canon$p05 | canon$median > canon$p95)) {
  stopf("Posterior median falls outside p05/p95 for at least one row.")
}

key <- paste(
  canon$EventID,
  canon$AnchorType,
  canon$Country,
  canon$ResponseVariable,
  canon$Horizon,
  sep = "||"
)

if (anyDuplicated(key)) {
  stopf("Duplicate event-anchor-country-variable-horizon rows in IRF summary.")
}

# -------------------------------------------------------------------------
# Transformation-aware variable interpretation
# -------------------------------------------------------------------------

rate_meaning <- if (RATE_MODE == "level") {
  "Interest-rate level increase"
} else {
  "Increase in quarterly interest-rate change (Delta r)"
}

rate_scale <- if (RATE_MODE == "level") {
  "DIRECT_LEVEL_DEVIATION"
} else {
  "FIRST_DIFFERENCE__CUMULATE_FOR_RATE_LEVEL_CHANGE"
}

eq_scale <- if (EQ_RETURN_MODE == "log_return") {
  "LOG_RETURN__DRAW_LEVEL_CUMULATION_NEEDED_FOR_EXACT_POSTERIOR_LEVEL_EFFECT"
} else {
  "RETURN__DO_NOT_TREAT_MEDIAN_SUM_AS_LOG_PRICE_EFFECT"
}

var_meaning <- c(
  r = rate_meaning,
  de = REER_POSITIVE_MEANING,
  deq = "Equity-return increase"
)

effect_scale <- c(
  r = rate_scale,
  de = "DLOG__CUMULATE_FOR_REER_LOG_LEVEL_EFFECT",
  deq = eq_scale
)

canon$PositiveDirectionMeaning <- unname(var_meaning[canon$ResponseVariable])
canon$EffectScale <- unname(effect_scale[canon$ResponseVariable])

canon$RequiresCumulationForLevelInterpretation <- ifelse(
  canon$ResponseVariable == "de",
  TRUE,
  ifelse(
    canon$ResponseVariable == "r" & RATE_MODE == "difference",
    TRUE,
    ifelse(
      canon$ResponseVariable == "deq" & EQ_RETURN_MODE == "log_return",
      TRUE,
      FALSE
    )
  )
)

scale_manifest <- data.frame(
  ResponseVariable = VARS,
  SourceField = unname(SOURCE_SUFFIX[VARS]),
  RuntimeMode = c(RATE_MODE, "fixed_upstream", EQ_RETURN_MODE),
  PositiveDirectionMeaning = unname(var_meaning[VARS]),
  EffectScale = unname(effect_scale[VARS]),
  RequiresCumulationForLevelInterpretation = c(
    RATE_MODE == "difference",
    TRUE,
    EQ_RETURN_MODE == "log_return"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scale_manifest,
  file.path(OUT, "00a_irf_variable_scale_manifest.csv"),
  row.names = FALSE
)

# 90% posterior credible classification.
canon$SignClass90 <- ifelse(
  canon$p05 > 0,
  "POSITIVE_CREDIBLE_90",
  ifelse(
    canon$p95 < 0,
    "NEGATIVE_CREDIBLE_90",
    "CROSSES_ZERO_90"
  )
)

canon$MedianSign <- ifelse(
  canon$median > 0,
  "POSITIVE",
  ifelse(canon$median < 0, "NEGATIVE", "ZERO")
)

canon$Credible90 <- canon$SignClass90 != "CROSSES_ZERO_90"
canon$BandWidth90 <- canon$p95 - canon$p05
canon$AbsoluteMedian <- abs(canon$median)

write.csv(
  canon,
  file.path(OUT, "00_irf_posterior_summary_audited.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Path-level audit
# -------------------------------------------------------------------------

group_key <- interaction(
  canon$EventID,
  canon$AnchorType,
  canon$Country,
  canon$ResponseVariable,
  drop = TRUE,
  lex.order = TRUE
)

split_groups <- split(canon, group_key)

audit_one <- function(x) {
  x <- x[order(x$Horizon), , drop = FALSE]

  h0 <- x[x$Horizon == 0, , drop = FALSE]
  if (!nrow(h0)) {
    stopf(
      "Missing horizon 0 for %s / %s / %s / %s",
      x$EventID[1], x$AnchorType[1], x$Country[1], x$ResponseVariable[1]
    )
  }

  peak_i <- which.max(abs(x$median))
  peak <- x[peak_i, , drop = FALSE]

  short <- x[x$Horizon >= 0 & x$Horizon <= 4, , drop = FALSE]
  medium <- x[x$Horizon >= 0 & x$Horizon <= 8, , drop = FALSE]

  cred <- x[x$Credible90, , drop = FALSE]
  credible_signs <- if (nrow(cred)) sign(cred$median) else numeric(0)
  credible_signs <- credible_signs[credible_signs != 0]

  sign_reversal <- FALSE
  if (length(credible_signs) >= 2L) {
    sign_reversal <- any(diff(credible_signs) != 0)
  }

  med_signs <- sign(x$median)
  med_signs <- med_signs[med_signs != 0]
  median_sign_reversal <- FALSE
  if (length(med_signs) >= 2L) {
    median_sign_reversal <- any(diff(med_signs) != 0)
  }

  # Path sums are retained for comparability, but they are NOT automatically
  # interpreted as cumulative level effects.
  pathsum04 <- sum(short$median, na.rm = TRUE)
  pathsum08 <- sum(medium$median, na.rm = TRUE)

  vv <- x$ResponseVariable[1]

  cumulative_level04 <- NA_real_
  cumulative_level08 <- NA_real_
  cumulative_level_interpretation <- "NOT_AVAILABLE_FROM_SUMMARY_MEDIANS"

  if (vv == "de") {
    cumulative_level04 <- pathsum04
    cumulative_level08 <- pathsum08
    cumulative_level_interpretation <- "SUM_OF_REER_DLOG_MEDIANS__APPROX_LOG_LEVEL_EFFECT"
  } else if (vv == "r" && RATE_MODE == "difference") {
    cumulative_level04 <- pathsum04
    cumulative_level08 <- pathsum08
    cumulative_level_interpretation <- "SUM_OF_DELTA_R_MEDIANS__RATE_LEVEL_CHANGE_APPROX"
  } else if (vv == "r" && RATE_MODE == "level") {
    cumulative_level_interpretation <- "RAW_IRF_ALREADY_RATE_LEVEL_DEVIATION__DO_NOT_CUMULATE_FOR_LEVEL"
  } else if (vv == "deq" && EQ_RETURN_MODE == "log_return") {
    cumulative_level04 <- pathsum04
    cumulative_level08 <- pathsum08
    cumulative_level_interpretation <- "LOG_RETURN_MODE__SUMMARY_MEDIAN_SUM_ONLY_APPROXIMATE__PREFER_DRAW_LEVEL_CUMULATION"
  } else if (vv == "deq") {
    cumulative_level_interpretation <- "SIMPLE_RETURN_MODE__DO_NOT_CALL_PATH_SUM_A_LOG_PRICE_EFFECT"
  }

  prob_pos_h0 <- h0$prob_positive[1]
  posterior_direction_h0 <- NA_character_

  if (is.finite(prob_pos_h0)) {
    posterior_direction_h0 <- ifelse(
      prob_pos_h0 >= 0.95,
      "P_POS_GE_0.95",
      ifelse(
        prob_pos_h0 <= 0.05,
        "P_NEG_GE_0.95",
        "DIRECTION_UNCERTAIN"
      )
    )
  }

  data.frame(
    EventID = x$EventID[1],
    EventSet = x$EventSet[1],
    EventLabel = x$EventLabel[1],
    ShockFamily = x$ShockFamily[1],
    AnchorType = x$AnchorType[1],
    AnchorQuarter = x$AnchorQuarter[1],
    Country = x$Country[1],
    ResponseVariable = vv,
    PositiveDirectionMeaning = x$PositiveDirectionMeaning[1],
    EffectScale = x$EffectScale[1],

    ImpactMedian = h0$median[1],
    ImpactP05 = h0$p05[1],
    ImpactP95 = h0$p95[1],
    ImpactMedianSign = h0$MedianSign[1],
    ImpactSignClass90 = h0$SignClass90[1],
    ImpactCredible90 = h0$Credible90[1],
    ImpactProbPositive = prob_pos_h0,
    ImpactPosteriorDirection = posterior_direction_h0,

    PeakMedian = peak$median[1],
    PeakAbsMedian = abs(peak$median[1]),
    PeakHorizon = peak$Horizon[1],
    PeakSignClass90 = peak$SignClass90[1],
    PeakCredible90 = peak$Credible90[1],

    SignificantCount_h0_4 = sum(short$Credible90),
    SignificantShare_h0_4 = mean(short$Credible90),
    SignificantCount_h0_8 = sum(medium$Credible90),
    SignificantShare_h0_8 = mean(medium$Credible90),

    # Legacy names retained; they are path sums, not automatically level effects.
    CumulativeMedian_h0_4 = pathsum04,
    CumulativeMedian_h0_8 = pathsum08,
    PathSumMedian_h0_4 = pathsum04,
    PathSumMedian_h0_8 = pathsum08,

    CumulativeLevelApprox_h0_4 = cumulative_level04,
    CumulativeLevelApprox_h0_8 = cumulative_level08,
    CumulativeLevelInterpretation = cumulative_level_interpretation,

    CredibleSignReversal = sign_reversal,
    MedianPathSignReversal = median_sign_reversal,
    EverCredible90 = any(x$Credible90),
    AllHorizonsCrossZero90 = !any(x$Credible90),

    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, lapply(split_groups, audit_one))
row.names(audit) <- NULL

# -------------------------------------------------------------------------
# Diagnostic flags
# -------------------------------------------------------------------------

audit$Flag_NoCredibleResponse <- audit$AllHorizonsCrossZero90
audit$Flag_ImpactNotCredible <- !audit$ImpactCredible90
audit$Flag_CredibleSignReversal <- audit$CredibleSignReversal
audit$Flag_PeakNotCredible <- !audit$PeakCredible90

# No universal "opposite to theory" flag is imposed.
audit$AnyDiagnosticFlag <- (
  audit$Flag_NoCredibleResponse |
  audit$Flag_CredibleSignReversal |
  audit$Flag_PeakNotCredible
)

write.csv(
  audit,
  file.path(OUT, "01_irf_path_significance_audit.csv"),
  row.names = FALSE
)

flagged <- audit[audit$AnyDiagnosticFlag, , drop = FALSE]
write.csv(
  flagged,
  file.path(OUT, "02_irf_flagged_cases.csv"),
  row.names = FALSE
)

impact_table <- audit[, c(
  "EventSet", "EventID", "AnchorType", "AnchorQuarter",
  "Country", "ResponseVariable", "PositiveDirectionMeaning", "EffectScale",
  "ImpactMedian", "ImpactP05", "ImpactP95",
  "ImpactSignClass90", "ImpactProbPositive",
  "PeakMedian", "PeakHorizon", "PeakSignClass90",
  "SignificantCount_h0_4",
  "PathSumMedian_h0_4",
  "CumulativeLevelApprox_h0_4",
  "CumulativeLevelInterpretation",
  "CredibleSignReversal"
)]

write.csv(
  impact_table,
  file.path(OUT, "03_irf_impact_and_peak_table.csv"),
  row.names = FALSE
)

var_summary <- do.call(rbind, lapply(VARS, function(v) {
  z <- audit[audit$ResponseVariable == v, , drop = FALSE]
  data.frame(
    ResponseVariable = v,
    PositiveDirectionMeaning = unname(var_meaning[v]),
    EffectScale = unname(effect_scale[v]),
    Paths = nrow(z),
    ImpactPositiveCredible90 = sum(z$ImpactSignClass90 == "POSITIVE_CREDIBLE_90"),
    ImpactNegativeCredible90 = sum(z$ImpactSignClass90 == "NEGATIVE_CREDIBLE_90"),
    ImpactCrossesZero90 = sum(z$ImpactSignClass90 == "CROSSES_ZERO_90"),
    EverCredible90 = sum(z$EverCredible90),
    CredibleSignReversal = sum(z$CredibleSignReversal),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  var_summary,
  file.path(OUT, "04_irf_variable_direction_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Integrity gate: numeric/data integrity only, not desired economic signs.
# -------------------------------------------------------------------------

audit_integrity <- data.frame(
  Check = c(
    "InputRowsFinite",
    "VariablesExactly_r_de_deq",
    "Horizon0AvailableEveryPath",
    "CredibleBandsOrdered",
    "NoDuplicateIRFRows",
    "RateModeValid",
    "EquityReturnModeValid",
    "TransformationAuditReadyWhenRequired"
  ),
  Pass = c(
    TRUE,
    setequal(sort(unique(canon$ResponseVariable)), sort(VARS)),
    all(vapply(split_groups, function(x) any(x$Horizon == 0), logical(1))),
    all(canon$p05 <= canon$median & canon$median <= canon$p95),
    !anyDuplicated(key),
    RATE_MODE %in% c("level", "difference"),
    EQ_RETURN_MODE %in% c("simple_return", "log_return"),
    if (REQUIRE_TRANSFORM_AUDIT) {
      exists("tg") && all(tg$Pass)
    } else {
      TRUE
    }
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit_integrity,
  file.path(OUT, "05_irf_audit_integrity_gate.csv"),
  row.names = FALSE
)

readme <- c(
  "IRF SIGNIFICANCE, DIRECTION, AND TRANSFORMATION-AWARE AUDIT",
  "===========================================================",
  "",
  sprintf("Input: %s", IRF_FILE),
  sprintf("Rate mode: %s", RATE_MODE),
  sprintf("Equity return mode: %s", EQ_RETURN_MODE),
  sprintf("Transformation audit required: %s", REQUIRE_TRANSFORM_AUDIT),
  "",
  "Repository response variables:",
  if (RATE_MODE == "level") {
    "- r   = RATE_LEVEL (direct rate-level response)"
  } else {
    "- r   = Delta RATE_LEVEL in robustness model (cumulate for a rate-level change)"
  },
  paste0("- de  = REER_DLOG (positive labelled: ", REER_POSITIVE_MEANING, ")"),
  if (EQ_RETURN_MODE == "log_return") {
    "- deq = EQ_RETURN, explicitly treated as log return by runtime setting"
  } else {
    "- deq = EQ_RETURN, treated conservatively as a generic/simple return"
  },
  "- CPI is NOT part of this 3-variable model.",
  "",
  "90% posterior credibility rule:",
  "- p05 > 0  -> POSITIVE_CREDIBLE_90",
  "- p95 < 0  -> NEGATIVE_CREDIBLE_90",
  "- otherwise -> CROSSES_ZERO_90",
  "",
  "Transformation discipline:",
  "- Raw IRF path sums are saved separately from cumulative level interpretations.",
  "- r level: raw IRF is already a level deviation; do not cumulate merely to obtain a rate-level response.",
  "- Delta r: cumulation is required to recover a rate-level change.",
  "- de = REER_DLOG: cumulation gives a cumulative REER log-level effect.",
  "- deq: do not interpret a sum as a log-price effect unless EQ_RETURN is explicitly verified as a log return.",
  "- Exact posterior cumulation should be done draw-by-draw before quantiles when draw-level IRFs are available.",
  "",
  "Interpretation discipline:",
  "- A positive/negative posterior median is not called credible if the 90% band crosses zero.",
  "- An opposite sign is not automatically treated as a coding error.",
  "- A credible sign reversal is flagged for substantive review.",
  "- No universal safe-haven sign is imposed on r, de, or deq.",
  "",
  "Key outputs:",
  "- 00_irf_posterior_summary_audited.csv",
  "- 00a_irf_variable_scale_manifest.csv",
  "- 01_irf_path_significance_audit.csv",
  "- 02_irf_flagged_cases.csv",
  "- 03_irf_impact_and_peak_table.csv",
  "- 04_irf_variable_direction_summary.csv",
  "- 05_irf_audit_integrity_gate.csv"
)

writeLines(readme, file.path(OUT, "README_irf_audit.txt"))

if (!all(audit_integrity$Pass)) {
  stopf(
    "IRF audit integrity gate failed. Inspect results/irf_audit/05_irf_audit_integrity_gate.csv"
  )
}

msg("IRF audit complete.")
msg("Input rows: %d", nrow(canon))
msg("IRF paths: %d", nrow(audit))
msg("Flagged paths: %d", nrow(flagged))
msg("Variables: %s", paste(sort(unique(canon$ResponseVariable)), collapse = ", "))
msg("Rate mode: %s", RATE_MODE)
msg("IRF audit integrity gate: READY")
