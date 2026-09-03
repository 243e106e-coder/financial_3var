#!/usr/bin/env Rscript

# =============================================================================
# 08n EVENT HETEROGENEITY + SAFE-HAVEN + TR ANOMALY AUDIT
# Financial 3-variable Delta-r TVP-GVAR
#
# IMPORTANT CONTRACT
# - This script does NOT re-estimate the model.
# - It reads only the accepted 08m artifact outputs.
# - Event-to-event comparisons are descriptive because 08m stores posterior
#   summaries, not paired draw-level event differences.
# - "TVP time variation" below is an IRF heterogeneity proxy across event
#   anchors; it is NOT a direct audit of latent TVP coefficient trajectories.
# =============================================================================

options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("FIN3_08M_ROOT", "source_08m")
OUT <- Sys.getenv(
  "FIN3_08N_OUT",
  file.path("results", "event_heterogeneity_safe_haven_tr_audit")
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PLOT_DIR <- file.path(OUT, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

SOURCE_RUN_ID <- Sys.getenv("FIN3_08M_RUN_ID", "")
SOURCE_HEAD_SHA <- Sys.getenv("FIN3_08M_HEAD_SHA", "")

EXPECTED_COUNTRIES <- c(
  "AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA"
)
EXPECTED_VARS <- c("r","de","deq")
EXPECTED_CORE_EVENTS <- c(
  "LEHMAN_2008",
  "SOVEREIGN_STRESS_2011",
  "BREXIT_2016",
  "COVID_2020",
  "RUSSIA_UKRAINE_2022",
  "IRAN_ISRAEL_2024"
)
EXPECTED_RAW_HORIZONS <- 0:12
SAFE_HAVEN_HORIZONS <- c(1L, 4L, 8L, 12L)
FOCUS_SAFE_HAVEN <- c("JP","CH","US")

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

read_required_csv <- function(path, required_cols = character()) {
  if (!file.exists(path)) stopf("Missing required 08m file: %s", path)
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  miss <- setdiff(required_cols, names(x))
  if (length(miss)) {
    stopf(
      "Malformed file %s; missing columns: %s",
      path, paste(miss, collapse = ", ")
    )
  }
  x
}

finite_or_na <- function(x) {
  is.na(x) | is.finite(x)
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

qfun <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = p, names = FALSE, type = 7))
}

med <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

sd_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x)
}

range_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  range(x)
}

formal_dir <- file.path(ROOT, "formal_irf_rate_diff")
audit_dir <- file.path(ROOT, "irf_audit")

gate_file <- file.path(formal_dir, "00_formal_rate_diff_irf_gate.csv")
raw_file <- file.path(formal_dir, "irf_posterior_summary.csv")
cum_reer_file <- file.path(formal_dir, "07_core_event_cumulative_reer.csv")
cum_rate_file <- file.path(formal_dir, "06_core_event_cumulative_rate_level.csv")
stability_file <- file.path(formal_dir, "01_irf_stability_by_anchor.csv")
scale_file <- file.path(formal_dir, "08_irf_scale_manifest.csv")
audit_gate_file <- file.path(audit_dir, "05_irf_audit_integrity_gate.csv")
direction_file <- file.path(audit_dir, "04_irf_variable_direction_summary.csv")

gate <- read_required_csv(
  gate_file,
  c(
    "Status","MainNetwork","PosteriorPanel","PosteriorFiles",
    "PosteriorDrawsTotal","EventAnchors","Horizon","RateMode",
    "WorstAnchorStableShare"
  )
)
if (nrow(gate) != 1L || gate$Status[1] != "READY_FOR_RATE_DIFF_IRF_AUDIT") {
  stopf("08m formal IRF gate is not READY_FOR_RATE_DIFF_IRF_AUDIT.")
}
if (tolower(as.character(gate$RateMode[1])) != "difference") {
  stopf("08n requires the accepted Delta-r 08m result.")
}

audit_gate <- read_required_csv(audit_gate_file, c("Check","Pass"))
if (!all(audit_gate$Pass)) {
  stopf("08m R55 integrity gate contains a failed check.")
}

raw <- read_required_csv(
  raw_file,
  c(
    "EventID","EventSet","EventLabel","ShockFamily","AnchorType",
    "AnchorQuarter","Country","ResponseVariable","ResponseScale","Horizon",
    "p05","p16","median","p84","p95","mean","prob_positive",
    "StablePosteriorDraws","TotalPosteriorDraws"
  )
)

cum_reer <- read_required_csv(
  cum_reer_file,
  c(
    "EventID","EventSet","EventLabel","ShockFamily","AnchorType",
    "AnchorQuarter","Country","ResponseVariable","ResponseScale","Horizon",
    "p05","p16","median","p84","p95","mean","prob_positive",
    "StablePosteriorDraws","TotalPosteriorDraws"
  )
)

cum_rate <- read_required_csv(
  cum_rate_file,
  c(
    "EventID","EventSet","EventLabel","ShockFamily","AnchorType",
    "AnchorQuarter","Country","ResponseVariable","ResponseScale","Horizon",
    "p05","p16","median","p84","p95","mean","prob_positive",
    "StablePosteriorDraws","TotalPosteriorDraws"
  )
)

stability <- read_required_csv(
  stability_file,
  c(
    "EventID","EventSet","AnchorType","AnchorQuarter","TotalDraws",
    "FiniteRhoShare","StableShare","G0OKShare","ValidIRFShare","RhoMedian"
  )
)

scale_manifest <- read_required_csv(
  scale_file,
  c("OutputVariable","Meaning","CumulationMethod","PublicationWarning")
)

direction_summary <- read_required_csv(
  direction_file,
  c(
    "ResponseVariable","Paths","ImpactPositiveCredible90",
    "ImpactNegativeCredible90","ImpactCrossesZero90","EverCredible90"
  )
)

# -----------------------------------------------------------------------------
# 00. Source integrity and shape checks
# -----------------------------------------------------------------------------

num_cols <- c("Horizon","p05","p16","median","p84","p95","mean","prob_positive")
for (cc in num_cols) {
  raw[[cc]] <- suppressWarnings(as.numeric(raw[[cc]]))
  if (!all(finite_or_na(raw[[cc]]))) stopf("Non-finite values in raw column %s", cc)
}
for (cc in num_cols) {
  cum_reer[[cc]] <- suppressWarnings(as.numeric(cum_reer[[cc]]))
  cum_rate[[cc]] <- suppressWarnings(as.numeric(cum_rate[[cc]]))
  if (!all(finite_or_na(cum_reer[[cc]]))) stopf("Non-finite values in cumulative REER column %s", cc)
  if (!all(finite_or_na(cum_rate[[cc]]))) stopf("Non-finite values in cumulative rate column %s", cc)
}

raw$Country <- toupper(trimws(raw$Country))
raw$ResponseVariable <- tolower(trimws(raw$ResponseVariable))
cum_reer$Country <- toupper(trimws(cum_reer$Country))
cum_rate$Country <- toupper(trimws(cum_rate$Country))

if (!setequal(unique(raw$Country), EXPECTED_COUNTRIES)) {
  stopf(
    "Country set mismatch. Found: %s",
    paste(sort(unique(raw$Country)), collapse = ",")
  )
}
if (!setequal(unique(raw$ResponseVariable), EXPECTED_VARS)) {
  stopf(
    "Response-variable set mismatch. Found: %s",
    paste(sort(unique(raw$ResponseVariable)), collapse = ",")
  )
}

core_t0 <- raw[
  raw$EventSet == "CORE" &
    raw$AnchorType == "EVENT_QUARTER_t0",
  ,
  drop = FALSE
]
if (!setequal(unique(core_t0$EventID), EXPECTED_CORE_EVENTS)) {
  stopf(
    "Core-event set mismatch. Found: %s",
    paste(sort(unique(core_t0$EventID)), collapse = ",")
  )
}
if (!setequal(sort(unique(core_t0$Horizon)), EXPECTED_RAW_HORIZONS)) {
  stopf("Raw core horizon grid is not exactly 0:12.")
}
expected_core_raw_rows <-
  length(EXPECTED_CORE_EVENTS) *
  length(EXPECTED_COUNTRIES) *
  length(EXPECTED_VARS) *
  length(EXPECTED_RAW_HORIZONS)
if (nrow(core_t0) != expected_core_raw_rows) {
  stopf(
    "Unexpected core t0 raw row count: found %d expected %d",
    nrow(core_t0), expected_core_raw_rows
  )
}

core_cum_reer <- cum_reer[
  cum_reer$EventSet == "CORE" &
    cum_reer$AnchorType == "EVENT_QUARTER_t0",
  ,
  drop = FALSE
]
if (!setequal(unique(core_cum_reer$EventID), EXPECTED_CORE_EVENTS)) {
  stopf("Cumulative REER core-event set mismatch.")
}
if (!setequal(sort(unique(core_cum_reer$Horizon)), c(0L, SAFE_HAVEN_HORIZONS))) {
  stopf("Cumulative REER selected horizon grid is not 0,1,4,8,12.")
}

source_manifest <- data.frame(
  Item = c(
    "Source08mRunID",
    "Source08mHeadSHA",
    "MainNetwork",
    "PosteriorPanel",
    "PosteriorFiles",
    "PosteriorDraws",
    "WorstAnchorStableShare",
    "RawCoreRows",
    "CoreEvents",
    "Countries",
    "Variables",
    "RawHorizons",
    "SafeHavenREERHorizons",
    "PairwiseEventDrawDifferencesAvailable",
    "TVPCoefficientTrajectoriesAvailable"
  ),
  Value = c(
    SOURCE_RUN_ID,
    SOURCE_HEAD_SHA,
    as.character(gate$MainNetwork[1]),
    as.character(gate$PosteriorPanel[1]),
    as.character(gate$PosteriorFiles[1]),
    as.character(gate$PosteriorDrawsTotal[1]),
    sprintf("%.8f", as.numeric(gate$WorstAnchorStableShare[1])),
    as.character(nrow(core_t0)),
    paste(EXPECTED_CORE_EVENTS, collapse = ";"),
    paste(EXPECTED_COUNTRIES, collapse = ";"),
    paste(EXPECTED_VARS, collapse = ";"),
    "0:12",
    paste(SAFE_HAVEN_HORIZONS, collapse = ";"),
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)
write.csv(source_manifest, file.path(OUT, "00_source_manifest.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 01. Core-event IRF profile similarity
# -----------------------------------------------------------------------------

pair_rows <- list()
kk <- 1L

event_pairs <- combn(EXPECTED_CORE_EVENTS, 2, simplify = FALSE)
scope_vars <- c("ALL", EXPECTED_VARS)

for (pp in event_pairs) {
  for (vv in scope_vars) {
    a <- core_t0[core_t0$EventID == pp[1], , drop = FALSE]
    b <- core_t0[core_t0$EventID == pp[2], , drop = FALSE]
    if (vv != "ALL") {
      a <- a[a$ResponseVariable == vv, , drop = FALSE]
      b <- b[b$ResponseVariable == vv, , drop = FALSE]
    }

    ma <- a[, c(
      "Country","ResponseVariable","Horizon",
      "median","p05","p95"
    )]
    mb <- b[, c(
      "Country","ResponseVariable","Horizon",
      "median","p05","p95"
    )]
    names(ma)[4:6] <- c("median_a","p05_a","p95_a")
    names(mb)[4:6] <- c("median_b","p05_b","p95_b")
    mm <- merge(
      ma, mb,
      by = c("Country","ResponseVariable","Horizon"),
      all = FALSE, sort = FALSE
    )
    if (!nrow(mm)) next

    absdiff <- abs(mm$median_a - mm$median_b)
    avg_band <- 0.5 * ((mm$p95_a - mm$p05_a) + (mm$p95_b - mm$p05_b))
    normdiff <- absdiff / pmax(avg_band, 1e-12)

    pair_rows[[kk]] <- data.frame(
      EventA = pp[1],
      EventB = pp[2],
      VariableScope = vv,
      Coordinates = nrow(mm),
      PearsonMedianProfile = safe_cor(mm$median_a, mm$median_b, "pearson"),
      SpearmanMedianProfile = safe_cor(mm$median_a, mm$median_b, "spearman"),
      MedianAbsoluteDifference = med(absdiff),
      P95AbsoluteDifference = qfun(absdiff, 0.95),
      MaxAbsoluteDifference = max(absdiff, na.rm = TRUE),
      MedianAbsDifferenceToAverageBandWidth90 = med(normdiff),
      P95AbsDifferenceToAverageBandWidth90 = qfun(normdiff, 0.95),
      InferenceLimit = "DESCRIPTIVE_ONLY__NO_PAIRED_EVENT_DRAW_DIFFERENCE_IN_08M",
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
}
pairwise_similarity <- do.call(rbind, pair_rows)
write.csv(
  pairwise_similarity,
  file.path(OUT, "01_core_event_profile_similarity.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 02. Event heterogeneity by coordinate and by country-variable
# -----------------------------------------------------------------------------

coord_key <- interaction(
  core_t0$Country,
  core_t0$ResponseVariable,
  core_t0$Horizon,
  drop = TRUE,
  lex.order = TRUE
)
coord_groups <- split(core_t0, coord_key)

coord_rows <- lapply(coord_groups, function(x) {
  rr <- range_finite(x$median)
  bw <- x$p95 - x$p05
  typical_abs <- med(abs(x$median))
  across_sd <- sd_finite(x$median)
  data.frame(
    Country = x$Country[1],
    ResponseVariable = x$ResponseVariable[1],
    Horizon = x$Horizon[1],
    CoreEvents = nrow(x),
    AcrossEventMeanMedianIRF = mean_finite(x$median),
    AcrossEventMedianIRF = med(x$median),
    AcrossEventSDMedianIRF = across_sd,
    AcrossEventMinMedianIRF = rr[1],
    AcrossEventMaxMedianIRF = rr[2],
    AcrossEventRangeMedianIRF = rr[2] - rr[1],
    TypicalAbsMedianIRF = typical_abs,
    MedianPosteriorBandWidth90 = med(bw),
    RelativeSDToTypicalAbsMedian =
      across_sd / pmax(typical_abs, 1e-12),
    RangeToMedianPosteriorBandWidth90 =
      (rr[2] - rr[1]) / pmax(med(bw), 1e-12),
    stringsAsFactors = FALSE
  )
})
hetero_coord <- do.call(rbind, coord_rows)
hetero_coord <- hetero_coord[
  order(hetero_coord$Country, hetero_coord$ResponseVariable, hetero_coord$Horizon),
  ,
  drop = FALSE
]
write.csv(
  hetero_coord,
  file.path(OUT, "02a_core_event_heterogeneity_by_coordinate.csv"),
  row.names = FALSE
)

cv_key <- interaction(
  hetero_coord$Country,
  hetero_coord$ResponseVariable,
  drop = TRUE,
  lex.order = TRUE
)
cv_groups <- split(hetero_coord, cv_key)
cv_rows <- lapply(cv_groups, function(x) {
  data.frame(
    Country = x$Country[1],
    ResponseVariable = x$ResponseVariable[1],
    Coordinates = nrow(x),
    MedianRelativeSDToTypicalAbsMedian =
      med(x$RelativeSDToTypicalAbsMedian),
    P95RelativeSDToTypicalAbsMedian =
      qfun(x$RelativeSDToTypicalAbsMedian, 0.95),
    MedianRangeToPosteriorBandWidth90 =
      med(x$RangeToMedianPosteriorBandWidth90),
    P95RangeToPosteriorBandWidth90 =
      qfun(x$RangeToMedianPosteriorBandWidth90, 0.95),
    MedianAcrossEventRange =
      med(x$AcrossEventRangeMedianIRF),
    MaxAcrossEventRange =
      max(x$AcrossEventRangeMedianIRF, na.rm = TRUE),
    Interpretation =
      "IRF_EVENT_HETEROGENEITY_PROXY__NOT_LATENT_TVP_COEFFICIENT_VARIATION",
    stringsAsFactors = FALSE
  )
})
hetero_cv <- do.call(rbind, cv_rows)
hetero_cv <- hetero_cv[
  order(hetero_cv$ResponseVariable, hetero_cv$Country),
  ,
  drop = FALSE
]
write.csv(
  hetero_cv,
  file.path(OUT, "02b_core_event_heterogeneity_by_country_variable.csv"),
  row.names = FALSE
)

# Overall proxy by response variable
proxy_rows <- lapply(EXPECTED_VARS, function(vv) {
  x <- hetero_coord[hetero_coord$ResponseVariable == vv, , drop = FALSE]
  data.frame(
    ResponseVariable = vv,
    Coordinates = nrow(x),
    MedianRelativeSDToTypicalAbsMedian =
      med(x$RelativeSDToTypicalAbsMedian),
    P95RelativeSDToTypicalAbsMedian =
      qfun(x$RelativeSDToTypicalAbsMedian, 0.95),
    MedianRangeToPosteriorBandWidth90 =
      med(x$RangeToMedianPosteriorBandWidth90),
    P95RangeToPosteriorBandWidth90 =
      qfun(x$RangeToMedianPosteriorBandWidth90, 0.95),
    DirectTVPCoefficientVariationIdentified = FALSE,
    Interpretation =
      "DESCRIPTIVE_IRF_TIME_HETEROGENEITY_PROXY_ONLY",
    stringsAsFactors = FALSE
  )
})
time_proxy <- do.call(rbind, proxy_rows)
write.csv(
  time_proxy,
  file.path(OUT, "11_irf_time_heterogeneity_proxy.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 03. t0 vs t+1 anchor similarity
# -----------------------------------------------------------------------------

anchor_rows <- list()
kk <- 1L
for (ev in EXPECTED_CORE_EVENTS) {
  for (vv in scope_vars) {
    a <- raw[
      raw$EventID == ev & raw$AnchorType == "EVENT_QUARTER_t0",
      ,
      drop = FALSE
    ]
    b <- raw[
      raw$EventID == ev & raw$AnchorType == "NEXT_QUARTER_t1",
      ,
      drop = FALSE
    ]
    if (vv != "ALL") {
      a <- a[a$ResponseVariable == vv, , drop = FALSE]
      b <- b[b$ResponseVariable == vv, , drop = FALSE]
    }

    ma <- a[, c("Country","ResponseVariable","Horizon","median","p05","p95")]
    mb <- b[, c("Country","ResponseVariable","Horizon","median","p05","p95")]
    names(ma)[4:6] <- c("median_t0","p05_t0","p95_t0")
    names(mb)[4:6] <- c("median_t1","p05_t1","p95_t1")
    mm <- merge(
      ma, mb,
      by = c("Country","ResponseVariable","Horizon"),
      all = FALSE, sort = FALSE
    )
    absdiff <- abs(mm$median_t1 - mm$median_t0)
    avg_band <- 0.5 * (
      (mm$p95_t0 - mm$p05_t0) +
      (mm$p95_t1 - mm$p05_t1)
    )
    anchor_rows[[kk]] <- data.frame(
      EventID = ev,
      VariableScope = vv,
      Coordinates = nrow(mm),
      PearsonMedianProfile_t0_t1 =
        safe_cor(mm$median_t0, mm$median_t1, "pearson"),
      SpearmanMedianProfile_t0_t1 =
        safe_cor(mm$median_t0, mm$median_t1, "spearman"),
      MedianAbsoluteDifference_t0_t1 =
        med(absdiff),
      P95AbsoluteDifference_t0_t1 =
        qfun(absdiff, 0.95),
      MedianAbsDifferenceToAverageBandWidth90 =
        med(absdiff / pmax(avg_band, 1e-12)),
      InferenceLimit =
        "DESCRIPTIVE_ANCHOR_COMPARISON__NO_PAIRED_DIFFERENCE_POSTERIOR",
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
}
anchor_similarity <- do.call(rbind, anchor_rows)
write.csv(
  anchor_similarity,
  file.path(OUT, "03_t0_t1_anchor_similarity.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 04. Safe-haven evidence from exact cumulative REER posterior summaries
# -----------------------------------------------------------------------------

safe_dat <- core_cum_reer[
  core_cum_reer$Horizon %in% SAFE_HAVEN_HORIZONS,
  ,
  drop = FALSE
]
safe_dat$PositiveCredible90 <- safe_dat$p05 > 0
safe_dat$NegativeCredible90 <- safe_dat$p95 < 0
safe_dat$CrossesZero90 <- !(safe_dat$PositiveCredible90 | safe_dat$NegativeCredible90)
safe_dat$MedianPositive <- safe_dat$median > 0

safe_groups <- split(safe_dat, safe_dat$Country)
safe_rows <- lapply(safe_groups, function(x) {
  eg <- split(x, x$EventID)
  events_all_pos_median <- sum(vapply(
    eg, function(z) all(z$median > 0), logical(1)
  ))
  events_any_pos_cred <- sum(vapply(
    eg, function(z) any(z$p05 > 0), logical(1)
  ))
  events_all_pos_cred <- sum(vapply(
    eg, function(z) all(z$p05 > 0), logical(1)
  ))

  pos_cred <- sum(x$PositiveCredible90)
  neg_cred <- sum(x$NegativeCredible90)
  med_prob <- med(x$prob_positive)

  label <- if (pos_cred > 0L) {
    "HAS_POSITIVE_REER_90_CREDIBLE_CELLS"
  } else if (is.finite(med_prob) && med_prob >= 0.75) {
    "DIRECTIONAL_POSITIVE_REER_ONLY"
  } else if (is.finite(med_prob) && med_prob > 0.50) {
    "WEAK_POSITIVE_REER_DIRECTION_ONLY"
  } else {
    "NO_POSITIVE_REER_DIRECTIONAL_EVIDENCE"
  }

  data.frame(
    Country = x$Country[1],
    Cells = nrow(x),
    PositiveMedianCells = sum(x$MedianPositive),
    PositiveMedianShare = mean(x$MedianPositive),
    PositiveCredible90Cells = pos_cred,
    NegativeCredible90Cells = neg_cred,
    CrossesZero90Cells = sum(x$CrossesZero90),
    MedianProbPositive = med_prob,
    MinProbPositive = min(x$prob_positive, na.rm = TRUE),
    MaxProbPositive = max(x$prob_positive, na.rm = TRUE),
    MedianCumulativeREER = med(x$median),
    EventsAllFourHorizonsPositiveMedian = events_all_pos_median,
    EventsAnyPositiveCredible90 = events_any_pos_cred,
    EventsAllFourHorizonsPositiveCredible90 = events_all_pos_cred,
    REEREvidenceLabel = label,
    SafeHavenIdentificationLimit =
      "REER_EVIDENCE_ALONE_IS_NOT_SUFFICIENT_FOR_SAFE_HAVEN_IDENTIFICATION",
    stringsAsFactors = FALSE
  )
})
safe_evidence <- do.call(rbind, safe_rows)
safe_evidence <- safe_evidence[
  order(
    -safe_evidence$PositiveCredible90Cells,
    -safe_evidence$MedianProbPositive,
    safe_evidence$Country
  ),
  ,
  drop = FALSE
]
safe_evidence$RankByREERPosteriorDirection <- seq_len(nrow(safe_evidence))
write.csv(
  safe_evidence,
  file.path(OUT, "04_safe_haven_reer_evidence_all_countries.csv"),
  row.names = FALSE
)

focus_safe <- safe_evidence[
  safe_evidence$Country %in% FOCUS_SAFE_HAVEN,
  ,
  drop = FALSE
]
focus_safe <- focus_safe[
  match(FOCUS_SAFE_HAVEN, focus_safe$Country),
  ,
  drop = FALSE
]
write.csv(
  focus_safe,
  file.path(OUT, "05_safe_haven_focus_jp_ch_us.csv"),
  row.names = FALSE
)

focus_financial <- core_t0[
  core_t0$Country %in% FOCUS_SAFE_HAVEN &
    core_t0$Horizon %in% c(0L, SAFE_HAVEN_HORIZONS),
  c(
    "EventID","EventLabel","ShockFamily","Country","ResponseVariable",
    "Horizon","median","p05","p95","prob_positive"
  ),
  drop = FALSE
]
focus_financial$SignClass90 <- ifelse(
  focus_financial$p05 > 0,
  "POSITIVE_CREDIBLE_90",
  ifelse(
    focus_financial$p95 < 0,
    "NEGATIVE_CREDIBLE_90",
    "CROSSES_ZERO_90"
  )
)
write.csv(
  focus_financial,
  file.path(OUT, "06_safe_haven_focus_raw_financial_irf.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 05. TR cumulative REER anomaly: concentration, cross-country rank, h=1 source
# -----------------------------------------------------------------------------

pos_cred_all <- safe_dat[safe_dat$PositiveCredible90, , drop = FALSE]
total_pos_cred <- nrow(pos_cred_all)
tr_pos_cred <- sum(pos_cred_all$Country == "TR")
tr_share_pos_cred <- if (total_pos_cred > 0L) {
  tr_pos_cred / total_pos_cred
} else {
  NA_real_
}

tr_cum_rows <- list()
kk <- 1L
for (ev in EXPECTED_CORE_EVENTS) {
  for (hh in SAFE_HAVEN_HORIZONS) {
    g <- safe_dat[
      safe_dat$EventID == ev & safe_dat$Horizon == hh,
      ,
      drop = FALSE
    ]
    tr <- g[g$Country == "TR", , drop = FALSE]
    if (nrow(tr) != 1L) stopf("Expected one TR cumulative REER row for %s h=%d", ev, hh)

    rr <- rank(g$median, ties.method = "average", na.last = "keep")
    tr_rank <- rr[match("TR", g$Country)]
    center <- med(g$median)
    scale <- mad(g$median, center = center, constant = 1.4826, na.rm = TRUE)
    robust_z <- if (is.finite(scale) && scale > 0) {
      (tr$median[1] - center) / scale
    } else {
      NA_real_
    }

    tr_cum_rows[[kk]] <- data.frame(
      EventID = ev,
      Horizon = hh,
      TRMedian = tr$median[1],
      TRP05 = tr$p05[1],
      TRP95 = tr$p95[1],
      TRProbPositive = tr$prob_positive[1],
      TRPositiveCredible90 = tr$p05[1] > 0,
      CrossCountryMedian = center,
      CrossCountryMADScale = scale,
      TRRobustZ_vs14 = robust_z,
      TRRankAmong14 = tr_rank,
      TRRankPercentileAmong14 = tr_rank / nrow(g),
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
}
tr_cum_anomaly <- do.call(rbind, tr_cum_rows)
write.csv(
  tr_cum_anomaly,
  file.path(OUT, "07_tr_cumulative_reer_anomaly_cells.csv"),
  row.names = FALSE
)

tr_raw <- core_t0[
  core_t0$Country == "TR" &
    core_t0$ResponseVariable == "de",
  c(
    "EventID","EventLabel","ShockFamily","Horizon",
    "median","p05","p95","prob_positive"
  ),
  drop = FALSE
]

# Add cross-country rank to every raw TR de horizon.
tr_rank_raw <- list()
kk <- 1L
for (ii in seq_len(nrow(tr_raw))) {
  ev <- tr_raw$EventID[ii]
  hh <- tr_raw$Horizon[ii]
  g <- core_t0[
    core_t0$EventID == ev &
      core_t0$ResponseVariable == "de" &
      core_t0$Horizon == hh,
    ,
    drop = FALSE
  ]
  rr <- rank(g$median, ties.method = "average", na.last = "keep")
  tr_rank <- rr[match("TR", g$Country)]
  center <- med(g$median)
  scale <- mad(g$median, center = center, constant = 1.4826, na.rm = TRUE)
  robust_z <- if (is.finite(scale) && scale > 0) {
    (tr_raw$median[ii] - center) / scale
  } else {
    NA_real_
  }
  tr_rank_raw[[kk]] <- data.frame(
    EventID = ev,
    Horizon = hh,
    CrossCountryMedian = center,
    CrossCountryMADScale = scale,
    TRRobustZ_vs14 = robust_z,
    TRRankAmong14 = tr_rank,
    TRRankPercentileAmong14 = tr_rank / nrow(g),
    stringsAsFactors = FALSE
  )
  kk <- kk + 1L
}
tr_rank_raw <- do.call(rbind, tr_rank_raw)
tr_raw <- merge(
  tr_raw,
  tr_rank_raw,
  by = c("EventID","Horizon"),
  all.x = TRUE,
  sort = FALSE
)
tr_raw <- tr_raw[order(match(tr_raw$EventID, EXPECTED_CORE_EVENTS), tr_raw$Horizon), ]
write.csv(
  tr_raw,
  file.path(OUT, "08_tr_raw_reer_path.csv"),
  row.names = FALSE
)

# Approximate median-path attribution.
# Exact cumulative quantiles are draw-level in 08m; summing posterior medians is
# not an exact posterior decomposition. Ratios below are therefore diagnostic.
dom_rows <- list()
kk <- 1L
for (ev in EXPECTED_CORE_EVENTS) {
  xr <- tr_raw[tr_raw$EventID == ev, , drop = FALSE]
  xc <- core_cum_reer[
    core_cum_reer$EventID == ev &
      core_cum_reer$Country == "TR",
    ,
    drop = FALSE
  ]
  raw_h0 <- xr$median[xr$Horizon == 0]
  raw_h1 <- xr$median[xr$Horizon == 1]
  cum_h1 <- xc$median[xc$Horizon == 1]
  cum_h12 <- xc$median[xc$Horizon == 12]

  if (length(raw_h0) != 1L || length(raw_h1) != 1L ||
      length(cum_h1) != 1L || length(cum_h12) != 1L) {
    stopf("TR h=1/h=12 decomposition inputs missing for %s", ev)
  }

  dom_rows[[kk]] <- data.frame(
    EventID = ev,
    RawREER_DLOG_Median_h0 = raw_h0,
    RawREER_DLOG_Median_h1 = raw_h1,
    ExactDrawLevelCumulativeREER_Median_h1 = cum_h1,
    ExactDrawLevelCumulativeREER_Median_h12 = cum_h12,
    SumOfRawMedianPath_h0_12 = sum(xr$median[xr$Horizon %in% 0:12]),
    SumAbsRawMedianPath_h0_12 = sum(abs(xr$median[xr$Horizon %in% 0:12])),
    ApproxH1MedianShareOfExactCumMedian_h1 =
      raw_h1 / pmax(abs(cum_h1), 1e-12),
    ApproxH1MedianShareOfExactCumMedian_h12 =
      raw_h1 / pmax(abs(cum_h12), 1e-12),
    H1ShareOfAbsoluteRawMedianPath =
      abs(raw_h1) / pmax(sum(abs(xr$median[xr$Horizon %in% 0:12])), 1e-12),
    AttributionWarning =
      "APPROXIMATE_MEDIAN_PATH_ATTRIBUTION__EXACT_CUMULATIVE_POSTERIOR_WAS_DRAW_LEVEL",
    stringsAsFactors = FALSE
  )
  kk <- kk + 1L
}
tr_dom <- do.call(rbind, dom_rows)
write.csv(
  tr_dom,
  file.path(OUT, "09_tr_h1_dominance_decomposition.csv"),
  row.names = FALSE
)

# Event invariance at raw REER h=1 for all countries.
h1_de <- core_t0[
  core_t0$ResponseVariable == "de" &
    core_t0$Horizon == 1,
  ,
  drop = FALSE
]
h1_groups <- split(h1_de, h1_de$Country)
h1_inv_rows <- lapply(h1_groups, function(x) {
  rr <- range_finite(x$median)
  data.frame(
    Country = x$Country[1],
    CoreEvents = nrow(x),
    MeanRawREER_DLOG_Median_h1 = mean_finite(x$median),
    SDRawREER_DLOG_Median_h1 = sd_finite(x$median),
    RangeRawREER_DLOG_Median_h1 = rr[2] - rr[1],
    RelativeSDToAbsMean_h1 =
      sd_finite(x$median) / pmax(abs(mean_finite(x$median)), 1e-12),
    MedianProbPositive_h1 = med(x$prob_positive),
    stringsAsFactors = FALSE
  )
})
h1_invariance <- do.call(rbind, h1_inv_rows)
h1_invariance <- h1_invariance[
  order(h1_invariance$SDRawREER_DLOG_Median_h1, h1_invariance$Country),
  ,
  drop = FALSE
]
h1_invariance$RankLowestEventSD <- seq_len(nrow(h1_invariance))
write.csv(
  h1_invariance,
  file.path(OUT, "10_de_h1_event_invariance_by_country.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 06. Claims guardrail and gate
# -----------------------------------------------------------------------------

all_pair_spearman <- pairwise_similarity$SpearmanMedianProfile[
  pairwise_similarity$VariableScope == "ALL"
]
min_pair_spearman <- min(all_pair_spearman, na.rm = TRUE)

tr_h1 <- tr_raw[tr_raw$Horizon == 1, , drop = FALSE]
tr_h1_top_all <- all(tr_h1$TRRankAmong14 == length(EXPECTED_COUNTRIES))
mean_h1_share_h12 <- mean_finite(tr_dom$ApproxH1MedianShareOfExactCumMedian_h12)
tr_concentration <- isTRUE(
  total_pos_cred > 0L &&
    tr_pos_cred == total_pos_cred &&
    tr_h1_top_all &&
    is.finite(mean_h1_share_h12) &&
    mean_h1_share_h12 > 0.80
)

focus_pos_cred <- setNames(
  focus_safe$PositiveCredible90Cells,
  focus_safe$Country
)

guardrails <- data.frame(
  Claim = c(
    "Core-event IRF posterior-median profiles are descriptively very similar",
    "A specific core event has a credibly different IRF from another core event",
    "JP has positive cumulative REER evidence at 90% central credibility",
    "CH has positive cumulative REER evidence at 90% central credibility",
    "US has positive cumulative REER evidence at 90% central credibility",
    "TR has positive cumulative REER cells at 90% central credibility",
    "TR should therefore be called a safe-haven currency",
    "TR cumulative REER anomaly is concentrated in the raw h=1 REER_DLOG response",
    "Latent TVP coefficients themselves show little time variation"
  ),
  SupportedBy08mArtifact = c(
    is.finite(min_pair_spearman) && min_pair_spearman >= 0.99,
    FALSE,
    focus_pos_cred["JP"] > 0,
    focus_pos_cred["CH"] > 0,
    focus_pos_cred["US"] > 0,
    tr_pos_cred > 0,
    FALSE,
    tr_concentration,
    FALSE
  ),
  EvidenceOrLimit = c(
    sprintf("Minimum all-variable pairwise Spearman across six core events = %.8f", min_pair_spearman),
    "08m stores event-specific posterior summaries, not paired draw-level event-difference posteriors.",
    sprintf("JP positive credible cumulative-REER cells = %d", focus_pos_cred["JP"]),
    sprintf("CH positive credible cumulative-REER cells = %d", focus_pos_cred["CH"]),
    sprintf("US positive credible cumulative-REER cells = %d", focus_pos_cred["US"]),
    sprintf(
      "TR positive credible cumulative-REER cells = %d of total positive credible cells = %d",
      tr_pos_cred, total_pos_cred
    ),
    "REER evidence alone is insufficient for safe-haven identification; the TR pattern is also an anomaly candidate.",
    sprintf(
      "TR is top-ranked at raw de h=1 in all core events = %s; mean approximate h1 share of cumulative h12 median = %.4f",
      tr_h1_top_all, mean_h1_share_h12
    ),
    "08m artifact does not contain latent coefficient trajectories; only event-anchor IRF heterogeneity can be proxied."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  guardrails,
  file.path(OUT, "12_claims_guardrail.csv"),
  row.names = FALSE
)

tr_status <- if (tr_concentration) {
  "TR_REER_H1_ANOMALY_CONCENTRATED"
} else if (tr_pos_cred > 0) {
  "TR_REER_ANOMALY_REQUIRES_REVIEW"
} else {
  "NO_TR_POSITIVE_CREDIBLE_REER_CONCENTRATION"
}

gate08n <- data.frame(
  Status = "DIAGNOSTIC_COMPLETE",
  Source08mRunID = SOURCE_RUN_ID,
  Source08mHeadSHA = SOURCE_HEAD_SHA,
  MainNetwork = as.character(gate$MainNetwork[1]),
  PosteriorDraws = as.integer(gate$PosteriorDrawsTotal[1]),
  Worst08mAnchorStableShare = as.numeric(gate$WorstAnchorStableShare[1]),
  CoreEvents = length(EXPECTED_CORE_EVENTS),
  Countries = length(EXPECTED_COUNTRIES),
  Variables = length(EXPECTED_VARS),
  MinCoreEventProfileSpearman = min_pair_spearman,
  TotalPositiveCredibleCumulativeREERCells = total_pos_cred,
  TRPositiveCredibleCumulativeREERCells = tr_pos_cred,
  TRShareOfPositiveCredibleCumulativeREERCells = tr_share_pos_cred,
  JPPositiveCredibleCumulativeREERCells = as.integer(focus_pos_cred["JP"]),
  CHPositiveCredibleCumulativeREERCells = as.integer(focus_pos_cred["CH"]),
  USPositiveCredibleCumulativeREERCells = as.integer(focus_pos_cred["US"]),
  TRRawREER_h1TopRankAllCoreEvents = tr_h1_top_all,
  MeanApproxTR_h1ShareOfCumREER_h12Median = mean_h1_share_h12,
  TRDiagnosticStatus = tr_status,
  PairwiseEventDifferencePosteriorAvailable = FALSE,
  DirectTVPCoefficientVariationIdentified = FALSE,
  Reason = paste(
    "08n is a source-only diagnostic of accepted 08m outputs.",
    "Event heterogeneity is descriptive, not a posterior test of event-to-event differences.",
    "Safe-haven evidence is reported conservatively from exact draw-level cumulative REER summaries.",
    "TR is audited as an anomaly candidate rather than relabelled as a safe haven."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  gate08n,
  file.path(OUT, "00_08n_gate.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 07. Lightweight diagnostic plots (base R only)
# -----------------------------------------------------------------------------

# Safe-haven REER directional ranking.
ord <- order(safe_evidence$MedianProbPositive, decreasing = TRUE)
png(
  file.path(PLOT_DIR, "01_safe_haven_reer_direction_rank.png"),
  width = 1400, height = 900, res = 140
)
op <- par(mar = c(9, 5, 4, 2) + 0.1)
barplot(
  safe_evidence$MedianProbPositive[ord],
  names.arg = safe_evidence$Country[ord],
  las = 2,
  ylim = c(0, 1),
  ylab = "Median posterior P(cumulative REER response > 0)",
  main = "Core-event cumulative REER directional evidence"
)
abline(h = 0.5, lty = 2)
par(op)
dev.off()

# TR raw REER path by event.
png(
  file.path(PLOT_DIR, "02_tr_raw_reer_path_by_event.png"),
  width = 1400, height = 900, res = 140
)
yr <- range(tr_raw$median, finite = TRUE)
plot(
  NA,
  xlim = range(EXPECTED_RAW_HORIZONS),
  ylim = yr,
  xlab = "Horizon (quarters)",
  ylab = "Posterior median raw REER_DLOG IRF",
  main = "TR raw REER response across six core event anchors"
)
for (ii in seq_along(EXPECTED_CORE_EVENTS)) {
  ev <- EXPECTED_CORE_EVENTS[ii]
  z <- tr_raw[tr_raw$EventID == ev, , drop = FALSE]
  z <- z[order(z$Horizon), , drop = FALSE]
  lines(z$Horizon, z$median, type = "o", pch = ii)
}
abline(h = 0, lty = 2)
legend(
  "topright",
  legend = EXPECTED_CORE_EVENTS,
  pch = seq_along(EXPECTED_CORE_EVENTS),
  lty = 1,
  cex = 0.8,
  bty = "n"
)
dev.off()

# Pairwise event-profile Spearman matrix.
all_pairs <- pairwise_similarity[
  pairwise_similarity$VariableScope == "ALL",
  ,
  drop = FALSE
]
cmat <- diag(1, length(EXPECTED_CORE_EVENTS))
rownames(cmat) <- colnames(cmat) <- EXPECTED_CORE_EVENTS
for (ii in seq_len(nrow(all_pairs))) {
  a <- all_pairs$EventA[ii]
  b <- all_pairs$EventB[ii]
  cmat[a, b] <- all_pairs$SpearmanMedianProfile[ii]
  cmat[b, a] <- all_pairs$SpearmanMedianProfile[ii]
}
write.csv(
  cmat,
  file.path(OUT, "01b_core_event_spearman_matrix.csv"),
  row.names = TRUE
)

png(
  file.path(PLOT_DIR, "03_core_event_profile_spearman.png"),
  width = 1400, height = 1200, res = 140
)
par(mar = c(11, 11, 4, 2) + 0.1)
image(
  x = seq_len(nrow(cmat)),
  y = seq_len(ncol(cmat)),
  z = cmat,
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Spearman similarity of core-event median IRF profiles"
)
axis(1, at = seq_len(nrow(cmat)), labels = rownames(cmat), las = 2)
axis(2, at = seq_len(ncol(cmat)), labels = colnames(cmat), las = 2)
box()
dev.off()

# -----------------------------------------------------------------------------
# 08. README
# -----------------------------------------------------------------------------

readme <- c(
  "08n EVENT HETEROGENEITY + SAFE-HAVEN + TR ANOMALY AUDIT",
  "=========================================================",
  "",
  sprintf("Source 08m run: %s", SOURCE_RUN_ID),
  sprintf("Source 08m head SHA: %s", SOURCE_HEAD_SHA),
  sprintf("Source model network: %s", gate$MainNetwork[1]),
  sprintf("Source posterior draws: %s", gate$PosteriorDrawsTotal[1]),
  sprintf("Worst 08m anchor stability share: %.6f", as.numeric(gate$WorstAnchorStableShare[1])),
  "",
  "Scope:",
  "- No MCMC or model re-estimation is performed.",
  "- Only accepted 08m artifact CSV outputs are read.",
  "- Event heterogeneity uses posterior medians and credible-band widths descriptively.",
  "- Exact cumulative REER summaries come from 08m draw-level cumulation.",
  "- JP/CH/US are the pre-specified safe-haven focus economies.",
  "- TR is treated as an anomaly candidate, not automatically as a safe haven.",
  "",
  "Critical inference limits:",
  "- 08m does NOT save paired draw-level event-difference posteriors.",
  "  Therefore 08n does NOT claim that event A differs credibly from event B.",
  "- 08m artifact does NOT save latent TVP coefficient trajectories.",
  "  Therefore 08n reports an IRF event-heterogeneity proxy, not direct TVP coefficient variation.",
  "- REER appreciation evidence alone is not sufficient to identify a safe-haven currency.",
  "",
  sprintf(
    "Minimum all-variable pairwise Spearman similarity across core events: %.8f",
    min_pair_spearman
  ),
  sprintf(
    "Positive 90%%-credible cumulative REER cells at h=1,4,8,12: total=%d, TR=%d, TR share=%s",
    total_pos_cred,
    tr_pos_cred,
    ifelse(is.finite(tr_share_pos_cred), sprintf("%.4f", tr_share_pos_cred), "NA")
  ),
  sprintf(
    "JP/CH/US positive 90%%-credible cumulative REER cells: JP=%d, CH=%d, US=%d",
    focus_pos_cred["JP"], focus_pos_cred["CH"], focus_pos_cred["US"]
  ),
  sprintf(
    "TR raw REER h=1 is top-ranked among 14 countries in every core event: %s",
    tr_h1_top_all
  ),
  sprintf(
    "Mean approximate TR h=1 median share of exact cumulative h=12 median: %.4f",
    mean_h1_share_h12
  ),
  sprintf("TR diagnostic status: %s", tr_status),
  "",
  "Key outputs:",
  "- 00_08n_gate.csv",
  "- 00_source_manifest.csv",
  "- 01_core_event_profile_similarity.csv",
  "- 01b_core_event_spearman_matrix.csv",
  "- 02a_core_event_heterogeneity_by_coordinate.csv",
  "- 02b_core_event_heterogeneity_by_country_variable.csv",
  "- 03_t0_t1_anchor_similarity.csv",
  "- 04_safe_haven_reer_evidence_all_countries.csv",
  "- 05_safe_haven_focus_jp_ch_us.csv",
  "- 06_safe_haven_focus_raw_financial_irf.csv",
  "- 07_tr_cumulative_reer_anomaly_cells.csv",
  "- 08_tr_raw_reer_path.csv",
  "- 09_tr_h1_dominance_decomposition.csv",
  "- 10_de_h1_event_invariance_by_country.csv",
  "- 11_irf_time_heterogeneity_proxy.csv",
  "- 12_claims_guardrail.csv",
  "- plots/*.png"
)
writeLines(readme, file.path(OUT, "README_08n.txt"))

cat("08N EVENT HETEROGENEITY + SAFE-HAVEN + TR ANOMALY AUDIT: DIAGNOSTIC_COMPLETE\n")
print(gate08n)
