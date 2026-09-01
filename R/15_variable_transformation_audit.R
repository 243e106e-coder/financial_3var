#!/usr/bin/env Rscript

# =============================================================================
# VARIABLE TRANSFORMATION + RATE LEVEL VS DIFFERENCE AUDIT
# Financial 3-variable TVP-GVAR
#
# Repository variables:
#   r   = RATE_LEVEL
#   de  = REER_DLOG
#   deq = EQ_RETURN
#
# Goals
# -----
# 1) Verify that REER and equity are already transformed upstream.
# 2) Diagnose stationarity/persistence without mechanically over-differencing.
# 3) Compare interest-rate LEVEL vs first difference on the SAME sample.
# 4) Recompute the static p=1, q=1 GVAR stability diagnostic for both rate forms.
# 5) Produce drop-in rate-difference robustness input without altering baseline data.
#
# IMPORTANT
# ---------
# - This script does NOT change the baseline panel.
# - It does NOT log the interest rate.
# - It does NOT difference REER_DLOG or EQ_RETURN again.
# - Stationarity tests are diagnostics, not an automatic model-selection rule.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("urca", quietly = TRUE)) {
  stopf("Package 'urca' is required for ADF/KPSS diagnostics.")
}

PANEL <- file.path(DERIVED_DIR, "panel_domestic_fin3.csv")
OUT <- file.path(RESULTS_DIR, "variable_transform")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL)) {
  stopf("Missing %s. Run R/10_build_financial_3var_input.R first.", PANEL)
}

d <- read.csv(PANEL, stringsAsFactors = FALSE, check.names = FALSE)

required <- c("Quarter", "Country", "r", "de", "deq", "gpr", "brent")
if (!all(required %in% names(d))) {
  stopf(
    "Panel is missing required columns: %s",
    paste(setdiff(required, names(d)), collapse = ", ")
  )
}

d$Country <- toupper(trimws(as.character(d$Country)))
d$Quarter <- toupper(trimws(as.character(d$Quarter)))
d$qid <- quarter_id(d$Quarter)

if (any(is.na(d$qid))) stopf("Unparseable Quarter value in baseline panel.")
if (any(!d$Country %in% COUNTRIES)) {
  stopf("Unknown country code(s) in baseline panel.")
}
if (!setequal(unique(d$Country), COUNTRIES)) {
  stopf("Baseline panel does not contain exactly the configured 14 economies.")
}

for (nm in c("r", "de", "deq", "gpr", "brent")) d[[nm]] <- num(d[[nm]])
if (any(!is.finite(as.matrix(d[, c("r","de","deq","gpr","brent")])))) {
  stopf("Non-finite values found in baseline model panel.")
}

quarters <- sort(unique(d$qid))
if (any(diff(quarters) != 1L)) stopf("Baseline panel is not continuous quarterly data.")

counts <- table(d$Country)
if (any(counts != length(quarters))) stopf("Baseline panel is not balanced.")

key <- paste(d$Country, d$qid, sep = "||")
if (anyDuplicated(key)) stopf("Duplicate country-quarter rows in baseline panel.")

d <- d[order(d$qid, match(d$Country, COUNTRIES)), , drop = FALSE]

# =============================================================================
# 1. Fixed variable-definition manifest
# =============================================================================

definitions <- data.frame(
  Variable = c("r", "de", "deq", "r_diff_robustness"),
  SourceField = c(
    unname(SOURCE_SUFFIX["r"]),
    unname(SOURCE_SUFFIX["de"]),
    unname(SOURCE_SUFFIX["deq"]),
    paste0("Delta(", unname(SOURCE_SUFFIX["r"]), ")")
  ),
  CurrentMeaning = c(
    "Short-term interest-rate level",
    "Quarterly log change in REER",
    "Quarterly equity return",
    "Quarter-to-quarter change in the short-term interest-rate level"
  ),
  BaselineUse = c(
    "KEEP_LEVEL",
    "KEEP_AS_IS",
    "KEEP_AS_IS",
    "ROBUSTNESS_ONLY"
  ),
  AdditionalLog = c("NO", "NO", "NO", "NO"),
  AdditionalDifference = c("NO", "NO", "NO", "ALREADY_FIRST_DIFFERENCE_OF_r"),
  Reason = c(
    "Levels preserve direct basis-point/percentage-point interpretation; negative or zero rates also make log(r) invalid.",
    "REER_DLOG is already a log difference; another difference would create a second difference.",
    "EQ_RETURN is already a return; another difference would create a change in returns rather than a return.",
    "Use only as a formal robustness specification for the persistence/stationarity concern in r."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  definitions,
  file.path(OUT, "00_variable_definitions.csv"),
  row.names = FALSE
)

# =============================================================================
# 2. Build robustness series: Delta r
#    This creates a separate panel; baseline is never overwritten.
# =============================================================================

d$dr <- NA_real_
for (cc in COUNTRIES) {
  ii <- which(d$Country == cc)
  ii <- ii[order(d$qid[ii])]
  d$dr[ii] <- c(NA_real_, diff(d$r[ii]))
}

# Exactly one initial missing dr per country is expected.
dr_missing_by_country <- tapply(!is.finite(d$dr), d$Country, sum)
if (any(dr_missing_by_country != 1L)) {
  stopf(
    "Expected exactly one initial missing Delta r per country. Found: %s",
    paste(names(dr_missing_by_country), dr_missing_by_country, sep = "=", collapse = ", ")
  )
}

common_qid <- quarters[-1L]
common_quarters <- quarter_label(common_qid)

level_common <- d[d$qid %in% common_qid, required, drop = FALSE]
rate_diff <- d[d$qid %in% common_qid, c(required, "dr"), drop = FALSE]
rate_diff$r <- rate_diff$dr
rate_diff$dr <- NULL

# Reorder and validate.
level_common <- level_common[
  order(quarter_id(level_common$Quarter), match(level_common$Country, COUNTRIES)),
  ,
  drop = FALSE
]
rate_diff <- rate_diff[
  order(quarter_id(rate_diff$Quarter), match(rate_diff$Country, COUNTRIES)),
  ,
  drop = FALSE
]

if (any(!is.finite(as.matrix(rate_diff[, c("r","de","deq","gpr","brent")])))) {
  stopf("Rate-difference robustness panel contains non-finite values.")
}

write.csv(
  level_common,
  file.path(DERIVED_DIR, "panel_domestic_fin3_level_common.csv"),
  row.names = FALSE
)
write.csv(
  rate_diff,
  file.path(DERIVED_DIR, "panel_domestic_fin3_rate_diff.csv"),
  row.names = FALSE
)

# =============================================================================
# 3. Descriptive statistics for r, Delta r, de, deq
# =============================================================================

qfun <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

series_rows <- list()
kk <- 0L

for (cc in COUNTRIES) {
  z <- d[d$Country == cc, , drop = FALSE]
  z <- z[order(z$qid), , drop = FALSE]

  series_list <- list(
    r = z$r,
    dr = z$dr,
    de = z$de,
    deq = z$deq
  )

  for (vv in names(series_list)) {
    x <- series_list[[vv]]
    x <- x[is.finite(x)]
    kk <- kk + 1L
    series_rows[[kk]] <- data.frame(
      Country = cc,
      Series = vv,
      N = length(x),
      Mean = mean(x),
      SD = stats::sd(x),
      Median = stats::median(x),
      P01 = qfun(x, 0.01),
      P05 = qfun(x, 0.05),
      P95 = qfun(x, 0.95),
      P99 = qfun(x, 0.99),
      Min = min(x),
      Max = max(x),
      MaxAbs = max(abs(x)),
      NonPositiveShare = mean(x <= 0),
      stringsAsFactors = FALSE
    )
  }
}

desc <- do.call(rbind, series_rows)
write.csv(
  desc,
  file.path(OUT, "01_descriptive_statistics.csv"),
  row.names = FALSE
)

# =============================================================================
# 4. ADF and KPSS diagnostics using urca
# =============================================================================

ADF_MAX_LAG <- as.integer(Sys.getenv("FIN3_ADF_MAX_LAG", "4"))
if (!is.finite(ADF_MAX_LAG) || ADF_MAX_LAG < 0L) ADF_MAX_LAG <- 4L

extract_adf <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20L || stats::sd(x) == 0) {
    return(data.frame(
      N = length(x),
      Deterministic = "drift",
      MaxLag = ADF_MAX_LAG,
      SelectedLag = NA_integer_,
      TestStatistic = NA_real_,
      CriticalValue5 = NA_real_,
      RejectUnitRoot5 = NA,
      stringsAsFactors = FALSE
    ))
  }

  fit <- urca::ur.df(
    x,
    type = "drift",
    lags = min(ADF_MAX_LAG, max(0L, floor(length(x) / 5L))),
    selectlags = "BIC"
  )

  tst <- fit@teststat
  cv <- fit@cval

  rn <- rownames(tst)
  if (is.null(rn)) rn <- rownames(cv)
  tau_idx <- grep("^tau", rn)
  if (!length(tau_idx)) tau_idx <- 1L
  tau_idx <- tau_idx[1]

  stat <- if (is.matrix(tst)) as.numeric(tst[tau_idx, 1]) else as.numeric(tst[tau_idx])

  cv5 <- NA_real_
  if (is.matrix(cv)) {
    cn <- colnames(cv)
    if (!is.null(cn) && "5pct" %in% cn) {
      cv5 <- as.numeric(cv[tau_idx, "5pct"])
    } else {
      cv5 <- as.numeric(cv[tau_idx, min(2L, ncol(cv))])
    }
  } else {
    nm <- names(cv)
    if (!is.null(nm) && any(nm == "5pct")) {
      cv5 <- as.numeric(cv[nm == "5pct"][1])
    } else {
      cv5 <- as.numeric(cv[min(2L, length(cv))])
    }
  }

  selected_lag <- if ("lags" %in% methods::slotNames(fit)) {
    as.integer(fit@lags)
  } else {
    NA_integer_
  }

  data.frame(
    N = length(x),
    Deterministic = "drift",
    MaxLag = ADF_MAX_LAG,
    SelectedLag = selected_lag,
    TestStatistic = stat,
    CriticalValue5 = cv5,
    RejectUnitRoot5 = is.finite(stat) && is.finite(cv5) && stat < cv5,
    stringsAsFactors = FALSE
  )
}

extract_kpss <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20L || stats::sd(x) == 0) {
    return(data.frame(
      N = length(x),
      Deterministic = "mu",
      LagRule = "short",
      TestStatistic = NA_real_,
      CriticalValue5 = NA_real_,
      RejectStationarity5 = NA,
      stringsAsFactors = FALSE
    ))
  }

  fit <- urca::ur.kpss(
    x,
    type = "mu",
    lags = "short"
  )

  stat <- as.numeric(fit@teststat)[1]
  cv <- fit@cval

  cv5 <- NA_real_
  if (is.matrix(cv)) {
    cn <- colnames(cv)
    if (!is.null(cn) && "5pct" %in% cn) {
      cv5 <- as.numeric(cv[1, "5pct"])
    } else {
      cv5 <- as.numeric(cv[1, min(2L, ncol(cv))])
    }
  } else {
    nm <- names(cv)
    if (!is.null(nm) && any(nm == "5pct")) {
      cv5 <- as.numeric(cv[nm == "5pct"][1])
    } else {
      # ur.kpss commonly stores 10%, 5%, 2.5%, 1% in that order.
      cv5 <- as.numeric(cv[min(2L, length(cv))])
    }
  }

  data.frame(
    N = length(x),
    Deterministic = "mu",
    LagRule = "short",
    TestStatistic = stat,
    CriticalValue5 = cv5,
    RejectStationarity5 = is.finite(stat) && is.finite(cv5) && stat > cv5,
    stringsAsFactors = FALSE
  )
}

adf_rows <- list()
kpss_rows <- list()
ar_rows <- list()
kk <- 0L

for (cc in COUNTRIES) {
  z <- d[d$Country == cc, , drop = FALSE]
  z <- z[order(z$qid), , drop = FALSE]

  series_list <- list(
    r = z$r,
    dr = z$dr,
    de = z$de,
    deq = z$deq
  )

  for (vv in names(series_list)) {
    x <- series_list[[vv]]

    kk <- kk + 1L

    aa <- extract_adf(x)
    aa$Country <- cc
    aa$Series <- vv
    adf_rows[[kk]] <- aa[, c("Country","Series", setdiff(names(aa), c("Country","Series")))]

    kp <- extract_kpss(x)
    kp$Country <- cc
    kp$Series <- vv
    kpss_rows[[kk]] <- kp[, c("Country","Series", setdiff(names(kp), c("Country","Series")))]

    xf <- x[is.finite(x)]
    ar1 <- NA_real_
    if (length(xf) >= 8L && stats::sd(xf) > 0) {
      fit_ar <- stats::lm(xf[-1] ~ xf[-length(xf)])
      ar1 <- unname(stats::coef(fit_ar)[2])
    }

    ar_rows[[kk]] <- data.frame(
      Country = cc,
      Series = vv,
      N = length(xf),
      AR1 = ar1,
      AbsAR1 = abs(ar1),
      stringsAsFactors = FALSE
    )
  }
}

adf <- do.call(rbind, adf_rows)
kpss <- do.call(rbind, kpss_rows)
persistence <- do.call(rbind, ar_rows)

write.csv(adf, file.path(OUT, "02_adf_tests.csv"), row.names = FALSE)
write.csv(kpss, file.path(OUT, "03_kpss_tests.csv"), row.names = FALSE)
write.csv(persistence, file.path(OUT, "04_persistence_ar1.csv"), row.names = FALSE)

joint <- merge(
  adf[, c("Country","Series","RejectUnitRoot5")],
  kpss[, c("Country","Series","RejectStationarity5")],
  by = c("Country","Series"),
  all = TRUE
)
joint <- merge(
  joint,
  persistence[, c("Country","Series","AR1")],
  by = c("Country","Series"),
  all = TRUE
)

joint$StationarityAssessment <- ifelse(
  joint$RejectUnitRoot5 & !joint$RejectStationarity5,
  "STATIONARY_SUPPORTED",
  ifelse(
    !joint$RejectUnitRoot5 & joint$RejectStationarity5,
    "UNIT_ROOT_CONCERN",
    ifelse(
      joint$RejectUnitRoot5 & joint$RejectStationarity5,
      "TEST_CONFLICT",
      "INCONCLUSIVE"
    )
  )
)

write.csv(
  joint,
  file.path(OUT, "05_stationarity_joint_assessment.csv"),
  row.names = FALSE
)

stationarity_summary <- do.call(
  rbind,
  lapply(c("r","dr","de","deq"), function(vv) {
    z <- joint[joint$Series == vv, , drop = FALSE]
    data.frame(
      Series = vv,
      Economies = nrow(z),
      ADFRejectUnitRoot5 = sum(z$RejectUnitRoot5, na.rm = TRUE),
      KPSSRejectStationarity5 = sum(z$RejectStationarity5, na.rm = TRUE),
      StationarySupported = sum(z$StationarityAssessment == "STATIONARY_SUPPORTED", na.rm = TRUE),
      UnitRootConcern = sum(z$StationarityAssessment == "UNIT_ROOT_CONCERN", na.rm = TRUE),
      TestConflict = sum(z$StationarityAssessment == "TEST_CONFLICT", na.rm = TRUE),
      Inconclusive = sum(z$StationarityAssessment == "INCONCLUSIVE", na.rm = TRUE),
      MedianAR1 = stats::median(z$AR1, na.rm = TRUE),
      MaxAbsAR1 = max(abs(z$AR1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  stationarity_summary,
  file.path(OUT, "06_stationarity_summary_by_series.csv"),
  row.names = FALSE
)

# =============================================================================
# 5. Static GVAR stability comparison:
#    BASELINE_FULL, BASELINE_COMMON, RATE_DIFF_COMMON
#
#    BASELINE_COMMON and RATE_DIFF_COMMON use exactly the same quarters.
# =============================================================================

lag1 <- function(x) c(NA_real_, x[-length(x)])

panel_to_inputs <- function(panel) {
  panel$Country <- toupper(trimws(panel$Country))
  panel$Quarter <- toupper(trimws(panel$Quarter))
  panel$qid <- quarter_id(panel$Quarter)

  qids <- sort(unique(panel$qid))
  qlabels <- quarter_label(qids)

  N <- length(COUNTRIES)
  K <- length(VARS)

  X <- array(
    NA_real_,
    c(length(qids), N, K),
    dimnames = list(qlabels, COUNTRIES, VARS)
  )

  for (i in seq_along(COUNTRIES)) {
    cc <- COUNTRIES[i]
    z <- panel[panel$Country == cc, , drop = FALSE]
    z <- z[match(qids, z$qid), , drop = FALSE]
    if (any(is.na(z$qid))) stopf("Unbalanced candidate panel for %s.", cc)
    X[, i, ] <- as.matrix(z[, VARS, drop = FALSE])
  }

  base <- panel[panel$Country == COUNTRIES[1], , drop = FALSE]
  base <- base[match(qids, base$qid), , drop = FALSE]

  list(
    X = X,
    qids = qids,
    qlabels = qlabels,
    gpr = num(base$gpr),
    brent = num(base$brent)
  )
}

make_star <- function(X, W) {
  N <- dim(X)[2]
  K <- dim(X)[3]
  S <- array(NA_real_, dim(X), dimnames = dimnames(X))
  for (i in seq_len(N)) {
    for (v in seq_len(K)) {
      S[, i, v] <- as.numeric(X[, , v] %*% W[i, ])
    }
  }
  S
}

selector <- function(i, K, N) {
  S <- matrix(0, K, N * K)
  S[, ((i - 1L) * K + 1L):(i * K)] <- diag(K)
  S
}

star_map <- function(i, W, K, N) {
  R <- matrix(0, K, N * K)
  for (j in seq_len(N)) {
    for (v in seq_len(K)) {
      R[v, (j - 1L) * K + v] <- W[i, j]
    }
  }
  R
}

fit_global_candidate <- function(panel, W, network_label, transform_label) {
  inp <- panel_to_inputs(panel)
  X <- inp$X
  gpr <- inp$gpr
  brent <- inp$brent

  if (!all(is.finite(X)) || !all(is.finite(gpr))) {
    stopf("Non-finite candidate inputs for %s.", transform_label)
  }

  include_brent <- all(is.finite(brent))

  N <- dim(X)[2]
  K <- dim(X)[3]
  S <- make_star(X, W)

  fit_country <- function(i) {
    Y <- X[, i, , drop = FALSE][, 1, ]
    Z <- S[, i, , drop = FALSE][, 1, ]

    rows <- 2:nrow(Y)
    D <- data.frame(const = rep(1, length(rows)))

    for (v in seq_len(K)) D[[paste0(VARS[v], "_L1")]] <- lag1(Y[, v])[rows]
    for (v in seq_len(K)) D[[paste0(VARS[v], "_star_0")]] <- Z[rows, v]
    for (v in seq_len(K)) D[[paste0(VARS[v], "_star_L1")]] <- lag1(Z[, v])[rows]

    D$gpr_0 <- gpr[rows]
    D$gpr_L1 <- lag1(gpr)[rows]

    if (include_brent) {
      D$brent_0 <- brent[rows]
      D$brent_L1 <- lag1(brent)[rows]
    }

    YY <- Y[rows, , drop = FALSE]
    ok <- complete.cases(D) & complete.cases(YY)
    D <- D[ok, , drop = FALSE]
    YY <- YY[ok, , drop = FALSE]
    Xm <- as.matrix(D)

    if (nrow(Xm) <= ncol(Xm) + 5L) {
      stopf("Too few observations for %s / %s.", COUNTRIES[i], transform_label)
    }

    B <- matrix(
      NA_real_,
      ncol(Xm),
      K,
      dimnames = list(colnames(Xm), VARS)
    )

    for (eq in seq_len(K)) {
      fit <- lm.fit(Xm, YY[, eq])
      if (any(!is.finite(fit$coefficients))) {
        stopf(
          "Rank deficiency: %s / %s / %s",
          transform_label, COUNTRIES[i], VARS[eq]
        )
      }
      B[, eq] <- fit$coefficients
    }

    A <- B0 <- B1 <- matrix(0, K, K, dimnames = list(VARS, VARS))
    for (eq in VARS) {
      for (v in VARS) {
        A[eq, v] <- B[paste0(v, "_L1"), eq]
        B0[eq, v] <- B[paste0(v, "_star_0"), eq]
        B1[eq, v] <- B[paste0(v, "_star_L1"), eq]
      }
    }

    rho <- max(Mod(eigen(A, only.values = TRUE)$values))

    list(
      A = A,
      B0 = B0,
      B1 = B1,
      n = nrow(YY),
      rho = rho,
      design_kappa = kappa(crossprod(Xm)),
      design_rcond = rcond(crossprod(Xm))
    )
  }

  fits <- lapply(seq_len(N), fit_country)

  G0 <- matrix(0, N * K, N * K)
  G1 <- matrix(0, N * K, N * K)

  for (i in seq_len(N)) {
    rr <- ((i - 1L) * K + 1L):(i * K)
    Si <- selector(i, K, N)
    Ri <- star_map(i, W, K, N)

    G0[rr, ] <- Si - fits[[i]]$B0 %*% Ri
    G1[rr, ] <- fits[[i]]$A %*% Si + fits[[i]]$B1 %*% Ri
  }

  rc <- rcond(G0)
  kap <- kappa(G0)

  F <- tryCatch(solve(G0, G1), error = function(e) NULL)
  ev <- if (is.null(F)) rep(NA_complex_, N * K) else eigen(F, only.values = TRUE)$values
  rho <- if (all(is.na(ev))) NA_real_ else max(Mod(ev), na.rm = TRUE)

  local <- data.frame(
    Transformation = transform_label,
    Network = network_label,
    Country = COUNTRIES,
    N = vapply(fits, `[[`, numeric(1), "n"),
    DomesticSpectralRadius = vapply(fits, `[[`, numeric(1), "rho"),
    DomesticStable = vapply(fits, function(z) z$rho < 1, logical(1)),
    DesignKappa = vapply(fits, `[[`, numeric(1), "design_kappa"),
    DesignRcond = vapply(fits, `[[`, numeric(1), "design_rcond"),
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    Transformation = transform_label,
    Network = network_label,
    SampleStart = inp$qlabels[1],
    SampleEnd = tail(inp$qlabels, 1),
    T = length(inp$qlabels),
    GlobalSpectralRadius = rho,
    Stable = is.finite(rho) && rho < 1,
    G0_rcond = rc,
    G0_kappa = kap,
    MaxLocalDomesticRadius = max(local$DomesticSpectralRadius),
    UnstableLocalCount = sum(!local$DomesticStable),
    stringsAsFactors = FALSE
  )

  list(summary = summary, local = local)
}

networks <- lapply(WEIGHT_FILES, read_weight_matrix)
if (file.exists(TRADE_WEIGHT_PATH)) {
  networks$legacy_trade <- read_weight_matrix(TRADE_WEIGHT_PATH)
}

candidates <- list(
  BASELINE_LEVEL_FULL = d[, required, drop = FALSE],
  BASELINE_LEVEL_COMMON = level_common,
  RATE_DIFF_COMMON = rate_diff
)

stability_rows <- list()
local_rows <- list()
ii <- 0L

for (tt in names(candidates)) {
  for (nn in names(networks)) {
    ii <- ii + 1L
    fit <- fit_global_candidate(
      candidates[[tt]],
      networks[[nn]],
      nn,
      tt
    )
    stability_rows[[ii]] <- fit$summary
    local_rows[[ii]] <- fit$local
  }
}

stability <- do.call(rbind, stability_rows)
local_stability <- do.call(rbind, local_rows)

write.csv(
  stability,
  file.path(OUT, "07_rate_level_vs_difference_stability.csv"),
  row.names = FALSE
)
write.csv(
  local_stability,
  file.path(OUT, "08_rate_level_vs_difference_local_stability.csv"),
  row.names = FALSE
)

# Same-sample direct comparison.
same_sample <- stability[
  stability$Transformation %in% c("BASELINE_LEVEL_COMMON", "RATE_DIFF_COMMON"),
  ,
  drop = FALSE
]

same_sample_wide <- merge(
  same_sample[same_sample$Transformation == "BASELINE_LEVEL_COMMON",
              c("Network","SampleStart","SampleEnd","T","GlobalSpectralRadius",
                "Stable","G0_rcond","G0_kappa","MaxLocalDomesticRadius",
                "UnstableLocalCount")],
  same_sample[same_sample$Transformation == "RATE_DIFF_COMMON",
              c("Network","GlobalSpectralRadius","Stable","G0_rcond","G0_kappa",
                "MaxLocalDomesticRadius","UnstableLocalCount")],
  by = "Network",
  suffixes = c("_RateLevel", "_RateDiff"),
  all = TRUE
)

write.csv(
  same_sample_wide,
  file.path(OUT, "09_same_sample_rate_spec_comparison.csv"),
  row.names = FALSE
)

# =============================================================================
# 6. Recommendation table
# =============================================================================

main_cmp <- same_sample_wide[same_sample_wide$Network == MAIN_NETWORK, , drop = FALSE]
if (nrow(main_cmp) != 1L) {
  stopf("Could not identify unique main-network same-sample comparison.")
}

rate_stationarity <- stationarity_summary[stationarity_summary$Series == "r", , drop = FALSE]
dr_stationarity <- stationarity_summary[stationarity_summary$Series == "dr", , drop = FALSE]
de_stationarity <- stationarity_summary[stationarity_summary$Series == "de", , drop = FALSE]
deq_stationarity <- stationarity_summary[stationarity_summary$Series == "deq", , drop = FALSE]

recommendation <- data.frame(
  Variable = c("r", "de", "deq"),
  Baseline = c("RATE_LEVEL", "REER_DLOG", "EQ_RETURN"),
  RecommendedAction = c(
    "KEEP_LEVEL_BASELINE__USE_DELTA_R_ROBUSTNESS",
    "KEEP_AS_IS__DO_NOT_DIFFERENCE_AGAIN",
    "KEEP_AS_IS__DO_NOT_DIFFERENCE_AGAIN"
  ),
  DiagnosticEvidence = c(
    sprintf(
      paste0(
        "Rate level: ADF rejects unit root in %d/%d economies; KPSS rejects stationarity in %d/%d; median AR1=%.4f. ",
        "Delta r: ADF rejects unit root in %d/%d; KPSS rejects stationarity in %d/%d; median AR1=%.4f. ",
        "Main-network same-sample spectral radius: level=%.6f, Delta r=%.6f."
      ),
      rate_stationarity$ADFRejectUnitRoot5,
      rate_stationarity$Economies,
      rate_stationarity$KPSSRejectStationarity5,
      rate_stationarity$Economies,
      rate_stationarity$MedianAR1,
      dr_stationarity$ADFRejectUnitRoot5,
      dr_stationarity$Economies,
      dr_stationarity$KPSSRejectStationarity5,
      dr_stationarity$Economies,
      dr_stationarity$MedianAR1,
      main_cmp$GlobalSpectralRadius_RateLevel,
      main_cmp$GlobalSpectralRadius_RateDiff
    ),
    sprintf(
      "REER_DLOG is already a log change. ADF rejects unit root in %d/%d; KPSS rejects stationarity in %d/%d; median AR1=%.4f.",
      de_stationarity$ADFRejectUnitRoot5,
      de_stationarity$Economies,
      de_stationarity$KPSSRejectStationarity5,
      de_stationarity$Economies,
      de_stationarity$MedianAR1
    ),
    sprintf(
      "EQ_RETURN is already a return. ADF rejects unit root in %d/%d; KPSS rejects stationarity in %d/%d; median AR1=%.4f.",
      deq_stationarity$ADFRejectUnitRoot5,
      deq_stationarity$Economies,
      deq_stationarity$KPSSRejectStationarity5,
      deq_stationarity$Economies,
      deq_stationarity$MedianAR1
    )
  ),
  Interpretation = c(
    "Baseline rate IRF is a direct rate-level deviation; Delta-r robustness IRF must be cumulated to recover a rate-level change.",
    "A de IRF is a response in quarterly REER log change; cumulative log-level interpretation requires cumulation.",
    "A deq IRF is a quarterly equity return response; do not call a sum a log-price effect unless upstream EQ_RETURN is explicitly verified as a log return."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  recommendation,
  file.path(OUT, "10_transformation_recommendation.csv"),
  row.names = FALSE
)

# =============================================================================
# 7. Integrity gate
#    Empirical stationarity outcomes do NOT make CI fail.
# =============================================================================

expected_suffix <- c(
  r = "RATE_LEVEL",
  de = "REER_DLOG",
  deq = "EQ_RETURN"
)

rate_nonpositive <- any(d$r <= 0)

integrity <- data.frame(
  Check = c(
    "ConfiguredSourceSuffixesMatchExpectedDefinitions",
    "BaselinePanelBalanced14Economies",
    "BaselinePanelContinuousQuarterly",
    "BaselinePanelFinite",
    "DeltaRHasExactlyOneInitialLossPerEconomy",
    "RateDifferencePanelFinite",
    "SameSampleComparisonUsesIdenticalDates",
    "MainNetworkComparisonAvailable",
    "REERIsAlreadyDLOG",
    "EquityIsAlreadyReturn",
    "RateLogTransformNotRequired"
  ),
  Pass = c(
    identical(unname(SOURCE_SUFFIX[c("r","de","deq")]), unname(expected_suffix)),
    all(table(d$Country) == length(quarters)) && length(unique(d$Country)) == length(COUNTRIES),
    all(diff(quarters) == 1L),
    all(is.finite(as.matrix(d[, c("r","de","deq","gpr","brent")]))),
    all(dr_missing_by_country == 1L),
    all(is.finite(as.matrix(rate_diff[, c("r","de","deq","gpr","brent")]))),
    identical(
      sort(unique(level_common$Quarter)),
      sort(unique(rate_diff$Quarter))
    ),
    nrow(main_cmp) == 1L,
    identical(unname(SOURCE_SUFFIX["de"]), "REER_DLOG"),
    identical(unname(SOURCE_SUFFIX["deq"]), "EQ_RETURN"),
    TRUE
  ),
  Detail = c(
    paste(names(SOURCE_SUFFIX), SOURCE_SUFFIX, sep = "=", collapse = "; "),
    sprintf("%d economies x %d quarters", length(unique(d$Country)), length(quarters)),
    sprintf("%s - %s", quarter_label(min(quarters)), quarter_label(max(quarters))),
    "All baseline r/de/deq/gpr/brent observations finite.",
    paste(names(dr_missing_by_country), dr_missing_by_country, sep = "=", collapse = "; "),
    sprintf("Rate-difference robustness sample: %s - %s", common_quarters[1], tail(common_quarters, 1)),
    sprintf("Common sample T=%d", length(common_qid)),
    sprintf(
      "Level rho=%.6f; Delta-r rho=%.6f",
      main_cmp$GlobalSpectralRadius_RateLevel,
      main_cmp$GlobalSpectralRadius_RateDiff
    ),
    "Do not apply Delta(REER_DLOG).",
    "Do not apply Delta(EQ_RETURN).",
    if (rate_nonpositive) {
      "At least one observed rate is <= 0, so log(r) would be mathematically invalid for the full panel."
    } else {
      "No non-positive rate occurs in this sample, but log(r) is still not required for the baseline economic specification."
    }
  ),
  stringsAsFactors = FALSE
)

write.csv(
  integrity,
  file.path(OUT, "11_transformation_integrity_gate.csv"),
  row.names = FALSE
)

# =============================================================================
# 8. Human-readable report
# =============================================================================

fmt_row <- function(tab, series) {
  z <- tab[tab$Series == series, , drop = FALSE]
  sprintf(
    "%s: ADF reject %d/%d; KPSS reject stationarity %d/%d; median AR1=%.4f",
    series,
    z$ADFRejectUnitRoot5,
    z$Economies,
    z$KPSSRejectStationarity5,
    z$Economies,
    z$MedianAR1
  )
}

report <- c(
  "FINANCIAL 3-VARIABLE TRANSFORMATION AUDIT",
  "==========================================",
  "",
  sprintf("Baseline sample: %s - %s (T=%d)",
          quarter_label(min(quarters)), quarter_label(max(quarters)), length(quarters)),
  sprintf("Same-sample rate robustness window: %s - %s (T=%d)",
          common_quarters[1], tail(common_quarters, 1), length(common_qid)),
  "",
  "BASELINE VARIABLE DEFINITIONS",
  "- r   = RATE_LEVEL",
  "- de  = REER_DLOG",
  "- deq = EQ_RETURN",
  "",
  "TRANSFORMATION DECISION",
  "- Keep r in levels for the baseline.",
  "- Use Delta r as a formal robustness specification.",
  "- Do NOT log r.",
  "- Do NOT difference REER_DLOG again.",
  "- Do NOT difference EQ_RETURN again.",
  "",
  "STATIONARITY / PERSISTENCE DIAGNOSTICS",
  fmt_row(stationarity_summary, "r"),
  fmt_row(stationarity_summary, "dr"),
  fmt_row(stationarity_summary, "de"),
  fmt_row(stationarity_summary, "deq"),
  "",
  "MAIN NETWORK SAME-SAMPLE STATIC STABILITY",
  sprintf("Rate level spectral radius: %.8f", main_cmp$GlobalSpectralRadius_RateLevel),
  sprintf("Rate level stable: %s", main_cmp$Stable_RateLevel),
  sprintf("Delta-r spectral radius: %.8f", main_cmp$GlobalSpectralRadius_RateDiff),
  sprintf("Delta-r stable: %s", main_cmp$Stable_RateDiff),
  "",
  "INTERPRETATION DISCIPLINE",
  "- r level IRF: direct interest-rate level deviation.",
  "- Delta-r IRF: cumulate across horizons to recover a rate-level change.",
  "- de IRF: quarterly REER log-change response; cumulation recovers a cumulative log-level effect.",
  "- deq IRF: quarterly equity-return response. Do not assume it is a log-return unless the upstream construction explicitly verifies that convention.",
  "",
  "OUTPUTS",
  "- results/variable_transform/00_variable_definitions.csv",
  "- results/variable_transform/01_descriptive_statistics.csv",
  "- results/variable_transform/02_adf_tests.csv",
  "- results/variable_transform/03_kpss_tests.csv",
  "- results/variable_transform/04_persistence_ar1.csv",
  "- results/variable_transform/05_stationarity_joint_assessment.csv",
  "- results/variable_transform/06_stationarity_summary_by_series.csv",
  "- results/variable_transform/07_rate_level_vs_difference_stability.csv",
  "- results/variable_transform/08_rate_level_vs_difference_local_stability.csv",
  "- results/variable_transform/09_same_sample_rate_spec_comparison.csv",
  "- results/variable_transform/10_transformation_recommendation.csv",
  "- results/variable_transform/11_transformation_integrity_gate.csv",
  "- data/derived/panel_domestic_fin3_level_common.csv",
  "- data/derived/panel_domestic_fin3_rate_diff.csv"
)

writeLines(report, file.path(OUT, "README_variable_transformation_audit.txt"))

gate_status <- if (all(integrity$Pass)) "READY" else "FAIL"
gate_text <- c(
  "VARIABLE TRANSFORMATION INTEGRITY GATE",
  "======================================",
  sprintf("Status: %s", gate_status),
  sprintf("Passed checks: %d/%d", sum(integrity$Pass), nrow(integrity)),
  "This gate checks data/transformation integrity only.",
  "Stationarity-test outcomes are reported, not forced into a mechanical model-selection rule."
)
writeLines(gate_text, file.path(OUT, "TRANSFORMATION_GATE.txt"))

cat(paste(report, collapse = "\n"), "\n")
cat("\n", paste(gate_text, collapse = "\n"), "\n", sep = "")

if (!all(integrity$Pass)) {
  stopf(
    "Variable-transformation integrity gate failed. Inspect %s",
    file.path(OUT, "11_transformation_integrity_gate.csv")
  )
}
