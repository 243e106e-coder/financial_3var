#!/usr/bin/env Rscript

# =============================================================================
# 20_build_financial_3var_input.R
# Version 2.0 — automatic common-sample detection
#
# Purpose
#   Build the 14-economy, 3-variable financial GVAR input WITHOUT hard-coding
#   2002Q2–2022Q4.
#
# Domestic variables
#   r   = short-term interest-rate level
#   de  = REER log change
#   deq = equity return
#
# Required common-sample series
#   14 economies × {r, de, deq}
#   Global GPR
#   Brent
#
# Default sample rule
#   Select the LONGEST CONTINUOUS quarterly interval for which every required
#   series is finite.
#
# Optional manual override
#   FIN3_SAMPLE_START=YYYYQn
#   FIN3_SAMPLE_END=YYYYQn
#
# If one or both are blank/unset, automatic sample detection is used.
#
# Inputs
#   source_repo/8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx
#   source_repo/8.12/gpr_quarterly_processed.csv
#   source_repo/8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv
#   data/GCAP_financial_W_2017.csv
#
# Outputs
#   data/derived/model_input_fin3.csv
#   data/derived/panel_fin3_long.csv
#   data/derived/financial_weights.csv
#   data/derived/build_summary.txt
#   data/derived/sample_audit/
#       01_source_range_summary.csv
#       02_series_coverage_audit.csv
#       03_quarter_completeness.csv
#       04_continuous_complete_runs.csv
#       05_final_sample_summary.csv
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("r","de","deq")
SOURCE_SUFFIX <- c(r="RATE_LEVEL", de="REER_DLOG", deq="EQ_RETURN")

SOURCE_DIR <- Sys.getenv("FIN3_SOURCE_DIR", "source_repo/8.12")

MACRO_PATH <- Sys.getenv(
  "FIN3_MACRO_XLSX",
  file.path(SOURCE_DIR, "TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx")
)

GPR_PATH <- Sys.getenv(
  "FIN3_GPR_CSV",
  file.path(SOURCE_DIR, "gpr_quarterly_processed.csv")
)

OIL_PATH <- Sys.getenv(
  "FIN3_OIL_CSV",
  file.path(SOURCE_DIR, "IMF_Brent_quarterly_log_2000Q1_2026Q2.csv")
)

WEIGHT_PATH <- Sys.getenv(
  "FIN3_WEIGHT_CSV",
  "data/GCAP_financial_W_2017.csv"
)

# Blank means AUTO.
SAMPLE_START_OVERRIDE <- trimws(Sys.getenv("FIN3_SAMPLE_START", ""))
SAMPLE_END_OVERRIDE   <- trimws(Sys.getenv("FIN3_SAMPLE_END", ""))

GPR_COLUMN <- Sys.getenv("FIN3_GPR_COLUMN", "LN_GPR_QMEAN")

OUT_DIR <- "data/derived"
AUDIT_DIR <- file.path(OUT_DIR, "sample_audit")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg <- function(...) cat(sprintf(...), "\n")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

quarter_id <- function(x) {
  sx <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))

  m <- regexec("^([12][0-9]{3})[^0-9]*Q([1-4])$", sx)
  mm <- regmatches(sx, m)
  ok <- lengths(mm) == 3L

  if (any(ok)) {
    yr <- as.integer(vapply(mm[ok], `[`, character(1), 2))
    qq <- as.integer(vapply(mm[ok], `[`, character(1), 3))
    out[ok] <- 4L * yr + qq
  }

  out
}

quarter_label <- function(qid) {
  yr <- (qid - 1L) %/% 4L
  qq <- qid - 4L * yr
  sprintf("%dQ%d", yr, qq)
}

range_text <- function(qid) {
  qid <- sort(unique(qid[is.finite(qid)]))
  if (!length(qid)) return(NA_character_)
  paste0(quarter_label(min(qid)), " - ", quarter_label(max(qid)))
}

for (p in c(MACRO_PATH, GPR_PATH, OIL_PATH, WEIGHT_PATH)) {
  if (!file.exists(p)) stopf("Required input not found: %s", p)
}

if (!requireNamespace("readxl", quietly = TRUE)) {
  stopf("Package 'readxl' is required.")
}

# =============================================================================
# 1. Read 14-economy macro workbook
# =============================================================================

sheets <- readxl::excel_sheets(MACRO_PATH)

required_cols <- c(
  "Quarter",
  unlist(
    lapply(
      COUNTRIES,
      function(cc) paste0(cc, "_", unname(SOURCE_SUFFIX))
    ),
    use.names = FALSE
  )
)

preferred_sheets <- c("MODEL_COMPLETE_14C", "MODEL_WIDE_ALL")
sheet_order <- c(
  intersect(preferred_sheets, sheets),
  setdiff(sheets, preferred_sheets)
)

macro <- NULL
macro_sheet <- NA_character_

for (ss in sheet_order) {

  d <- tryCatch(
    as.data.frame(
      readxl::read_excel(
        MACRO_PATH,
        sheet = ss,
        col_names = TRUE,
        .name_repair = "minimal",
        guess_max = 3000
      ),
      check.names = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(d)) next

  names(d) <- trimws(names(d))

  if (all(required_cols %in% names(d))) {
    macro <- d
    macro_sheet <- ss
    break
  }
}

if (is.null(macro)) {
  stopf(
    paste0(
      "No workbook sheet contains the exact required mapping for all 14 economies: ",
      "Quarter + RATE_LEVEL + REER_DLOG + EQ_RETURN."
    )
  )
}

macro_qid <- quarter_id(macro$Quarter)

if (mean(!is.na(macro_qid)) < 0.95) {
  stopf("Quarter parsing failed in macro sheet '%s'.", macro_sheet)
}

macro_wide <- data.frame(
  qid = macro_qid,
  Quarter = quarter_label(macro_qid),
  check.names = FALSE
)

for (cc in COUNTRIES) {
  for (v in VARS) {
    src <- paste0(cc, "_", SOURCE_SUFFIX[[v]])
    dst <- paste0(cc, "_", v)
    macro_wide[[dst]] <- num(macro[[src]])
  }
}

macro_wide <- macro_wide[!is.na(macro_wide$qid), , drop = FALSE]
macro_wide <- macro_wide[!duplicated(macro_wide$qid), , drop = FALSE]
macro_wide <- macro_wide[order(macro_wide$qid), , drop = FALSE]

if (!nrow(macro_wide)) stopf("No valid macro observations after parsing.")

# =============================================================================
# 2. Read GPR and Brent
# =============================================================================

read_qcsv <- function(path, exact = NULL, label = "series") {

  d <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  q_rates <- vapply(
    d,
    function(z) mean(!is.na(quarter_id(z))),
    numeric(1)
  )
  q_rates[!is.finite(q_rates)] <- 0

  date_col <- which.max(q_rates)

  if (!length(date_col) || q_rates[date_col] < 0.50) {
    stopf("Could not detect quarter column in %s.", label)
  }

  value_col <- NA_integer_

  if (!is.null(exact) && exact %in% names(d)) {
    value_col <- match(exact, names(d))
  }

  if (is.na(value_col)) {
    candidates <- setdiff(seq_along(d), date_col)

    valid_rate <- vapply(
      candidates,
      function(j) mean(is.finite(num(d[[j]]))),
      numeric(1)
    )

    if (!length(valid_rate) || max(valid_rate) < 0.50) {
      stopf("Could not detect numeric value column in %s.", label)
    }

    value_col <- candidates[which.max(valid_rate)]
  }

  z <- data.frame(
    qid = quarter_id(d[[date_col]]),
    value = num(d[[value_col]]),
    stringsAsFactors = FALSE
  )

  z <- z[!is.na(z$qid), , drop = FALSE]

  # If duplicate quarter rows exist, average only finite observations.
  q_unique <- sort(unique(z$qid))
  out <- data.frame(
    qid = q_unique,
    value = NA_real_
  )

  for (i in seq_along(q_unique)) {
    vals <- z$value[z$qid == q_unique[i]]
    vals <- vals[is.finite(vals)]
    if (length(vals)) out$value[i] <- mean(vals)
  }

  attr(out, "value_column") <- names(d)[value_col]
  out
}

gpr <- read_qcsv(GPR_PATH, GPR_COLUMN, "GPR")
oil <- read_qcsv(OIL_PATH, NULL, "Brent")

# =============================================================================
# 3. Build master quarterly grid over the union of all source ranges
# =============================================================================

all_qid <- sort(unique(c(
  macro_wide$qid,
  gpr$qid,
  oil$qid
)))

if (!length(all_qid)) stopf("No quarterly observations available.")

master_qid <- seq(min(all_qid), max(all_qid), by = 1L)

master <- data.frame(
  qid = master_qid,
  Quarter = quarter_label(master_qid),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

macro_match <- match(master$qid, macro_wide$qid)

macro_cols <- unlist(
  lapply(
    COUNTRIES,
    function(cc) paste0(cc, "_", VARS)
  ),
  use.names = FALSE
)

for (nm in macro_cols) {
  master[[nm]] <- macro_wide[[nm]][macro_match]
}

master$gpr <- gpr$value[match(master$qid, gpr$qid)]
master$brent <- oil$value[match(master$qid, oil$qid)]

# =============================================================================
# 4. Series-level coverage audit
# =============================================================================

audit_rows <- list()
k <- 0L

add_audit <- function(group, country, variable, x) {
  idx <- which(is.finite(x))
  k <<- k + 1L

  audit_rows[[k]] <<- data.frame(
    Group = group,
    Country = country,
    Variable = variable,
    FirstFiniteQuarter = if (length(idx)) master$Quarter[min(idx)] else NA_character_,
    LastFiniteQuarter = if (length(idx)) master$Quarter[max(idx)] else NA_character_,
    FiniteObservations = length(idx),
    TotalMasterQuarters = nrow(master),
    CoverageShare = length(idx) / nrow(master),
    InternalMissingCount = if (length(idx) >= 2L) {
      sum(!is.finite(x[min(idx):max(idx)]))
    } else {
      NA_integer_
    },
    stringsAsFactors = FALSE
  )
}

for (cc in COUNTRIES) {
  for (v in VARS) {
    add_audit(
      group = "Domestic",
      country = cc,
      variable = v,
      x = master[[paste0(cc, "_", v)]]
    )
  }
}

add_audit("Global", "GLOBAL", "gpr", master$gpr)
add_audit("Global", "GLOBAL", "brent", master$brent)

series_audit <- do.call(rbind, audit_rows)

write.csv(
  series_audit,
  file.path(AUDIT_DIR, "02_series_coverage_audit.csv"),
  row.names = FALSE,
  na = ""
)

# =============================================================================
# 5. Quarter-level completeness and automatic longest continuous run
# =============================================================================

required_series <- c(macro_cols, "gpr", "brent")

finite_matrix <- sapply(
  master[, required_series, drop = FALSE],
  is.finite
)

if (is.null(dim(finite_matrix))) {
  finite_matrix <- matrix(finite_matrix, ncol = 1L)
}

master$RequiredSeriesCount <- length(required_series)
master$FiniteRequiredSeriesCount <- rowSums(finite_matrix)
master$MissingRequiredSeriesCount <- length(required_series) - master$FiniteRequiredSeriesCount
master$CompleteForModel <- master$MissingRequiredSeriesCount == 0L

write.csv(
  master[, c(
    "Quarter",
    "RequiredSeriesCount",
    "FiniteRequiredSeriesCount",
    "MissingRequiredSeriesCount",
    "CompleteForModel"
  )],
  file.path(AUDIT_DIR, "03_quarter_completeness.csv"),
  row.names = FALSE
)

complete_idx <- which(master$CompleteForModel)

if (!length(complete_idx)) {
  stopf("No quarter has complete coverage for all 14×3 domestic variables + GPR + Brent.")
}

# Split complete observations into truly continuous runs.
run_break <- c(
  TRUE,
  diff(complete_idx) != 1L
)

run_id <- cumsum(run_break)
run_list <- split(complete_idx, run_id)

runs <- do.call(
  rbind,
  lapply(seq_along(run_list), function(i) {
    ii <- run_list[[i]]

    data.frame(
      Run = i,
      StartQuarter = master$Quarter[min(ii)],
      EndQuarter = master$Quarter[max(ii)],
      StartQID = master$qid[min(ii)],
      EndQID = master$qid[max(ii)],
      NQuarters = length(ii),
      stringsAsFactors = FALSE
    )
  })
)

runs <- runs[
  order(-runs$NQuarters, runs$StartQID),
  ,
  drop = FALSE
]

write.csv(
  runs,
  file.path(AUDIT_DIR, "04_continuous_complete_runs.csv"),
  row.names = FALSE
)

auto_start_qid <- runs$StartQID[1]
auto_end_qid <- runs$EndQID[1]

# =============================================================================
# 6. Optional manual sample override
# =============================================================================

manual_override <- nzchar(SAMPLE_START_OVERRIDE) || nzchar(SAMPLE_END_OVERRIDE)

if (manual_override) {

  if (!nzchar(SAMPLE_START_OVERRIDE) || !nzchar(SAMPLE_END_OVERRIDE)) {
    stopf(
      paste0(
        "Manual sample override requires BOTH FIN3_SAMPLE_START and ",
        "FIN3_SAMPLE_END."
      )
    )
  }

  override_start_qid <- quarter_id(SAMPLE_START_OVERRIDE)
  override_end_qid <- quarter_id(SAMPLE_END_OVERRIDE)

  if (
    length(override_start_qid) != 1L ||
    length(override_end_qid) != 1L ||
    is.na(override_start_qid) ||
    is.na(override_end_qid)
  ) {
    stopf("Invalid manual sample override.")
  }

  if (override_start_qid > override_end_qid) {
    stopf("FIN3_SAMPLE_START is after FIN3_SAMPLE_END.")
  }

  final_start_qid <- override_start_qid
  final_end_qid <- override_end_qid
  sample_rule <- "MANUAL_OVERRIDE"

} else {

  final_start_qid <- auto_start_qid
  final_end_qid <- auto_end_qid
  sample_rule <- "AUTO_LONGEST_CONTINUOUS_COMPLETE_RUN"
}

final_mask <- master$qid >= final_start_qid &
              master$qid <= final_end_qid

final <- master[final_mask, , drop = FALSE]

if (!nrow(final)) stopf("Final sample is empty.")

if (any(diff(final$qid) != 1L)) {
  stopf("Final sample is not continuous.")
}

bad_final <- !complete.cases(
  data.frame(
    lapply(
      final[, required_series, drop = FALSE],
      function(x) is.finite(x)
    )
  )
)

if (any(bad_final)) {
  bad_q <- final$Quarter[which(bad_final)[1]]

  stopf(
    paste0(
      "Final sample contains missing/non-finite required observations. ",
      "First bad quarter: %s"
    ),
    bad_q
  )
}

# =============================================================================
# 7. Source-range summary
# =============================================================================

source_summary <- data.frame(
  Source = c(
    "Macro workbook raw quarter column",
    "3-variable macro panel union",
    "GPR",
    "Brent",
    "Automatic complete common run",
    "Final estimation sample"
  ),
  Range = c(
    range_text(macro_wide$qid),
    range_text(
      master$qid[
        rowSums(
          sapply(master[, macro_cols, drop = FALSE], is.finite)
        ) > 0
      ]
    ),
    range_text(gpr$qid[is.finite(gpr$value)]),
    range_text(oil$qid[is.finite(oil$value)]),
    paste0(
      quarter_label(auto_start_qid),
      " - ",
      quarter_label(auto_end_qid)
    ),
    paste0(
      quarter_label(final_start_qid),
      " - ",
      quarter_label(final_end_qid)
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  source_summary,
  file.path(AUDIT_DIR, "01_source_range_summary.csv"),
  row.names = FALSE
)

sample_summary <- data.frame(
  Item = c(
    "SampleRule",
    "SourceWorkbookSheet",
    "SourceWorkbookObservedRange",
    "AutomaticLongestCompleteStart",
    "AutomaticLongestCompleteEnd",
    "AutomaticLongestCompleteN",
    "FinalSampleStart",
    "FinalSampleEnd",
    "FinalSampleN",
    "ManualOverrideUsed",
    "RequiredDomesticSeries",
    "RequiredGlobalSeries"
  ),
  Value = c(
    sample_rule,
    macro_sheet,
    range_text(macro_wide$qid),
    quarter_label(auto_start_qid),
    quarter_label(auto_end_qid),
    auto_end_qid - auto_start_qid + 1L,
    quarter_label(final_start_qid),
    quarter_label(final_end_qid),
    nrow(final),
    manual_override,
    paste0(length(COUNTRIES), " countries × ", length(VARS), " variables = ", length(macro_cols)),
    "GPR + Brent"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_summary,
  file.path(AUDIT_DIR, "05_final_sample_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# 8. Read and validate GCAP financial W
# =============================================================================

wtab <- read.csv(
  WEIGHT_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

names(wtab) <- trimws(names(wtab))

row_ids <- toupper(trimws(as.character(wtab[[1]])))
col_ids <- toupper(trimws(names(wtab)[-1]))

if (!all(COUNTRIES %in% row_ids)) {
  stopf(
    "Missing weight rows: %s",
    paste(setdiff(COUNTRIES, row_ids), collapse = ", ")
  )
}

if (!all(COUNTRIES %in% col_ids)) {
  stopf(
    "Missing weight columns: %s",
    paste(setdiff(COUNTRIES, col_ids), collapse = ", ")
  )
}

W <- matrix(
  NA_real_,
  nrow = length(COUNTRIES),
  ncol = length(COUNTRIES),
  dimnames = list(COUNTRIES, COUNTRIES)
)

for (i in COUNTRIES) {
  rr <- match(i, row_ids)

  for (j in COUNTRIES) {
    cc <- match(j, col_ids) + 1L
    W[i, j] <- num(wtab[[cc]][rr])
  }
}

if (any(!is.finite(W))) stopf("Financial W contains non-finite values.")
if (any(W < -1e-12)) stopf("Financial W contains negative values.")

if (max(abs(diag(W))) > 1e-8) {
  stopf("Financial-weight diagonal is not zero.")
}

if (max(abs(rowSums(W) - 1)) > 1e-6) {
  stopf("Financial weights are not row normalized.")
}

diag(W) <- 0
W <- W / rowSums(W)

# =============================================================================
# 9. Construct foreign variables using FINAL detected sample
# =============================================================================

for (v in VARS) {

  Xv <- do.call(
    cbind,
    lapply(
      COUNTRIES,
      function(cc) final[[paste0(cc, "_", v)]]
    )
  )

  colnames(Xv) <- COUNTRIES

  for (i in COUNTRIES) {
    final[[paste0(i, "_", v, "_star")]] <-
      as.numeric(Xv %*% W[i, ])
  }
}

# =============================================================================
# 10. Save model inputs
# =============================================================================

out_cols <- c(
  "Quarter",
  unlist(
    lapply(
      COUNTRIES,
      function(cc) {
        c(
          paste0(cc, "_", VARS),
          paste0(cc, "_", VARS, "_star")
        )
      }
    ),
    use.names = FALSE
  ),
  "gpr",
  "brent"
)

model <- final[, out_cols, drop = FALSE]

long <- do.call(
  rbind,
  lapply(
    COUNTRIES,
    function(cc) {
      data.frame(
        Quarter = final$Quarter,
        Country = cc,

        r = final[[paste0(cc, "_r")]],
        de = final[[paste0(cc, "_de")]],
        deq = final[[paste0(cc, "_deq")]],

        r_star = final[[paste0(cc, "_r_star")]],
        de_star = final[[paste0(cc, "_de_star")]],
        deq_star = final[[paste0(cc, "_deq_star")]],

        gpr = final$gpr,
        brent = final$brent,

        stringsAsFactors = FALSE
      )
    }
  )
)

if (any(!complete.cases(model))) {
  stopf("Model input unexpectedly contains missing values.")
}

if (any(!complete.cases(long))) {
  stopf("Long panel unexpectedly contains missing values.")
}

write.csv(
  model,
  file.path(OUT_DIR, "model_input_fin3.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  long,
  file.path(OUT_DIR, "panel_fin3_long.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  data.frame(
    Country = COUNTRIES,
    W,
    check.names = FALSE
  ),
  file.path(OUT_DIR, "financial_weights.csv"),
  row.names = FALSE,
  na = ""
)

# =============================================================================
# 11. Human-readable build summary
# =============================================================================

summary_lines <- c(
  "14-economy 3-variable financial GVAR input — AUTO SAMPLE VERSION",
  "",
  paste0("Macro source: ", MACRO_PATH),
  paste0("Macro sheet: ", macro_sheet),
  paste0("Macro observed range: ", range_text(macro_wide$qid)),
  paste0("GPR finite range: ", range_text(gpr$qid[is.finite(gpr$value)])),
  paste0("Brent finite range: ", range_text(oil$qid[is.finite(oil$value)])),
  "",
  paste0(
    "Automatic longest complete common sample: ",
    quarter_label(auto_start_qid),
    " - ",
    quarter_label(auto_end_qid),
    " (",
    auto_end_qid - auto_start_qid + 1L,
    " quarters)"
  ),
  paste0("Sample rule used: ", sample_rule),
  paste0(
    "FINAL ESTIMATION SAMPLE: ",
    quarter_label(final_start_qid),
    " - ",
    quarter_label(final_end_qid),
    " (T = ",
    nrow(final),
    ")"
  ),
  "",
  "Required domestic variables: r, de, deq",
  "Required global variables: GPR, Brent",
  "Foreign variables: r*, de*, deq* using GCAP 2017 financial W",
  paste0("GPR source column: ", attr(gpr, "value_column")),
  paste0("Brent source column: ", attr(oil, "value_column")),
  "",
  paste0(
    "Max |diag(W)| = ",
    format(max(abs(diag(W))), scientific = TRUE)
  ),
  paste0(
    "Max |rowsum(W)-1| = ",
    format(max(abs(rowSums(W) - 1)), scientific = TRUE)
  ),
  "",
  "No new interpolation or imputation is performed by this script.",
  "The final sample is selected only from actually finite observed/processed series."
)

writeLines(
  summary_lines,
  file.path(OUT_DIR, "build_summary.txt")
)

cat(paste(summary_lines, collapse = "\n"), "\n")
