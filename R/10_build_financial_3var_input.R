#!/usr/bin/env Rscript

source("R/00_config.R")

AUDIT_DIR <- file.path(DERIVED_DIR, "sample_audit")
dir.create(DERIVED_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(AUDIT_DIR, recursive=TRUE, showWarnings=FALSE)

for (p in c(MACRO_PATH, GPR_PATH, OIL_PATH)) {
  if (!file.exists(p)) stopf("Required source input not found: %s", p)
}
if (!requireNamespace("readxl", quietly=TRUE)) stopf("Package 'readxl' is required.")

required_macro_columns <- c(
  "Quarter",
  unlist(lapply(COUNTRIES, function(cc) paste0(cc, "_", unname(SOURCE_SUFFIX))),
         use.names=FALSE)
)

sheets <- readxl::excel_sheets(MACRO_PATH)
preferred <- c("MODEL_COMPLETE_14C", "MODEL_WIDE_ALL")
sheet_order <- c(intersect(preferred, sheets), setdiff(sheets, preferred))
macro <- NULL
macro_sheet <- NA_character_

for (ss in sheet_order) {
  z <- tryCatch(
    as.data.frame(readxl::read_excel(
      MACRO_PATH, sheet=ss, col_names=TRUE, .name_repair="minimal", guess_max=3000
    ), check.names=FALSE),
    error=function(e) NULL
  )
  if (is.null(z)) next
  names(z) <- trimws(names(z))
  if (all(required_macro_columns %in% names(z))) {
    macro <- z
    macro_sheet <- ss
    break
  }
}
if (is.null(macro)) {
  stopf("No macro sheet contains Quarter plus all 14 x RATE_LEVEL/REER_DLOG/EQ_RETURN columns.")
}

macro_qid <- quarter_id(macro$Quarter)
if (mean(!is.na(macro_qid)) < 0.95) stopf("Quarter parsing failed in macro sheet %s", macro_sheet)

macro_wide <- data.frame(qid=macro_qid, Quarter=quarter_label(macro_qid),
                         stringsAsFactors=FALSE, check.names=FALSE)
for (cc in COUNTRIES) {
  for (v in VARS) {
    src <- paste0(cc, "_", SOURCE_SUFFIX[[v]])
    macro_wide[[paste0(cc, "_", v)]] <- num(macro[[src]])
  }
}
macro_wide <- macro_wide[!is.na(macro_wide$qid), , drop=FALSE]
macro_wide <- macro_wide[!duplicated(macro_wide$qid), , drop=FALSE]
macro_wide <- macro_wide[order(macro_wide$qid), , drop=FALSE]

read_quarterly_csv <- function(path, preferred_value=NULL, label="series") {
  d <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  q_rate <- vapply(d, function(x) mean(!is.na(quarter_id(x))), numeric(1))
  q_rate[!is.finite(q_rate)] <- 0
  date_col <- which.max(q_rate)
  if (!length(date_col) || q_rate[date_col] < 0.5) stopf("Cannot detect quarter column in %s", label)

  value_col <- if (!is.null(preferred_value) && preferred_value %in% names(d)) {
    match(preferred_value, names(d))
  } else {
    candidates <- setdiff(seq_along(d), date_col)
    finite_rate <- vapply(candidates, function(j) mean(is.finite(num(d[[j]]))), numeric(1))
    if (!length(finite_rate) || max(finite_rate) < 0.5) stopf("Cannot detect value column in %s", label)
    candidates[which.max(finite_rate)]
  }

  z <- data.frame(qid=quarter_id(d[[date_col]]), value=num(d[[value_col]]))
  z <- z[!is.na(z$qid), , drop=FALSE]
  ids <- sort(unique(z$qid))
  out <- data.frame(qid=ids, value=NA_real_)
  for (i in seq_along(ids)) {
    vals <- z$value[z$qid == ids[i] & is.finite(z$value)]
    if (length(vals)) out$value[i] <- mean(vals)
  }
  attr(out, "value_column") <- names(d)[value_col]
  out
}

gpr <- read_quarterly_csv(GPR_PATH, GPR_COLUMN, "GPR")
oil <- read_quarterly_csv(OIL_PATH, NULL, "Brent")

all_qid <- sort(unique(c(macro_wide$qid, gpr$qid, oil$qid)))
if (!length(all_qid)) stopf("No quarterly observations found.")
master <- data.frame(
  qid=seq(min(all_qid), max(all_qid), by=1L),
  stringsAsFactors=FALSE,
  check.names=FALSE
)
master$Quarter <- quarter_label(master$qid)

macro_cols <- unlist(lapply(COUNTRIES, function(cc) paste0(cc, "_", VARS)),
                     use.names=FALSE)
idx <- match(master$qid, macro_wide$qid)
for (nm in macro_cols) master[[nm]] <- macro_wide[[nm]][idx]
master$gpr <- gpr$value[match(master$qid, gpr$qid)]
master$brent <- oil$value[match(master$qid, oil$qid)]

coverage <- list()
k <- 0L
add_coverage <- function(group, country, variable, x) {
  good <- which(is.finite(x))
  k <<- k + 1L
  coverage[[k]] <<- data.frame(
    Group=group, Country=country, Variable=variable,
    FirstFiniteQuarter=if (length(good)) master$Quarter[min(good)] else NA_character_,
    LastFiniteQuarter=if (length(good)) master$Quarter[max(good)] else NA_character_,
    FiniteObservations=length(good), TotalMasterQuarters=nrow(master),
    CoverageShare=length(good)/nrow(master),
    InternalMissingCount=if (length(good)>=2L) sum(!is.finite(x[min(good):max(good)])) else NA_integer_,
    stringsAsFactors=FALSE
  )
}
for (cc in COUNTRIES) for (v in VARS) {
  add_coverage("Domestic", cc, v, master[[paste0(cc,"_",v)]])
}
add_coverage("Global", "GLOBAL", "gpr", master$gpr)
add_coverage("Global", "GLOBAL", "brent", master$brent)
coverage_df <- do.call(rbind, coverage)
write.csv(coverage_df, file.path(AUDIT_DIR, "01_series_coverage.csv"), row.names=FALSE, na="")

required_series <- c(macro_cols, "gpr", "brent")
finite_matrix <- sapply(master[,required_series,drop=FALSE], is.finite)
master$RequiredSeriesCount <- length(required_series)
master$FiniteRequiredSeriesCount <- rowSums(finite_matrix)
master$MissingRequiredSeriesCount <- length(required_series) - master$FiniteRequiredSeriesCount
master$CompleteForModel <- master$MissingRequiredSeriesCount == 0L
write.csv(master[,c("Quarter","RequiredSeriesCount","FiniteRequiredSeriesCount",
                    "MissingRequiredSeriesCount","CompleteForModel")],
          file.path(AUDIT_DIR,"02_quarter_completeness.csv"), row.names=FALSE)

complete_idx <- which(master$CompleteForModel)
if (!length(complete_idx)) stopf("No complete quarter exists for 14 x 3 variables plus GPR and Brent.")
run_id <- cumsum(c(TRUE, diff(complete_idx) != 1L))
run_list <- split(complete_idx, run_id)
runs <- do.call(rbind, lapply(seq_along(run_list), function(i) {
  ii <- run_list[[i]]
  data.frame(Run=i, StartQuarter=master$Quarter[min(ii)], EndQuarter=master$Quarter[max(ii)],
             StartQID=master$qid[min(ii)], EndQID=master$qid[max(ii)], NQuarters=length(ii))
}))
runs <- runs[order(-runs$NQuarters, runs$StartQID), , drop=FALSE]
write.csv(runs, file.path(AUDIT_DIR,"03_continuous_complete_runs.csv"), row.names=FALSE)

auto_start <- runs$StartQID[1]
auto_end <- runs$EndQID[1]
manual <- nzchar(SAMPLE_START_OVERRIDE) || nzchar(SAMPLE_END_OVERRIDE)
if (manual) {
  if (!nzchar(SAMPLE_START_OVERRIDE) || !nzchar(SAMPLE_END_OVERRIDE)) {
    stopf("Manual sample override requires both FIN3_SAMPLE_START and FIN3_SAMPLE_END.")
  }
  final_start <- quarter_id(SAMPLE_START_OVERRIDE)
  final_end <- quarter_id(SAMPLE_END_OVERRIDE)
  if (length(final_start)!=1L || length(final_end)!=1L || is.na(final_start) ||
      is.na(final_end) || final_start > final_end) stopf("Invalid manual sample override.")
  sample_rule <- "MANUAL_OVERRIDE"
} else {
  final_start <- auto_start
  final_end <- auto_end
  sample_rule <- "AUTO_LONGEST_CONTINUOUS_COMPLETE_RUN"
}

final <- master[master$qid >= final_start & master$qid <= final_end, , drop=FALSE]
if (!nrow(final) || any(diff(final$qid) != 1L)) stopf("Final sample is empty or discontinuous.")
if (any(!sapply(final[,required_series,drop=FALSE], function(x) all(is.finite(x))))) {
  stopf("Final sample contains non-finite required observations.")
}

panel <- do.call(rbind, lapply(COUNTRIES, function(cc) {
  data.frame(
    Quarter=final$Quarter, Country=cc,
    r=final[[paste0(cc,"_r")]],
    de=final[[paste0(cc,"_de")]],
    deq=final[[paste0(cc,"_deq")]],
    gpr=final$gpr, brent=final$brent,
    stringsAsFactors=FALSE
  )
}))
if (any(!complete.cases(panel))) stopf("Final panel contains missing values.")
write.csv(panel, file.path(DERIVED_DIR,"panel_domestic_fin3.csv"), row.names=FALSE)

sample_summary <- data.frame(
  Item=c("SampleRule","MacroSheet","AutomaticStart","AutomaticEnd","AutomaticN",
         "FinalStart","FinalEnd","FinalN","ManualOverride","Variables","GlobalControls"),
  Value=c(sample_rule,macro_sheet,quarter_label(auto_start),quarter_label(auto_end),
          auto_end-auto_start+1L,quarter_label(final_start),quarter_label(final_end),
          nrow(final),manual,paste(VARS,collapse=", "),"gpr, brent")
)
write.csv(sample_summary, file.path(AUDIT_DIR,"04_final_sample_summary.csv"), row.names=FALSE)

summary_lines <- c(
  "CLEAN 14-ECONOMY FINANCIAL 3-VARIABLE INPUT",
  "================================================",
  sprintf("Macro source: %s", MACRO_PATH),
  sprintf("Macro sheet: %s", macro_sheet),
  sprintf("GPR value column: %s", attr(gpr,"value_column")),
  sprintf("Brent value column: %s", attr(oil,"value_column")),
  sprintf("Sample rule: %s", sample_rule),
  sprintf("Final sample: %s - %s (T=%d)", quarter_label(final_start),
          quarter_label(final_end), nrow(final)),
  "Variables: r, de, deq",
  "No project-level interpolation or imputation is performed."
)
writeLines(summary_lines, file.path(DERIVED_DIR,"build_summary.txt"))
cat(paste(summary_lines, collapse="\n"), "\n")
