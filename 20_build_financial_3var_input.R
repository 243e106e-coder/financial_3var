#!/usr/bin/env Rscript

# ============================================================
# 20_build_financial_3var_input.R
# Build the 14-economy, 3-variable financial-market GVAR input.
#
# Domestic variables:
#   r   = short-term interest rate (level)
#   de  = change in REER
#   deq = equity return
#
# Cross-country aggregation:
#   GCAP 2017 bilateral portfolio-investment weights.
#
# Global controls retained from the existing project:
#   global GPR (LN_GPR_QMEAN)
#   Brent oil price (if available; retained to avoid changing several
#   dimensions of the old specification at once)
#
# Outputs:
#   financial_3var/data/model_input_fin3.csv
#   financial_3var/data/panel_fin3_long.csv
#   financial_3var/data/financial_weights.csv
#   financial_3var/data/build_summary.txt
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("r", "de", "deq")
SOURCE_SUFFIX <- c(r = "RATE_LEVEL", de = "REER_DLOG", deq = "EQ_RETURN")

SAMPLE_START <- Sys.getenv("FIN3_SAMPLE_START", "2002Q2")
SAMPLE_END   <- Sys.getenv("FIN3_SAMPLE_END",   "2022Q4")
GPR_COLUMN   <- Sys.getenv("FIN3_GPR_COLUMN", "LN_GPR_QMEAN")

OUT_DIR <- "financial_3var/data"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

resolve_file <- function(env_name, candidates, required = TRUE) {
  x <- Sys.getenv(env_name, "")
  if (nzchar(x)) {
    if (!file.exists(x)) stopf("%s points to missing file: %s", env_name, x)
    return(x)
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[[1]])
  if (required) stopf("Could not locate %s. Tried: %s", env_name, paste(candidates, collapse = ", "))
  NA_character_
}

macro_file <- resolve_file(
  "FIN3_MACRO_XLSX",
  c(
    "8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx",
    ".github/workflows/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx",
    "TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx"
  )
)
weight_file <- resolve_file(
  "FIN3_WEIGHT_CSV",
  c("financial_3var/GCAP_financial_W_2017.csv")
)
gpr_file <- resolve_file(
  "FIN3_GPR_CSV",
  c("8.12/gpr_quarterly_processed.csv", "data/gpr_quarterly_processed.csv", "GPR_quarterly_processed_2000Q1_2026Q2.csv")
)
oil_file <- resolve_file(
  "FIN3_OIL_CSV",
  c("8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv", "IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"),
  required = FALSE
)

if (!requireNamespace("readxl", quietly = TRUE)) stopf("Package readxl is required")

# ---------- 1. Read the validated macro workbook ----------
sheets <- readxl::excel_sheets(macro_file)
preferred <- c("MODEL_COMPLETE_14C", "MODEL_WIDE_ALL")
macro_sheet_env <- Sys.getenv("FIN3_MACRO_SHEET", "")
if (nzchar(macro_sheet_env)) {
  if (!macro_sheet_env %in% sheets) stopf("Requested macro sheet does not exist: %s", macro_sheet_env)
  candidate_sheets <- macro_sheet_env
} else {
  candidate_sheets <- c(intersect(preferred, sheets), setdiff(sheets, preferred))
}

required_macro_cols <- c(
  "Quarter",
  unlist(lapply(COUNTRIES, function(cc) paste0(cc, "_", unname(SOURCE_SUFFIX))), use.names = FALSE)
)

macro <- NULL
macro_sheet <- NA_character_
for (ss in candidate_sheets) {
  tmp <- tryCatch(
    as.data.frame(readxl::read_excel(macro_file, sheet = ss, .name_repair = "minimal"), check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tmp)) next
  names(tmp) <- trimws(names(tmp))
  if (all(required_macro_cols %in% names(tmp))) {
    macro <- tmp
    macro_sheet <- ss
    break
  }
}
if (is.null(macro)) {
  stopf("No sheet contains the exact 14-country RATE_LEVEL / REER_DLOG / EQ_RETURN mapping")
}

qid <- quarter_id(macro$Quarter)
if (mean(!is.na(qid)) < 0.95) stopf("Quarter column could not be parsed reliably in %s", macro_sheet)

wide <- data.frame(qid = qid, Quarter = quarter_label(qid), check.names = FALSE)
for (cc in COUNTRIES) {
  for (v in VARS) {
    src <- paste0(cc, "_", SOURCE_SUFFIX[[v]])
    wide[[paste0(cc, "_", v)]] <- num(macro[[src]])
  }
}
wide <- wide[!is.na(wide$qid), , drop = FALSE]
wide <- wide[!duplicated(wide$qid), , drop = FALSE]
wide <- wide[order(wide$qid), , drop = FALSE]

q0 <- quarter_id(SAMPLE_START)
q1 <- quarter_id(SAMPLE_END)
if (is.na(q0) || is.na(q1) || q0 > q1) stopf("Invalid sample bounds: %s to %s", SAMPLE_START, SAMPLE_END)
wide <- wide[wide$qid >= q0 & wide$qid <= q1, , drop = FALSE]
if (!nrow(wide)) stopf("No macro observations inside requested sample")
if (any(diff(wide$qid) != 1L)) stopf("Macro sample has missing quarters inside %s-%s", SAMPLE_START, SAMPLE_END)

macro_cols <- unlist(lapply(COUNTRIES, function(cc) paste0(cc, "_", VARS)), use.names = FALSE)
if (any(!complete.cases(wide[, macro_cols, drop = FALSE]))) {
  miss <- which(!complete.cases(wide[, macro_cols, drop = FALSE]))
  stopf("NA values remain in 3-variable macro sample; first affected quarter: %s", wide$Quarter[miss[1]])
}

# ---------- 2. Read and validate the GCAP financial matrix ----------
wtab <- read.csv(weight_file, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(wtab) <- trimws(names(wtab))
if (ncol(wtab) < 15L) stopf("Financial weight table is unexpectedly narrow")
row_ids <- toupper(trimws(as.character(wtab[[1]])))
col_ids <- toupper(trimws(names(wtab)[-1]))
if (!all(COUNTRIES %in% row_ids)) stopf("Missing financial-weight rows: %s", paste(setdiff(COUNTRIES, row_ids), collapse = ", "))
if (!all(COUNTRIES %in% col_ids)) stopf("Missing financial-weight columns: %s", paste(setdiff(COUNTRIES, col_ids), collapse = ", "))
if (anyDuplicated(row_ids)) stopf("Duplicate financial-weight reporter rows")

W <- matrix(NA_real_, length(COUNTRIES), length(COUNTRIES), dimnames = list(COUNTRIES, COUNTRIES))
for (i in COUNTRIES) {
  rr <- match(i, row_ids)
  for (j in COUNTRIES) {
    cc <- match(j, col_ids) + 1L
    W[i, j] <- num(wtab[[cc]][rr])
  }
}
if (any(!is.finite(W)) || any(W < -1e-12)) stopf("Financial matrix contains invalid values")
if (max(abs(diag(W))) > 1e-8) stopf("Financial matrix diagonal is not zero")
rs <- rowSums(W)
if (any(rs <= 0)) stopf("Financial matrix has non-positive row sum")
if (max(abs(rs - 1)) > 1e-6) stopf("Financial matrix is not row-normalized; max deviation %.8g", max(abs(rs - 1)))
diag(W) <- 0
W <- W / rowSums(W)

# ---------- 3. Read global GPR and Brent ----------
read_quarterly_csv <- function(path, value_col = NULL, label = "series") {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  q_rates <- vapply(d, function(z) mean(!is.na(quarter_id(z))), numeric(1))
  q_rates[!is.finite(q_rates)] <- 0
  dc <- which.max(q_rates)
  if (!length(dc) || q_rates[dc] < 0.50) stopf("Could not detect quarter column in %s", label)
  vc <- NA_integer_
  if (!is.null(value_col) && value_col %in% names(d)) vc <- match(value_col, names(d))
  if (is.na(vc)) {
    cand <- setdiff(seq_along(d), dc)
    nr <- vapply(cand, function(j) mean(is.finite(num(d[[j]]))), numeric(1))
    if (!length(nr) || max(nr) < 0.50) stopf("Could not detect numeric value column in %s", label)
    vc <- cand[which.max(nr)]
  }
  z <- data.frame(qid = quarter_id(d[[dc]]), value = num(d[[vc]]))
  z <- z[!is.na(z$qid) & is.finite(z$value), , drop = FALSE]
  aggregate(value ~ qid, z, mean)
}

gpr <- read_quarterly_csv(gpr_file, GPR_COLUMN, "GPR")
wide$gpr <- gpr$value[match(wide$qid, gpr$qid)]
if (any(!is.finite(wide$gpr))) stopf("GPR does not fully cover the requested sample")

if (!is.na(oil_file)) {
  oil <- read_quarterly_csv(oil_file, NULL, "Brent")
  wide$brent <- oil$value[match(wide$qid, oil$qid)]
  if (any(!is.finite(wide$brent))) stopf("Brent does not fully cover the requested sample")
} else {
  wide$brent <- NA_real_
  warning("Brent file was not found. Input will be built without a usable Brent control.")
}

# ---------- 4. Construct foreign (star) variables ----------
for (v in VARS) {
  domestic_matrix <- do.call(cbind, lapply(COUNTRIES, function(cc) wide[[paste0(cc, "_", v)]]))
  colnames(domestic_matrix) <- COUNTRIES
  for (i in COUNTRIES) {
    wide[[paste0(i, "_", v, "_star")]] <- as.numeric(domestic_matrix %*% W[i, ])
  }
}

# Deterministic output ordering.
out_cols <- c(
  "Quarter",
  unlist(lapply(COUNTRIES, function(cc) c(paste0(cc, "_", VARS), paste0(cc, "_", VARS, "_star"))), use.names = FALSE),
  "gpr",
  "brent"
)
model <- wide[, out_cols, drop = FALSE]

# ---------- 5. Long panel for diagnostics ----------
long_rows <- vector("list", length(COUNTRIES))
for (k in seq_along(COUNTRIES)) {
  cc <- COUNTRIES[k]
  long_rows[[k]] <- data.frame(
    Quarter = wide$Quarter,
    Country = cc,
    r = wide[[paste0(cc, "_r")]],
    de = wide[[paste0(cc, "_de")]],
    deq = wide[[paste0(cc, "_deq")]],
    r_star = wide[[paste0(cc, "_r_star")]],
    de_star = wide[[paste0(cc, "_de_star")]],
    deq_star = wide[[paste0(cc, "_deq_star")]],
    gpr = wide$gpr,
    brent = wide$brent,
    stringsAsFactors = FALSE
  )
}
panel_long <- do.call(rbind, long_rows)

write.csv(model, file.path(OUT_DIR, "model_input_fin3.csv"), row.names = FALSE, na = "")
write.csv(panel_long, file.path(OUT_DIR, "panel_fin3_long.csv"), row.names = FALSE, na = "")
write.csv(data.frame(Country = COUNTRIES, W, check.names = FALSE), file.path(OUT_DIR, "financial_weights.csv"), row.names = FALSE, na = "")

summary_lines <- c(
  "14-economy / 3-variable financial GVAR input",
  paste0("Macro workbook: ", macro_file),
  paste0("Macro sheet: ", macro_sheet),
  paste0("Financial weights: ", weight_file),
  "Financial-weight concept: GCAP 2017 Position Residency portfolio weights",
  paste0("Sample: ", model$Quarter[1], " - ", tail(model$Quarter, 1)),
  paste0("Observations: ", nrow(model), " quarters"),
  paste0("Countries: ", paste(COUNTRIES, collapse = ", ")),
  "Domestic variables: r, de, deq",
  "Foreign variables: financial-weighted r*, de*, deq*",
  paste0("GPR column: ", GPR_COLUMN),
  paste0("Brent included: ", ifelse(all(is.finite(model$brent)), "YES", "NO")),
  paste0("Max financial-weight row-sum deviation: ", format(max(abs(rowSums(W) - 1)), scientific = TRUE)),
  paste0("Max financial-weight diagonal: ", format(max(abs(diag(W))), scientific = TRUE)),
  "No interpolation or imputation is performed by this script."
)
writeLines(summary_lines, file.path(OUT_DIR, "build_summary.txt"))

msg("Built financial 3-variable input: %s to %s (%d quarters)", model$Quarter[1], tail(model$Quarter, 1), nrow(model))
msg("Financial W validation: max diag = %.3e; max row-sum deviation = %.3e", max(abs(diag(W))), max(abs(rowSums(W) - 1)))
