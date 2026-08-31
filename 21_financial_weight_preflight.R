#!/usr/bin/env Rscript

# ============================================================
# 21_financial_weight_preflight.R
# Validate and describe the GCAP 2017 financial-weight matrix.
# Also compare it with the legacy trade matrix when available.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
FIN_W <- "financial_3var/data/financial_weights.csv"
OUT_DIR <- "financial_3var/results/preflight"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

read_weight_matrix <- function(path) {
  if (!file.exists(path)) stopf("Weight file not found: %s", path)
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(d) <- trimws(names(d))
  rows <- toupper(trimws(as.character(d[[1]])))
  cols <- toupper(trimws(names(d)[-1]))
  if (!all(COUNTRIES %in% rows) || !all(COUNTRIES %in% cols)) stopf("Incomplete 14-country matrix in %s", path)
  W <- matrix(NA_real_, 14, 14, dimnames = list(COUNTRIES, COUNTRIES))
  for (i in COUNTRIES) {
    rr <- match(i, rows)
    for (j in COUNTRIES) W[i,j] <- num(d[[match(j, cols) + 1L]][rr])
  }
  W
}

W <- read_weight_matrix(FIN_W)
validation <- data.frame(
  metric = c("n_rows", "n_cols", "nonfinite_count", "negative_count", "max_abs_diagonal", "max_abs_row_sum_minus_1", "min_weight", "max_weight"),
  value = c(nrow(W), ncol(W), sum(!is.finite(W)), sum(W < -1e-12, na.rm = TRUE), max(abs(diag(W))), max(abs(rowSums(W) - 1)), min(W), max(W)),
  stringsAsFactors = FALSE
)
write.csv(validation, file.path(OUT_DIR, "01_weight_validation.csv"), row.names = FALSE)

if (any(!is.finite(W))) stopf("Financial W contains non-finite values")
if (any(W < -1e-12)) stopf("Financial W contains negative values")
if (max(abs(diag(W))) > 1e-8) stopf("Financial W diagonal is not zero")
if (max(abs(rowSums(W) - 1)) > 1e-6) stopf("Financial W is not row-normalized")

# Top 5 counterparties by reporter.
top_rows <- list()
for (i in COUNTRIES) {
  z <- sort(W[i, setdiff(COUNTRIES, i)], decreasing = TRUE)
  kk <- seq_len(min(5L, length(z)))
  top_rows[[i]] <- data.frame(
    Country = i,
    Rank = kk,
    Counterpart = names(z)[kk],
    Weight = as.numeric(z[kk]),
    stringsAsFactors = FALSE
  )
}
top_df <- do.call(rbind, top_rows)
write.csv(top_df, file.path(OUT_DIR, "02_top_counterparts.csv"), row.names = FALSE)

concentration <- data.frame(
  Country = COUNTRIES,
  HHI = vapply(COUNTRIES, function(i) sum(W[i,]^2), numeric(1)),
  MaxWeight = vapply(COUNTRIES, function(i) max(W[i,]), numeric(1)),
  EffectivePartners = vapply(COUNTRIES, function(i) 1 / sum(W[i,]^2), numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(concentration, file.path(OUT_DIR, "03_weight_concentration.csv"), row.names = FALSE)

trade_candidates <- c(
  "8.12/Trade_Weights_14_Economies_2000_2014.csv",
  ".github/workflows/Trade_Weights_14_Economies_2000_2014.csv",
  "Trade_Weights_14_Economies_2000_2014.csv"
)
trade_path <- trade_candidates[file.exists(trade_candidates)]
comparison <- data.frame()
if (length(trade_path)) {
  WT <- read_weight_matrix(trade_path[[1]])
  diag(WT) <- 0
  WT <- WT / rowSums(WT)
  comparison <- do.call(rbind, lapply(COUNTRIES, function(i) {
    a <- W[i, setdiff(COUNTRIES, i)]
    b <- WT[i, setdiff(COUNTRIES, i)]
    data.frame(
      Country = i,
      Correlation = suppressWarnings(cor(a, b)),
      L1_Distance = sum(abs(a - b)),
      MaxAbsDifference = max(abs(a - b)),
      FinancialTop = names(which.max(W[i,])),
      TradeTop = names(which.max(WT[i,])),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(comparison, file.path(OUT_DIR, "04_financial_vs_trade_weight_shift.csv"), row.names = FALSE)
}

readme <- c(
  "Financial-weight preflight: PASS",
  paste0("Countries: ", paste(COUNTRIES, collapse = ", ")),
  paste0("Max |diagonal|: ", format(max(abs(diag(W))), scientific = TRUE)),
  paste0("Max |row sum - 1|: ", format(max(abs(rowSums(W) - 1)), scientific = TRUE)),
  paste0("Minimum weight: ", signif(min(W), 8)),
  paste0("Maximum weight: ", signif(max(W), 8)),
  paste0("Trade comparison available: ", ifelse(nrow(comparison) > 0, "YES", "NO")),
  "This step performs no model estimation and no imputation."
)
writeLines(readme, file.path(OUT_DIR, "README_preflight.txt"))
cat(paste(readme, collapse = "\n"), "\n")
