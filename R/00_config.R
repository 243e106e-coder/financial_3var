options(stringsAsFactors = FALSE, warn = 1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("r","de","deq")
SOURCE_SUFFIX <- c(r="RATE_LEVEL", de="REER_DLOG", deq="EQ_RETURN")

MAIN_NETWORK <- "main_restated_exdom_2017"
WEIGHT_FILES <- c(
  main_restated_exdom_2017 = "data/weights/W_main_restated_exdom_2017.csv",
  residency_2017 = "data/weights/W_residency_2017.csv",
  equity_restated_2017 = "data/weights/W_equity_restated_2017.csv",
  bonds_restated_2017 = "data/weights/W_bonds_restated_2017.csv"
)

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
TRADE_WEIGHT_PATH <- Sys.getenv(
  "FIN3_TRADE_WEIGHT_CSV",
  file.path(SOURCE_DIR, "Trade_Weights_14_Economies_2000_2014.csv")
)
GPR_COLUMN <- Sys.getenv("FIN3_GPR_COLUMN", "LN_GPR_QMEAN")
SAMPLE_START_OVERRIDE <- trimws(Sys.getenv("FIN3_SAMPLE_START", ""))
SAMPLE_END_OVERRIDE <- trimws(Sys.getenv("FIN3_SAMPLE_END", ""))

DERIVED_DIR <- "data/derived"
RESULTS_DIR <- "results"

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

read_weight_matrix <- function(path, normalize = TRUE) {
  if (!file.exists(path)) stopf("Weight file not found: %s", path)
  d <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE,
                fileEncoding="UTF-8-BOM")
  names(d) <- trimws(names(d))
  rr <- toupper(trimws(as.character(d[[1]])))
  cc <- toupper(trimws(names(d)[-1]))
  if (!all(COUNTRIES %in% rr) || !all(COUNTRIES %in% cc)) {
    stopf("Weight matrix does not contain the exact 14-economy sample: %s", path)
  }
  W <- matrix(NA_real_, length(COUNTRIES), length(COUNTRIES),
              dimnames=list(COUNTRIES, COUNTRIES))
  for (i in COUNTRIES) {
    ri <- match(i, rr)
    for (j in COUNTRIES) W[i,j] <- num(d[[match(j,cc)+1L]][ri])
  }
  if (any(!is.finite(W))) stopf("Non-finite financial weight in %s", path)
  if (any(W < -1e-12)) stopf("Negative financial weight in %s", path)
  if (max(abs(diag(W))) > 1e-8) stopf("Non-zero diagonal in %s", path)
  if (any(rowSums(W) <= 0)) stopf("Non-positive row sum in %s", path)
  if (normalize) {
    diag(W) <- 0
    W <- W / rowSums(W)
  }
  W
}
