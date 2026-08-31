#!/usr/bin/env Rscript

# ============================================================
# 20_build_financial_3var_input.R
# Standalone repository version.
#
# Domestic variables:
#   r   = short-term interest-rate level
#   de  = REER log change
#   deq = equity return
#
# Financial network:
#   GCAP 2017 residence-based bilateral portfolio weights
#
# Source macro/GPR/Brent files are read from SOURCE_DIR, which
# the GitHub Actions workflow populates from the legacy project.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("r","de","deq")
SOURCE_SUFFIX <- c(r="RATE_LEVEL", de="REER_DLOG", deq="EQ_RETURN")

SOURCE_DIR   <- Sys.getenv("FIN3_SOURCE_DIR", "source_repo/8.12")
MACRO_PATH   <- Sys.getenv("FIN3_MACRO_XLSX",
                           file.path(SOURCE_DIR, "TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx"))
GPR_PATH     <- Sys.getenv("FIN3_GPR_CSV",
                           file.path(SOURCE_DIR, "gpr_quarterly_processed.csv"))
OIL_PATH     <- Sys.getenv("FIN3_OIL_CSV",
                           file.path(SOURCE_DIR, "IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"))
WEIGHT_PATH  <- Sys.getenv("FIN3_WEIGHT_CSV", "data/GCAP_financial_W_2017.csv")

SAMPLE_START <- Sys.getenv("FIN3_SAMPLE_START", "2002Q2")
SAMPLE_END   <- Sys.getenv("FIN3_SAMPLE_END", "2022Q4")
GPR_COLUMN   <- Sys.getenv("FIN3_GPR_COLUMN", "LN_GPR_QMEAN")

OUT_DIR <- "data/derived"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...), call.=FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
msg <- function(...) cat(sprintf(...), "\n")

quarter_id <- function(x) {
  sx <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))
  m <- regexec("^([12][0-9]{3})[^0-9]*Q([1-4])$", sx)
  mm <- regmatches(sx, m)
  ok <- lengths(mm) == 3L
  if (any(ok)) {
    yr <- as.integer(vapply(mm[ok], `[`, character(1), 2))
    qq <- as.integer(vapply(mm[ok], `[`, character(1), 3))
    out[ok] <- 4L*yr + qq
  }
  out
}
quarter_label <- function(qid) {
  yr <- (qid-1L) %/% 4L
  qq <- qid - 4L*yr
  sprintf("%dQ%d", yr, qq)
}

for (p in c(MACRO_PATH, GPR_PATH, WEIGHT_PATH)) {
  if (!file.exists(p)) stopf("Required input not found: %s", p)
}
if (!requireNamespace("readxl", quietly=TRUE)) stopf("Package readxl is required")

# ---- 1. Macro data ----
sheets <- readxl::excel_sheets(MACRO_PATH)
required_cols <- c(
  "Quarter",
  unlist(lapply(COUNTRIES, function(cc) paste0(cc, "_", unname(SOURCE_SUFFIX))), use.names=FALSE)
)

preferred <- c("MODEL_COMPLETE_14C","MODEL_WIDE_ALL")
sheet_order <- c(intersect(preferred, sheets), setdiff(sheets, preferred))
macro <- NULL
macro_sheet <- NA_character_

for (ss in sheet_order) {
  d <- tryCatch(
    as.data.frame(readxl::read_excel(MACRO_PATH, sheet=ss, .name_repair="minimal"),
                  check.names=FALSE),
    error=function(e) NULL
  )
  if (is.null(d)) next
  names(d) <- trimws(names(d))
  if (all(required_cols %in% names(d))) {
    macro <- d
    macro_sheet <- ss
    break
  }
}
if (is.null(macro)) stopf("No sheet contains exact RATE_LEVEL / REER_DLOG / EQ_RETURN mapping")

qid <- quarter_id(macro$Quarter)
if (mean(!is.na(qid)) < .95) stopf("Quarter parsing failed in sheet %s", macro_sheet)

wide <- data.frame(qid=qid, Quarter=quarter_label(qid), check.names=FALSE)
for (cc in COUNTRIES) {
  for (v in VARS) {
    wide[[paste0(cc,"_",v)]] <- num(macro[[paste0(cc,"_",SOURCE_SUFFIX[[v]])]])
  }
}
wide <- wide[!is.na(wide$qid), , drop=FALSE]
wide <- wide[!duplicated(wide$qid), , drop=FALSE]
wide <- wide[order(wide$qid), , drop=FALSE]

q0 <- quarter_id(SAMPLE_START)
q1 <- quarter_id(SAMPLE_END)
wide <- wide[wide$qid >= q0 & wide$qid <= q1, , drop=FALSE]
if (!nrow(wide)) stopf("No observations in requested sample")
if (any(diff(wide$qid) != 1L)) stopf("Requested sample contains missing quarters")

macro_cols <- unlist(lapply(COUNTRIES, function(cc) paste0(cc,"_",VARS)), use.names=FALSE)
bad <- which(!complete.cases(wide[,macro_cols,drop=FALSE]))
if (length(bad)) stopf("3-variable sample still has NA; first affected quarter = %s", wide$Quarter[bad[1]])

# ---- 2. GCAP W ----
wtab <- read.csv(WEIGHT_PATH, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
names(wtab) <- trimws(names(wtab))
rows <- toupper(trimws(as.character(wtab[[1]])))
cols <- toupper(trimws(names(wtab)[-1]))

if (!all(COUNTRIES %in% rows)) stopf("Missing weight rows: %s",
                                     paste(setdiff(COUNTRIES,rows),collapse=", "))
if (!all(COUNTRIES %in% cols)) stopf("Missing weight cols: %s",
                                     paste(setdiff(COUNTRIES,cols),collapse=", "))

W <- matrix(NA_real_, length(COUNTRIES), length(COUNTRIES),
            dimnames=list(COUNTRIES,COUNTRIES))
for (i in COUNTRIES) {
  rr <- match(i, rows)
  for (j in COUNTRIES) W[i,j] <- num(wtab[[match(j,cols)+1L]][rr])
}

if (any(!is.finite(W)) || any(W < -1e-12)) stopf("Invalid financial weights")
if (max(abs(diag(W))) > 1e-8) stopf("Financial-weight diagonal is not zero")
if (max(abs(rowSums(W)-1)) > 1e-6) stopf("Financial weights are not row normalized")
diag(W) <- 0
W <- W / rowSums(W)

# ---- 3. Global quarterly series ----
read_qcsv <- function(path, exact=NULL, label="series") {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  qr <- vapply(d, function(z) mean(!is.na(quarter_id(z))), numeric(1))
  qr[!is.finite(qr)] <- 0
  dc <- which.max(qr)
  if (!length(dc) || qr[dc] < .5) stopf("Cannot detect quarter column in %s", label)

  vc <- NA_integer_
  if (!is.null(exact) && exact %in% names(d)) vc <- match(exact,names(d))
  if (is.na(vc)) {
    cand <- setdiff(seq_along(d),dc)
    rate <- vapply(cand, function(j) mean(is.finite(num(d[[j]]))), numeric(1))
    if (!length(rate) || max(rate) < .5) stopf("Cannot detect value column in %s", label)
    vc <- cand[which.max(rate)]
  }
  z <- data.frame(qid=quarter_id(d[[dc]]), value=num(d[[vc]]))
  z <- z[!is.na(z$qid) & is.finite(z$value),,drop=FALSE]
  aggregate(value ~ qid, z, mean)
}

gpr <- read_qcsv(GPR_PATH, GPR_COLUMN, "GPR")
wide$gpr <- gpr$value[match(wide$qid,gpr$qid)]
if (any(!is.finite(wide$gpr))) stopf("GPR does not fully cover sample")

oil <- read_qcsv(OIL_PATH, NULL, "Brent")
if (is.null(oil)) {
  wide$brent <- NA_real_
  warning("Brent input not found; continuing without usable Brent")
} else {
  wide$brent <- oil$value[match(wide$qid,oil$qid)]
  if (any(!is.finite(wide$brent))) stopf("Brent does not fully cover sample")
}

# ---- 4. Foreign variables ----
for (v in VARS) {
  Xv <- do.call(cbind, lapply(COUNTRIES, function(cc) wide[[paste0(cc,"_",v)]]))
  colnames(Xv) <- COUNTRIES
  for (i in COUNTRIES) {
    wide[[paste0(i,"_",v,"_star")]] <- as.numeric(Xv %*% W[i,])
  }
}

out_cols <- c(
  "Quarter",
  unlist(lapply(COUNTRIES, function(cc)
    c(paste0(cc,"_",VARS),paste0(cc,"_",VARS,"_star"))), use.names=FALSE),
  "gpr","brent"
)
model <- wide[,out_cols,drop=FALSE]

long <- do.call(rbind,lapply(COUNTRIES,function(cc) {
  data.frame(
    Quarter=wide$Quarter, Country=cc,
    r=wide[[paste0(cc,"_r")]],
    de=wide[[paste0(cc,"_de")]],
    deq=wide[[paste0(cc,"_deq")]],
    r_star=wide[[paste0(cc,"_r_star")]],
    de_star=wide[[paste0(cc,"_de_star")]],
    deq_star=wide[[paste0(cc,"_deq_star")]],
    gpr=wide$gpr, brent=wide$brent,
    stringsAsFactors=FALSE
  )
}))

write.csv(model,file.path(OUT_DIR,"model_input_fin3.csv"),row.names=FALSE,na="")
write.csv(long,file.path(OUT_DIR,"panel_fin3_long.csv"),row.names=FALSE,na="")
write.csv(data.frame(Country=COUNTRIES,W,check.names=FALSE),
          file.path(OUT_DIR,"financial_weights.csv"),row.names=FALSE,na="")

summary <- c(
  "14-economy 3-variable financial GVAR input",
  paste0("Macro source: ",MACRO_PATH),
  paste0("Macro sheet: ",macro_sheet),
  paste0("Weight source: ",WEIGHT_PATH),
  paste0("Sample: ",model$Quarter[1]," - ",tail(model$Quarter,1)),
  paste0("T = ",nrow(model)),
  "Variables: r, de, deq",
  "Foreign variables: r*, de*, deq* using GCAP 2017 W",
  paste0("GPR: ",GPR_COLUMN),
  paste0("Brent included: ",all(is.finite(model$brent))),
  paste0("Max |diag(W)| = ",format(max(abs(diag(W))),scientific=TRUE)),
  paste0("Max |rowsum(W)-1| = ",format(max(abs(rowSums(W)-1)),scientific=TRUE)),
  "No new interpolation/imputation is performed by this script."
)
writeLines(summary,file.path(OUT_DIR,"build_summary.txt"))
cat(paste(summary,collapse="\n"),"\n")
