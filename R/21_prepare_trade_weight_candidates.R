#!/usr/bin/env Rscript

# =============================================================================
# 21_prepare_trade_weight_candidates.R
#
# Reconstruct and audit the two legacy trade-network candidates used by this
# project:
#
#   trade_2000_2012
#       Preferred trade baseline candidate. It ends before the documented
#       2013-2014 mirror-completion step in the later 2000-2014 construction.
#
#   trade_2000_2014_mirror
#       Robustness network. The legacy source code explicitly documents that
#       missing 2013-2014 bilateral trade observations were mirror-completed
#       before this final matrix was constructed.
#
# IMPORTANT
# ---------
# This script does NOT claim that every underlying 2000-2012 raw trade flow is
# directly observed. It establishes a narrower, auditable claim:
# the 2000-2012 matrix excludes the documented 2013-2014 mirror-completion step.
#
# It also records, but never uses as model input, the new IMF PIP / GCAP audit
# workbook already uploaded to data/weights.
# =============================================================================

suppressPackageStartupMessages(library(readxl))

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")

SOURCE_ROOT <- trimws(Sys.getenv("FIN3_TRADE_SOURCE_ROOT", "source_repo/8.12"))
OUT_WEIGHT <- "data/weights"
OUT <- "results/trade_weight_audit"
dir.create(OUT_WEIGHT, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...), call.=FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

find_one <- function(candidates, label) {
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stopf("Cannot locate %s. Tried: %s", label, paste(candidates, collapse=", "))
  hit[1]
}

src2012 <- find_one(
  file.path(
    SOURCE_ROOT,
    c(
      "Trade_Weights_14_Economies_2000_2012(2).xlsx",
      "Trade_Weights_14_Economies_2000_2012(2)(2).xlsx",
      "Trade_Weights_14_Economies_2000_2012.xlsx"
    )
  ),
  "2000-2012 trade-weight workbook"
)

src2014 <- find_one(
  file.path(
    SOURCE_ROOT,
    c(
      "Trade_Weights_14_Economies_2000_2014.csv",
      "Trade_Weights_14_Economies_2000_2014.xlsx"
    )
  ),
  "2000-2014 trade-weight matrix"
)

read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    d <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE,
                  fileEncoding="UTF-8-BOM")
    sheet_used <- NA_character_
  } else if (ext %in% c("xlsx","xls")) {
    ss <- readxl::excel_sheets(path)
    sheet_used <- if ("Trade_Weights" %in% ss) "Trade_Weights" else ss[1]
    d <- as.data.frame(
      suppressWarnings(readxl::read_excel(path, sheet=sheet_used, .name_repair="minimal")),
      check.names=FALSE, stringsAsFactors=FALSE
    )
  } else {
    stopf("Unsupported trade-weight source: %s", path)
  }
  names(d) <- trimws(sub("^\ufeff","",names(d)))
  attr(d, "sheet_used") <- sheet_used
  d
}

extract_matrix <- function(path, label) {
  d <- read_table(path)
  if (ncol(d) < 15L) stopf("%s source is too narrow.", label)

  rr <- toupper(trimws(as.character(d[[1]])))
  cc <- toupper(trimws(names(d)[-1]))

  if (anyDuplicated(rr)) stopf("%s has duplicate reporter rows.", label)
  if (anyDuplicated(cc)) stopf("%s has duplicate counterpart columns.", label)
  if (!all(COUNTRIES %in% rr) || !all(COUNTRIES %in% cc)) {
    stopf("%s does not contain the exact 14-economy matrix.", label)
  }

  W <- matrix(
    NA_real_, length(COUNTRIES), length(COUNTRIES),
    dimnames=list(COUNTRIES, COUNTRIES)
  )
  for (i in COUNTRIES) {
    ri <- match(i, rr)
    for (j in COUNTRIES) W[i,j] <- num(d[[match(j, cc)+1L]][ri])
  }

  if (any(!is.finite(W))) stopf("%s contains non-finite cells.", label)
  if (any(W < -1e-12)) stopf("%s contains negative weights.", label)

  maxdiag <- max(abs(diag(W)))
  rs <- rowSums(W)
  maxdev <- max(abs(rs - 1))

  if (maxdiag > 1e-8) stopf("%s has non-zero diagonal.", label)
  if (any(rs <= 0)) stopf("%s has non-positive row sum.", label)
  if (maxdev > 1e-6) {
    stopf("%s is not already a row-normalized final matrix; max row-sum deviation=%.8g", label, maxdev)
  }

  # Remove floating-point drift only.
  diag(W) <- 0
  W <- W / rowSums(W)

  list(
    W=W,
    sheet=attr(d, "sheet_used"),
    md5=unname(tools::md5sum(path)),
    size=file.info(path)$size,
    maxdiag=maxdiag,
    max_row_sum_deviation=maxdev
  )
}

a12 <- extract_matrix(src2012, "trade_2000_2012")
a14 <- extract_matrix(src2014, "trade_2000_2014_mirror")

write_matrix <- function(W, path) {
  write.csv(
    data.frame(Country=rownames(W), W, check.names=FALSE),
    path, row.names=FALSE, quote=FALSE
  )
}

p12 <- file.path(OUT_WEIGHT, "W_trade_2000_2012.csv")
p14 <- file.path(OUT_WEIGHT, "W_trade_2000_2014_mirror.csv")
write_matrix(a12$W, p12)
write_matrix(a14$W, p14)

# Existing financial network comparison.
fin_path <- file.path(OUT_WEIGHT, "W_main_restated_exdom_2017.csv")
if (!file.exists(fin_path)) stopf("Missing committed GCAP main matrix: %s", fin_path)
fin <- extract_matrix(fin_path, "main_restated_exdom_2017")$W

matrix_distance <- function(A, B, a, b) {
  D <- A-B
  ij <- which(abs(D)==max(abs(D)), arr.ind=TRUE)[1,]
  data.frame(
    NetworkA=a,
    NetworkB=b,
    MaxAbsDifference=max(abs(D)),
    MeanAbsDifference=mean(abs(D)),
    FrobeniusNorm=sqrt(sum(D^2)),
    LargestDifferenceFrom=rownames(D)[ij[1]],
    LargestDifferenceTo=colnames(D)[ij[2]],
    stringsAsFactors=FALSE
  )
}

dist <- rbind(
  matrix_distance(a12$W, fin, "trade_2000_2012", "main_restated_exdom_2017"),
  matrix_distance(a14$W, fin, "trade_2000_2014_mirror", "main_restated_exdom_2017"),
  matrix_distance(a12$W, a14$W, "trade_2000_2012", "trade_2000_2014_mirror")
)
write.csv(dist, file.path(OUT, "03_network_distance.csv"), row.names=FALSE)

validation <- rbind(
  data.frame(
    Network="trade_2000_2012",
    Source=src2012,
    Sheet=ifelse(is.na(a12$sheet),"",a12$sheet),
    Rows=14, Columns=14,
    MaxAbsDiagonal=a12$maxdiag,
    MaxAbsRowSumMinus1=a12$max_row_sum_deviation,
    Nonfinite=0, Negative=0,
    Status="OK",
    stringsAsFactors=FALSE
  ),
  data.frame(
    Network="trade_2000_2014_mirror",
    Source=src2014,
    Sheet=ifelse(is.na(a14$sheet),"",a14$sheet),
    Rows=14, Columns=14,
    MaxAbsDiagonal=a14$maxdiag,
    MaxAbsRowSumMinus1=a14$max_row_sum_deviation,
    Nonfinite=0, Negative=0,
    Status="OK",
    stringsAsFactors=FALSE
  )
)
write.csv(validation, file.path(OUT, "01_trade_matrix_validation.csv"), row.names=FALSE)

lineage <- data.frame(
  Network=c("trade_2000_2012","trade_2000_2014_mirror"),
  SourceFile=c(src2012,src2014),
  SourceMD5=c(a12$md5,a14$md5),
  AveragingWindow=c("2000-2012","2000-2014"),
  Documented2013_2014MirrorCompletionIncluded=c(FALSE,TRUE),
  ProjectImputationPerformedByThisScript=c(FALSE,FALSE),
  PermittedRole=c(
    "BASELINE_CANDIDATE__EXCLUDES_DOCUMENTED_2013_2014_MIRROR_STEP",
    "ROBUSTNESS_ONLY__DOCUMENTED_2013_2014_MIRROR_COMPLETION"
  ),
  ImportantLimitation=c(
    "This audit does not infer that every underlying pre-2013 raw trade observation was directly observed; it verifies the final 14x14 matrix and excludes the documented 2013-2014 mirror-completion step.",
    "The legacy construction explicitly states that missing 2013-2014 bilateral trade observations were mirror-completed before the final matrix was built."
  ),
  stringsAsFactors=FALSE
)
write.csv(lineage, file.path(OUT_WEIGHT, "trade_weight_lineage.csv"), row.names=FALSE)
write.csv(lineage, file.path(OUT, "02_trade_weight_lineage.csv"), row.names=FALSE)

# Register the user's new PIP/GCAP audit evidence without using it as model input.
audit_xlsx <- list.files(
  OUT_WEIGHT,
  pattern="^IMF_PIP_GCAP_financial_weight_audit_2015_2024.*[.]xlsx$",
  full.names=TRUE
)
audit_rows <- data.frame()
if (length(audit_xlsx)) {
  audit_rows <- do.call(rbind, lapply(audit_xlsx, function(f) {
    data.frame(
      File=f,
      SizeBytes=file.info(f)$size,
      MD5=unname(tools::md5sum(f)),
      Role="AUDIT_EVIDENCE_ONLY__NOT_MODEL_WEIGHT_INPUT",
      stringsAsFactors=FALSE
    )
  }))
}
write.csv(
  audit_rows,
  file.path(OUT, "04_pip_gcap_audit_workbook_registry.csv"),
  row.names=FALSE
)

# Validate newly uploaded GCAP validation evidence if present.
gcap_val <- file.path(OUT_WEIGHT, "precomputed_validation_checks.csv")
gcap_validation_ok <- NA
if (file.exists(gcap_val)) {
  gv <- read.csv(gcap_val, stringsAsFactors=FALSE, check.names=FALSE)
  if (!"Status" %in% names(gv)) stopf("Malformed precomputed_validation_checks.csv")
  gcap_validation_ok <- all(trimws(gv$Status)=="OK")
  if (!gcap_validation_ok) stopf("At least one uploaded GCAP validation row is not OK.")
}

gate <- data.frame(
  Status="READY_FOR_TRADE_2000_2012_MODEL",
  PreferredTradeNetwork="trade_2000_2012",
  TradeRobustnessNetwork="trade_2000_2014_mirror",
  PreferredWeightFile=p12,
  RobustnessWeightFile=p14,
  PreferredMatrixFinite=all(is.finite(a12$W)),
  PreferredMatrixNonnegative=all(a12$W>=-1e-12),
  PreferredMatrixZeroDiagonal=max(abs(diag(a12$W)))<=1e-8,
  PreferredMatrixRowNormalized=max(abs(rowSums(a12$W)-1))<=1e-8,
  Documented2013_2014MirrorCompletionIncludedInPreferred=FALSE,
  ProjectImputationPerformed=FALSE,
  UploadedPIPGCAPAuditWorkbookRegistered=length(audit_xlsx)>=1L,
  UploadedGCAPValidationEvidencePass=gcap_validation_ok,
  stringsAsFactors=FALSE
)
write.csv(gate, file.path(OUT, "00_trade_weight_gate.csv"), row.names=FALSE)

txt <- c(
  "TRADE WEIGHT PREPARATION: READY_FOR_TRADE_2000_2012_MODEL",
  "==========================================================",
  sprintf("Preferred: %s", p12),
  sprintf("Robustness: %s", p14),
  "",
  "The preferred 2000-2012 matrix ends before the documented 2013-2014 mirror-completion step.",
  "No missing trade cell is created or imputed by this script.",
  "The PIP/GCAP audit workbook is registered as audit evidence only and is not fed into the model.",
  "",
  sprintf("Trade-2012 vs GCAP max absolute cell difference: %.6f",
          dist$MaxAbsDifference[dist$NetworkA=="trade_2000_2012" & dist$NetworkB=="main_restated_exdom_2017"])
)
writeLines(txt, file.path(OUT, "README_trade_weight_audit.txt"))
cat(paste(txt, collapse="\n"), "\n")
