#!/usr/bin/env Rscript

# =============================================================================
# 23_gcap_multiyear_coverage_audit.R
#
# Purpose
#   1) Download/read the official GCAP Restated Bilateral External Portfolios.
#   2) Audit investor-year coverage for the 14-economy financial GVAR sample.
#   3) Identify STRICT common years that are methodologically comparable with
#      the current 2017 financial-weight reconstruction:
#        - United States: Fund Holdings; E + BC + BG + BSF
#        - All other investors: Issuance; EF + B
#   4) Rebuild yearly 14x14 financial-weight matrices for every strict common year.
#   5) If >= 2 strict common years exist, construct a multi-year average W.
#      If only one common year exists, DO NOT manufacture a multi-year matrix.
#   6) Compare rebuilt 2017 W with data/GCAP_financial_W_2017.csv if present.
#
# Outputs
#   results/gcap_coverage/01_investor_year_coverage.csv
#   results/gcap_coverage/02_common_years.csv
#   results/gcap_coverage/03_pair_coverage_by_year.csv
#   results/gcap_coverage/04_rebuilt_weight_validation.csv
#   results/gcap_coverage/05_2017_rebuild_comparison.csv   (if reference exists)
#   results/gcap_coverage/README_gcap_coverage.txt
#   results/gcap_coverage/yearly_weights/W_GCAP_<year>.csv
#   data/GCAP_financial_W_multiyear_average.csv             (only if >=2 common years)
#
# Reproducibility
#   Official source:
#   https://globalcapitalallocation.s3.us-east-2.amazonaws.com/
#     Restated_Bilateral_External_Portfolios.dta
#
# Notes
#   - No new interpolation or extrapolation is performed here.
#   - Position_Residency is used to remain consistent with the current 2017 W.
#   - EA counterpart is a fixed EA19 issuer aggregation, matching the current file.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# ---------- helpers -----------------------------------------------------------

msg <- function(...) cat(sprintf(...), "\n")

clean_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  make.unique(x, sep = "_")
}

pick_col <- function(nms, candidates, label) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0L) {
    stop(
      sprintf(
        "Cannot find required column '%s'. Available columns: %s",
        label, paste(nms, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  hit[1L]
}

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

row_normalize <- function(Fmat) {
  W <- Fmat
  diag(W) <- 0
  rs <- rowSums(W, na.rm = TRUE)
  if (any(!is.finite(rs) | rs <= 0)) {
    bad <- rownames(W)[!is.finite(rs) | rs <= 0]
    stop("Non-positive/non-finite row totals for: ", paste(bad, collapse = ", "))
  }
  W <- W / rs
  diag(W) <- 0
  W
}

write_matrix_csv <- function(M, path) {
  out <- data.frame(Country = rownames(M), M, check.names = FALSE)
  utils::write.csv(out, path, row.names = FALSE, na = "")
}

read_weight_csv <- function(path, countries) {
  z <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  first <- names(z)[1]
  rn <- as.character(z[[first]])
  z[[first]] <- NULL
  M <- as.matrix(z)
  storage.mode(M) <- "double"
  rownames(M) <- rn

  if (!all(countries %in% rownames(M)) || !all(countries %in% colnames(M))) {
    stop("Reference weight matrix does not contain the full 14-country sample.")
  }
  M[countries, countries, drop = FALSE]
}

longest_consecutive_run <- function(years) {
  years <- sort(unique(as.integer(years)))
  if (length(years) == 0L) return(integer(0))
  groups <- cumsum(c(TRUE, diff(years) != 1L))
  runs <- split(years, groups)
  runs[[which.max(lengths(runs))]]
}

# ---------- sample definition -------------------------------------------------

countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")

investor_map <- c(
  AU="AUS", BR="BRA", CA="CAN", CH="CHE", CN="CHN", EA="EMU",
  UK="GBR", JP="JPN", KR="KOR", NO="NOR", SG="SGP", TR="TUR",
  US="USA", ZA="ZAF"
)

# Fixed EA19 issuer aggregation, matching the existing 2017 reconstruction.
ea19 <- c(
  "AUT","BEL","CYP","EST","FIN","FRA","DEU","GRC","IRL","ITA",
  "LVA","LTU","LUX","MLT","NLD","PRT","SVK","SVN","ESP"
)

issuer_map <- list(
  AU="AUS", BR="BRA", CA="CAN", CH="CHE", CN="CHN", EA=ea19,
  UK="GBR", JP="JPN", KR="KOR", NO="NOR", SG="SGP", TR="TUR",
  US="USA", ZA="ZAF"
)

# Strict methodology chosen to reproduce the current 2017 construction.
strict_method <- setNames(rep("Issuance", length(countries)), countries)
strict_method["US"] <- "Fund Holdings"

required_assets <- lapply(countries, function(cc) {
  if (cc == "US") c("E","BC","BG","BSF") else c("EF","B")
})
names(required_assets) <- countries

# ---------- directories -------------------------------------------------------

out_dir <- Sys.getenv("OUT_DIR", unset = "results/gcap_coverage")
yearly_dir <- file.path(out_dir, "yearly_weights")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(yearly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

# ---------- input discovery / official download -------------------------------

official_url <- paste0(
  "https://globalcapitalallocation.s3.us-east-2.amazonaws.com/",
  "Restated_Bilateral_External_Portfolios.dta"
)

gcap_file <- Sys.getenv("GCAP_FILE", unset = "")

candidates <- c(
  gcap_file,
  "data/raw/Restated_Bilateral_External_Portfolios.dta",
  "data/Restated_Bilateral_External_Portfolios.dta",
  "Restated_Bilateral_External_Portfolios.dta",
  "data/raw/Restated_Bilateral_External_Portfolios.xlsx",
  "data/Restated_Bilateral_External_Portfolios.xlsx",
  "Restated_Bilateral_External_Portfolios.xlsx"
)
candidates <- unique(candidates[nzchar(candidates)])
existing <- candidates[file.exists(candidates)]

downloaded_temp <- FALSE
if (length(existing) > 0L) {
  input_file <- existing[1L]
  msg("Using local GCAP source: %s", input_file)
} else {
  input_file <- tempfile(fileext = ".dta")
  msg("Local GCAP source not found. Downloading official DTA...")
  msg("Source: %s", official_url)
  ok <- tryCatch({
    utils::download.file(official_url, input_file, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) {
    msg("Download failed: %s", conditionMessage(e))
    FALSE
  })

  if (!ok || !file.exists(input_file) || file.info(input_file)$size <= 0) {
    stop(
      paste0(
        "Could not obtain the GCAP raw file.\n",
        "Either commit Restated_Bilateral_External_Portfolios.dta to data/raw/\n",
        "or set GCAP_FILE to its path."
      ),
      call. = FALSE
    )
  }
  downloaded_temp <- TRUE
}

# ---------- read source -------------------------------------------------------

ext <- tolower(tools::file_ext(input_file))

if (ext == "dta") {
  if (!requireNamespace("haven", quietly = TRUE)) {
    stop("Package 'haven' is required for .dta input.", call. = FALSE)
  }
  raw <- as.data.frame(haven::read_dta(input_file))
} else if (ext %in% c("xlsx","xls")) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required for Excel input.", call. = FALSE)
  }
  sheets <- readxl::excel_sheets(input_file)
  if (length(sheets) == 0L) stop("Excel file has no sheets.", call. = FALSE)

  # Read the largest-looking sheet among the first few sheets.
  tmp <- lapply(sheets, function(sh) {
    tryCatch(as.data.frame(readxl::read_excel(input_file, sheet = sh)),
             error = function(e) NULL)
  })
  sizes <- vapply(tmp, function(x) if (is.null(x)) -1 else nrow(x), numeric(1))
  raw <- tmp[[which.max(sizes)]]
} else if (ext == "csv") {
  raw <- utils::read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  stop("Unsupported GCAP input format: ", ext, call. = FALSE)
}

if (nrow(raw) == 0L) stop("GCAP input has zero rows.", call. = FALSE)

names(raw) <- clean_names(names(raw))
nms <- names(raw)

col_method <- pick_col(nms, c("methodology"), "Methodology")
col_year   <- pick_col(nms, c("year"), "Year")
col_inv    <- pick_col(nms, c("investor"), "Investor")
col_asset  <- pick_col(nms, c("asset_class_code","assetclasscode"), "Asset Class Code")
col_issuer <- pick_col(nms, c("issuer"), "Issuer")
col_pos    <- pick_col(
  nms,
  c("position_residency","position_residence","position_residency_"),
  "Position Residency"
)

d <- data.frame(
  methodology = trimws(as.character(raw[[col_method]])),
  year        = as.integer(safe_numeric(raw[[col_year]])),
  investor    = toupper(trimws(as.character(raw[[col_inv]]))),
  asset_code  = toupper(trimws(as.character(raw[[col_asset]]))),
  issuer      = toupper(trimws(as.character(raw[[col_issuer]]))),
  position    = safe_numeric(raw[[col_pos]]),
  stringsAsFactors = FALSE
)

d <- d[
  !is.na(d$year) &
  nzchar(d$investor) &
  nzchar(d$asset_code) &
  nzchar(d$issuer),
  , drop = FALSE
]

model_from_iso <- setNames(names(investor_map), unname(investor_map))
d$model_country <- unname(model_from_iso[d$investor])
d_sample <- d[!is.na(d$model_country), , drop = FALSE]

if (nrow(d_sample) == 0L) {
  stop("None of the 14 model investor codes were found in the GCAP file.", call. = FALSE)
}

all_years <- sort(unique(d_sample$year))
msg("Raw GCAP years found for sample investors: %s",
    paste(range(all_years), collapse = " - "))

# ---------- strict investor-year coverage ------------------------------------

coverage_rows <- list()
k <- 0L

for (cc in countries) {
  iso <- investor_map[[cc]]
  meth <- strict_method[[cc]]
  req <- required_assets[[cc]]

  for (yy in all_years) {
    z <- d_sample[
      d_sample$investor == iso &
      d_sample$year == yy &
      d_sample$methodology == meth,
      , drop = FALSE
    ]

    assets_present <- sort(unique(z$asset_code))
    complete <- all(req %in% assets_present)

    # Loose coverage: same required asset scheme, any methodology.
    za <- d_sample[
      d_sample$investor == iso &
      d_sample$year == yy,
      , drop = FALSE
    ]
    loose_complete <- all(req %in% sort(unique(za$asset_code)))

    k <- k + 1L
    coverage_rows[[k]] <- data.frame(
      Country = cc,
      InvestorISO = iso,
      Year = yy,
      StrictMethodology = meth,
      RequiredAssets = paste(req, collapse = "+"),
      StrictRows = nrow(z),
      StrictIssuerCount = length(unique(z$issuer)),
      StrictPositionNA = sum(is.na(z$position)),
      StrictComplete = complete,
      AnyMethodRows = nrow(za),
      AnyMethodologies = paste(sort(unique(za$methodology)), collapse = " | "),
      LooseAssetComplete = loose_complete,
      stringsAsFactors = FALSE
    )
  }
}

coverage <- do.call(rbind, coverage_rows)
utils::write.csv(
  coverage,
  file.path(out_dir, "01_investor_year_coverage.csv"),
  row.names = FALSE,
  na = ""
)

common_summary <- do.call(rbind, lapply(all_years, function(yy) {
  z <- coverage[coverage$Year == yy, , drop = FALSE]
  data.frame(
    Year = yy,
    StrictCompleteInvestors = sum(z$StrictComplete),
    LooseCompleteInvestors = sum(z$LooseAssetComplete),
    StrictAll14 = all(z$StrictComplete),
    LooseAll14 = all(z$LooseAssetComplete),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(
  common_summary,
  file.path(out_dir, "02_common_years.csv"),
  row.names = FALSE
)

strict_common_years <- common_summary$Year[common_summary$StrictAll14]
loose_common_years  <- common_summary$Year[common_summary$LooseAll14]
strict_run <- longest_consecutive_run(strict_common_years)

# ---------- matrix builder ----------------------------------------------------

build_position_matrix <- function(yy) {
  F <- matrix(
    0,
    nrow = length(countries),
    ncol = length(countries),
    dimnames = list(countries, countries)
  )

  diag_na_counts <- 0L
  used_rows <- 0L

  for (i in countries) {
    inv_iso <- investor_map[[i]]
    meth <- strict_method[[i]]
    req <- required_assets[[i]]

    zi <- d_sample[
      d_sample$investor == inv_iso &
      d_sample$year == yy &
      d_sample$methodology == meth &
      d_sample$asset_code %in% req,
      , drop = FALSE
    ]

    for (j in countries) {
      if (i == j) {
        F[i, j] <- 0
        next
      }

      issuer_codes <- issuer_map[[j]]
      zij <- zi[zi$issuer %in% issuer_codes, , drop = FALSE]
      used_rows <- used_rows + nrow(zij)

      if (nrow(zij) == 0L) {
        F[i, j] <- 0
      } else {
        diag_na_counts <- diag_na_counts + sum(is.na(zij$position))
        vals <- zij$position[is.finite(zij$position)]
        F[i, j] <- if (length(vals) == 0L) 0 else sum(vals)
      }
    }
  }

  list(
    F = F,
    W = row_normalize(F),
    n_position_na_in_sample_pairs = diag_na_counts,
    used_rows = used_rows
  )
}

# ---------- build all strict common-year matrices -----------------------------

validation <- list()
pair_coverage <- list()
weights_by_year <- list()

if (length(strict_common_years) == 0L) {
  warning("No strict common year across all 14 economies.")
} else {
  for (yy in strict_common_years) {
    built <- build_position_matrix(yy)
    F <- built$F
    W <- built$W
    weights_by_year[[as.character(yy)]] <- W

    write_matrix_csv(
      W,
      file.path(yearly_dir, sprintf("W_GCAP_%d.csv", yy))
    )

    offdiag <- row(F) != col(F)
    pair_coverage[[as.character(yy)]] <- data.frame(
      Year = yy,
      SampleDirectedPairs = sum(offdiag),
      PositivePositionPairs = sum(F[offdiag] > 0, na.rm = TRUE),
      ZeroPositionPairs = sum(F[offdiag] == 0, na.rm = TRUE),
      NonfinitePositionPairs = sum(!is.finite(F[offdiag])),
      PositionNAEncounteredInSourceRows = built$n_position_na_in_sample_pairs,
      stringsAsFactors = FALSE
    )

    validation[[as.character(yy)]] <- data.frame(
      Year = yy,
      Rows = nrow(W),
      Columns = ncol(W),
      Nonfinite = sum(!is.finite(W)),
      Negative = sum(W < -1e-14, na.rm = TRUE),
      MaxAbsDiagonal = max(abs(diag(W))),
      MaxAbsRowSumMinus1 = max(abs(rowSums(W) - 1)),
      MinWeight = min(W),
      MaxWeight = max(W),
      stringsAsFactors = FALSE
    )
  }
}

pair_cov_df <- if (length(pair_coverage)) do.call(rbind, pair_coverage) else data.frame()
val_df <- if (length(validation)) do.call(rbind, validation) else data.frame()

utils::write.csv(
  pair_cov_df,
  file.path(out_dir, "03_pair_coverage_by_year.csv"),
  row.names = FALSE
)
utils::write.csv(
  val_df,
  file.path(out_dir, "04_rebuilt_weight_validation.csv"),
  row.names = FALSE
)

# ---------- compare rebuilt 2017 with current repo matrix ---------------------

ref_path <- "data/GCAP_financial_W_2017.csv"
comparison_2017 <- NULL

if ("2017" %in% names(weights_by_year) && file.exists(ref_path)) {
  W_new <- weights_by_year[["2017"]]
  W_ref <- read_weight_csv(ref_path, countries)
  D <- W_new - W_ref

  comparison_2017 <- data.frame(
    Metric = c(
      "MaxAbsDifference",
      "MeanAbsDifference",
      "RMSE",
      "CorrelationOffDiagonal",
      "RowSumMaxAbsDifference"
    ),
    Value = c(
      max(abs(D)),
      mean(abs(D)),
      sqrt(mean(D^2)),
      suppressWarnings(cor(W_new[row(W_new) != col(W_new)],
                           W_ref[row(W_ref) != col(W_ref)])),
      max(abs(rowSums(W_new) - rowSums(W_ref)))
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    comparison_2017,
    file.path(out_dir, "05_2017_rebuild_comparison.csv"),
    row.names = FALSE
  )
}

# ---------- construct multi-year average only if defensible -------------------

multiyear_written <- FALSE
multiyear_path <- "data/GCAP_financial_W_multiyear_average.csv"
average_method <- NA_character_

if (length(strict_common_years) >= 2L) {
  # Preferred: average annual normalized weights, then renormalize to remove
  # tiny numerical drift. This avoids mechanically overweighting later years
  # simply because the global stock of financial assets is larger.
  A <- Reduce("+", weights_by_year) / length(weights_by_year)
  diag(A) <- 0
  A <- A / rowSums(A)

  write_matrix_csv(A, multiyear_path)
  multiyear_written <- TRUE
  average_method <- "Mean of annual row-normalized W across strict common years"
} else if (file.exists(multiyear_path)) {
  # Safety: do not leave a stale average matrix from a prior run.
  file.remove(multiyear_path)
}

# ---------- README / decision -------------------------------------------------

strict_text <- if (length(strict_common_years)) {
  paste(strict_common_years, collapse = ", ")
} else {
  "NONE"
}

loose_text <- if (length(loose_common_years)) {
  paste(loose_common_years, collapse = ", ")
} else {
  "NONE"
}

run_text <- if (length(strict_run)) {
  sprintf("%d-%d (%d years)", min(strict_run), max(strict_run), length(strict_run))
} else {
  "NONE"
}

decision <- if (length(strict_common_years) >= 3L) {
  paste0(
    "PASS FOR MULTI-YEAR BASELINE: at least 3 strict common years exist. ",
    "Use data/GCAP_financial_W_multiyear_average.csv as the candidate baseline, ",
    "then rerun weight preflight and GVAR stability."
  )
} else if (length(strict_common_years) == 2L) {
  paste0(
    "LIMITED MULTI-YEAR OPTION: only 2 strict common years exist. ",
    "Treat the 2-year average as robustness, not automatically as baseline."
  )
} else if (length(strict_common_years) == 1L) {
  paste0(
    "KEEP 2017 FIXED FINANCIAL W: only one strict common year exists across all 14 economies. ",
    "Do not manufacture a 2007-2017 average from unbalanced investor coverage."
  )
} else {
  paste0(
    "NO STRICT COMMON YEAR FOUND. Check source version, methodology labels, and asset codes ",
    "before using any financial-weight matrix."
  )
}

readme <- c(
  "GCAP MULTI-YEAR COVERAGE AUDIT",
  "================================",
  "",
  sprintf("Input: %s", if (downloaded_temp) official_url else input_file),
  sprintf("Raw sample-investor year range: %s", paste(range(all_years), collapse = " - ")),
  "",
  "Strict construction rule (matches current 2017 design):",
  "  US: Fund Holdings; E + BC + BG + BSF",
  "  Other 13 investors: Issuance; EF + B",
  "  Position concept: Position_Residency",
  "  EA counterpart: fixed EA19 issuer aggregation",
  "  Diagonal: zero",
  "  Normalization: within 14-economy sample, row sum = 1",
  "",
  sprintf("Strict all-14 common years: %s", strict_text),
  sprintf("Loose all-14 common years (any methodology): %s", loose_text),
  sprintf("Longest strict consecutive run: %s", run_text),
  "",
  sprintf("Multi-year matrix written: %s", if (multiyear_written) "YES" else "NO"),
  sprintf("Multi-year path: %s", if (multiyear_written) multiyear_path else "N/A"),
  sprintf("Average method: %s", if (!is.na(average_method)) average_method else "N/A"),
  "",
  paste0("DECISION: ", decision),
  "",
  "Interpretation:",
  "  A strict common year requires all 14 investor economies to have the required",
  "  asset classes under the same methodology rule used for the current 2017 matrix.",
  "  This is stricter than merely having some GCAP observations in a year.",
  "",
  "Important source caveat:",
  "  This script performs no new interpolation/extrapolation. GCAP itself may harmonize",
  "  intermittent Position_Residency gaps as documented by the dataset authors.",
  ""
)

if (!is.null(comparison_2017)) {
  readme <- c(
    readme,
    "2017 rebuild cross-check:",
    sprintf("  Max absolute difference vs current data/GCAP_financial_W_2017.csv: %.12g",
            comparison_2017$Value[comparison_2017$Metric == "MaxAbsDifference"]),
    ""
  )
}

writeLines(readme, file.path(out_dir, "README_gcap_coverage.txt"))

msg("")
msg("=== GCAP COVERAGE AUDIT COMPLETE ===")
msg("Strict common years: %s", strict_text)
msg("Loose common years:  %s", loose_text)
msg("Decision: %s", decision)
msg("Outputs: %s", out_dir)

if (multiyear_written) {
  msg("Candidate multi-year W: %s", multiyear_path)
}
