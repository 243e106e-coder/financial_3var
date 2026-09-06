#!/usr/bin/env Rscript

# =============================================================================
# 71_pip_multiyear_financial_weight_coverage_audit.R
#
# 09a: STRICT IMF PIP/CPIS multi-year financial-weight coverage audit
#
# DATA SOURCE
# -----------
# IMF Portfolio Investment Positions by Counterpart Economy (PIP; formerly CPIS)
# Public SDMX endpoint:
#   providerId = IMF_DATA
#   flowRef    = IMF.STA,PIP
#
# STRICT PRINCIPLES
# -----------------
# 1) Use published ASSET positions only. No mirror / liability reconstruction.
# 2) Use the published TOTAL portfolio-investment position indicator:
#      P_TOTINV_P_USD
#    with total holder sector S1 and total counterpart sector S1.
# 3) Annual end-period observations only.
# 4) Missing / suppressed / ambiguous observations are NEVER converted to zero.
# 5) Published zero is accepted as zero.
# 6) Negative published positions are retained in the audit, but a negative
#    final directed sample-pair total is NOT allowed into a nonnegative GVAR W.
# 7) Euro Area is held fixed at EA19 for comparability with the current project:
#      AUT BEL CYP DEU ESP EST FIN FRA GRC IRL ITA LTU LUX LVA
#      MLT NLD PRT SVK SVN
#    We do NOT use a dynamic "Euro Area" aggregate that changes membership.
# 8) No imputation, interpolation, carry-forward, mirror filling, or smoothing.
# 9) A year is STRICT_MATRIX_READY only if all 14 x 13 = 182 directed
#    off-diagonal sample pairs are published/constructible and nonnegative.
#
# If strict years exist, per-year W matrices are produced. If >=3 strict years
# exist, the baseline candidate is formed by averaging RAW bilateral positions
# across strict years first and row-normalizing once afterward.
# =============================================================================

source("R/00_config.R")

if (!requireNamespace("rsdmx", quietly = TRUE)) {
  stopf("Package 'rsdmx' is required.")
}

get_env_int <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(as.integer(default))
  x <- suppressWarnings(as.integer(z))
  if (is.na(x)) stopf("Environment variable %s is not an integer.", name)
  x
}
get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}

START_YEAR <- get_env_int("FIN3_PIP_START_YEAR", 2015L)
END_YEAR   <- get_env_int("FIN3_PIP_END_YEAR", 2024L)
OUT <- get_env_chr(
  "FIN3_PIP_AUDIT_OUT",
  file.path(RESULTS_DIR, "pip_multiyear_coverage_audit")
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MATRIX_DIR <- file.path(OUT, "strict_year_matrices")
dir.create(MATRIX_DIR, recursive = TRUE, showWarnings = FALSE)

if (START_YEAR > END_YEAR) stopf("START_YEAR > END_YEAR.")
YEARS <- START_YEAR:END_YEAR

FLOW_REF <- "IMF.STA,PIP"
INDICATOR <- "P_TOTINV_P_USD"
ACCOUNTING_ENTRY <- "A"
HOLDER_SECTOR <- "S1"
COUNTERPART_SECTOR <- "S1"
FREQUENCY <- "A"

headers <- list(
  "User-Agent" = "financial_3var-research-audit/1.0",
  "Accept" = "application/vnd.sdmx.genericdata+xml;version=2.1"
)

# Project 14-economy mapping.
outside_map <- c(
  AU="AUS", BR="BRA", CA="CAN", CH="CHE", CN="CHN",
  UK="GBR", JP="JPN", KR="KOR", NO="NOR", SG="SGP",
  TR="TUR", US="USA", ZA="ZAF"
)
ea19 <- c(
  "AUT","BEL","CYP","DEU","ESP","EST","FIN","FRA","GRC",
  "IRL","ITA","LTU","LUX","LVA","MLT","NLD","PRT","SVK","SVN"
)
if (!identical(setdiff(COUNTRIES, "EA"), names(outside_map))) {
  stopf("Project country mapping is not aligned with R/00_config.R.")
}
all_underlying <- unique(c(unname(outside_map), ea19))

# -----------------------------------------------------------------------------
# 0. DSD / schema audit
# -----------------------------------------------------------------------------

msg("Fetching current IMF PIP DSD ...")
dsd <- rsdmx::readSDMX(
  providerId = "IMF_DATA",
  resource = "datastructure",
  resourceId = "DSD_PIP",
  headers = headers,
  references = "all",
  verbose = TRUE
)

ds <- methods::slot(dsd, "datastructures")@datastructures[[1]]
dims <- methods::slot(ds, "Components")@Dimensions
dim_names <- vapply(
  dims,
  function(x) methods::slot(x, "conceptRef"),
  FUN.VALUE = ""
)
write.csv(
  data.frame(Order=seq_along(dim_names), Dimension=dim_names),
  file.path(OUT, "01_pip_dsd_dimensions.csv"),
  row.names = FALSE
)

expected_key_dims <- c(
  "COUNTRY","ACCOUNTING_ENTRY","INDICATOR","SECTOR",
  "COUNTERPART_SECTOR","COUNTERPART_COUNTRY","FREQUENCY"
)
if (length(dim_names) < length(expected_key_dims) ||
    !identical(dim_names[seq_along(expected_key_dims)], expected_key_dims)) {
  stopf(
    "PIP key dimension order changed. Found leading dimensions: %s",
    paste(head(dim_names, 10), collapse = ",")
  )
}

# Extract labels for key codes when codelists are available. This is an audit
# aid only; the query itself uses the documented SDMX codes.
get_codelist_df <- function(dsd_obj, dimension_id) {
  cls <- methods::slot(dsd_obj, "codelists")@codelists
  hits <- Filter(
    function(cl) endsWith(methods::slot(cl, "id"), dimension_id),
    cls
  )
  if (!length(hits)) {
    return(data.frame(Code=character(), Label=character()))
  }
  cl <- hits[[1]]
  codes <- methods::slot(cl, "Code")
  data.frame(
    Code = vapply(codes, function(cd) methods::slot(cd, "id"), ""),
    Label = vapply(codes, function(cd) {
      lab <- methods::slot(cd, "label")
      vals <- unname(unlist(lab))
      if (!length(vals)) "" else as.character(vals[1])
    }, ""),
    stringsAsFactors = FALSE
  )
}

schema_codes <- list(
  COUNTRY = unique(c(all_underlying)),
  ACCOUNTING_ENTRY = ACCOUNTING_ENTRY,
  INDICATOR = INDICATOR,
  SECTOR = HOLDER_SECTOR,
  COUNTERPART_SECTOR = COUNTERPART_SECTOR,
  FREQUENCY = FREQUENCY
)
schema_rows <- list()
for (nm in names(schema_codes)) {
  cl <- get_codelist_df(dsd, nm)
  req <- schema_codes[[nm]]
  for (code in req) {
    hit <- cl[cl$Code == code, , drop=FALSE]
    schema_rows[[length(schema_rows)+1L]] <- data.frame(
      Dimension = nm,
      Code = code,
      FoundInCurrentDSD = nrow(hit) == 1L,
      Label = if (nrow(hit) == 1L) hit$Label[1] else NA_character_,
      stringsAsFactors = FALSE
    )
  }
}
schema_audit <- do.call(rbind, schema_rows)
write.csv(
  schema_audit,
  file.path(OUT, "02_requested_pip_codes_schema_audit.csv"),
  row.names = FALSE
)
if (any(!schema_audit$FoundInCurrentDSD)) {
  bad <- schema_audit[!schema_audit$FoundInCurrentDSD,]
  stopf(
    "At least one requested PIP code is absent from the current DSD: %s",
    paste(paste0(bad$Dimension, "=", bad$Code), collapse = "; ")
  )
}

# -----------------------------------------------------------------------------
# 1. Download one strict target extract from the official IMF endpoint
# -----------------------------------------------------------------------------

# One batched request:
# COUNTRY . ACCOUNTING_ENTRY . INDICATOR . SECTOR .
# COUNTERPART_SECTOR . COUNTERPART_COUNTRY . FREQUENCY
key <- list(
  all_underlying,
  ACCOUNTING_ENTRY,
  INDICATOR,
  HOLDER_SECTOR,
  COUNTERPART_SECTOR,
  all_underlying,
  FREQUENCY
)

msg(
  "Downloading official PIP target extract: %d reporters x %d counterparts, %d-%d ...",
  length(all_underlying), length(all_underlying), START_YEAR, END_YEAR
)

sdmx_data <- rsdmx::readSDMX(
  providerId = "IMF_DATA",
  resource = "data",
  flowRef = FLOW_REF,
  key = key,
  start = START_YEAR,
  end = END_YEAR,
  dsd = TRUE,
  headers = headers,
  verbose = TRUE
)
raw <- as.data.frame(sdmx_data, labels = FALSE)

required <- c(
  "COUNTRY","ACCOUNTING_ENTRY","INDICATOR","SECTOR",
  "COUNTERPART_SECTOR","COUNTERPART_COUNTRY","FREQUENCY",
  "TIME_PERIOD","OBS_VALUE"
)
if (!all(required %in% names(raw))) {
  stopf(
    "Downloaded PIP extract is missing columns: %s",
    paste(setdiff(required, names(raw)), collapse = ",")
  )
}

raw$COUNTRY <- toupper(trimws(as.character(raw$COUNTRY)))
raw$COUNTERPART_COUNTRY <- toupper(trimws(as.character(raw$COUNTERPART_COUNTRY)))
raw$TIME_PERIOD <- trimws(as.character(raw$TIME_PERIOD))
raw$OBS_VALUE_NUM <- suppressWarnings(as.numeric(as.character(raw$OBS_VALUE)))

# Verify that the API did not return a broader category than requested.
strict_filter_ok <- (
  raw$ACCOUNTING_ENTRY == ACCOUNTING_ENTRY &
  raw$INDICATOR == INDICATOR &
  raw$SECTOR == HOLDER_SECTOR &
  raw$COUNTERPART_SECTOR == COUNTERPART_SECTOR &
  raw$FREQUENCY == FREQUENCY
)
if (any(!strict_filter_ok, na.rm = TRUE)) {
  stopf("IMF API returned observations outside the requested strict PIP filter.")
}

raw <- raw[
  raw$COUNTRY %in% all_underlying &
  raw$COUNTERPART_COUNTRY %in% all_underlying &
  grepl("^[12][0-9]{3}$", raw$TIME_PERIOD),
  ,
  drop = FALSE
]
raw$Year <- as.integer(raw$TIME_PERIOD)
raw <- raw[raw$Year %in% YEARS, , drop=FALSE]

# Keep relevant publication-status metadata if present.
status_cols <- intersect(
  c(
    "OBS_STATUS","STATUS","CONF_STATUS","CONFIDENTIALITY_STATUS",
    "DV_TYPE","DERIVATION_TYPE","SCALE","UNIT"
  ),
  names(raw)
)

# P_TOTINV_P_USD is a single indicator. SCALE must not silently vary across
# observations. Relative W is invariant to a common scale, but mixed scale would
# invalidate raw-position averaging.
if ("SCALE" %in% names(raw)) {
  sc <- unique(raw$SCALE[!is.na(raw$SCALE) & trimws(as.character(raw$SCALE)) != ""])
  if (length(sc) > 1L) {
    stopf("PIP target extract contains multiple SCALE values: %s", paste(sc, collapse=","))
  }
}

write.csv(
  raw[, unique(c(required, "Year", "OBS_VALUE_NUM", status_cols)), drop=FALSE],
  file.path(OUT, "03_raw_official_pip_target_extract.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 2. Underlying reporter-counterpart-year publication status
# -----------------------------------------------------------------------------

status_text <- function(z) {
  if (!length(status_cols) || !nrow(z)) return("")
  vals <- unlist(lapply(status_cols, function(nm) {
    v <- unique(trimws(as.character(z[[nm]])))
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) return(character())
    paste0(nm, "=", paste(v, collapse="|"))
  }))
  paste(vals, collapse=";")
}

cell_rows <- vector(
  "list",
  length(all_underlying) * length(all_underlying) * length(YEARS)
)
kk <- 0L

for (rep in all_underlying) {
  for (cp in all_underlying) {
    for (yr in YEARS) {
      kk <- kk + 1L
      z <- raw[
        raw$COUNTRY == rep &
        raw$COUNTERPART_COUNTRY == cp &
        raw$Year == yr,
        ,
        drop = FALSE
      ]

      if (nrow(z) == 0L) {
        cls <- "NO_PUBLISHED_OBSERVATION"
        val <- NA_real_
      } else if (nrow(z) > 1L) {
        # Never silently select one duplicate.
        cls <- "AMBIGUOUS_MULTIPLE_OBSERVATIONS"
        val <- NA_real_
      } else if (!is.finite(z$OBS_VALUE_NUM[1])) {
        cls <- "PUBLISHED_NONNUMERIC_OR_SUPPRESSED"
        val <- NA_real_
      } else {
        val <- z$OBS_VALUE_NUM[1]
        if (val > 0) {
          cls <- "OBSERVED_POSITIVE"
        } else if (val == 0) {
          cls <- "OBSERVED_ZERO"
        } else {
          cls <- "OBSERVED_NEGATIVE"
        }
      }

      cell_rows[[kk]] <- data.frame(
        Reporter = rep,
        Counterpart = cp,
        Year = yr,
        PublicationClass = cls,
        PublishedValue = val,
        PublicationMetadata = status_text(z),
        DirectNumericPublished = cls %in% c(
          "OBSERVED_POSITIVE","OBSERVED_ZERO","OBSERVED_NEGATIVE"
        ),
        NonnegativePublished = cls %in% c(
          "OBSERVED_POSITIVE","OBSERVED_ZERO"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
underlying <- do.call(rbind, cell_rows)
write.csv(
  underlying,
  file.path(OUT, "04_underlying_cell_publication_status.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Build exact 14-economy directed-pair coverage with fixed EA19
# -----------------------------------------------------------------------------

lookup_cell <- function(rep, cp, yr) {
  z <- underlying[
    underlying$Reporter == rep &
    underlying$Counterpart == cp &
    underlying$Year == yr,
    ,
    drop = FALSE
  ]
  if (nrow(z) != 1L) stopf("Internal underlying-cell lookup error.")
  z
}

pair_rows <- list()
kk <- 0L

for (yr in YEARS) {
  for (i in COUNTRIES) {
    for (j in COUNTRIES) {
      if (i == j) next

      kk <- kk + 1L

      if (i != "EA" && j != "EA") {
        constituents <- lookup_cell(outside_map[[i]], outside_map[[j]], yr)
        pair_type <- "STANDALONE_TO_STANDALONE"
      } else if (i != "EA" && j == "EA") {
        constituents <- do.call(rbind, lapply(
          ea19,
          function(cp) lookup_cell(outside_map[[i]], cp, yr)
        ))
        pair_type <- "STANDALONE_TO_FIXED_EA19"
      } else if (i == "EA" && j != "EA") {
        constituents <- do.call(rbind, lapply(
          ea19,
          function(rep) lookup_cell(rep, outside_map[[j]], yr)
        ))
        pair_type <- "FIXED_EA19_TO_STANDALONE"
      } else {
        stopf("Unexpected EA/EA off-diagonal branch.")
      }

      all_numeric <- all(constituents$DirectNumericPublished)
      missing_n <- sum(!constituents$DirectNumericPublished)
      negative_n <- sum(
        constituents$PublicationClass == "OBSERVED_NEGATIVE",
        na.rm = TRUE
      )

      if (all_numeric) {
        pair_value <- sum(constituents$PublishedValue)
        pair_nonnegative <- is.finite(pair_value) && pair_value >= 0
        pair_class <- if (!pair_nonnegative) {
          "COMPLETE_BUT_NEGATIVE_PAIR_TOTAL"
        } else if (pair_value == 0) {
          "STRICT_COMPLETE_ZERO"
        } else {
          "STRICT_COMPLETE_POSITIVE"
        }
      } else {
        pair_value <- NA_real_
        pair_nonnegative <- FALSE
        pair_class <- "INCOMPLETE_UNDERLYING_PUBLICATION"
      }

      pair_rows[[kk]] <- data.frame(
        Year = yr,
        From = i,
        To = j,
        PairType = pair_type,
        UnderlyingCellsRequired = nrow(constituents),
        UnderlyingCellsDirectNumeric = sum(constituents$DirectNumericPublished),
        UnderlyingCellsMissingOrAmbiguous = missing_n,
        UnderlyingNegativePublishedCells = negative_n,
        PairPublishedComplete = all_numeric,
        PairValue = pair_value,
        PairNonnegative = pair_nonnegative,
        StrictWeightReady = all_numeric && pair_nonnegative,
        PairClass = pair_class,
        stringsAsFactors = FALSE
      )
    }
  }
}
pairs <- do.call(rbind, pair_rows)
if (nrow(pairs) != length(YEARS) * 14L * 13L) {
  stopf("Expected %d sample directed-pair rows; found %d.",
        length(YEARS)*182L, nrow(pairs))
}
write.csv(
  pairs,
  file.path(OUT, "05_sample_pair_coverage_strict.csv"),
  row.names = FALSE
)

# Detailed missing constituents for only the sample-relevant pairs.
missing_detail <- list()
mm <- 0L
for (yr in YEARS) {
  for (i in COUNTRIES) {
    for (j in COUNTRIES) {
      if (i == j) next

      if (i != "EA" && j != "EA") {
        z <- lookup_cell(outside_map[[i]], outside_map[[j]], yr)
      } else if (i != "EA" && j == "EA") {
        z <- do.call(rbind, lapply(
          ea19, function(cp) lookup_cell(outside_map[[i]], cp, yr)
        ))
      } else if (i == "EA" && j != "EA") {
        z <- do.call(rbind, lapply(
          ea19, function(rep) lookup_cell(rep, outside_map[[j]], yr)
        ))
      } else next

      bad <- z[
        !z$DirectNumericPublished |
        z$PublicationClass == "OBSERVED_NEGATIVE",
        ,
        drop=FALSE
      ]
      if (nrow(bad)) {
        mm <- mm + 1L
        bad$SampleFrom <- i
        bad$SampleTo <- j
        missing_detail[[mm]] <- bad
      }
    }
  }
}
missing_detail_df <- if (length(missing_detail)) {
  do.call(rbind, missing_detail)
} else {
  data.frame(
    Reporter=character(), Counterpart=character(), Year=integer(),
    PublicationClass=character(), PublishedValue=double(),
    PublicationMetadata=character(), DirectNumericPublished=logical(),
    NonnegativePublished=logical(), SampleFrom=character(),
    SampleTo=character(), stringsAsFactors=FALSE
  )
}
write.csv(
  missing_detail_df,
  file.path(OUT, "06_missing_ambiguous_negative_underlying_cells.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Year-level and country-row coverage gates
# -----------------------------------------------------------------------------

year_rows <- lapply(YEARS, function(yr) {
  z <- pairs[pairs$Year == yr, , drop=FALSE]
  stand <- z[z$PairType == "STANDALONE_TO_STANDALONE",]
  ea <- z[z$PairType != "STANDALONE_TO_STANDALONE",]

  ready <- all(z$StrictWeightReady)
  F <- matrix(
    NA_real_, length(COUNTRIES), length(COUNTRIES),
    dimnames=list(COUNTRIES, COUNTRIES)
  )
  diag(F) <- 0
  for (rr in seq_len(nrow(z))) {
    F[z$From[rr], z$To[rr]] <- z$PairValue[rr]
  }
  row_positive <- ready && all(rowSums(F) > 0)

  data.frame(
    Year = yr,
    DirectedPairsExpected = 182L,
    StrictWeightReadyPairs = sum(z$StrictWeightReady),
    IncompletePairs = sum(!z$PairPublishedComplete),
    NegativePairTotals = sum(z$PairPublishedComplete & !z$PairNonnegative),
    StandalonePairsExpected = nrow(stand),
    StandalonePairsReady = sum(stand$StrictWeightReady),
    EACompositePairsExpected = nrow(ea),
    EACompositePairsReady = sum(ea$StrictWeightReady),
    All182PairsReady = ready,
    AllRowPositionSumsPositive = row_positive,
    MatrixReadyStrict = ready && row_positive,
    stringsAsFactors = FALSE
  )
})
year_summary <- do.call(rbind, year_rows)
write.csv(
  year_summary,
  file.path(OUT, "07_year_coverage_summary.csv"),
  row.names = FALSE
)

row_rows <- list()
kk <- 0L
for (yr in YEARS) {
  for (i in COUNTRIES) {
    kk <- kk + 1L
    z <- pairs[pairs$Year == yr & pairs$From == i, , drop=FALSE]
    row_rows[[kk]] <- data.frame(
      Year=yr,
      Country=i,
      CounterpartsExpected=13L,
      CounterpartsReady=sum(z$StrictWeightReady),
      RowComplete=all(z$StrictWeightReady),
      RowPublishedPositionSum=if (all(z$StrictWeightReady)) sum(z$PairValue) else NA_real_,
      stringsAsFactors=FALSE
    )
  }
}
row_summary <- do.call(rbind, row_rows)
write.csv(
  row_summary,
  file.path(OUT, "08_country_year_row_coverage.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Build matrices ONLY for strictly ready years
# -----------------------------------------------------------------------------

strict_years <- year_summary$Year[year_summary$MatrixReadyStrict]
F_list <- list()

write_matrix <- function(A, path) {
  write.csv(
    data.frame(Country=rownames(A), A, check.names=FALSE),
    path, row.names=FALSE, quote=FALSE
  )
}

for (yr in strict_years) {
  z <- pairs[pairs$Year == yr, , drop=FALSE]
  F <- matrix(
    0, length(COUNTRIES), length(COUNTRIES),
    dimnames=list(COUNTRIES, COUNTRIES)
  )
  for (rr in seq_len(nrow(z))) F[z$From[rr], z$To[rr]] <- z$PairValue[rr]
  diag(F) <- 0

  if (any(!is.finite(F)) || any(F < 0) || any(rowSums(F) <= 0)) {
    stopf("Internal strict-year matrix integrity failure for %d.", yr)
  }

  W <- F / rowSums(F)
  F_list[[as.character(yr)]] <- F
  write_matrix(F, file.path(MATRIX_DIR, sprintf("F_pip_total_position_%d.csv", yr)))
  write_matrix(W, file.path(MATRIX_DIR, sprintf("W_pip_total_position_%d.csv", yr)))
}

multiyear_file <- NA_character_
if (length(strict_years) >= 3L) {
  Favg <- Reduce("+", F_list) / length(F_list)
  diag(Favg) <- 0
  if (any(Favg < 0) || any(rowSums(Favg) <= 0)) {
    stopf("Multi-year average raw position matrix is invalid.")
  }
  Wavg <- Favg / rowSums(Favg)
  span <- paste0(min(strict_years), "_", max(strict_years))
  multiyear_file <- sprintf("W_pip_multiyear_position_average_%s.csv", span)
  write_matrix(
    Favg,
    file.path(OUT, sprintf("F_pip_multiyear_position_average_%s.csv", span))
  )
  write_matrix(Wavg, file.path(OUT, multiyear_file))
}

status <- if (length(strict_years) >= 3L) {
  "AUDIT_COMPLETE_MULTIYEAR_BASELINE_AVAILABLE"
} else if (length(strict_years) >= 1L) {
  "AUDIT_COMPLETE_STRICT_YEARS_FOUND_BUT_LT3"
} else {
  "AUDIT_COMPLETE_NO_STRICT_COMPLETE_YEAR"
}

gate <- data.frame(
  Status = status,
  Source = "IMF PIP/CPIS official SDMX",
  FlowRef = FLOW_REF,
  Indicator = INDICATOR,
  AccountingEntry = ACCOUNTING_ENTRY,
  HolderSector = HOLDER_SECTOR,
  CounterpartSector = COUNTERPART_SECTOR,
  Frequency = FREQUENCY,
  StartYear = START_YEAR,
  EndYear = END_YEAR,
  CandidateYearsAudited = length(YEARS),
  StrictMatrixReadyYears = length(strict_years),
  StrictYears = if (length(strict_years)) paste(strict_years, collapse=";") else "",
  FixedEuroAreaDefinition = "EA19",
  ExpectedDirectedPairsPerYear = 182L,
  MissingTreatedAsZero = FALSE,
  SuppressedTreatedAsZero = FALSE,
  NegativePairClippedToZero = FALSE,
  MirrorImputationUsed = FALSE,
  InterpolationUsed = FALSE,
  CarryForwardUsed = FALSE,
  MultiyearBaselineCandidateFile =
    if (is.na(multiyear_file)) "" else multiyear_file,
  stringsAsFactors = FALSE
)
write.csv(gate, file.path(OUT, "00_pip_coverage_gate.csv"), row.names = FALSE)

readme <- c(
  sprintf("IMF PIP STRICT MULTI-YEAR COVERAGE AUDIT: %s", status),
  "============================================================",
  "",
  sprintf("Official flow: %s", FLOW_REF),
  sprintf("Indicator: %s", INDICATOR),
  sprintf("Filter: Assets / S1 holder / S1 counterpart / annual"),
  sprintf("Candidate years: %d-%d", START_YEAR, END_YEAR),
  sprintf("Project sample: %s", paste(COUNTRIES, collapse=", ")),
  "Euro Area: fixed EA19, matching the project's existing GCAP convention.",
  "",
  "STRICT RULE:",
  "- every off-diagonal sample pair must be directly published/constructible;",
  "- missing/suppressed/ambiguous observations are NOT zero;",
  "- no mirror filling, interpolation, carry-forward, or imputation;",
  "- negative final sample-pair totals are not admitted to a nonnegative GVAR W;",
  "- 14 x 13 = 182 directed pairs must pass in the same year.",
  "",
  sprintf(
    "Strict matrix-ready years (%d): %s",
    length(strict_years),
    if (length(strict_years)) paste(strict_years, collapse=", ") else "NONE"
  ),
  "",
  if (length(strict_years) >= 3L) {
    paste0(
      "MULTI-YEAR BASELINE CANDIDATE CREATED: ",
      multiyear_file,
      " (raw positions averaged across strict years, then row-normalized)."
    )
  } else {
    "No multi-year baseline is created unless at least 3 strict years pass."
  },
  "",
  "Inspect first:",
  "- 07_year_coverage_summary.csv",
  "- 08_country_year_row_coverage.csv",
  "- 06_missing_ambiguous_negative_underlying_cells.csv",
  "- 05_sample_pair_coverage_strict.csv",
  "",
  "No missing value has been fabricated."
)
writeLines(readme, file.path(OUT, "README_pip_coverage_audit.txt"))
cat(paste(readme, collapse="\n"), "\n")
