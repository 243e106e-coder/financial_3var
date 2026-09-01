# Reproducible builder for the 14-economy financial-weight matrices.
# Input: Restated_Bilateral_External_Portfolios.xlsx (GCAP v1.8, Dec. 2022)
# No missing value is converted to zero. The script stops if an off-diagonal
# selected cell is missing. GCAP's own estimates/interpolation remain publisher
# adjustments and are disclosed in the output documentation.

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with install.packages('readxl').")
}

input_file <- "data/raw/Restated_Bilateral_External_Portfolios.xlsx"
output_dir <- "data/processed/financial_weights_14"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reference_year <- 2017L
country_codes <- c(AU="AUS", BR="BRA", CA="CAN", CH="CHE", CN="CHN",
                   EA="EMU", UK="GBR", JP="JPN", KR="KOR", NO="NOR",
                   SG="SGP", TR="TUR", US="USA", ZA="ZAF")
ea_2017 <- c("AUT","BEL","CYP","DEU","ESP","EST","FIN","FRA","GRC",
             "IRL","ITA","LTU","LUX","LVA","MLT","NLD","PRT","SVK","SVN")
fund_investors <- c("AUS","CAN","CHE","EMU","GBR","NOR","USA")

raw <- as.data.frame(readxl::read_excel(input_file, sheet = "Sheet1"))
required <- c("Methodology","Year","Investor","Asset_Class_Code","Issuer",
              "Position_Residency","Restatement_Ex_Domestic")
if (length(setdiff(required, names(raw))) > 0L) {
  stop("Missing input columns: ", paste(setdiff(required, names(raw)), collapse=", "))
}

sum_na_safe <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm=TRUE)
pair_rows <- list()

for (i in names(country_codes)) {
  investor_iso <- unname(country_codes[i])
  method <- if (investor_iso %in% fund_investors) "Fund Holdings" else "Issuance"
  assets <- if (investor_iso == "USA" && method == "Fund Holdings") {
    c("E","BC","BG","BSF")
  } else {
    c("EF","B")
  }
  z <- raw[raw$Year == reference_year & raw$Investor == investor_iso &
             raw$Methodology == method & raw$Asset_Class_Code %in% assets, ]
  if (nrow(z) == 0L) stop("No source rows for ", i)

  for (j in names(country_codes)) {
    issuer_iso <- unname(country_codes[j])
    take <- if (j == "EA") z$Issuer %in% ea_2017 else z$Issuer == issuer_iso
    q <- z[take, , drop=FALSE]
    is_equity <- q$Asset_Class_Code %in% c("E","EF")
    pair_rows[[length(pair_rows)+1L]] <- data.frame(
      Investor_Code=i, Issuer_Code=j, Year=reference_year,
      Methodology=method, Asset_Class_Codes=paste(assets, collapse="+"),
      Position_Residency_USD_mn=sum_na_safe(q$Position_Residency),
      Restated_Ex_Domestic_USD_mn=sum_na_safe(q$Restatement_Ex_Domestic),
      Restated_Equity_USD_mn=sum_na_safe(q$Restatement_Ex_Domestic[is_equity]),
      Restated_Bonds_USD_mn=sum_na_safe(q$Restatement_Ex_Domestic[!is_equity]),
      Source_Rows=nrow(q), Source_NA=sum(is.na(q$Restatement_Ex_Domestic)),
      stringsAsFactors=FALSE
    )
  }
}

lineage <- do.call(rbind, pair_rows)
if (any(is.na(lineage$Restated_Ex_Domestic_USD_mn[lineage$Investor_Code != lineage$Issuer_Code]))) {
  stop("Off-diagonal source values are missing; no imputation was performed.")
}

make_weight <- function(value_col) {
  codes <- names(country_codes)
  mat <- matrix(NA_real_, length(codes), length(codes), dimnames=list(codes,codes))
  for (r in seq_len(nrow(lineage))) {
    mat[lineage$Investor_Code[r], lineage$Issuer_Code[r]] <- lineage[[value_col]][r]
  }
  diag(mat) <- 0
  rs <- rowSums(mat)
  if (any(!is.finite(rs) | rs <= 0)) stop("Invalid row total in ", value_col)
  mat / rs
}

outputs <- list(
  `01_W_main_restated_2017.csv`=make_weight("Restated_Ex_Domestic_USD_mn"),
  `02_W_residency_input_2017.csv`=make_weight("Position_Residency_USD_mn"),
  `03_W_equity_restated_2017.csv`=make_weight("Restated_Equity_USD_mn"),
  `04_W_bonds_restated_2017.csv`=make_weight("Restated_Bonds_USD_mn")
)
for (f in names(outputs)) {
  write.csv(outputs[[f]], file.path(output_dir, f), row.names=TRUE, quote=FALSE)
}
write.csv(lineage, file.path(output_dir, "05_cell_lineage_core.csv"), row.names=FALSE)

cat("Generated four 14x14 matrices. No project-level imputation was used.\n")
cat("All row sums:", paste(round(rowSums(outputs[[1]]), 12), collapse=", "), "\n")
