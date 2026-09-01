#!/usr/bin/env Rscript

source("R/00_config.R")

YEAR <- 2017L
OUT <- file.path(RESULTS_DIR,"weight_rebuild")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

official_url <- paste0(
  "https://globalcapitalallocation.s3.us-east-2.amazonaws.com/",
  "Restated_Bilateral_External_Portfolios.dta"
)

input_file <- trimws(Sys.getenv("GCAP_FILE",""))
if (!nzchar(input_file)) {
  input_file <- tempfile(fileext=".dta")
  msg("Downloading published GCAP source: %s",official_url)
  ok <- tryCatch({
    utils::download.file(official_url,input_file,mode="wb",quiet=FALSE)
    TRUE
  },error=function(e) {msg("Download failed: %s",conditionMessage(e)); FALSE})
  if (!ok || !file.exists(input_file) || file.info(input_file)$size<=0) {
    stopf("Could not download the published GCAP source. Set GCAP_FILE to a local copy.")
  }
}
if (!file.exists(input_file)) stopf("GCAP input not found: %s",input_file)

ext <- tolower(tools::file_ext(input_file))
if (ext=="dta") {
  if (!requireNamespace("haven",quietly=TRUE)) stopf("Package 'haven' is required.")
  raw <- as.data.frame(haven::read_dta(input_file))
} else if (ext %in% c("xlsx","xls")) {
  if (!requireNamespace("readxl",quietly=TRUE)) stopf("Package 'readxl' is required.")
  raw <- as.data.frame(readxl::read_excel(input_file,sheet=1))
} else if (ext=="csv") {
  raw <- read.csv(input_file,check.names=FALSE,stringsAsFactors=FALSE)
} else stopf("Unsupported GCAP input format: %s",ext)
if (!nrow(raw)) stopf("GCAP input contains no rows.")

norm_name <- function(x) gsub("[^a-z0-9]","",tolower(x))
pick <- function(candidates,label) {
  z <- match(norm_name(candidates),norm_name(names(raw)))
  z <- z[!is.na(z)]
  if (!length(z)) stopf("Cannot identify GCAP column: %s",label)
  names(raw)[z[1]]
}

col_method <- pick(c("Methodology"),"Methodology")
col_year <- pick(c("Year"),"Year")
col_investor <- pick(c("Investor"),"Investor")
col_asset <- pick(c("Asset_Class_Code","AssetClassCode"),"Asset Class Code")
col_issuer <- pick(c("Issuer"),"Issuer")
col_residency <- pick(c("Position_Residency","PositionResidence"),"Position Residency")
col_restated <- pick(c("Restatement_Ex_Domestic","RestatementExDomestic"),"Restatement Ex Domestic")

d <- data.frame(
  methodology=trimws(as.character(raw[[col_method]])),
  year=as.integer(raw[[col_year]]),
  investor=toupper(trimws(as.character(raw[[col_investor]]))),
  asset=toupper(trimws(as.character(raw[[col_asset]]))),
  issuer=toupper(trimws(as.character(raw[[col_issuer]]))),
  residency=num(raw[[col_residency]]),
  restated=num(raw[[col_restated]]),
  stringsAsFactors=FALSE
)

investor_map <- c(AU="AUS",BR="BRA",CA="CAN",CH="CHE",CN="CHN",EA="EMU",
                  UK="GBR",JP="JPN",KR="KOR",NO="NOR",SG="SGP",TR="TUR",
                  US="USA",ZA="ZAF")
ea19 <- c("AUT","BEL","CYP","DEU","ESP","EST","FIN","FRA","GRC",
          "IRL","ITA","LTU","LUX","LVA","MLT","NLD","PRT","SVK","SVN")
issuer_map <- lapply(investor_map,function(x)x)
issuer_map$EA <- ea19
fund_investors <- c("AUS","CAN","CHE","EMU","GBR","NOR","USA")

methods <- setNames(ifelse(investor_map %in% fund_investors,"Fund Holdings","Issuance"),
                    names(investor_map))
assets <- lapply(names(investor_map),function(cc) {
  if (cc=="US") c("E","BC","BG","BSF") else c("EF","B")
})
names(assets) <- names(investor_map)

positions <- list(
  main="restated", residency="residency", equity="restated", bonds="restated"
)
F <- lapply(positions,function(x) matrix(NA_real_,length(COUNTRIES),length(COUNTRIES),
                                        dimnames=list(COUNTRIES,COUNTRIES)))
coverage <- list()

for (i in COUNTRIES) {
  iso <- investor_map[[i]]; method <- methods[[i]]; req <- assets[[i]]
  zi <- d[d$year==YEAR & d$investor==iso & d$methodology==method & d$asset %in% req,,drop=FALSE]
  if (!nrow(zi)) stopf("No GCAP rows for %s/%s/%d",i,method,YEAR)
  if (!setequal(unique(zi$asset),req)) stopf("Asset-class coverage is incomplete for %s",i)
  coverage[[i]] <- data.frame(
    Country=i,InvestorISO=iso,Year=YEAR,Methodology=method,
    AssetClasses=paste(req,collapse="+"),SourceRows=nrow(zi),
    SourceIssuers=length(unique(zi$issuer)),RestatedNA=sum(is.na(zi$restated)),
    ResidencyNA=sum(is.na(zi$residency)),OurImputation="none",
    stringsAsFactors=FALSE
  )

  for (j in COUNTRIES) {
    if (i==j) {
      for (nm in names(F)) F[[nm]][i,j] <- 0
      next
    }
    zij <- zi[zi$issuer %in% issuer_map[[j]],,drop=FALSE]
    if (!nrow(zij)) stopf("No source rows for directed pair %s -> %s",i,j)
    if (any(!is.finite(zij$restated)) || any(!is.finite(zij$residency))) {
      stopf("Non-finite published source value for directed pair %s -> %s",i,j)
    }
    equity <- zij$asset %in% c("E","EF")
    F$main[i,j] <- sum(zij$restated)
    F$residency[i,j] <- sum(zij$residency)
    F$equity[i,j] <- sum(zij$restated[equity])
    F$bonds[i,j] <- sum(zij$restated[!equity])
  }
}

normalize <- function(A,label) {
  if (any(!is.finite(A))) stopf("Non-finite cell in rebuilt %s matrix",label)
  if (any(A < -1e-12)) stopf("Negative cell in rebuilt %s matrix",label)
  diag(A) <- 0
  if (any(rowSums(A)<=0)) stopf("Non-positive row sum in rebuilt %s matrix",label)
  A/rowSums(A)
}
W <- Map(normalize,F,names(F))
names(W) <- names(F)

out_names <- c(
  main="W_main_restated_exdom_2017.csv",
  residency="W_residency_2017.csv",
  equity="W_equity_restated_2017.csv",
  bonds="W_bonds_restated_2017.csv"
)
committed <- c(
  main=WEIGHT_FILES[["main_restated_exdom_2017"]],
  residency=WEIGHT_FILES[["residency_2017"]],
  equity=WEIGHT_FILES[["equity_restated_2017"]],
  bonds=WEIGHT_FILES[["bonds_restated_2017"]]
)

comparison <- list()
for (nm in names(W)) {
  write.csv(data.frame(Country=COUNTRIES,W[[nm]],check.names=FALSE),
            file.path(OUT,out_names[[nm]]),row.names=FALSE,quote=FALSE)
  ref <- read_weight_matrix(committed[[nm]])
  delta <- W[[nm]]-ref
  comparison[[nm]] <- data.frame(
    Matrix=nm,MaxAbsDifference=max(abs(delta)),MeanAbsDifference=mean(abs(delta)),
    MaxAbsDiagonal=max(abs(diag(W[[nm]]))),
    MaxAbsRowSumMinus1=max(abs(rowSums(W[[nm]])-1)),
    Status=if(max(abs(delta))<=1e-10)"MATCH" else "MISMATCH",
    stringsAsFactors=FALSE
  )
}
comparison_df <- do.call(rbind,comparison)
write.csv(do.call(rbind,coverage),file.path(OUT,"01_country_methodology_coverage.csv"),row.names=FALSE)
write.csv(comparison_df,file.path(OUT,"02_committed_matrix_comparison.csv"),row.names=FALSE)

if (any(comparison_df$Status!="MATCH")) stopf("At least one committed matrix does not match the published-source reconstruction.")

txt <- c(
  "GCAP PUBLISHED-SOURCE REBUILD: PASS",
  "===================================",
  sprintf("Reference year: %d",YEAR),
  "Main field: Restatement_Ex_Domestic",
  "Methodology: Fund Holdings for AU/CA/CH/EA/UK/NO/US; Issuance for BR/CN/JP/KR/SG/TR/ZA",
  "Project-level imputation: none",
  sprintf("Largest committed-vs-rebuilt difference: %.12g",max(comparison_df$MaxAbsDifference)),
  "All four committed matrices match the published-source reconstruction."
)
writeLines(txt,file.path(OUT,"README_weight_rebuild.txt"))
cat(paste(txt,collapse="\n"),"\n")
