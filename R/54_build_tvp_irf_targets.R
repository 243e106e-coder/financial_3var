#!/usr/bin/env Rscript

# Build TVP-IRF target grids from the audited event calendar.
#
# This script does not estimate the TVP-GVAR. It defines exactly which
# event-quarter anchors, economies, variables, and horizons the formal
# Bayesian TVP-IRF code must produce.

source("R/00_config.R")

EVENT_DIR <- file.path(RESULTS_DIR, "events")
EVENT_FILE <- file.path(EVENT_DIR, "00_event_calendar.csv")
OUT <- EVENT_DIR

if (!file.exists(EVENT_FILE)) {
  stopf("Missing %s. Run R/53_event_calendar.R first.", EVENT_FILE)
}

events <- read.csv(EVENT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

need <- c(
  "EventID", "EventSet", "EventLabel", "ShockFamily",
  "EventQuarter", "AnchorQuarter_t0", "AnchorQuarter_t1",
  "EventWindowReady"
)
if (!all(need %in% names(events))) {
  stopf("Event calendar is missing required columns.")
}

if (any(!events$EventWindowReady)) {
  bad <- paste(events$EventID[!events$EventWindowReady], collapse = ", ")
  stopf("Cannot build IRF targets; event windows not ready: %s", bad)
}

# Formal IRF horizons. 12 quarters captures immediate and medium-run
# financial-market transmission while remaining compact for a quarterly model.
IRF_HORIZONS <- 0:12

# Every event is evaluated at t0 and t+1.
anchor_rows <- do.call(rbind, lapply(seq_len(nrow(events)), function(i) {
  e <- events[i, , drop = FALSE]
  data.frame(
    EventID = e$EventID,
    EventSet = e$EventSet,
    EventLabel = e$EventLabel,
    ShockFamily = e$ShockFamily,
    EventQuarter = e$EventQuarter,
    AnchorType = c("EVENT_QUARTER_t0", "NEXT_QUARTER_t1"),
    AnchorQuarter = c(e$AnchorQuarter_t0, e$AnchorQuarter_t1),
    stringsAsFactors = FALSE
  )
}))

# Full formal target grid:
# event x t0/t1 anchor x economy x response variable x horizon.
targets <- merge(
  anchor_rows,
  data.frame(Country = COUNTRIES, stringsAsFactors = FALSE),
  all = TRUE
)
targets <- merge(
  targets,
  data.frame(ResponseVariable = VARS, stringsAsFactors = FALSE),
  all = TRUE
)
targets <- merge(
  targets,
  data.frame(Horizon = IRF_HORIZONS, stringsAsFactors = FALSE),
  all = TRUE
)

# Preserve meaningful order.
event_order <- events$EventID
targets$EventOrder <- match(targets$EventID, event_order)
targets$CountryOrder <- match(targets$Country, COUNTRIES)
targets$VariableOrder <- match(targets$ResponseVariable, VARS)
targets$AnchorOrder <- match(
  targets$AnchorType,
  c("EVENT_QUARTER_t0", "NEXT_QUARTER_t1")
)

targets <- targets[
  order(
    targets$EventOrder,
    targets$AnchorOrder,
    targets$CountryOrder,
    targets$VariableOrder,
    targets$Horizon
  ),
  ,
  drop = FALSE
]

targets$TargetID <- sprintf(
  "%s__%s__%s__%s__h%02d",
  targets$EventID,
  ifelse(targets$AnchorType == "EVENT_QUARTER_t0", "t0", "t1"),
  targets$Country,
  targets$ResponseVariable,
  targets$Horizon
)

keep <- c(
  "TargetID",
  "EventSet",
  "EventID",
  "EventLabel",
  "ShockFamily",
  "EventQuarter",
  "AnchorType",
  "AnchorQuarter",
  "Country",
  "ResponseVariable",
  "Horizon"
)
targets <- targets[, keep, drop = FALSE]

write.csv(
  targets,
  file.path(OUT, "05_formal_tvp_irf_targets.csv"),
  row.names = FALSE
)

write.csv(
  targets[targets$EventSet == "CORE", , drop = FALSE],
  file.path(OUT, "06_core_tvp_irf_targets.csv"),
  row.names = FALSE
)

write.csv(
  targets[targets$EventSet == "APPENDIX", , drop = FALSE],
  file.path(OUT, "07_appendix_tvp_irf_targets.csv"),
  row.names = FALSE
)

# Event-level t0 vs t1 comparison map. This is useful later when the formal
# IRF script computes shock-dependent safe-haven behavior.
pair_map <- do.call(rbind, lapply(seq_len(nrow(events)), function(i) {
  e <- events[i, , drop = FALSE]
  data.frame(
    EventSet = e$EventSet,
    EventID = e$EventID,
    EventLabel = e$EventLabel,
    ShockFamily = e$ShockFamily,
    EventQuarter = e$AnchorQuarter_t0,
    NextQuarter = e$AnchorQuarter_t1,
    Comparison = sprintf(
      "IRF(%s) versus IRF(%s)",
      e$AnchorQuarter_t0,
      e$AnchorQuarter_t1
    ),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  pair_map,
  file.path(OUT, "08_event_t0_t1_comparison_map.csv"),
  row.names = FALSE
)

# Safe-haven output contract for the later formal model.
# This file specifies what will be classified; it does NOT impose a sign rule
# prematurely because r, REER and equity responses require different economic
# interpretations and uncertainty criteria.
safe_haven_contract <- expand.grid(
  EventID = events$EventID,
  AnchorType = c("EVENT_QUARTER_t0", "NEXT_QUARTER_t1"),
  Country = COUNTRIES,
  ResponseVariable = VARS,
  stringsAsFactors = FALSE
)
safe_haven_contract$RequiredPosteriorOutputs <- paste(
  "median",
  "p16",
  "p84",
  "p05",
  "p95",
  "sign_probability",
  sep = ";"
)
safe_haven_contract$ClassificationStatus <- "TO_BE_CLASSIFIED_AFTER_FORMAL_IRF"

write.csv(
  safe_haven_contract,
  file.path(OUT, "09_safe_haven_classification_contract.csv"),
  row.names = FALSE
)

summary <- data.frame(
  Metric = c(
    "EventsTotal",
    "CoreEvents",
    "AppendixEvents",
    "AnchorsPerEvent",
    "Countries",
    "Variables",
    "IRFHorizons",
    "CoreIRFTargetRows",
    "AppendixIRFTargetRows",
    "AllIRFTargetRows"
  ),
  Value = c(
    nrow(events),
    sum(events$EventSet == "CORE"),
    sum(events$EventSet == "APPENDIX"),
    2L,
    length(COUNTRIES),
    length(VARS),
    length(IRF_HORIZONS),
    sum(targets$EventSet == "CORE"),
    sum(targets$EventSet == "APPENDIX"),
    nrow(targets)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary,
  file.path(OUT, "10_irf_target_summary.csv"),
  row.names = FALSE
)

msg("Formal TVP-IRF target grid created.")
msg("Events: %d (%d core, %d appendix)",
    nrow(events),
    sum(events$EventSet == "CORE"),
    sum(events$EventSet == "APPENDIX"))
msg("Countries: %d; variables: %d; horizons: %d",
    length(COUNTRIES), length(VARS), length(IRF_HORIZONS))
msg("Total target rows: %d", nrow(targets))
