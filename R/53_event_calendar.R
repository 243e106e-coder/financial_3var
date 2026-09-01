#!/usr/bin/env Rscript

# Event calendar and sample audit for the financial 3-variable TVP-GVAR.
#
# Raw panel sample expected from R/10_build_financial_3var_input.R:
#   2000Q2 - 2025Q3
#
# Effective TVP sample under p = 1:
#   one quarter later than the raw sample start (normally 2000Q3).
#
# This script does NOT estimate IRFs. It creates a single audited event
# registry that later formal TVP-GVAR / TVP-IRF code can consume.

source("R/00_config.R")

PANEL_PATH <- file.path(DERIVED_DIR, "panel_domestic_fin3.csv")
OUT <- file.path(RESULTS_DIR, "events")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL_PATH)) {
  stopf("Missing %s. Run R/10_build_financial_3var_input.R first.", PANEL_PATH)
}

panel <- read.csv(PANEL_PATH, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("Quarter", "Country") %in% names(panel))) {
  stopf("Panel must contain Quarter and Country columns.")
}

panel$Quarter <- toupper(trimws(as.character(panel$Quarter)))
qid <- quarter_id(panel$Quarter)
if (any(!is.finite(qid))) stopf("Unparseable quarter labels in panel.")

raw_qid <- sort(unique(qid))
raw_quarters <- quarter_label(raw_qid)

if (length(raw_qid) < 2L) stopf("Need at least two quarters in the raw panel.")

raw_start_id <- min(raw_qid)
raw_end_id <- max(raw_qid)

# With p=1, the first quarter cannot enter the dynamic design because lag 1
# is unavailable. The effective sample therefore starts at raw_start + 1.
effective_start_id <- raw_start_id + 1L
effective_end_id <- raw_end_id

RAW_SAMPLE <- sprintf("%s - %s",
                      quarter_label(raw_start_id),
                      quarter_label(raw_end_id))
EFFECTIVE_SAMPLE <- sprintf("%s - %s",
                            quarter_label(effective_start_id),
                            quarter_label(effective_end_id))

# -------------------------------------------------------------------------
# Event registry
# -------------------------------------------------------------------------
# Core events are the main-paper comparison set.
# Appendix events expand the historical/geopolitical/financial coverage.
#
# "Timing" is used only to motivate the t0 / t0+1 robustness window.
# All events are evaluated at both t0 and t0+1 for consistency, while t-1
# is retained as a pre-event comparison quarter.

events <- data.frame(
  EventID = c(
    "LEHMAN_2008",
    "SOVEREIGN_STRESS_2011",
    "BREXIT_2016",
    "COVID_2020",
    "RUSSIA_UKRAINE_2022",
    "IRAN_ISRAEL_2024",
    "SEPT11_2001",
    "IRAQ_2003",
    "BANKING_STRESS_2023",
    "ISRAEL_HAMAS_2023"
  ),
  EventSet = c(
    rep("CORE", 6),
    rep("APPENDIX", 4)
  ),
  EventQuarter = c(
    "2008Q3",
    "2011Q3",
    "2016Q2",
    "2020Q1",
    "2022Q1",
    "2024Q2",
    "2001Q3",
    "2003Q1",
    "2023Q1",
    "2023Q4"
  ),
  EventLabel = c(
    "Lehman / global financial crisis escalation",
    "Euro-area sovereign stress and US sovereign-rating shock",
    "Brexit referendum",
    "COVID-19 global market stress",
    "Russia-Ukraine war outbreak",
    "Iran-Israel direct military escalation",
    "September 11 attacks",
    "Iraq war outbreak",
    "SVB / Credit Suisse banking stress",
    "Israel-Hamas war outbreak"
  ),
  ShockFamily = c(
    "GLOBAL_FINANCIAL",
    "SOVEREIGN_CREDIT",
    "POLITICAL_INSTITUTIONAL",
    "GLOBAL_LIQUIDITY",
    "GEOPOLITICAL_WAR_SANCTIONS",
    "GEOPOLITICAL_MILITARY",
    "GEOPOLITICAL_TERROR",
    "GEOPOLITICAL_WAR",
    "BANKING_RATE_RISK",
    "GEOPOLITICAL_WAR"
  ),
  Timing = c(
    "LATE_QUARTER",
    "MID_TO_LATE_QUARTER",
    "LATE_QUARTER",
    "LATE_QUARTER",
    "MID_QUARTER",
    "EARLY_QUARTER",
    "LATE_QUARTER",
    "LATE_QUARTER",
    "LATE_QUARTER",
    "EARLY_QUARTER"
  ),
  MainPaper = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    FALSE, FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

events$EventQID <- quarter_id(events$EventQuarter)
if (any(!is.finite(events$EventQID))) stopf("Event registry contains invalid quarter labels.")

events$PreQuarter <- quarter_label(events$EventQID - 1L)
events$AnchorQuarter_t0 <- quarter_label(events$EventQID)
events$AnchorQuarter_t1 <- quarter_label(events$EventQID + 1L)

events$InRawSample_t0 <- (
  events$EventQID >= raw_start_id &
  events$EventQID <= raw_end_id
)
events$InEffectiveSample_t0 <- (
  events$EventQID >= effective_start_id &
  events$EventQID <= effective_end_id
)
events$InEffectiveSample_t1 <- (
  events$EventQID + 1L >= effective_start_id &
  events$EventQID + 1L <= effective_end_id
)
events$PreQuarterAvailable <- (
  events$EventQID - 1L >= effective_start_id &
  events$EventQID - 1L <= effective_end_id
)

events$EventWindowReady <- (
  events$InEffectiveSample_t0 &
  events$InEffectiveSample_t1
)

if (any(!events$InRawSample_t0)) {
  bad <- paste(events$EventID[!events$InRawSample_t0], collapse = ", ")
  stopf("Event t0 falls outside raw panel sample: %s", bad)
}

if (any(!events$EventWindowReady)) {
  bad <- paste(events$EventID[!events$EventWindowReady], collapse = ", ")
  stopf("Event t0/t0+1 window falls outside effective TVP sample: %s", bad)
}

# Save audited master registry.
write.csv(
  events,
  file.path(OUT, "00_event_calendar.csv"),
  row.names = FALSE
)

# Long event-window representation.
window_rows <- lapply(seq_len(nrow(events)), function(i) {
  e <- events[i, , drop = FALSE]
  data.frame(
    EventID = e$EventID,
    EventSet = e$EventSet,
    EventLabel = e$EventLabel,
    ShockFamily = e$ShockFamily,
    EventQuarter = e$EventQuarter,
    WindowPosition = c("t_minus_1", "t0", "t_plus_1"),
    Quarter = c(e$PreQuarter, e$AnchorQuarter_t0, e$AnchorQuarter_t1),
    IsPrimaryAnchor = c(FALSE, TRUE, FALSE),
    IsSecondaryAnchor = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
})
event_windows <- do.call(rbind, window_rows)

event_windows$QuarterID <- quarter_id(event_windows$Quarter)
event_windows$InEffectiveSample <- (
  event_windows$QuarterID >= effective_start_id &
  event_windows$QuarterID <= effective_end_id
)

write.csv(
  event_windows[event_windows$EventSet == "CORE", , drop = FALSE],
  file.path(OUT, "01_core_event_windows.csv"),
  row.names = FALSE
)

write.csv(
  event_windows[event_windows$EventSet == "APPENDIX", , drop = FALSE],
  file.path(OUT, "02_appendix_event_windows.csv"),
  row.names = FALSE
)

sample_audit <- data.frame(
  Item = c(
    "RawSampleStart",
    "RawSampleEnd",
    "EffectiveTVPSampleStart_p1",
    "EffectiveTVPSampleEnd",
    "RawQuarterCount",
    "EffectiveQuarterCount_p1",
    "CoreEventCount",
    "AppendixEventCount",
    "AllEvent_t0_InRawSample",
    "AllEvent_t0_InEffectiveSample",
    "AllEvent_t1_InEffectiveSample",
    "AllCorePreQuartersAvailable"
  ),
  Value = c(
    quarter_label(raw_start_id),
    quarter_label(raw_end_id),
    quarter_label(effective_start_id),
    quarter_label(effective_end_id),
    as.character(length(raw_qid)),
    as.character(length(raw_qid) - 1L),
    as.character(sum(events$EventSet == "CORE")),
    as.character(sum(events$EventSet == "APPENDIX")),
    as.character(all(events$InRawSample_t0)),
    as.character(all(events$InEffectiveSample_t0)),
    as.character(all(events$InEffectiveSample_t1)),
    as.character(all(events$PreQuarterAvailable[events$EventSet == "CORE"]))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_audit,
  file.path(OUT, "03_event_sample_audit.csv"),
  row.names = FALSE
)

# Compact paper-ready event table.
paper_table <- events[, c(
  "EventSet",
  "EventID",
  "EventQuarter",
  "AnchorQuarter_t1",
  "EventLabel",
  "ShockFamily",
  "Timing",
  "EventWindowReady"
)]
names(paper_table)[names(paper_table) == "AnchorQuarter_t1"] <- "NextQuarter"

write.csv(
  paper_table,
  file.path(OUT, "04_paper_event_table.csv"),
  row.names = FALSE
)

readme <- c(
  "FINANCIAL 3-VARIABLE TVP-GVAR EVENT CALENDAR",
  "=============================================",
  sprintf("Raw panel sample: %s", RAW_SAMPLE),
  sprintf("Effective TVP sample (p=1): %s", EFFECTIVE_SAMPLE),
  "",
  "Core main-paper events:",
  paste0(
    "- ",
    events$EventQuarter[events$EventSet == "CORE"],
    ": ",
    events$EventLabel[events$EventSet == "CORE"]
  ),
  "",
  "Appendix events:",
  paste0(
    "- ",
    events$EventQuarter[events$EventSet == "APPENDIX"],
    ": ",
    events$EventLabel[events$EventSet == "APPENDIX"]
  ),
  "",
  "Window convention:",
  "- t-1 = pre-event comparison quarter",
  "- t0  = event quarter / primary TVP-IRF anchor",
  "- t+1 = next-quarter robustness anchor",
  "",
  "Important:",
  "- Event selection is kept separate from estimation code.",
  "- All core and appendix t0/t+1 anchors must be inside the effective TVP sample.",
  "- 2001Q3 is retained only as an appendix event because it lies near the sample start.",
  "- The formal TVP-IRF stage should consume these CSV files rather than hard-code quarters."
)

writeLines(readme, file.path(OUT, "README_event_calendar.txt"))

msg("Event calendar created.")
msg("Raw sample: %s", RAW_SAMPLE)
msg("Effective TVP sample (p=1): %s", EFFECTIVE_SAMPLE)
msg("Core events: %d", sum(events$EventSet == "CORE"))
msg("Appendix events: %d", sum(events$EventSet == "APPENDIX"))
msg("All t0/t+1 event windows ready: %s", all(events$EventWindowReady))
