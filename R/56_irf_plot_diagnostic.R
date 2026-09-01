#!/usr/bin/env Rscript

# =============================================================================
# TRANSFORMATION-AWARE DIAGNOSTIC IRF FIGURES
# Financial 3-variable TVP-GVAR
#
# Reads only audited numerical IRFs from R/55_irf_significance_audit.R.
# Significance is never inferred from the image.
# =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

AUDIT_DIR <- file.path(RESULTS_DIR, "irf_audit")
DETAIL_FILE <- file.path(AUDIT_DIR, "00_irf_posterior_summary_audited.csv")
SCALE_FILE <- file.path(AUDIT_DIR, "00a_irf_variable_scale_manifest.csv")
PATH_FILE <- file.path(AUDIT_DIR, "01_irf_path_significance_audit.csv")
FLAG_FILE <- file.path(AUDIT_DIR, "02_irf_flagged_cases.csv")

PLOT_DIR <- file.path(AUDIT_DIR, "plots")
CORE_DIR <- file.path(PLOT_DIR, "core")
APP_DIR <- file.path(PLOT_DIR, "appendix")
FLAG_DIR <- file.path(PLOT_DIR, "flagged")

dir.create(CORE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(APP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FLAG_DIR, recursive = TRUE, showWarnings = FALSE)

MAX_FLAGGED <- as.integer(Sys.getenv("FIN3_MAX_FLAGGED_PLOTS", "100"))
if (!is.finite(MAX_FLAGGED) || MAX_FLAGGED < 0) MAX_FLAGGED <- 100L

RATE_MODE <- tolower(trimws(Sys.getenv("FIN3_RATE_MODE", "level")))
if (!RATE_MODE %in% c("level", "difference")) {
  stopf("FIN3_RATE_MODE must be 'level' or 'difference'.")
}

EQ_RETURN_MODE <- tolower(trimws(Sys.getenv("FIN3_EQ_RETURN_MODE", "simple_return")))
if (!EQ_RETURN_MODE %in% c("simple_return", "log_return")) {
  stopf("FIN3_EQ_RETURN_MODE must be 'simple_return' or 'log_return'.")
}

for (f in c(DETAIL_FILE, SCALE_FILE, PATH_FILE, FLAG_FILE)) {
  if (!file.exists(f)) {
    stopf("Missing %s. Run R/55_irf_significance_audit.R first.", f)
  }
}

d <- read.csv(DETAIL_FILE, stringsAsFactors = FALSE, check.names = FALSE)
scales <- read.csv(SCALE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
paths <- read.csv(PATH_FILE, stringsAsFactors = FALSE, check.names = FALSE)
flagged <- read.csv(FLAG_FILE, stringsAsFactors = FALSE, check.names = FALSE)

rate_title <- if (RATE_MODE == "level") {
  "Interest-rate level"
} else {
  "Interest-rate change (Delta r)"
}

rate_axis <- if (RATE_MODE == "level") {
  "Interest-rate level response"
} else {
  "Delta-interest-rate response"
}

var_title <- c(
  r = rate_title,
  de = "REER log change",
  deq = "Equity return"
)

var_axis <- c(
  r = rate_axis,
  de = "REER log-change response",
  deq = "Equity-return response"
)

sanitize <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

event_label_for <- function(z) {
  lab <- unique(z$EventLabel[!is.na(z$EventLabel) & nzchar(z$EventLabel)])
  if (length(lab)) lab[1] else unique(z$EventID)[1]
}

anchor_label <- function(z) {
  aq <- unique(z$AnchorQuarter[!is.na(z$AnchorQuarter) & nzchar(z$AnchorQuarter)])
  at <- unique(z$AnchorType)
  if (length(aq)) sprintf("%s (%s)", at[1], aq[1]) else at[1]
}

subtitle_for_var <- function(z, v) {
  base <- "Median with 90% posterior credible band; crossing zero = not 90% credible"

  if (v == "r" && RATE_MODE == "level") {
    return(paste0(
      base,
      "; r is a level response, so the plotted path is already the interest-rate deviation"
    ))
  }

  if (v == "r" && RATE_MODE == "difference") {
    return(paste0(
      base,
      "; plotted r is Delta r; cumulate draw-level IRFs for a rate-level change"
    ))
  }

  if (v == "de") {
    meaning <- unique(z$PositiveDirectionMeaning)[1]
    return(paste0(
      base,
      "; de = REER_DLOG; positive de = ", meaning,
      "; cumulate draw-level IRFs for a REER log-level effect"
    ))
  }

  if (v == "deq" && EQ_RETURN_MODE == "log_return") {
    return(paste0(
      base,
      "; deq treated as log return by runtime setting; use draw-level cumulation for price-level effects"
    ))
  }

  paste0(
    base,
    "; deq = equity return; no log-price cumulation is assumed"
  )
}

# -------------------------------------------------------------------------
# Main diagnostic figure: event x anchor x variable, 14 country facets
# -------------------------------------------------------------------------

main_keys <- unique(d[, c("EventSet", "EventID", "AnchorType", "ResponseVariable")])

manifest <- list()
midx <- 0L

for (i in seq_len(nrow(main_keys))) {
  k <- main_keys[i, , drop = FALSE]

  z <- d[
    d$EventID == k$EventID &
      d$AnchorType == k$AnchorType &
      d$ResponseVariable == k$ResponseVariable,
    ,
    drop = FALSE
  ]

  z$Country <- factor(z$Country, levels = COUNTRIES)

  v <- k$ResponseVariable
  set_name <- ifelse(is.na(k$EventSet) || !nzchar(k$EventSet), "UNKNOWN", k$EventSet)
  target_dir <- if (set_name == "CORE") CORE_DIR else APP_DIR

  ttl <- sprintf(
    "%s | %s | %s response",
    event_label_for(z),
    anchor_label(z),
    unname(var_title[v])
  )

  p <- ggplot(z, aes(x = Horizon, y = median)) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.20) +
    geom_line(linewidth = 0.55) +
    facet_wrap(~Country, scales = "free_y", ncol = 4) +
    labs(
      title = ttl,
      subtitle = subtitle_for_var(z, v),
      x = "Horizon (quarters)",
      y = unname(var_axis[v])
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8.5),
      strip.text = element_text(face = "bold")
    )

  fn <- sprintf(
    "%s__%s__%s.png",
    sanitize(k$EventID),
    ifelse(k$AnchorType == "EVENT_QUARTER_t0", "t0", "t1"),
    sanitize(v)
  )

  path <- file.path(target_dir, fn)

  ggsave(
    filename = path,
    plot = p,
    width = 12,
    height = 9,
    dpi = 180
  )

  midx <- midx + 1L
  manifest[[midx]] <- data.frame(
    PlotType = "FACET_14_COUNTRIES",
    EventSet = set_name,
    EventID = k$EventID,
    AnchorType = k$AnchorType,
    Country = "ALL_14",
    ResponseVariable = v,
    RuntimeRateMode = RATE_MODE,
    RuntimeEquityReturnMode = EQ_RETURN_MODE,
    File = path,
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------------
# Individual diagnostic figures for flagged paths
# -------------------------------------------------------------------------

if (nrow(flagged) && MAX_FLAGGED > 0L) {
  flagged$Priority <- ifelse(
    flagged$Flag_CredibleSignReversal,
    1L,
    ifelse(flagged$Flag_NoCredibleResponse, 2L, 3L)
  )

  flagged <- flagged[
    order(flagged$Priority, flagged$EventID, flagged$Country, flagged$ResponseVariable),
    ,
    drop = FALSE
  ]

  flagged <- head(flagged, MAX_FLAGGED)

  for (i in seq_len(nrow(flagged))) {
    f <- flagged[i, , drop = FALSE]

    z <- d[
      d$EventID == f$EventID &
        d$AnchorType == f$AnchorType &
        d$Country == f$Country &
        d$ResponseVariable == f$ResponseVariable,
      ,
      drop = FALSE
    ]

    if (!nrow(z)) next

    v <- f$ResponseVariable

    status_bits <- c(
      sprintf("impact: %s", f$ImpactSignClass90),
      sprintf("peak h=%s: %s", f$PeakHorizon, f$PeakSignClass90)
    )

    if (isTRUE(f$CredibleSignReversal)) {
      status_bits <- c(status_bits, "credible sign reversal")
    }
    if (isTRUE(f$AllHorizonsCrossZero90)) {
      status_bits <- c(status_bits, "all horizons cross zero")
    }

    ttl <- sprintf(
      "%s | %s | %s | %s",
      event_label_for(z),
      f$Country,
      anchor_label(z),
      unname(var_title[v])
    )

    p <- ggplot(z, aes(x = Horizon, y = median)) +
      geom_hline(yintercept = 0, linewidth = 0.40) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.22) +
      geom_line(linewidth = 0.75) +
      geom_point(size = 1.3) +
      labs(
        title = ttl,
        subtitle = paste(
          paste(status_bits, collapse = " | "),
          subtitle_for_var(z, v),
          sep = "\n"
        ),
        x = "Horizon (quarters)",
        y = unname(var_axis[v])
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.5)
      )

    fn <- sprintf(
      "%s__%s__%s__%s.png",
      sanitize(f$EventID),
      ifelse(f$AnchorType == "EVENT_QUARTER_t0", "t0", "t1"),
      sanitize(f$Country),
      sanitize(v)
    )

    path <- file.path(FLAG_DIR, fn)

    ggsave(
      filename = path,
      plot = p,
      width = 8,
      height = 5.5,
      dpi = 180
    )

    midx <- midx + 1L
    manifest[[midx]] <- data.frame(
      PlotType = "FLAGGED_INDIVIDUAL",
      EventSet = f$EventSet,
      EventID = f$EventID,
      AnchorType = f$AnchorType,
      Country = f$Country,
      ResponseVariable = v,
      RuntimeRateMode = RATE_MODE,
      RuntimeEquityReturnMode = EQ_RETURN_MODE,
      File = path,
      stringsAsFactors = FALSE
    )
  }
}

manifest_df <- if (length(manifest)) do.call(rbind, manifest) else data.frame()

write.csv(
  manifest_df,
  file.path(AUDIT_DIR, "06_irf_plot_manifest.csv"),
  row.names = FALSE
)

readme <- c(
  "TRANSFORMATION-AWARE IRF DIAGNOSTIC FIGURES",
  "============================================",
  "",
  "Figures are generated from audited numerical posterior summaries.",
  "No significance decision is made from image appearance alone.",
  "",
  sprintf("Rate mode: %s", RATE_MODE),
  sprintf("Equity return mode: %s", EQ_RETURN_MODE),
  "",
  "Variable meanings:",
  if (RATE_MODE == "level") {
    "- r   = interest-rate level"
  } else {
    "- r   = Delta interest rate"
  },
  "- de  = REER_DLOG",
  "- deq = equity return",
  "- CPI is not present.",
  "",
  "Main plots:",
  "- results/irf_audit/plots/core/",
  "- results/irf_audit/plots/appendix/",
  "",
  "Flagged individual plots:",
  "- results/irf_audit/plots/flagged/",
  "",
  "Every main figure shows:",
  "- posterior median",
  "- 90% posterior credible interval (p05-p95)",
  "- zero line",
  "- 14 country panels",
  "",
  "Interpretation:",
  "- interval fully above zero = positive credible response at 90%",
  "- interval fully below zero = negative credible response at 90%",
  "- interval crossing zero = direction is not 90% credible",
  "- r level IRFs are already level deviations",
  "- Delta-r IRFs require cumulation for rate-level changes",
  "- de is already REER_DLOG",
  "- deq is already an equity return"
)

writeLines(readme, file.path(AUDIT_DIR, "README_irf_plots.txt"))

msg("IRF diagnostic plots complete.")
msg("Manifest rows: %d", nrow(manifest_df))
msg("Rate mode: %s", RATE_MODE)
msg("Core plot directory: %s", CORE_DIR)
msg("Appendix plot directory: %s", APP_DIR)
msg("Flagged plot directory: %s", FLAG_DIR)
