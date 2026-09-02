#!/usr/bin/env Rscript

# =============================================================================
# 64_domestic_lag_source_diagnostic.R
#
# 08i — Domestic-lag source diagnostic for the accepted formal TVP-GVAR posterior.
#
# Purpose
# -------
# 08h showed that removing the complete domestic lag block A_i produces by far
# the largest improvement in posterior global stability. 08i therefore drills
# into A_i without re-estimating the model.
#
# For each event anchor and posterior draw, the accepted R62 system is
#
#   G0 = I - C0
#   G1 = D + L
#   F  = solve(G0) %*% G1
#
# where
#   D = blockdiag(A_1, ..., A_N)  is the domestic lag contribution,
#   C0                            is the contemporaneous foreign contribution,
#   L                             is the lagged foreign contribution.
#
# 08i keeps G0 and L unchanged and modifies only entries of D.
#
# Diagnostic families
# -------------------
# GLOBAL:
#   BASELINE
#   NO_DOMESTIC_LAG
#   ZERO_ALL_OWN_LAGS
#   ZERO_ALL_CROSS_LAGS
#
# COUNTRY_GROUP:
#   ZERO_COUNTRY_OWN_<CC>
#   ZERO_COUNTRY_CROSS_<CC>
#
# EQUATION_ROW:
#   ZERO_EQ_<CC>_<response>
#
# PREDICTOR_COLUMN:
#   ZERO_PRED_<CC>_<predictor>
#
# INDIVIDUAL_COEFFICIENT:
#   ZERO_COEF_<CC>_<response>_FROM_<predictor>
#   (only for FIN3_08I_TARGET_COUNTRIES; default ZA,BR,JP)
#
# IMPORTANT
# ---------
# - Diagnostic only: no re-estimation.
# - The formal R62 >= 90% stability gate is NOT relaxed.
# - A zeroed coefficient/block is an attribution counterfactual, not a valid
#   alternative publication specification.
# - Publication changes must be separately estimated and revalidated.
# =============================================================================

source("R/00_config.R")

get_env_num <- function(name, default) {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) return(default)
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stopf("Environment variable %s is not numeric: %s", name, z)
  out
}

get_env_chr <- function(name, default = "") {
  z <- trimws(Sys.getenv(name, ""))
  if (!nzchar(z)) default else z
}

PARTS_ROOT <- get_env_chr("FIN3_PARTS_ROOT", "posterior_parts")
CONV_DIR <- get_env_chr("FIN3_CONVERGENCE_DIR", file.path(RESULTS_DIR, "formal_tvp"))
REF08H_DIR <- get_env_chr("FIN3_08H_REFERENCE_DIR", "reference_08h")
NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 2000))
DIAG_DRAWS_PER_CHAIN <- as.integer(get_env_num("FIN3_DIAG_DRAWS_PER_CHAIN", 250))
MIN_G0_RCOND <- get_env_num("FIN3_IRF_MIN_G0_RCOND", 1e-10)
DIAG_CORES <- as.integer(get_env_num("FIN3_DIAG_CORES", 2))
TARGET_RAW <- get_env_chr("FIN3_08I_TARGET_COUNTRIES", "ZA,BR,JP")

if (NCHAINS != 4L) stopf("08i expects the accepted formal 4-chain posterior.")
if (EXPECTED_STORED < 200L) stopf("Expected stored draws are implausibly small.")
if (DIAG_DRAWS_PER_CHAIN < 50L) stopf("Use at least 50 diagnostic draws per chain.")
if (DIAG_DRAWS_PER_CHAIN > EXPECTED_STORED) DIAG_DRAWS_PER_CHAIN <- EXPECTED_STORED
if (DIAG_CORES < 1L) DIAG_CORES <- 1L

TARGET_COUNTRIES <- unique(toupper(trimws(strsplit(TARGET_RAW, ",", fixed = TRUE)[[1]])))
TARGET_COUNTRIES <- TARGET_COUNTRIES[nzchar(TARGET_COUNTRIES)]
bad_target <- setdiff(TARGET_COUNTRIES, COUNTRIES)
if (length(bad_target)) {
  stopf("Unknown FIN3_08I_TARGET_COUNTRIES: %s", paste(bad_target, collapse = ", "))
}
if (!length(TARGET_COUNTRIES)) stopf("At least one target country is required.")

CONV_GATE <- file.path(CONV_DIR, "00_formal_mcmc_gate.csv")
if (!file.exists(CONV_GATE)) stopf("Formal MCMC convergence gate missing: %s", CONV_GATE)
cg <- read.csv(CONV_GATE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(cg) != 1L || !"Status" %in% names(cg)) stopf("Malformed formal MCMC gate.")
if (!identical(trimws(cg$Status[1]), "READY_FOR_FORMAL_IRF")) {
  stopf("Formal MCMC gate is not READY_FOR_FORMAL_IRF.")
}

OUT <- file.path(RESULTS_DIR, "domestic_lag_source")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load and validate exact accepted 14-country x 4-chain posterior grid
# =============================================================================

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

expected_parts <- length(COUNTRIES) * NCHAINS
if (length(files) != expected_parts) {
  stopf("Expected exactly %d posterior RDS files; found %d.", expected_parts, length(files))
}

parts <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)
parts_index <- vector("list", expected_parts)
ip <- 0L

for (cc in COUNTRIES) {
  parts[[cc]] <- vector("list", NCHAINS)
  for (ch in seq_len(NCHAINS)) {
    candidates <- files[
      grepl(sprintf("formal_tvp_%s_chain%d\\.rds$", cc, ch), files)
    ]
    if (length(candidates) != 1L) {
      stopf("Expected one posterior part for %s chain %d; found %d.", cc, ch, length(candidates))
    }

    z <- readRDS(candidates)
    if (!identical(as.character(z$meta$country), cc)) stopf("Country metadata mismatch: %s", candidates)
    if (as.integer(z$meta$chain) != ch) stopf("Chain metadata mismatch: %s", candidates)
    if (as.integer(z$meta$stored_draws) != EXPECTED_STORED) {
      stopf("Stored draw mismatch: %s chain %d.", cc, ch)
    }
    if (is.null(z$event_coef) || any(!is.finite(z$event_coef))) {
      stopf("Missing/non-finite event_coef: %s chain %d.", cc, ch)
    }

    parts[[cc]][[ch]] <- z
    ip <- ip + 1L
    parts_index[[ip]] <- data.frame(
      Country = cc,
      Chain = ch,
      Burn = as.integer(z$meta$burn),
      InternalPostBurn = as.integer(z$meta$keep_internal),
      StoredDraws = as.integer(z$meta$stored_draws),
      File = candidates,
      stringsAsFactors = FALSE
    )
  }
}

parts_index <- do.call(rbind, parts_index)
write.csv(parts_index, file.path(OUT, "00_parts_index.csv"), row.names = FALSE)

ref <- parts[[COUNTRIES[1]]][[1]]
anchors <- ref$anchors
KEY_TERMS <- ref$meta$key_terms

required_key_terms <- c(
  paste0(VARS, "_L1"),
  paste0(VARS, "_star_0"),
  paste0(VARS, "_star_L1"),
  "gpr_0", "gpr_L1"
)
if (!all(required_key_terms %in% KEY_TERMS)) {
  stopf("Posterior parts do not contain the coefficient blocks required by R62/08h.")
}

anchor_key <- paste(anchors$EventID, anchors$AnchorType, anchors$AnchorQuarter)
for (cc in COUNTRIES) {
  for (ch in seq_len(NCHAINS)) {
    z <- parts[[cc]][[ch]]
    if (!identical(z$meta$key_terms, KEY_TERMS)) {
      stopf("Key-term ordering mismatch: %s chain %d.", cc, ch)
    }
    zk <- paste(z$anchors$EventID, z$anchors$AnchorType, z$anchors$AnchorQuarter)
    if (!identical(zk, anchor_key)) {
      stopf("Event-anchor ordering mismatch: %s chain %d.", cc, ch)
    }
  }
}

draw_idx <- unique(as.integer(round(seq(
  1, EXPECTED_STORED,
  length.out = DIAG_DRAWS_PER_CHAIN
))))
if (length(draw_idx) < 50L) stopf("Too few unique diagnostic draw indices after thinning.")
ND_PER_CHAIN <- length(draw_idx)
ND_TOTAL <- ND_PER_CHAIN * NCHAINS

# =============================================================================
# 2. Exact R62/08h matrix mapping
# =============================================================================

N <- length(COUNTRIES)
K <- length(VARS)
NK <- N * K

if (K != 3L) stopf("08i is designed for the current 3-variable financial system.")

W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])

term_dom <- paste0(VARS, "_L1")
term_star0 <- paste0(VARS, "_star_0")
term_star1 <- paste0(VARS, "_star_L1")

idx_dom <- match(term_dom, KEY_TERMS)
idx_star0 <- match(term_star0, KEY_TERMS)
idx_star1 <- match(term_star1, KEY_TERMS)
if (any(is.na(c(idx_dom, idx_star0, idx_star1)))) {
  stopf("Failed to map domestic/foreign blocks.")
}

selectors <- vector("list", N)
star_maps <- vector("list", N)

for (i in seq_len(N)) {
  S <- matrix(0, K, NK)
  S[, ((i - 1L) * K + 1L):(i * K)] <- diag(K)
  selectors[[i]] <- S

  R <- matrix(0, K, NK)
  for (j in seq_len(N)) {
    for (v in seq_len(K)) {
      R[v, (j - 1L) * K + v] <- W[i, j]
    }
  }
  star_maps[[i]] <- R
}

global_pos <- function(i, v) (i - 1L) * K + v

extract_blocks <- function(ch, dd, a) {
  A <- vector("list", N)
  B0 <- vector("list", N)
  B1 <- vector("list", N)

  for (i in seq_len(N)) {
    cc <- COUNTRIES[i]
    raw <- parts[[cc]][[ch]]$event_coef[dd, a, , , drop = FALSE]
    coef <- matrix(
      raw,
      nrow = K,
      ncol = length(KEY_TERMS),
      dimnames = list(VARS, KEY_TERMS)
    )
    A[[i]] <- coef[, idx_dom, drop = FALSE]
    B0[[i]] <- coef[, idx_star0, drop = FALSE]
    B1[[i]] <- coef[, idx_star1, drop = FALSE]
  }

  list(A = A, B0 = B0, B1 = B1)
}

build_components <- function(blocks) {
  C0 <- matrix(0, NK, NK)
  D  <- matrix(0, NK, NK)
  L  <- matrix(0, NK, NK)

  for (i in seq_len(N)) {
    rr <- ((i - 1L) * K + 1L):(i * K)
    C0[rr, ] <- blocks$B0[[i]] %*% star_maps[[i]]
    D[rr, ]  <- blocks$A[[i]] %*% selectors[[i]]
    L[rr, ]  <- blocks$B1[[i]] %*% star_maps[[i]]
  }

  list(C0 = C0, D = D, L = L)
}

safe_impact <- function(G0) {
  rc <- tryCatch(rcond(G0), error = function(e) NA_real_)
  if (!is.finite(rc) || rc < MIN_G0_RCOND) {
    return(list(rcond = rc, impact = NULL))
  }
  imp <- tryCatch(solve(G0), error = function(e) NULL)
  if (is.null(imp) || any(!is.finite(imp))) imp <- NULL
  list(rcond = rc, impact = imp)
}

safe_rho <- function(impact, G1) {
  if (is.null(impact)) return(NA_real_)
  Fmat <- impact %*% G1
  tryCatch(
    max(Mod(eigen(Fmat, only.values = TRUE)$values)),
    error = function(e) NA_real_
  )
}

# =============================================================================
# 3. Scenario contract
# =============================================================================

scenario_rows <- list()
sid <- 0L

add_scenario <- function(scenario, family, country = NA_character_,
                         response = NA_character_, predictor = NA_character_,
                         definition = "") {
  sid <<- sid + 1L
  scenario_rows[[sid]] <<- data.frame(
    Scenario = scenario,
    ScenarioFamily = family,
    TargetCountry = country,
    Response = response,
    Predictor = predictor,
    Definition = definition,
    stringsAsFactors = FALSE
  )
}

add_scenario(
  "BASELINE", "GLOBAL",
  definition = "Accepted R62 system: G0=I-C0; G1=D+L"
)
add_scenario(
  "NO_DOMESTIC_LAG", "GLOBAL",
  definition = "Set all domestic lag entries D=0; retain G0 and L"
)
add_scenario(
  "ZERO_ALL_OWN_LAGS", "GLOBAL",
  definition = "Set every country A_i diagonal element to zero; retain domestic cross-lags"
)
add_scenario(
  "ZERO_ALL_CROSS_LAGS", "GLOBAL",
  definition = "Set every country A_i off-diagonal element to zero; retain domestic own-lags"
)

for (cc in COUNTRIES) {
  add_scenario(
    paste0("ZERO_COUNTRY_OWN_", cc), "COUNTRY_GROUP", country = cc,
    definition = paste0("Set the three domestic own-lags in ", cc, " to zero")
  )
  add_scenario(
    paste0("ZERO_COUNTRY_CROSS_", cc), "COUNTRY_GROUP", country = cc,
    definition = paste0("Set the six domestic cross-lags in ", cc, " to zero")
  )
}

for (cc in COUNTRIES) {
  for (resp in VARS) {
    add_scenario(
      paste0("ZERO_EQ_", cc, "_", resp),
      "EQUATION_ROW", country = cc, response = resp,
      definition = paste0("Zero the entire domestic lag row in ", cc, " response equation ", resp)
    )
  }
}

for (cc in COUNTRIES) {
  for (pred in VARS) {
    add_scenario(
      paste0("ZERO_PRED_", cc, "_", pred),
      "PREDICTOR_COLUMN", country = cc, predictor = pred,
      definition = paste0("Zero domestic lag predictor ", pred, "_L1 across all ", cc, " equations")
    )
  }
}

for (cc in TARGET_COUNTRIES) {
  for (resp in VARS) {
    for (pred in VARS) {
      add_scenario(
        paste0("ZERO_COEF_", cc, "_", resp, "_FROM_", pred),
        "INDIVIDUAL_COEFFICIENT",
        country = cc,
        response = resp,
        predictor = pred,
        definition = paste0("Zero only ", cc, " A[", resp, ",", pred, "]")
      )
    }
  }
}

scenario_contract <- do.call(rbind, scenario_rows)
write.csv(
  scenario_contract,
  file.path(OUT, "00_scenario_contract.csv"),
  row.names = FALSE
)

SCENARIOS <- scenario_contract$Scenario
NS <- length(SCENARIOS)
scenario_pos <- setNames(seq_len(NS), SCENARIOS)

# Precompute matrix-index operations for every scenario.
ops <- vector("list", NS)

for (s in seq_len(NS)) {
  sc <- scenario_contract[s, ]

  if (sc$Scenario == "BASELINE") {
    ops[[s]] <- list(type = "NONE")
    next
  }

  if (sc$Scenario == "NO_DOMESTIC_LAG") {
    ops[[s]] <- list(type = "ZERO_ALL")
    next
  }

  if (sc$Scenario == "ZERO_ALL_OWN_LAGS") {
    rows <- cols <- integer(0)
    for (i in seq_len(N)) {
      for (v in seq_len(K)) {
        rows <- c(rows, global_pos(i, v))
        cols <- c(cols, global_pos(i, v))
      }
    }
    ops[[s]] <- list(type = "CELLS", rows = rows, cols = cols)
    next
  }

  if (sc$Scenario == "ZERO_ALL_CROSS_LAGS") {
    rows <- cols <- integer(0)
    for (i in seq_len(N)) {
      for (resp in seq_len(K)) {
        for (pred in seq_len(K)) {
          if (resp != pred) {
            rows <- c(rows, global_pos(i, resp))
            cols <- c(cols, global_pos(i, pred))
          }
        }
      }
    }
    ops[[s]] <- list(type = "CELLS", rows = rows, cols = cols)
    next
  }

  cc <- sc$TargetCountry
  ci <- match(cc, COUNTRIES)

  if (sc$ScenarioFamily == "COUNTRY_GROUP") {
    rows <- cols <- integer(0)
    if (grepl("^ZERO_COUNTRY_OWN_", sc$Scenario)) {
      for (v in seq_len(K)) {
        rows <- c(rows, global_pos(ci, v))
        cols <- c(cols, global_pos(ci, v))
      }
    } else {
      for (resp in seq_len(K)) {
        for (pred in seq_len(K)) {
          if (resp != pred) {
            rows <- c(rows, global_pos(ci, resp))
            cols <- c(cols, global_pos(ci, pred))
          }
        }
      }
    }
    ops[[s]] <- list(type = "CELLS", rows = rows, cols = cols)
    next
  }

  if (sc$ScenarioFamily == "EQUATION_ROW") {
    ri <- match(sc$Response, VARS)
    rows <- rep(global_pos(ci, ri), K)
    cols <- vapply(seq_len(K), function(pred) global_pos(ci, pred), integer(1))
    ops[[s]] <- list(type = "CELLS", rows = rows, cols = cols)
    next
  }

  if (sc$ScenarioFamily == "PREDICTOR_COLUMN") {
    pi <- match(sc$Predictor, VARS)
    rows <- vapply(seq_len(K), function(resp) global_pos(ci, resp), integer(1))
    cols <- rep(global_pos(ci, pi), K)
    ops[[s]] <- list(type = "CELLS", rows = rows, cols = cols)
    next
  }

  if (sc$ScenarioFamily == "INDIVIDUAL_COEFFICIENT") {
    ri <- match(sc$Response, VARS)
    pi <- match(sc$Predictor, VARS)
    ops[[s]] <- list(
      type = "CELLS",
      rows = global_pos(ci, ri),
      cols = global_pos(ci, pi)
    )
    next
  }

  stopf("Unhandled scenario: %s", sc$Scenario)
}

apply_domestic_operation <- function(D, op) {
  if (op$type == "NONE") return(D)
  if (op$type == "ZERO_ALL") return(D * 0)
  out <- D
  out[cbind(op$rows, op$cols)] <- 0
  out
}

# =============================================================================
# 4. Anchor-level posterior stability under each domestic-lag counterfactual
# =============================================================================

summarize_rho <- function(rho, rc) {
  finite <- is.finite(rho)
  g0ok <- is.finite(rc) & rc >= MIN_G0_RCOND
  stable <- finite & rho < 1 & g0ok
  z <- rho[finite]

  data.frame(
    DiagnosticDraws = length(rho),
    FiniteRhoShare = mean(finite),
    G0OKShare = mean(g0ok),
    StableShare = mean(stable),
    RhoMedian = if (length(z)) median(z) else NA_real_,
    RhoP90 = if (length(z)) unname(quantile(z, .90)) else NA_real_,
    RhoP95 = if (length(z)) unname(quantile(z, .95)) else NA_real_,
    RhoMax = if (length(z)) max(z) else NA_real_,
    stringsAsFactors = FALSE
  )
}

process_anchor <- function(a) {
  meta <- anchors[a, , drop = FALSE]

  rho_mat <- matrix(NA_real_, ND_TOTAL, NS, dimnames = list(NULL, SCENARIOS))
  rc_vec <- rep(NA_real_, ND_TOTAL)

  # Store actual domestic coefficients for posterior magnitude summaries.
  coef_arr <- array(
    NA_real_,
    dim = c(ND_TOTAL, N, K, K),
    dimnames = list(NULL, COUNTRIES, VARS, VARS)
  )

  gd <- 0L

  for (ch in seq_len(NCHAINS)) {
    for (dd in draw_idx) {
      gd <- gd + 1L

      blocks <- extract_blocks(ch, dd, a)
      comp <- build_components(blocks)

      G0 <- diag(NK) - comp$C0
      imp <- safe_impact(G0)
      rc_vec[gd] <- imp$rcond

      for (i in seq_len(N)) {
        coef_arr[gd, i, , ] <- blocks$A[[i]]
      }

      if (is.null(imp$impact)) next

      for (s in seq_len(NS)) {
        Dcf <- apply_domestic_operation(comp$D, ops[[s]])
        rho_mat[gd, s] <- safe_rho(imp$impact, Dcf + comp$L)
      }
    }
  }

  scenario_rows <- vector("list", NS)
  for (s in seq_len(NS)) {
    sm <- summarize_rho(rho_mat[, s], rc_vec)
    scenario_rows[[s]] <- cbind(
      data.frame(
        EventID = meta$EventID,
        EventSet = meta$EventSet,
        EventLabel = meta$EventLabel,
        ShockFamily = meta$ShockFamily,
        EventQuarter = meta$EventQuarter,
        AnchorType = meta$AnchorType,
        AnchorQuarter = meta$AnchorQuarter,
        Scenario = scenario_contract$Scenario[s],
        ScenarioFamily = scenario_contract$ScenarioFamily[s],
        TargetCountry = scenario_contract$TargetCountry[s],
        Response = scenario_contract$Response[s],
        Predictor = scenario_contract$Predictor[s],
        stringsAsFactors = FALSE
      ),
      sm
    )
  }

  coef_rows <- vector("list", N * K * K)
  kk <- 0L
  for (i in seq_len(N)) {
    for (resp in seq_len(K)) {
      for (pred in seq_len(K)) {
        kk <- kk + 1L
        x <- coef_arr[, i, resp, pred]
        x <- x[is.finite(x)]
        coef_rows[[kk]] <- data.frame(
          EventID = meta$EventID,
          EventSet = meta$EventSet,
          AnchorType = meta$AnchorType,
          AnchorQuarter = meta$AnchorQuarter,
          Country = COUNTRIES[i],
          Response = VARS[resp],
          Predictor = VARS[pred],
          IsOwnLag = resp == pred,
          N = length(x),
          Mean = if (length(x)) mean(x) else NA_real_,
          Median = if (length(x)) median(x) else NA_real_,
          MeanAbs = if (length(x)) mean(abs(x)) else NA_real_,
          P05 = if (length(x)) unname(quantile(x, .05)) else NA_real_,
          P95 = if (length(x)) unname(quantile(x, .95)) else NA_real_,
          ProbPositive = if (length(x)) mean(x > 0) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    scenario = do.call(rbind, scenario_rows),
    coef = do.call(rbind, coef_rows)
  )
}

anchor_ids <- seq_len(nrow(anchors))
if (.Platform$OS.type == "unix" && DIAG_CORES > 1L) {
  ans <- parallel::mclapply(
    anchor_ids,
    process_anchor,
    mc.cores = min(DIAG_CORES, length(anchor_ids)),
    mc.preschedule = TRUE
  )
} else {
  ans <- lapply(anchor_ids, process_anchor)
}

scenario_anchor <- do.call(rbind, lapply(ans, `[[`, "scenario"))
coef_anchor <- do.call(rbind, lapply(ans, `[[`, "coef"))

# Baseline comparison within each anchor.
baseline_lookup <- scenario_anchor[
  scenario_anchor$Scenario == "BASELINE",
  c("EventID", "AnchorType", "AnchorQuarter", "StableShare", "RhoMedian")
]
names(baseline_lookup)[4:5] <- c("BaselineStableShare", "BaselineRhoMedian")

scenario_anchor <- merge(
  scenario_anchor,
  baseline_lookup,
  by = c("EventID", "AnchorType", "AnchorQuarter"),
  all.x = TRUE,
  sort = FALSE
)
scenario_anchor$DeltaStableShareVsBaseline <-
  scenario_anchor$StableShare - scenario_anchor$BaselineStableShare
scenario_anchor$DeltaRhoMedianVsBaseline <-
  scenario_anchor$RhoMedian - scenario_anchor$BaselineRhoMedian

write.csv(
  scenario_anchor,
  file.path(OUT, "01_scenario_stability_by_anchor.csv"),
  row.names = FALSE
)

# =============================================================================
# 5. Overall scenario rankings
# =============================================================================

overall <- do.call(rbind, lapply(SCENARIOS, function(sc) {
  d <- scenario_anchor[scenario_anchor$Scenario == sc, , drop = FALSE]
  meta <- scenario_contract[scenario_contract$Scenario == sc, , drop = FALSE]
  data.frame(
    Scenario = sc,
    ScenarioFamily = meta$ScenarioFamily,
    TargetCountry = meta$TargetCountry,
    Response = meta$Response,
    Predictor = meta$Predictor,
    Anchors = nrow(d),
    MeanStableShare = mean(d$StableShare, na.rm = TRUE),
    MinStableShare = min(d$StableShare, na.rm = TRUE),
    MedianStableShare = median(d$StableShare, na.rm = TRUE),
    MeanRhoMedian = mean(d$RhoMedian, na.rm = TRUE),
    MaxRhoP95 = max(d$RhoP95, na.rm = TRUE),
    MeanDeltaStableShareVsBaseline = mean(d$DeltaStableShareVsBaseline, na.rm = TRUE),
    MeanDeltaRhoMedianVsBaseline = mean(d$DeltaRhoMedianVsBaseline, na.rm = TRUE),
    AnchorsStableShareGE90 = sum(d$StableShare >= .90, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

overall <- overall[
  order(-overall$MeanDeltaStableShareVsBaseline, overall$MeanRhoMedian),
  ,
  drop = FALSE
]
write.csv(overall, file.path(OUT, "02_scenario_overall_summary.csv"), row.names = FALSE)

global_summary <- overall[overall$ScenarioFamily == "GLOBAL", , drop = FALSE]
country_group <- overall[overall$ScenarioFamily == "COUNTRY_GROUP", , drop = FALSE]
eq_summary <- overall[overall$ScenarioFamily == "EQUATION_ROW", , drop = FALSE]
pred_summary <- overall[overall$ScenarioFamily == "PREDICTOR_COLUMN", , drop = FALSE]
coef_scenario_summary <- overall[
  overall$ScenarioFamily == "INDIVIDUAL_COEFFICIENT",
  ,
  drop = FALSE
]

country_group <- country_group[
  order(-country_group$MeanDeltaStableShareVsBaseline),
  ,
  drop = FALSE
]
eq_summary <- eq_summary[
  order(-eq_summary$MeanDeltaStableShareVsBaseline),
  ,
  drop = FALSE
]
pred_summary <- pred_summary[
  order(-pred_summary$MeanDeltaStableShareVsBaseline),
  ,
  drop = FALSE
]
coef_scenario_summary <- coef_scenario_summary[
  order(-coef_scenario_summary$MeanDeltaStableShareVsBaseline),
  ,
  drop = FALSE
]

write.csv(global_summary, file.path(OUT, "03_global_own_vs_cross_summary.csv"), row.names = FALSE)
write.csv(country_group, file.path(OUT, "04_country_own_cross_attribution.csv"), row.names = FALSE)
write.csv(eq_summary, file.path(OUT, "05_equation_row_attribution.csv"), row.names = FALSE)
write.csv(pred_summary, file.path(OUT, "06_predictor_column_attribution.csv"), row.names = FALSE)
write.csv(
  coef_scenario_summary,
  file.path(OUT, "07_target_individual_coefficient_attribution.csv"),
  row.names = FALSE
)

# Posterior magnitude summary of actual A coefficients, all countries.
write.csv(
  coef_anchor,
  file.path(OUT, "08_domestic_coefficient_posterior_by_anchor.csv"),
  row.names = FALSE
)

coef_overall <- do.call(
  rbind,
  lapply(split(
    coef_anchor,
    interaction(
      coef_anchor$Country,
      coef_anchor$Response,
      coef_anchor$Predictor,
      drop = TRUE
    )
  ), function(d) {
    data.frame(
      Country = d$Country[1],
      Response = d$Response[1],
      Predictor = d$Predictor[1],
      IsOwnLag = d$IsOwnLag[1],
      Anchors = nrow(d),
      MeanPosteriorMedian = mean(d$Median, na.rm = TRUE),
      MeanAbsCoefficient = mean(d$MeanAbs, na.rm = TRUE),
      MeanProbPositive = mean(d$ProbPositive, na.rm = TRUE),
      MeanP05 = mean(d$P05, na.rm = TRUE),
      MeanP95 = mean(d$P95, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
coef_overall <- coef_overall[
  order(-coef_overall$MeanAbsCoefficient),
  ,
  drop = FALSE
]
write.csv(
  coef_overall,
  file.path(OUT, "09_domestic_coefficient_posterior_overall.csv"),
  row.names = FALSE
)

# =============================================================================
# 6. Integrity comparison to accepted 08h baseline
# =============================================================================

baseline_row <- overall[overall$Scenario == "BASELINE", , drop = FALSE]
nodom_row <- overall[overall$Scenario == "NO_DOMESTIC_LAG", , drop = FALSE]
if (nrow(baseline_row) != 1L || nrow(nodom_row) != 1L) {
  stopf("Malformed baseline/global domestic-lag summaries.")
}

ref08h_file <- file.path(REF08H_DIR, "00_dynamic_stability_source_gate.csv")
ref_match_status <- "REFERENCE_NOT_AVAILABLE"
ref08h_baseline <- NA_real_
baseline_gap <- NA_real_

if (file.exists(ref08h_file)) {
  rg <- read.csv(ref08h_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (
    nrow(rg) == 1L &&
    all(c("DiagnosticDrawsPerChain", "BaselineMeanStableShare") %in% names(rg))
  ) {
    ref08h_baseline <- as.numeric(rg$BaselineMeanStableShare[1])
    if (
      is.finite(ref08h_baseline) &&
      as.integer(rg$DiagnosticDrawsPerChain[1]) == ND_PER_CHAIN
    ) {
      baseline_gap <- abs(baseline_row$MeanStableShare - ref08h_baseline)
      ref_match_status <- if (baseline_gap <= 1e-10) {
        "EXACT_08H_BASELINE_MATCH"
      } else {
        "08H_BASELINE_MISMATCH"
      }
    } else {
      ref_match_status <- "REFERENCE_DRAW_COUNT_DIFFERS"
    }
  } else {
    ref_match_status <- "REFERENCE_MALFORMED"
  }
}

integrity_ok <- all(scenario_anchor$FiniteRhoShare >= .99) &&
  all(scenario_anchor$G0OKShare >= .99) &&
  nrow(parts_index) == expected_parts &&
  !identical(ref_match_status, "08H_BASELINE_MISMATCH")

family_best <- function(fam) {
  d <- overall[overall$ScenarioFamily == fam, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(-d$MeanDeltaStableShareVsBaseline), , drop = FALSE]
  d[1, , drop = FALSE]
}

best_country_group <- family_best("COUNTRY_GROUP")
best_eq <- family_best("EQUATION_ROW")
best_pred <- family_best("PREDICTOR_COLUMN")
best_coef <- family_best("INDIVIDUAL_COEFFICIENT")

status <- if (integrity_ok) "DIAGNOSTIC_COMPLETE" else "DIAGNOSTIC_INTEGRITY_WARNING"

gate <- data.frame(
  Status = status,
  MainNetwork = MAIN_NETWORK,
  Countries = N,
  Chains = NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  DiagnosticDrawsPerChain = ND_PER_CHAIN,
  DiagnosticDrawsTotal = ND_TOTAL,
  EventAnchors = nrow(anchors),
  Scenarios = NS,
  TargetCountries = paste(TARGET_COUNTRIES, collapse = ","),
  BaselineMeanStableShare = baseline_row$MeanStableShare,
  NoDomesticLagMeanStableShare = nodom_row$MeanStableShare,
  NoDomesticLagDeltaStableShare = nodom_row$MeanDeltaStableShareVsBaseline,
  BestCountryGroupScenario = if (!is.null(best_country_group)) best_country_group$Scenario else NA_character_,
  BestCountryGroupDelta = if (!is.null(best_country_group)) best_country_group$MeanDeltaStableShareVsBaseline else NA_real_,
  BestEquationScenario = if (!is.null(best_eq)) best_eq$Scenario else NA_character_,
  BestEquationDelta = if (!is.null(best_eq)) best_eq$MeanDeltaStableShareVsBaseline else NA_real_,
  BestPredictorScenario = if (!is.null(best_pred)) best_pred$Scenario else NA_character_,
  BestPredictorDelta = if (!is.null(best_pred)) best_pred$MeanDeltaStableShareVsBaseline else NA_real_,
  BestIndividualCoefficientScenario = if (!is.null(best_coef)) best_coef$Scenario else NA_character_,
  BestIndividualCoefficientDelta = if (!is.null(best_coef)) best_coef$MeanDeltaStableShareVsBaseline else NA_real_,
  Reference08hBaseline = ref08h_baseline,
  Reference08hBaselineAbsGap = baseline_gap,
  Reference08hStatus = ref_match_status,
  FormalIRFGateWasRelaxed = FALSE,
  ReestimationPerformed = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  gate,
  file.path(OUT, "00_domestic_lag_source_gate.csv"),
  row.names = FALSE
)

contract <- data.frame(
  Item = c(
    "Purpose",
    "Posterior",
    "MatrixDefinition",
    "Thinning",
    "FormalGate",
    "GlobalOwnLag",
    "GlobalCrossLag",
    "CountryGroups",
    "EquationRows",
    "PredictorColumns",
    "IndividualCoefficients",
    "PublicationUse"
  ),
  Value = c(
    "Drill into the domestic lag block A_i after 08h identified it as the leading source of instability",
    "Exact accepted 14-country x 4-chain posterior grid (52 original non-AU + 4 rescued AU)",
    "Same G0/G1 construction as R62 and 08h; only D=blockdiag(A_i) is modified",
    sprintf("Deterministic balanced thinning: %d of %d stored draws per chain", ND_PER_CHAIN, EXPECTED_STORED),
    "Unchanged; no stability threshold is lowered",
    "Zero all diagonal A_i entries across all countries",
    "Zero all off-diagonal A_i entries across all countries",
    "For every country separately: zero own-lag group or cross-lag group",
    "For every country x response equation: zero all three domestic lag coefficients in that equation",
    "For every country x predictor variable: zero that domestic lag predictor across all three equations",
    paste0("For ", paste(TARGET_COUNTRIES, collapse = ","), ": zero each of 9 A_i entries individually"),
    "All zeroing exercises are attribution diagnostics; any model change requires re-estimation"
  ),
  stringsAsFactors = FALSE
)
write.csv(contract, file.path(OUT, "10_diagnostic_contract.csv"), row.names = FALSE)

fmt_best <- function(x, label) {
  if (is.null(x)) return(sprintf("- %s: unavailable", label))
  sprintf(
    "- %s: %s (stable-share delta %.6f)",
    label, x$Scenario[1], x$MeanDeltaStableShareVsBaseline[1]
  )
}

readme <- c(
  sprintf("08i DOMESTIC LAG SOURCE DIAGNOSTIC: %s", status),
  "============================================================",
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf("Countries: %d", N),
  sprintf("Chains: %d", NCHAINS),
  sprintf("Stored draws per chain: %d", EXPECTED_STORED),
  sprintf("Diagnostic draws per chain: %d", ND_PER_CHAIN),
  sprintf("Total diagnostic draws per anchor: %d", ND_TOTAL),
  sprintf("Event anchors: %d", nrow(anchors)),
  sprintf("Scenarios: %d", NS),
  sprintf("Individual-coefficient target countries: %s", paste(TARGET_COUNTRIES, collapse = ",")),
  sprintf("Baseline mean stable share: %.6f", baseline_row$MeanStableShare),
  sprintf("NO_DOMESTIC_LAG mean stable share: %.6f", nodom_row$MeanStableShare),
  sprintf("NO_DOMESTIC_LAG improvement: %.6f", nodom_row$MeanDeltaStableShareVsBaseline),
  sprintf("08h reference status: %s", ref_match_status),
  "",
  "Leading 08i attributions:",
  fmt_best(best_country_group, "Best country own/cross group"),
  fmt_best(best_eq, "Best country-response equation row"),
  fmt_best(best_pred, "Best country-predictor column"),
  fmt_best(best_coef, "Best target individual domestic coefficient"),
  "",
  "Interpretation:",
  "- ZERO_ALL_OWN_LAGS vs ZERO_ALL_CROSS_LAGS separates domestic persistence from domestic cross-variable feedback.",
  "- COUNTRY_GROUP identifies whether a country's own-lags or cross-lags matter more.",
  "- EQUATION_ROW identifies the response equation that contributes most.",
  "- PREDICTOR_COLUMN identifies which lagged domestic variable contributes most.",
  "- INDIVIDUAL_COEFFICIENT drills to one A_i[response,predictor] entry for ZA,BR,JP by default.",
  "- A large improvement is attribution evidence, not permission to delete that coefficient from the publication model.",
  "- Do not lower the formal R62 >=90% posterior stability requirement.",
  "",
  "Outputs:",
  "- 00_domestic_lag_source_gate.csv",
  "- 00_parts_index.csv",
  "- 00_scenario_contract.csv",
  "- 01_scenario_stability_by_anchor.csv",
  "- 02_scenario_overall_summary.csv",
  "- 03_global_own_vs_cross_summary.csv",
  "- 04_country_own_cross_attribution.csv",
  "- 05_equation_row_attribution.csv",
  "- 06_predictor_column_attribution.csv",
  "- 07_target_individual_coefficient_attribution.csv",
  "- 08_domestic_coefficient_posterior_by_anchor.csv",
  "- 09_domestic_coefficient_posterior_overall.csv",
  "- 10_diagnostic_contract.csv"
)

writeLines(readme, file.path(OUT, "README_08i_domestic_lag_source.txt"))
cat(paste(readme, collapse = "\n"), "\n")
