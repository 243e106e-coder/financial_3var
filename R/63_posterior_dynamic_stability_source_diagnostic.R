#!/usr/bin/env Rscript

# =============================================================================
# 63_posterior_dynamic_stability_source_diagnostic.R
#
# Diagnostic-only decomposition of posterior GVAR dynamic instability.
#
# IMPORTANT:
# - Does NOT re-estimate the TVP-GVAR.
# - Does NOT relax the formal R62 stability gate.
# - Uses the same posterior event coefficients and the same G0/G1 definitions
#   as R/62_structural_gpr_tvp_irf.R.
# - Counterfactual scenarios are diagnostics, not alternative publication
#   specifications unless separately estimated and validated.
#
# Baseline at each event anchor t and posterior draw d:
#   G0 = I - contemporaneous foreign block
#   G1 = domestic lag block + foreign lag block
#   F  = solve(G0) %*% G1
#
# Diagnostics:
#   1) BASELINE
#   2) NO_FOREIGN_LAG       : B1_i = 0
#   3) NO_FOREIGN_CONTEMP   : B0_i = 0
#   4) NO_DOMESTIC_LAG      : A_i  = 0
#   5) DOMESTIC_ONLY        : B0_i = B1_i = 0
#   6) ZERO_DYNAMIC_ROW_<CC>: set country CC's G1 row block to zero while
#                             retaining full G0. Pure attribution diagnostic.
#   7) LOO_<CC>             : remove CC from the system, subset and renormalize
#                             W among remaining countries, without re-estimation.
#
# Default posterior thinning is deterministic and balanced across all 4 chains.
# This script is for source diagnosis, not posterior inference.
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
NCHAINS <- as.integer(get_env_num("FIN3_NCHAINS", 4))
EXPECTED_STORED <- as.integer(get_env_num("FIN3_STORED_PER_CHAIN", 2000))
DIAG_DRAWS_PER_CHAIN <- as.integer(get_env_num("FIN3_DIAG_DRAWS_PER_CHAIN", 250))
MIN_G0_RCOND <- get_env_num("FIN3_IRF_MIN_G0_RCOND", 1e-10)
DIAG_CORES <- as.integer(get_env_num("FIN3_DIAG_CORES", 2))
LOO_COUNTRIES_RAW <- get_env_chr("FIN3_LOO_COUNTRIES", "ZA,JP,BR,US,SG")

if (NCHAINS != 4L) stopf("08h expects the formal 4-chain posterior.")
if (EXPECTED_STORED < 200L) stopf("Expected stored draws are implausibly small.")
if (DIAG_DRAWS_PER_CHAIN < 50L) stopf("Use at least 50 diagnostic draws per chain.")
if (DIAG_DRAWS_PER_CHAIN > EXPECTED_STORED) {
  DIAG_DRAWS_PER_CHAIN <- EXPECTED_STORED
}
if (DIAG_CORES < 1L) DIAG_CORES <- 1L

LOO_COUNTRIES <- unique(toupper(trimws(strsplit(LOO_COUNTRIES_RAW, ",", fixed = TRUE)[[1]])))
LOO_COUNTRIES <- LOO_COUNTRIES[nzchar(LOO_COUNTRIES)]
bad_loo <- setdiff(LOO_COUNTRIES, COUNTRIES)
if (length(bad_loo)) stopf("Unknown FIN3_LOO_COUNTRIES: %s", paste(bad_loo, collapse = ", "))

CONV_GATE <- file.path(CONV_DIR, "00_formal_mcmc_gate.csv")
if (!file.exists(CONV_GATE)) stopf("Formal MCMC convergence gate missing: %s", CONV_GATE)
cg <- read.csv(CONV_GATE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(cg) != 1L || !"Status" %in% names(cg)) stopf("Malformed formal MCMC gate.")
if (!identical(trimws(cg$Status[1]), "READY_FOR_FORMAL_IRF")) {
  stopf("Formal MCMC gate is not READY_FOR_FORMAL_IRF.")
}

OUT <- file.path(RESULTS_DIR, "dynamic_stability_source")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load the exact 14-country x 4-chain posterior grid
# =============================================================================

files <- list.files(
  PARTS_ROOT,
  pattern = "^formal_tvp_[A-Z]{2}_chain[0-9]+\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) != length(COUNTRIES) * NCHAINS) {
  stopf(
    "Expected exactly %d posterior RDS files; found %d.",
    length(COUNTRIES) * NCHAINS, length(files)
  )
}

parts <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)
parts_index <- vector("list", length(COUNTRIES) * NCHAINS)
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
    if (!identical(as.character(z$meta$country), cc)) {
      stopf("Country metadata mismatch in %s.", candidates)
    }
    if (as.integer(z$meta$chain) != ch) {
      stopf("Chain metadata mismatch in %s.", candidates)
    }
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
  stopf("Posterior parts do not contain the coefficient blocks required by R62.")
}

anchor_key <- paste(anchors$EventID, anchors$AnchorType, anchors$AnchorQuarter)
for (cc in COUNTRIES) {
  for (ch in seq_len(NCHAINS)) {
    z <- parts[[cc]][[ch]]
    if (!identical(z$meta$key_terms, KEY_TERMS)) {
      stopf("Key-term ordering mismatch: %s chain %d.", cc, ch)
    }
    if (!identical(
      paste(z$anchors$EventID, z$anchors$AnchorType, z$anchors$AnchorQuarter),
      anchor_key
    )) {
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
# 2. Matrix definitions exactly aligned with R62
# =============================================================================

N <- length(COUNTRIES)
K <- length(VARS)
NK <- N * K
W <- read_weight_matrix(WEIGHT_FILES[[MAIN_NETWORK]])

term_dom <- paste0(VARS, "_L1")
term_star0 <- paste0(VARS, "_star_0")
term_star1 <- paste0(VARS, "_star_L1")

idx_dom <- match(term_dom, KEY_TERMS)
idx_star0 <- match(term_star0, KEY_TERMS)
idx_star1 <- match(term_star1, KEY_TERMS)
if (any(is.na(c(idx_dom, idx_star0, idx_star1)))) {
  stopf("Failed to map domestic/foreign lag blocks.")
}

make_structure <- function(orig_idx, W_sub) {
  n <- length(orig_idx)
  nk <- n * K

  selectors <- vector("list", n)
  star_maps <- vector("list", n)

  for (ii in seq_len(n)) {
    S <- matrix(0, K, nk)
    S[, ((ii - 1L) * K + 1L):(ii * K)] <- diag(K)
    selectors[[ii]] <- S

    R <- matrix(0, K, nk)
    for (jj in seq_len(n)) {
      for (v in seq_len(K)) {
        R[v, (jj - 1L) * K + v] <- W_sub[ii, jj]
      }
    }
    star_maps[[ii]] <- R
  }

  list(
    orig_idx = orig_idx,
    countries = COUNTRIES[orig_idx],
    W = W_sub,
    n = n,
    nk = nk,
    selectors = selectors,
    star_maps = star_maps
  )
}

FULL_STRUCT <- make_structure(seq_len(N), W)

LOO_STRUCTS <- setNames(vector("list", length(LOO_COUNTRIES)), LOO_COUNTRIES)
for (cc in LOO_COUNTRIES) {
  drop_i <- match(cc, COUNTRIES)
  keep <- setdiff(seq_len(N), drop_i)
  Wsub <- W[keep, keep, drop = FALSE]
  rs <- rowSums(Wsub)
  if (any(!is.finite(rs)) || any(rs <= 0)) {
    stopf("LOO_%s leaves a non-positive financial-weight row sum.", cc)
  }
  Wsub <- Wsub / rs
  diag(Wsub) <- 0
  LOO_STRUCTS[[cc]] <- make_structure(keep, Wsub)
}

extract_blocks <- function(ch, dd, a) {
  A <- vector("list", N)
  B0 <- vector("list", N)
  B1 <- vector("list", N)
  local_rho <- rep(NA_real_, N)

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

    local_rho[i] <- tryCatch(
      max(Mod(eigen(A[[i]], only.values = TRUE)$values)),
      error = function(e) NA_real_
    )
  }

  list(A = A, B0 = B0, B1 = B1, local_rho = local_rho)
}

build_components <- function(blocks, st) {
  C0 <- matrix(0, st$nk, st$nk)  # contemporaneous foreign contribution
  D  <- matrix(0, st$nk, st$nk)  # domestic lag contribution
  L  <- matrix(0, st$nk, st$nk)  # foreign lag contribution

  for (ii in seq_len(st$n)) {
    oi <- st$orig_idx[ii]
    rr <- ((ii - 1L) * K + 1L):(ii * K)
    S <- st$selectors[[ii]]
    R <- st$star_maps[[ii]]

    C0[rr, ] <- blocks$B0[[oi]] %*% R
    D[rr, ]  <- blocks$A[[oi]] %*% S
    L[rr, ]  <- blocks$B1[[oi]] %*% R
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

block_scenarios <- data.frame(
  Scenario = c(
    "BASELINE",
    "NO_FOREIGN_LAG",
    "NO_FOREIGN_CONTEMP",
    "NO_DOMESTIC_LAG",
    "DOMESTIC_ONLY"
  ),
  ScenarioFamily = "BLOCK",
  TargetCountry = NA_character_,
  Definition = c(
    "R62 baseline: G0=I-C0; G1=D+L",
    "Counterfactual: foreign lag block L=0; retain baseline G0",
    "Counterfactual: contemporaneous foreign block C0=0; retain baseline G1",
    "Counterfactual: domestic lag block D=0; retain baseline G0",
    "Counterfactual: C0=0 and L=0, so F=D"
  ),
  stringsAsFactors = FALSE
)

row_scenarios <- data.frame(
  Scenario = paste0("ZERO_DYNAMIC_ROW_", COUNTRIES),
  ScenarioFamily = "COUNTRY_DYNAMIC_ROW",
  TargetCountry = COUNTRIES,
  Definition = paste0(
    "Attribution only: zero ", COUNTRIES,
    "'s rows in full-system G1 while retaining baseline G0"
  ),
  stringsAsFactors = FALSE
)

loo_scenarios <- data.frame(
  Scenario = paste0("LOO_", LOO_COUNTRIES),
  ScenarioFamily = "LEAVE_ONE_COUNTRY_OUT",
  TargetCountry = LOO_COUNTRIES,
  Definition = paste0(
    "Reduced-system counterfactual: remove ", LOO_COUNTRIES,
    ", subset and renormalize W; no coefficient re-estimation"
  ),
  stringsAsFactors = FALSE
)

scenario_contract <- rbind(block_scenarios, row_scenarios, loo_scenarios)
write.csv(
  scenario_contract,
  file.path(OUT, "00_scenario_contract.csv"),
  row.names = FALSE
)

SCENARIOS <- scenario_contract$Scenario
NS <- length(SCENARIOS)
scenario_pos <- setNames(seq_len(NS), SCENARIOS)

# =============================================================================
# 4. Process one anchor
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
  rc_mat <- matrix(NA_real_, ND_TOTAL, NS, dimnames = list(NULL, SCENARIOS))
  local_rho_mat <- matrix(NA_real_, ND_TOTAL, N, dimnames = list(NULL, COUNTRIES))
  chain_vec <- integer(ND_TOTAL)
  stored_draw_vec <- integer(ND_TOTAL)

  gd <- 0L

  for (ch in seq_len(NCHAINS)) {
    for (dd in draw_idx) {
      gd <- gd + 1L
      chain_vec[gd] <- ch
      stored_draw_vec[gd] <- dd

      blocks <- extract_blocks(ch, dd, a)
      local_rho_mat[gd, ] <- blocks$local_rho

      comp <- build_components(blocks, FULL_STRUCT)
      I_full <- diag(NK)

      # Baseline G0 is shared by BASELINE, NO_FOREIGN_LAG,
      # NO_DOMESTIC_LAG, and all ZERO_DYNAMIC_ROW scenarios.
      base_imp <- safe_impact(I_full - comp$C0)
      rc_base <- base_imp$rcond

      # BASELINE
      s <- scenario_pos[["BASELINE"]]
      rc_mat[gd, s] <- rc_base
      rho_mat[gd, s] <- safe_rho(base_imp$impact, comp$D + comp$L)

      # NO_FOREIGN_LAG
      s <- scenario_pos[["NO_FOREIGN_LAG"]]
      rc_mat[gd, s] <- rc_base
      rho_mat[gd, s] <- safe_rho(base_imp$impact, comp$D)

      # NO_DOMESTIC_LAG
      s <- scenario_pos[["NO_DOMESTIC_LAG"]]
      rc_mat[gd, s] <- rc_base
      rho_mat[gd, s] <- safe_rho(base_imp$impact, comp$L)

      # NO_FOREIGN_CONTEMP => G0 = I
      s <- scenario_pos[["NO_FOREIGN_CONTEMP"]]
      rc_mat[gd, s] <- 1
      rho_mat[gd, s] <- safe_rho(I_full, comp$D + comp$L)

      # DOMESTIC_ONLY => G0 = I, G1 = D
      s <- scenario_pos[["DOMESTIC_ONLY"]]
      rc_mat[gd, s] <- 1
      rho_mat[gd, s] <- safe_rho(I_full, comp$D)

      # Country dynamic-row attribution. This is intentionally not a
      # re-estimated model; it isolates each country's row contribution to G1.
      G1_base <- comp$D + comp$L
      for (i in seq_len(N)) {
        sc <- paste0("ZERO_DYNAMIC_ROW_", COUNTRIES[i])
        s <- scenario_pos[[sc]]
        G1_cf <- G1_base
        rr <- ((i - 1L) * K + 1L):(i * K)
        G1_cf[rr, ] <- 0
        rc_mat[gd, s] <- rc_base
        rho_mat[gd, s] <- safe_rho(base_imp$impact, G1_cf)
      }

      # Reduced-system leave-one-country-out diagnostics.
      for (cc in LOO_COUNTRIES) {
        sc <- paste0("LOO_", cc)
        s <- scenario_pos[[sc]]
        st <- LOO_STRUCTS[[cc]]
        comp_sub <- build_components(blocks, st)
        I_sub <- diag(st$nk)
        imp_sub <- safe_impact(I_sub - comp_sub$C0)
        rc_mat[gd, s] <- imp_sub$rcond
        rho_mat[gd, s] <- safe_rho(imp_sub$impact, comp_sub$D + comp_sub$L)
      }
    }
  }

  # Scenario-level anchor summaries.
  scenario_rows <- vector("list", NS)
  for (s in seq_len(NS)) {
    sm <- summarize_rho(rho_mat[, s], rc_mat[, s])
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
        stringsAsFactors = FALSE
      ),
      sm
    )
  }
  scenario_anchor <- do.call(rbind, scenario_rows)

  bpos <- scenario_pos[["BASELINE"]]
  baseline_rho <- rho_mat[, bpos]
  baseline_rc <- rc_mat[, bpos]
  baseline_stable <- is.finite(baseline_rho) & baseline_rho < 1 &
    is.finite(baseline_rc) & baseline_rc >= MIN_G0_RCOND
  baseline_unstable <- is.finite(baseline_rho) & baseline_rho >= 1 &
    is.finite(baseline_rc) & baseline_rc >= MIN_G0_RCOND

  # Local-country diagnostics conditional on the same baseline global draws.
  local_rows <- vector("list", N)
  for (i in seq_len(N)) {
    lr <- local_rho_mat[, i]
    ok <- is.finite(lr) & is.finite(baseline_rho)

    spear <- NA_real_
    if (sum(ok) >= 10L && stats::sd(lr[ok]) > 0 && stats::sd(baseline_rho[ok]) > 0) {
      spear <- suppressWarnings(stats::cor(
        lr[ok], baseline_rho[ok],
        method = "spearman",
        use = "complete.obs"
      ))
    }

    local_rows[[i]] <- data.frame(
      EventID = meta$EventID,
      EventSet = meta$EventSet,
      AnchorType = meta$AnchorType,
      AnchorQuarter = meta$AnchorQuarter,
      Country = COUNTRIES[i],
      DiagnosticDraws = ND_TOTAL,
      LocalStableShare = if (any(is.finite(lr))) mean(lr[is.finite(lr)] < 1) else NA_real_,
      LocalRhoMedian = if (any(is.finite(lr))) median(lr[is.finite(lr)]) else NA_real_,
      LocalRhoP95 = if (any(is.finite(lr))) unname(quantile(lr[is.finite(lr)], .95)) else NA_real_,
      SpearmanLocalVsGlobalRho = spear,
      LocalUnstableGivenGlobalUnstable = if (sum(baseline_unstable) > 0) {
        mean(lr[baseline_unstable] >= 1, na.rm = TRUE)
      } else NA_real_,
      MeanLocalRho_GlobalStable = if (sum(baseline_stable) > 0) {
        mean(lr[baseline_stable], na.rm = TRUE)
      } else NA_real_,
      MeanLocalRho_GlobalUnstable = if (sum(baseline_unstable) > 0) {
        mean(lr[baseline_unstable], na.rm = TRUE)
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  list(
    scenario_anchor = scenario_anchor,
    local_driver = do.call(rbind, local_rows)
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

scenario_anchor <- do.call(rbind, lapply(ans, `[[`, "scenario_anchor"))
local_driver <- do.call(rbind, lapply(ans, `[[`, "local_driver"))

# Add within-anchor changes versus baseline.
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
# 5. Overall source-attribution summaries
# =============================================================================

overall_rows <- lapply(SCENARIOS, function(sc) {
  d <- scenario_anchor[scenario_anchor$Scenario == sc, , drop = FALSE]
  meta <- scenario_contract[scenario_contract$Scenario == sc, , drop = FALSE]
  data.frame(
    Scenario = sc,
    ScenarioFamily = meta$ScenarioFamily,
    TargetCountry = meta$TargetCountry,
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
})
overall <- do.call(rbind, overall_rows)
overall <- overall[order(-overall$MeanDeltaStableShareVsBaseline, overall$MeanRhoMedian), ]
write.csv(
  overall,
  file.path(OUT, "02_scenario_overall_summary.csv"),
  row.names = FALSE
)

block_summary <- overall[overall$ScenarioFamily == "BLOCK", , drop = FALSE]
write.csv(
  block_summary,
  file.path(OUT, "03_block_source_summary.csv"),
  row.names = FALSE
)

row_attr <- overall[overall$ScenarioFamily == "COUNTRY_DYNAMIC_ROW", , drop = FALSE]
row_attr <- row_attr[order(-row_attr$MeanDeltaStableShareVsBaseline), ]
write.csv(
  row_attr,
  file.path(OUT, "04_country_dynamic_row_attribution.csv"),
  row.names = FALSE
)

loo_summary <- overall[overall$ScenarioFamily == "LEAVE_ONE_COUNTRY_OUT", , drop = FALSE]
loo_summary <- loo_summary[order(-loo_summary$MeanDeltaStableShareVsBaseline), ]
write.csv(
  loo_summary,
  file.path(OUT, "05_leave_one_country_out_summary.csv"),
  row.names = FALSE
)

write.csv(
  local_driver,
  file.path(OUT, "06_baseline_local_driver_by_anchor_country.csv"),
  row.names = FALSE
)

local_overall <- do.call(rbind, lapply(COUNTRIES, function(cc) {
  d <- local_driver[local_driver$Country == cc, , drop = FALSE]
  data.frame(
    Country = cc,
    Anchors = nrow(d),
    MeanLocalStableShare = mean(d$LocalStableShare, na.rm = TRUE),
    MeanLocalRhoMedian = mean(d$LocalRhoMedian, na.rm = TRUE),
    MeanLocalRhoP95 = mean(d$LocalRhoP95, na.rm = TRUE),
    MeanSpearmanLocalVsGlobalRho = mean(d$SpearmanLocalVsGlobalRho, na.rm = TRUE),
    MeanLocalUnstableGivenGlobalUnstable = mean(d$LocalUnstableGivenGlobalUnstable, na.rm = TRUE),
    MeanLocalRho_GlobalStable = mean(d$MeanLocalRho_GlobalStable, na.rm = TRUE),
    MeanLocalRho_GlobalUnstable = mean(d$MeanLocalRho_GlobalUnstable, na.rm = TRUE),
    MeanLocalRhoGap_UnstableMinusStable = mean(
      d$MeanLocalRho_GlobalUnstable - d$MeanLocalRho_GlobalStable,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}))
local_overall <- local_overall[
  order(-local_overall$MeanSpearmanLocalVsGlobalRho,
        local_overall$MeanLocalStableShare),
]
write.csv(
  local_overall,
  file.path(OUT, "07_baseline_local_driver_overall.csv"),
  row.names = FALSE
)

# =============================================================================
# 6. Diagnostic gate and README
# =============================================================================

baseline_overall <- overall[overall$Scenario == "BASELINE", , drop = FALSE]
if (nrow(baseline_overall) != 1L) stopf("Baseline overall summary is malformed.")

best_block_nonbaseline <- block_summary[
  block_summary$Scenario != "BASELINE",
  ,
  drop = FALSE
]
best_block_nonbaseline <- best_block_nonbaseline[
  order(-best_block_nonbaseline$MeanDeltaStableShareVsBaseline),
  ,
  drop = FALSE
]

top_row <- if (nrow(row_attr)) row_attr[1, , drop = FALSE] else NULL
top_loo <- if (nrow(loo_summary)) loo_summary[1, , drop = FALSE] else NULL

integrity_ok <- all(scenario_anchor$FiniteRhoShare >= .99) &&
  all(scenario_anchor$G0OKShare >= .99) &&
  nrow(parts_index) == length(COUNTRIES) * NCHAINS

status <- if (integrity_ok) "DIAGNOSTIC_COMPLETE" else "DIAGNOSTIC_INTEGRITY_WARNING"

gate <- data.frame(
  Status = status,
  MainNetwork = MAIN_NETWORK,
  Countries = length(COUNTRIES),
  Chains = NCHAINS,
  StoredDrawsPerChain = EXPECTED_STORED,
  DiagnosticDrawsPerChain = ND_PER_CHAIN,
  DiagnosticDrawsTotal = ND_TOTAL,
  EventAnchors = nrow(anchors),
  Scenarios = NS,
  BaselineMeanStableShare = baseline_overall$MeanStableShare,
  BaselineMinStableShare = baseline_overall$MinStableShare,
  BestBlockScenario = if (nrow(best_block_nonbaseline)) best_block_nonbaseline$Scenario[1] else NA_character_,
  BestBlockDeltaStableShare = if (nrow(best_block_nonbaseline)) {
    best_block_nonbaseline$MeanDeltaStableShareVsBaseline[1]
  } else NA_real_,
  TopDynamicRowCountry = if (!is.null(top_row)) top_row$TargetCountry[1] else NA_character_,
  TopDynamicRowDeltaStableShare = if (!is.null(top_row)) {
    top_row$MeanDeltaStableShareVsBaseline[1]
  } else NA_real_,
  TopLOOCountry = if (!is.null(top_loo)) top_loo$TargetCountry[1] else NA_character_,
  TopLOODeltaStableShare = if (!is.null(top_loo)) {
    top_loo$MeanDeltaStableShareVsBaseline[1]
  } else NA_real_,
  FormalIRFGateWasRelaxed = FALSE,
  ReestimationPerformed = FALSE,
  stringsAsFactors = FALSE
)
write.csv(
  gate,
  file.path(OUT, "00_dynamic_stability_source_gate.csv"),
  row.names = FALSE
)

contract <- data.frame(
  Item = c(
    "Purpose",
    "Posterior",
    "MatrixDefinition",
    "Thinning",
    "FormalGate",
    "CountryRowCounterfactual",
    "LOOCounterfactual"
  ),
  Value = c(
    "Diagnose why posterior F=G0^{-1}G1 is unstable; not a replacement IRF",
    "Exact accepted 14-country x 4-chain posterior grid",
    "Same G0/G1 construction as R/62_structural_gpr_tvp_irf.R",
    sprintf("Deterministic balanced thinning: %d of %d stored draws per chain", ND_PER_CHAIN, EXPECTED_STORED),
    "Unchanged; no stable-share threshold is lowered",
    "Zeros one country's G1 row block while retaining full G0; attribution only",
    paste0("Removes country and renormalizes W among remaining economies; no re-estimation: ",
           paste(LOO_COUNTRIES, collapse = ", "))
  ),
  stringsAsFactors = FALSE
)
write.csv(
  contract,
  file.path(OUT, "08_diagnostic_contract.csv"),
  row.names = FALSE
)

readme <- c(
  sprintf("08h POSTERIOR DYNAMIC STABILITY SOURCE DIAGNOSTIC: %s", status),
  "=================================================================",
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf("Countries: %d", length(COUNTRIES)),
  sprintf("Chains: %d", NCHAINS),
  sprintf("Stored draws per chain: %d", EXPECTED_STORED),
  sprintf("Diagnostic draws per chain: %d", ND_PER_CHAIN),
  sprintf("Total diagnostic draws per anchor: %d", ND_TOTAL),
  sprintf("Event anchors: %d", nrow(anchors)),
  sprintf("Scenarios: %d", NS),
  sprintf("Baseline mean stable share: %.6f", gate$BaselineMeanStableShare),
  sprintf("Baseline minimum anchor stable share: %.6f", gate$BaselineMinStableShare),
  "",
  "Leading diagnostics (screening only):",
  sprintf(
    "- Best block counterfactual by stable-share improvement: %s (delta %.6f)",
    gate$BestBlockScenario, gate$BestBlockDeltaStableShare
  ),
  sprintf(
    "- Largest dynamic-row attribution: %s (delta %.6f)",
    gate$TopDynamicRowCountry, gate$TopDynamicRowDeltaStableShare
  ),
  sprintf(
    "- Largest leave-one-country-out improvement: %s (delta %.6f)",
    gate$TopLOOCountry, gate$TopLOODeltaStableShare
  ),
  "",
  "Interpretation rules:",
  "- BASELINE reproduces the R62 dynamic-stability object on a deterministic thinned posterior.",
  "- NO_FOREIGN_LAG isolates the role of foreign lag coefficients B1.",
  "- NO_FOREIGN_CONTEMP isolates the role of contemporaneous foreign coefficients B0.",
  "- NO_DOMESTIC_LAG isolates the role of domestic lag coefficients A.",
  "- DOMESTIC_ONLY removes both contemporaneous and lagged foreign propagation.",
  "- ZERO_DYNAMIC_ROW_* is a matrix attribution diagnostic, NOT a valid re-estimated model.",
  "- LOO_* is a reduced-system counterfactual with W renormalized, NOT a re-estimated model.",
  "- Do not use counterfactual IRFs for publication from this script.",
  "- Do not lower the formal R62 >=90% posterior stability requirement because of this diagnostic.",
  "",
  "Outputs:",
  "- 00_dynamic_stability_source_gate.csv",
  "- 00_parts_index.csv",
  "- 00_scenario_contract.csv",
  "- 01_scenario_stability_by_anchor.csv",
  "- 02_scenario_overall_summary.csv",
  "- 03_block_source_summary.csv",
  "- 04_country_dynamic_row_attribution.csv",
  "- 05_leave_one_country_out_summary.csv",
  "- 06_baseline_local_driver_by_anchor_country.csv",
  "- 07_baseline_local_driver_overall.csv",
  "- 08_diagnostic_contract.csv"
)
writeLines(readme, file.path(OUT, "README_08h_dynamic_stability_source.txt"))
cat(paste(readme, collapse = "\n"), "\n")
