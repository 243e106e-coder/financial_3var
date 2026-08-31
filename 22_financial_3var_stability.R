#!/usr/bin/env Rscript

# ============================================================
# 22_financial_3var_stability.R
# First-stage stability test for the 3-variable financial GVAR.
#
# Specification held deliberately simple and comparable:
#   x_i,t = a_i + A_i x_i,t-1
#                 + B_i0 x*_i,t + B_i1 x*_i,t-1
#                 + global controls (GPR, Brent; current + lag 1)
#                 + u_i,t
#   x_i = (r, de, deq)'
#
# p = 1, q = 1.
# GPR/Brent are treated as global controls for this stability diagnostic.
# The endogenous global system contains 14 * 3 = 42 country variables.
#
# It estimates the same 3-variable model twice when possible:
#   (1) GCAP financial W  [main diagnostic]
#   (2) legacy trade W    [network-only comparison]
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("r", "de", "deq")
K <- length(VARS)
N <- length(COUNTRIES)
P <- 1L
Q <- 1L

PANEL_PATH <- "financial_3var/data/panel_fin3_long.csv"
FIN_W_PATH <- "financial_3var/data/financial_weights.csv"
OUT_DIR <- "financial_3var/results/stability"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
lag1 <- function(x) c(NA_real_, x[-length(x)])

read_weight_matrix <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(d) <- trimws(names(d))
  rows <- toupper(trimws(as.character(d[[1]])))
  cols <- toupper(trimws(names(d)[-1]))
  if (!all(COUNTRIES %in% rows) || !all(COUNTRIES %in% cols)) stopf("Incomplete 14-country matrix: %s", path)
  W <- matrix(NA_real_, N, N, dimnames = list(COUNTRIES, COUNTRIES))
  for (i in COUNTRIES) {
    rr <- match(i, rows)
    for (j in COUNTRIES) W[i,j] <- num(d[[match(j, cols) + 1L]][rr])
  }
  if (any(!is.finite(W)) || any(W < -1e-12)) stopf("Invalid weights in %s", path)
  diag(W) <- 0
  rs <- rowSums(W)
  if (any(rs <= 0)) stopf("Non-positive row sum in %s", path)
  W / rs
}

if (!file.exists(PANEL_PATH)) stopf("Run 20_build_financial_3var_input.R first")
panel <- read.csv(PANEL_PATH, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("Quarter", "Country", VARS, paste0(VARS, "_star"), "gpr", "brent")
if (!all(need %in% names(panel))) stopf("Panel input is missing required columns")

# Keep one global sequence and one domestic sequence per country.
quarters <- unique(panel$Quarter)
if (length(quarters) < 40L) stopf("Sample is unexpectedly short")
if (any(table(panel$Country) != length(quarters))) stopf("Panel is not balanced across countries")

FIN_W <- read_weight_matrix(FIN_W_PATH)

trade_candidates <- c(
  "8.12/Trade_Weights_14_Economies_2000_2014.csv",
  ".github/workflows/Trade_Weights_14_Economies_2000_2014.csv",
  "Trade_Weights_14_Economies_2000_2014.csv"
)
trade_hit <- trade_candidates[file.exists(trade_candidates)]
TRADE_W <- if (length(trade_hit)) read_weight_matrix(trade_hit[[1]]) else NULL

# Build domestic arrays from the long panel. Stars are recomputed for each W
# so the trade-vs-financial comparison changes only the network.
X <- array(NA_real_, dim = c(length(quarters), N, K), dimnames = list(quarters, COUNTRIES, VARS))
for (i in seq_along(COUNTRIES)) {
  cc <- COUNTRIES[i]
  z <- panel[panel$Country == cc, , drop = FALSE]
  z <- z[match(quarters, z$Quarter), , drop = FALSE]
  X[, i, ] <- as.matrix(z[, VARS, drop = FALSE])
}
base <- panel[panel$Country == COUNTRIES[1], , drop = FALSE]
base <- base[match(quarters, base$Quarter), , drop = FALSE]
gpr <- num(base$gpr)
brent <- num(base$brent)
include_brent <- all(is.finite(brent))
if (!all(is.finite(gpr))) stopf("GPR contains missing values")

make_stars <- function(W) {
  Xs <- array(NA_real_, dim = dim(X), dimnames = dimnames(X))
  for (i in seq_along(COUNTRIES)) {
    for (v in seq_len(K)) Xs[, i, v] <- as.numeric(X[, , v] %*% W[i, ])
  }
  Xs
}

fit_country <- function(i, Xs) {
  Y <- X[, i, , drop = FALSE][,1,]
  S <- Xs[, i, , drop = FALSE][,1,]
  rows <- 2:nrow(Y)

  D <- data.frame(const = rep(1, length(rows)))
  for (v in seq_len(K)) D[[paste0(VARS[v], "_L1")]] <- lag1(Y[,v])[rows]
  for (v in seq_len(K)) D[[paste0(VARS[v], "_star_0")]] <- S[rows,v]
  for (v in seq_len(K)) D[[paste0(VARS[v], "_star_L1")]] <- lag1(S[,v])[rows]
  D$gpr_0 <- gpr[rows]
  D$gpr_L1 <- lag1(gpr)[rows]
  if (include_brent) {
    D$brent_0 <- brent[rows]
    D$brent_L1 <- lag1(brent)[rows]
  }

  YY <- Y[rows,,drop=FALSE]
  ok <- complete.cases(D) & complete.cases(YY)
  D <- D[ok,,drop=FALSE]
  YY <- YY[ok,,drop=FALSE]
  Xm <- as.matrix(D)
  if (nrow(Xm) <= ncol(Xm) + 5L) stopf("Too few effective observations for %s", COUNTRIES[i])

  B <- matrix(NA_real_, nrow = ncol(Xm), ncol = K, dimnames = list(colnames(Xm), VARS))
  residuals <- matrix(NA_real_, nrow = nrow(YY), ncol = K, dimnames = list(NULL, VARS))
  for (eq in seq_len(K)) {
    fit <- lm.fit(Xm, YY[,eq])
    if (any(!is.finite(fit$coefficients))) stopf("Rank deficiency for %s equation %s", COUNTRIES[i], VARS[eq])
    B[,eq] <- fit$coefficients
    residuals[,eq] <- fit$residuals
  }

  A <- matrix(0, K, K, dimnames = list(VARS, VARS))
  B0 <- matrix(0, K, K, dimnames = list(VARS, VARS))
  B1 <- matrix(0, K, K, dimnames = list(VARS, VARS))
  for (eq in VARS) {
    for (v in VARS) {
      A[eq,v]  <- B[paste0(v, "_L1"), eq]
      B0[eq,v] <- B[paste0(v, "_star_0"), eq]
      B1[eq,v] <- B[paste0(v, "_star_L1"), eq]
    }
  }

  evA <- eigen(A, only.values = TRUE)$values
  list(
    country = COUNTRIES[i], A = A, B0 = B0, B1 = B1, B = B,
    n = nrow(YY), domestic_rho = max(Mod(evA)),
    design_kappa = kappa(crossprod(Xm)), design_rcond = rcond(crossprod(Xm)),
    sigma = crossprod(residuals) / nrow(residuals)
  )
}

selection_matrix <- function(i) {
  S <- matrix(0, K, N*K)
  cols <- ((i-1L)*K+1L):(i*K)
  S[, cols] <- diag(K)
  S
}

star_matrix <- function(i, W) {
  R <- matrix(0, K, N*K)
  for (j in seq_len(N)) {
    for (v in seq_len(K)) R[v, (j-1L)*K + v] <- W[i,j]
  }
  R
}

fit_global <- function(W, label) {
  Xs <- make_stars(W)
  fits <- lapply(seq_len(N), fit_country, Xs = Xs)

  G0 <- matrix(0, N*K, N*K)
  G1 <- matrix(0, N*K, N*K)
  for (i in seq_len(N)) {
    rr <- ((i-1L)*K+1L):(i*K)
    Si <- selection_matrix(i)
    Ri <- star_matrix(i, W)
    G0[rr,] <- Si - fits[[i]]$B0 %*% Ri
    G1[rr,] <- fits[[i]]$A %*% Si + fits[[i]]$B1 %*% Ri
  }

  rc <- rcond(G0)
  kap <- kappa(G0)
  F <- tryCatch(solve(G0, G1), error = function(e) NULL)
  if (is.null(F)) {
    ev <- rep(NA_complex_, N*K)
    rho <- NA_real_
  } else {
    ev <- eigen(F, only.values = TRUE)$values
    rho <- max(Mod(ev))
  }

  local <- data.frame(
    Network = label,
    Country = COUNTRIES,
    N = vapply(fits, `[[`, numeric(1), "n"),
    DomesticSpectralRadius = vapply(fits, `[[`, numeric(1), "domestic_rho"),
    DesignKappa = vapply(fits, `[[`, numeric(1), "design_kappa"),
    DesignRcond = vapply(fits, `[[`, numeric(1), "design_rcond"),
    DomesticStable = vapply(fits, function(z) z$domestic_rho < 1, logical(1)),
    stringsAsFactors = FALSE
  )

  list(
    label = label, fits = fits, G0 = G0, G1 = G1, F = F, eigen = ev,
    rho = rho, rcond = rc, kappa = kap, local = local
  )
}

fin <- fit_global(FIN_W, "GCAP_financial_2017")
trade <- if (!is.null(TRADE_W)) fit_global(TRADE_W, "legacy_trade") else NULL

summary_rows <- list(data.frame(
  Network = fin$label,
  GlobalSpectralRadius = fin$rho,
  Stable = is.finite(fin$rho) && fin$rho < 1,
  G0_rcond = fin$rcond,
  G0_kappa = fin$kappa,
  MaxLocalDomesticRadius = max(fin$local$DomesticSpectralRadius),
  UnstableLocalCount = sum(!fin$local$DomesticStable),
  stringsAsFactors = FALSE
))
if (!is.null(trade)) summary_rows[[2]] <- data.frame(
  Network = trade$label,
  GlobalSpectralRadius = trade$rho,
  Stable = is.finite(trade$rho) && trade$rho < 1,
  G0_rcond = trade$rcond,
  G0_kappa = trade$kappa,
  MaxLocalDomesticRadius = max(trade$local$DomesticSpectralRadius),
  UnstableLocalCount = sum(!trade$local$DomesticStable),
  stringsAsFactors = FALSE
)
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "01_global_stability_comparison.csv"), row.names = FALSE)
write.csv(fin$local, file.path(OUT_DIR, "02_local_financial_stability.csv"), row.names = FALSE)
if (!is.null(trade)) write.csv(trade$local, file.path(OUT_DIR, "03_local_trade_stability.csv"), row.names = FALSE)

write_eigs <- function(obj, path) {
  ev <- obj$eigen
  d <- data.frame(
    Rank = seq_along(ev),
    Real = Re(ev),
    Imag = Im(ev),
    Modulus = Mod(ev),
    stringsAsFactors = FALSE
  )
  d <- d[order(d$Modulus, decreasing = TRUE), , drop = FALSE]
  d$Rank <- seq_len(nrow(d))
  write.csv(d, path, row.names = FALSE)
}
write_eigs(fin, file.path(OUT_DIR, "04_financial_global_eigenvalues.csv"))
if (!is.null(trade)) write_eigs(trade, file.path(OUT_DIR, "05_trade_global_eigenvalues.csv"))

spec <- data.frame(
  Item = c("Countries", "Domestic variables", "Domestic lag p", "Foreign lag q", "Contemporaneous foreign variables", "Global GPR", "Brent", "Sample", "Financial W"),
  Value = c(
    paste(COUNTRIES, collapse = ", "), paste(VARS, collapse = ", "), P, Q, "YES",
    "current + lag 1", ifelse(include_brent, "current + lag 1", "not available"),
    paste0(quarters[1], " - ", tail(quarters, 1)), "GCAP 2017 Position Residency"
  ),
  stringsAsFactors = FALSE
)
write.csv(spec, file.path(OUT_DIR, "06_model_specification.csv"), row.names = FALSE)

readme <- c(
  "3-variable financial GVAR stability diagnostic",
  paste0("Sample: ", quarters[1], " - ", tail(quarters, 1)),
  "Variables: r, de, deq",
  "Lag structure: p=1, q=1; contemporaneous foreign variables included",
  paste0("Brent included: ", ifelse(include_brent, "YES", "NO")),
  paste0("Financial global spectral radius: ", signif(fin$rho, 8)),
  paste0("Financial stable (<1): ", is.finite(fin$rho) && fin$rho < 1),
  paste0("Financial G0 rcond: ", signif(fin$rcond, 8)),
  paste0("Financial G0 kappa: ", signif(fin$kappa, 8)),
  if (!is.null(trade)) paste0("Trade global spectral radius under the SAME 3-variable specification: ", signif(trade$rho, 8)) else "Trade comparison: unavailable",
  "This is a diagnostic OLS GVAR, not the final Bayesian TVP estimation.",
  "Do not proceed to the expensive TVP stage until this diagnostic has been reviewed."
)
writeLines(readme, file.path(OUT_DIR, "README_stability.txt"))
cat(paste(readme, collapse = "\n"), "\n")
