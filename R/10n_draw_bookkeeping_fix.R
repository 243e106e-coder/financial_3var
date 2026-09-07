# =============================================================================
# 10n draw-bookkeeping fix
# Replace the block in R/69_ttvp_recovery_chain.R beginning at:
#     post <- fit$posterior
# and ending after:
#     svp <- post$svparms[keep_draw,,drop=FALSE]
#
# IMPORTANT:
# - MODE == "static": raw draw 1 is the static-block burn-boundary/storage row.
#   Discard it, then retain exactly STORED strictly post-burn draws.
# - MODE == "unrestricted": fit$posterior is already post-burn storage.
#   Do NOT require or discard an extra boundary row.
# - This changes storage bookkeeping only. It does NOT change the statistical
#   model, priors, burn-in, internal post-burn iterations, TVS, SV, or gate.
# =============================================================================

post <- fit$posterior

if (is.null(post$A) || length(dim(post$A)) != 3L) {
  stopf("Unexpected A structure.")
}

nstored_raw <- dim(post$A)[1]

# Return the first dimension / row count for posterior objects that are expected
# to share the same draw index.
posterior_draw_count <- function(x, name) {
  if (is.null(x)) stopf("Missing posterior object: %s", name)

  dx <- dim(x)
  if (is.null(dx)) {
    stopf("Posterior object %s has no draw dimension.", name)
  }

  n <- dx[1]
  if (!is.finite(n) || n < 1L) {
    stopf("Invalid raw draw count for %s: %s", name, as.character(n))
  }
  as.integer(n)
}

raw_counts <- c(
  A          = posterior_draw_count(post$A,          "A"),
  D_dyn      = posterior_draw_count(post$D_dyn,      "D_dyn"),
  thresholds = posterior_draw_count(post$thresholds, "thresholds"),
  omega      = posterior_draw_count(post$omega,      "omega"),
  V0         = posterior_draw_count(post$V0,         "V0"),
  H          = posterior_draw_count(post$H,          "H"),
  svparms    = posterior_draw_count(post$svparms,    "svparms")
)

if (length(unique(raw_counts)) != 1L) {
  stopf(
    "Raw posterior draw-count mismatch: %s",
    paste(sprintf("%s=%d", names(raw_counts), raw_counts), collapse = ", ")
  )
}

if (raw_counts[["A"]] != nstored_raw) {
  stopf("Internal raw-draw bookkeeping mismatch for A.")
}

# -----------------------------------------------------------------------------
# MODE-SPECIFIC storage contract
# -----------------------------------------------------------------------------
# static:
#   tested static-block extension can expose one leading burn-boundary/storage
#   row.  That row is not a strictly post-burn posterior draw and must be
#   discarded.
#
# unrestricted:
#   original pinned threshtvp branch returns posterior storage directly.  If it
#   returns exactly STORED rows, all STORED rows are valid posterior draws.
#   There is no additional boundary row to remove.
# -----------------------------------------------------------------------------

if (MODE == "static") {
  if (nstored_raw < STORED + 1L) {
    stopf(
      paste0(
        "Too few raw draws in static mode: need at least %d ",
        "(1 boundary + %d post-burn), found %d."
      ),
      STORED + 1L, STORED, nstored_raw
    )
  }

  candidates <- seq.int(2L, nstored_raw)
  boundary_dropped <- 1L

} else if (MODE == "unrestricted") {
  if (nstored_raw < STORED) {
    stopf(
      paste0(
        "Too few raw draws in unrestricted mode: need at least %d ",
        "post-burn draws, found %d."
      ),
      STORED, nstored_raw
    )
  }

  candidates <- seq_len(nstored_raw)
  boundary_dropped <- 0L

} else {
  stopf("Unknown sampler MODE in draw bookkeeping: %s", MODE)
}

n_candidates <- length(candidates)
if (n_candidates < STORED) {
  stopf(
    "Too few eligible posterior candidates in %s mode: need %d, found %d.",
    MODE, STORED, n_candidates
  )
}

# Select exactly STORED draws.  If the sampler exposes more eligible storage
# rows than requested (e.g. due to floating-point construction of the thinning
# grid), choose a deterministic approximately-even subset.  This does not alter
# the MCMC target.
if (n_candidates == STORED) {
  select_pos <- seq_len(STORED)
} else {
  select_pos <- as.integer(round(
    seq(1, n_candidates, length.out = STORED)
  ))
}

if (length(select_pos) != STORED ||
    length(unique(select_pos)) != STORED ||
    min(select_pos) < 1L ||
    max(select_pos) > n_candidates) {
  stopf(
    paste0(
      "Deterministic posterior selection failed in %s mode: ",
      "candidates=%d selected=%d unique=%d."
    ),
    MODE, n_candidates, length(select_pos), length(unique(select_pos))
  )
}

keep_draw <- candidates[select_pos]

if (length(keep_draw) != STORED || any(diff(keep_draw) <= 0L)) {
  stopf("Posterior draw contract failed after selection in %s mode.", MODE)
}

if (MODE == "static" && any(keep_draw <= 1L)) {
  stopf("Static-mode posterior selection retained the boundary draw.")
}
if (MODE == "unrestricted" && any(keep_draw < 1L)) {
  stopf("Unrestricted-mode posterior selection produced an invalid index.")
}

msg(
  paste0(
    "10n draw bookkeeping: mode=%s; raw=%d; boundary_dropped=%d; ",
    "eligible_candidates=%d; selected=%d; first_raw_index=%d; ",
    "last_raw_index=%d"
  ),
  MODE,
  nstored_raw,
  boundary_dropped,
  n_candidates,
  length(keep_draw),
  min(keep_draw),
  max(keep_draw)
)

Astd       <- post$A[keep_draw,,,drop=FALSE]
Ddyn       <- post$D_dyn[keep_draw,,,drop=FALSE]
thresholds <- post$thresholds[keep_draw,,drop=FALSE]
omega      <- post$omega[keep_draw,,drop=FALSE]
V0         <- post$V0[keep_draw,,drop=FALSE]
H          <- post$H[keep_draw,,drop=FALSE]
svp        <- post$svparms[keep_draw,,drop=FALSE]

# Hard machine check: every posterior object used downstream must contain
# exactly STORED draws after mode-specific selection.
selected_counts <- c(
  A          = dim(Astd)[1],
  D_dyn      = dim(Ddyn)[1],
  thresholds = nrow(thresholds),
  omega      = nrow(omega),
  V0         = nrow(V0),
  H          = nrow(H),
  svparms    = nrow(svp)
)

if (any(selected_counts != STORED)) {
  stopf(
    "Selected posterior draw-count mismatch: %s; expected each=%d.",
    paste(
      sprintf("%s=%d", names(selected_counts), selected_counts),
      collapse = ", "
    ),
    STORED
  )
}

msg(
  "10n posterior draw contract PASS: mode=%s; all 7 posterior objects=%d draws.",
  MODE, STORED
)
