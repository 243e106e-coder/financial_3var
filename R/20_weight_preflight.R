#!/usr/bin/env Rscript

source("R/00_config.R")

OUT <- file.path(RESULTS_DIR, "weights")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

matrices <- lapply(WEIGHT_FILES, read_weight_matrix, normalize=FALSE)

validation <- do.call(rbind, lapply(names(matrices), function(nm) {
  W <- matrices[[nm]]
  data.frame(
    Network=nm,
    Rows=nrow(W), Columns=ncol(W),
    Nonfinite=sum(!is.finite(W)),
    Negative=sum(W < -1e-12, na.rm=TRUE),
    MaxAbsDiagonal=max(abs(diag(W))),
    MaxAbsRowSumMinus1=max(abs(rowSums(W)-1)),
    MinWeight=min(W), MaxWeight=max(W),
    Status=if (all(is.finite(W)) && all(W >= -1e-12) &&
              max(abs(diag(W))) <= 1e-8 &&
              max(abs(rowSums(W)-1)) <= 1e-8) "OK" else "CHECK",
    stringsAsFactors=FALSE
  )
}))
write.csv(validation, file.path(OUT,"01_weight_validation.csv"), row.names=FALSE)

if (any(validation$Status != "OK")) stopf("At least one committed financial matrix failed validation.")

top_rows <- list(); k <- 0L
for (nm in names(matrices)) {
  W <- matrices[[nm]]
  for (i in COUNTRIES) {
    z <- sort(W[i,setdiff(COUNTRIES,i)], decreasing=TRUE)
    for (rank in seq_len(min(5L,length(z)))) {
      k <- k+1L
      top_rows[[k]] <- data.frame(
        Network=nm, Country=i, Rank=rank,
        Counterpart=names(z)[rank], Weight=as.numeric(z[rank])
      )
    }
  }
}
write.csv(do.call(rbind,top_rows), file.path(OUT,"02_top_counterparts.csv"), row.names=FALSE)

concentration <- do.call(rbind, lapply(names(matrices), function(nm) {
  W <- matrices[[nm]]
  data.frame(
    Network=nm, Country=COUNTRIES,
    HHI=rowSums(W^2),
    EffectivePartners=1/rowSums(W^2),
    MaxWeight=apply(W,1,max),
    stringsAsFactors=FALSE
  )
}))
write.csv(concentration, file.path(OUT,"03_weight_concentration.csv"), row.names=FALSE)

pairs <- combn(names(matrices), 2, simplify=FALSE)
distance <- do.call(rbind, lapply(pairs, function(p) {
  A <- matrices[[p[1]]]; B <- matrices[[p[2]]]
  D <- A-B
  data.frame(
    NetworkA=p[1], NetworkB=p[2],
    MaxAbsDifference=max(abs(D)),
    MeanAbsDifference=mean(abs(D)),
    FrobeniusNorm=sqrt(sum(D^2)),
    LargestDifferenceFrom=rownames(D)[which(abs(D)==max(abs(D)),arr.ind=TRUE)[1,1]],
    LargestDifferenceTo=colnames(D)[which(abs(D)==max(abs(D)),arr.ind=TRUE)[1,2]],
    stringsAsFactors=FALSE
  )
}))
write.csv(distance, file.path(OUT,"04_network_distance.csv"), row.names=FALSE)

txt <- c(
  "FINANCIAL WEIGHT PREFLIGHT: PASS",
  "================================",
  sprintf("Networks checked: %s", paste(names(matrices),collapse=", ")),
  sprintf("Main network: %s", MAIN_NETWORK),
  sprintf("Maximum committed row-sum error: %.12g", max(validation$MaxAbsRowSumMinus1)),
  sprintf("Largest main-vs-residency cell difference: %.12g",
          distance$MaxAbsDifference[distance$NetworkA==MAIN_NETWORK & distance$NetworkB=="residency_2017"]),
  "No missing, negative, non-normalized, or non-zero-diagonal weights detected."
)
writeLines(txt, file.path(OUT,"README_weight_preflight.txt"))
cat(paste(txt,collapse="\n"),"\n")
