# =====================================================================
# beta_diversity.R -- betalinkr pairwise interaction beta-diversity
# =====================================================================

#' Pairwise interaction beta-diversity across all time points
#'
#' Partitions `bipartite::betalinkr()` S / OS / WN / ST components for every
#' pair of time points.
#'
#' @param networks A named list of matrices, ordered by time point.
#' @param binary Passed to `betalinkr()` (TRUE = presence/absence).
#' @param partitioning Passed to `betalinkr()` (default `"commondenom"`).
#' @return A data frame, one row per pair of time points.
#' @export
pairwise_betadiversity <- function(networks, binary = TRUE,
                                   partitioning = "commondenom") {
  .bd_require("bipartite", "beta diversity")
  labs <- names(networks); idx <- utils::combn(seq_along(networks), 2)
  out <- data.frame()
  for (k in seq_len(ncol(idx))) {
    i <- idx[1, k]; j <- idx[2, k]
    pair <- stats::setNames(list(networks[[i]], networks[[j]]), c(labs[i], labs[j]))
    b <- bipartite::betalinkr(pair, partitioning = partitioning, binary = binary)
    out <- rbind(out, data.frame(comparison = paste(labs[i], "vs", labs[j]),
                                 S = b["S"], OS = b["OS"], WN = b["WN"],
                                 ST = if ("ST" %in% names(b)) b["ST"] else NA,
                                 row.names = NULL, stringsAsFactors = FALSE))
  }
  out
}
