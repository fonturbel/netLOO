# =====================================================================
# metrics.R -- per-network structural metrics
# =====================================================================

#' Connectance-corrected NODF (Song et al. 2017)
#'
#' Needs the `maxnodf` and `bipartite` packages.
#'
#' @param web A bipartite matrix.
#' @param quality Passed to `maxnodf::maxnodf()`.
#' @return A list with `NODF_raw`, `max_NODF`, `NODF_n`, `connectance`,
#'   `S_geometric`, `NODF_c`.
#' @export
calc_NODFc <- function(web, quality = 2) {
  .bd_require(c("bipartite", "vegan", "maxnodf"), "NODFc")
  web_bin <- (web > 0) * 1
  nr <- nrow(web_bin); nc <- ncol(web_bin)
  nodf_raw <- vegan::nestednodf(web_bin)$statistic["NODF"]
  invisible(utils::capture.output(
    max_val <- suppressWarnings(maxnodf::maxnodf(web_bin, quality = quality)$max_nodf)))
  nodf_n <- nodf_raw / max_val
  C <- sum(web_bin) / (nr * nc)
  S <- sqrt(nr * nc)
  list(NODF_raw = as.numeric(nodf_raw), max_NODF = max_val,
       NODF_n = as.numeric(nodf_n), connectance = C, S_geometric = S,
       NODF_c = as.numeric(nodf_n / (C * log(S))))
}

#' H2' specialization with the bipartite 2.23 bug fix
#'
#' Manually recomputes H2 from `H2fun()`'s uncorrected/min/max values,
#' working around a normalization bug present in bipartite 2.23.
#'
#' @param web A bipartite matrix.
#' @param weight_scale Multiplier applied before calling `H2fun()`.
#' @param ... Passed to `bipartite::H2fun()`.
#' @return The named numeric vector returned by `H2fun()`, with `H2` fixed.
#' @export
H2_fixed <- function(web, weight_scale = 5, ...) {
  .bd_require("bipartite", "H2'")
  res <- bipartite::H2fun(web * weight_scale, ...)
  res["H2"] <- (res["H2uncorr"] - res["H2min"]) / (res["H2max"] - res["H2min"])
  res
}

#' Robustness of the lower level to removal of higher-level species
#'
#' @param web A bipartite matrix.
#' @param method `"random"` (averaged over `reps` random removal orders) or
#'   `"degree"` (most-connected higher-level species removed first).
#' @param reps Replicates for `method = "random"`.
#' @return A numeric vector: fraction of lower-level species retained after
#'   each successive removal.
#' @export
calc_robustness <- function(web, method = "random", reps = 100) {
  one_run <- function(order_idx) {
    wc <- web
    frac <- numeric(ncol(web))
    for (i in seq_along(order_idx)) {
      wc[, order_idx[i]] <- 0
      frac[i] <- sum(rowSums(wc) > 0) / nrow(web)
    }
    frac
  }
  if (method == "degree") {
    one_run(order(colSums(web), decreasing = TRUE))
  } else {
    rowMeans(replicate(reps, one_run(sample(seq_len(ncol(web))))))
  }
}

#' All single-network metrics for one web at one time point
#'
#' @param web A bipartite matrix.
#' @param label Time-point label.
#' @param taxonomy Optional taxonomy lookup (see [load_taxonomy()]).
#' @param do_nodfc Compute connectance-corrected NODF (needs `maxnodf`).
#' @param robustness_reps Random-removal replicates for the robustness curve.
#' @return A list of per-network metrics.
#' @export
calculate_network_metrics <- function(web, label, taxonomy = NULL,
                                       do_nodfc = TRUE, robustness_reps = 100) {
  .bd_require("bipartite", "network metrics")
  n_lower <- nrow(web); n_higher <- ncol(web)
  n_links <- sum(web > 0)
  connectance <- n_links / (n_lower * n_higher)

  nestedness <- tryCatch(bipartite::networklevel(web, index = "NODF"),
                         error = function(e) NA)
  modularity <- tryCatch(bipartite::computeModules(web)@likelihood,
                         error = function(e) NA)
  other <- tryCatch(
    bipartite::networklevel(web, index = c("web asymmetry", "links per species",
                                           "H2", "generality", "vulnerability")),
    error = function(e) c(`web asymmetry` = NA, `links per species` = NA,
                          H2 = NA, generality = NA, vulnerability = NA))
  nodfc <- if (do_nodfc) tryCatch(calc_NODFc(web)$NODF_c,
                                  error = function(e) NA) else NA

  robustness <- mean(calc_robustness(web, "random", robustness_reps))

  tax_richness <- NULL
  if (!is.null(taxonomy)) {
    grp <- classify_species(colnames(web), taxonomy)
    tax_richness <- as.data.frame(table(Group = grp),
                                  responseName = "Species_Count")
  }

  list(label = label, n_lower = n_lower, n_higher = n_higher,
       n_links = n_links, connectance = connectance,
       nestedness = if (is.list(nestedness)) nestedness$NODF else nestedness,
       NODFc = nodfc, modularity = modularity,
       web_asymmetry = unname(other["web asymmetry"]),
       links_per_species = unname(other["links per species"]),
       H2 = unname(other["H2"]),
       generality = unname(other["generality"]),
       vulnerability = unname(other["vulnerability"]),
       robustness = robustness, taxonomic_richness = tax_richness)
}

#' Tidy per-network metrics into one data frame
#'
#' @param temporal_metrics A named list of metric lists, one per time point
#'   (as produced inside [run_network_pipeline()]).
#' @return A data frame, one row per period.
#' @export
metrics_to_df <- function(temporal_metrics) {
  do.call(rbind, lapply(temporal_metrics, function(m) data.frame(
    label = m$label, n_lower = m$n_lower, n_higher = m$n_higher,
    n_links = m$n_links, connectance = m$connectance,
    nestedness = m$nestedness, NODFc = m$NODFc, modularity = m$modularity,
    web_asymmetry = m$web_asymmetry, links_per_species = m$links_per_species,
    H2 = m$H2, generality = m$generality, vulnerability = m$vulnerability,
    robustness = m$robustness, stringsAsFactors = FALSE)))
}
