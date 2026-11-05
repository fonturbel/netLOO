# =====================================================================
# metrics.R -- per-network structural metrics
# =====================================================================

#' Connectance-corrected NODF (Song et al. 2017)
#'
#' NODF_c = NODF_n / (C * log10(S)), where NODF_n is NODF divided by the
#' maximum NODF attainable for the same dimensions and fill (both on a 0-1
#' scale), C is connectance and S the geometric mean of the two richnesses.
#' This is the same formula as `maxnodf::NODFc()`. Versions <= 0.0.3 mixed a
#' 0-100 and a 0-1 scale and used the natural log, inflating NODF_c ~43x.
#'
#' Needs the `vegan` and `maxnodf` packages.
#'
#' @param web A bipartite matrix.
#' @param quality Passed to `maxnodf::maxnodf()` (0 = greedy, 1 = + hill
#'   climbing, 2 = + simulated annealing).
#' @return A list with `NODF_raw` (0-1), `max_NODF`, `NODF_n`,
#'   `connectance`, `S_geometric`, `NODF_c`.
#' @export
calc_NODFc <- function(web, quality = 1) {
  .bd_require(c("vegan", "maxnodf"), "NODFc")
  web_bin <- (web > 0) * 1
  nr <- nrow(web_bin); nc <- ncol(web_bin)
  nodf_raw <- unname(vegan::nestednodf(web_bin)$statistic["NODF"]) / 100
  invisible(utils::capture.output(
    max_val <- suppressWarnings(maxnodf::maxnodf(web_bin, quality = quality)$max_nodf)))
  if (!is.finite(max_val) || max_val <= 0)
    stop("Maximum NODF is 0 for this matrix; NODF_c undefined.")
  nodf_n <- nodf_raw / max_val
  C <- sum(web_bin) / (nr * nc)
  S <- sqrt(nr * nc)
  list(NODF_raw = nodf_raw, max_NODF = max_val, NODF_n = nodf_n,
       connectance = C, S_geometric = S,
       NODF_c = nodf_n / (C * log10(S)))
}

#' Network-level specialisation H2' (Bluethgen et al. 2006)
#'
#' H2' needs interaction frequencies. In a binary web the observed entropy
#' always equals its maximum, so H2' carries no information; this returns NA
#' with a warning instead of a misleading number. For weighted (integer) webs
#' it returns `bipartite::H2fun()` unchanged, where
#' H2' = (H2max - H2uncorr) / (H2max - H2min).
#'
#' @param web A bipartite matrix.
#' @param ... Passed to `bipartite::H2fun()`.
#' @return The named numeric vector from `H2fun()` (`H2`, `H2min`, `H2max`,
#'   `H2uncorr`), all NA for binary webs.
#' @export
calc_H2 <- function(web, ...) {
  .bd_require("bipartite", "H2'")
  if (is_binary_web(web)) {
    warning("H2' is not defined for binary (presence/absence) webs; returning NA.",
            call. = FALSE)
    return(c(H2 = NA_real_, H2min = NA_real_, H2max = NA_real_, H2uncorr = NA_real_))
  }
  bipartite::H2fun(web, ...)
}

#' Deprecated: use calc_H2()
#'
#' Versions <= 0.0.3 multiplied the web by `weight_scale` and returned
#' (H2uncorr - H2min) / (H2max - H2min), which is 1 - H2'. `bipartite::H2fun()`
#' was correct all along (H2' = 0 is the expected value for binary webs).
#' This now calls [calc_H2()] and ignores `weight_scale`.
#'
#' @param web A bipartite matrix.
#' @param weight_scale Ignored.
#' @param ... Passed to [calc_H2()].
#' @return See [calc_H2()].
#' @export
H2_fixed <- function(web, weight_scale = NULL, ...) {
  .Deprecated("calc_H2", msg = paste(
    "H2_fixed() returned 1 - H2' on a rescaled binary web and is deprecated.",
    "Use calc_H2() (NA for binary webs)."))
  calc_H2(web, ...)
}

# Pull one value from a networklevel() result whose names may carry a
# .HL / .LL suffix (e.g. "generality.HL", "vulnerability.LL").
.nl_value <- function(x, name) {
  hit <- names(x)[names(x) == name | startsWith(names(x), paste0(name, "."))]
  if (length(hit)) unname(x[hit[1]]) else NA_real_
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
#' @param nodfc_quality Passed to [calc_NODFc()].
#' @return A list of per-network metrics. `H2` is NA for binary webs.
#' @export
calculate_network_metrics <- function(web, label, taxonomy = NULL,
                                       do_nodfc = TRUE, robustness_reps = 100,
                                       nodfc_quality = 1) {
  .bd_require("bipartite", "network metrics")
  n_lower <- nrow(web); n_higher <- ncol(web)
  n_links <- sum(web > 0)
  connectance <- n_links / (n_lower * n_higher)

  nestedness <- tryCatch(unname(bipartite::networklevel(web, index = "NODF")),
                         error = function(e) NA)
  modularity <- tryCatch(bipartite::computeModules(web)@likelihood,
                         error = function(e) NA)
  other <- tryCatch(
    bipartite::networklevel(web, index = c("web asymmetry", "links per species",
                                           "generality", "vulnerability")),
    error = function(e) c(`web asymmetry` = NA, `links per species` = NA,
                          generality = NA, vulnerability = NA))
  h2 <- if (is_binary_web(web)) NA_real_ else
    tryCatch(unname(calc_H2(web)["H2"]), error = function(e) NA_real_)
  nodfc <- if (do_nodfc) tryCatch(calc_NODFc(web, nodfc_quality)$NODF_c,
                                  error = function(e) {
                                    warning(label, ": NODFc failed: ",
                                            conditionMessage(e), call. = FALSE)
                                    NA }) else NA

  robustness <- mean(calc_robustness(web, "random", robustness_reps))

  tax_richness <- NULL
  if (!is.null(taxonomy)) {
    grp <- classify_species(colnames(web), taxonomy)
    tax_richness <- as.data.frame(table(Group = grp),
                                  responseName = "Species_Count")
  }

  list(label = label, n_lower = n_lower, n_higher = n_higher,
       n_links = n_links, connectance = connectance,
       nestedness = nestedness, NODFc = nodfc, modularity = modularity,
       web_asymmetry = .nl_value(other, "web asymmetry"),
       links_per_species = .nl_value(other, "links per species"),
       H2 = h2,
       generality = .nl_value(other, "generality"),
       vulnerability = .nl_value(other, "vulnerability"),
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
