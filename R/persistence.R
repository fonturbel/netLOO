# =====================================================================
# persistence.R -- species persistence, turnover, core species, trends
# =====================================================================

#' Species persistence and turnover across time points
#'
#' @param networks A named list of matrices, ordered by time point.
#' @return A list with `lower_presence`, `higher_presence` (binary
#'   species x period matrices), `lower_persistence`, `higher_persistence`
#'   (fraction of periods each species is present in), and `turnover`
#'   (a data frame of Jaccard similarity and gains/losses for each
#'   consecutive pair of periods).
#' @export
analyze_species_persistence <- function(networks) {
  all_lower  <- unique(unlist(lapply(networks, rownames)))
  all_higher <- unique(unlist(lapply(networks, colnames)))
  L <- matrix(0, length(all_lower), length(networks),
              dimnames = list(all_lower, names(networks)))
  H <- matrix(0, length(all_higher), length(networks),
              dimnames = list(all_higher, names(networks)))
  for (yr in names(networks)) {
    L[rownames(networks[[yr]]), yr] <- 1
    H[colnames(networks[[yr]]), yr] <- 1
  }
  lower_persist  <- rowSums(L) / ncol(L)
  higher_persist <- rowSums(H) / ncol(H)

  labels <- names(networks); turnover <- data.frame()
  for (i in seq_len(length(labels) - 1)) {
    a <- labels[i]; b <- labels[i + 1]
    la <- rownames(networks[[a]]); lb <- rownames(networks[[b]])
    ha <- colnames(networks[[a]]); hb <- colnames(networks[[b]])
    turnover <- rbind(turnover, data.frame(
      from = a, to = b,
      lower_jaccard  = length(intersect(la, lb)) / length(union(la, lb)),
      higher_jaccard = length(intersect(ha, hb)) / length(union(ha, hb)),
      lower_gains  = length(setdiff(lb, la)), lower_losses  = length(setdiff(la, lb)),
      higher_gains = length(setdiff(hb, ha)), higher_losses = length(setdiff(ha, hb)),
      stringsAsFactors = FALSE))
  }
  list(lower_presence = L, higher_presence = H,
       lower_persistence = lower_persist, higher_persistence = higher_persist,
       turnover = turnover)
}

#' Identify core species (present in a threshold fraction of periods)
#'
#' @param networks A named list of matrices, ordered by time point.
#' @param threshold Minimum persistence fraction to call a species "core".
#' @param persistence Optional precomputed result of
#'   [analyze_species_persistence()] (recomputed if NULL).
#' @return A list with `core_lower`, `core_higher` (species names) and
#'   `core_metrics` (per-period core-subnetwork stats).
#' @export
identify_core_species <- function(networks, threshold = 0.5,
                                  persistence = NULL) {
  if (is.null(persistence)) persistence <- analyze_species_persistence(networks)
  core_lower  <- names(persistence$lower_persistence[persistence$lower_persistence  >= threshold])
  core_higher <- names(persistence$higher_persistence[persistence$higher_persistence >= threshold])
  core_metrics <- list()
  for (yr in names(networks)) {
    web <- networks[[yr]]
    cl <- intersect(core_lower, rownames(web))
    ch <- intersect(core_higher, colnames(web))
    if (length(cl) && length(ch)) {
      cw <- web[cl, ch, drop = FALSE]
      core_metrics[[yr]] <- list(
        n_core_lower = length(cl), n_core_higher = length(ch),
        core_links = sum(cw > 0),
        core_connectance = sum(cw > 0) / (nrow(cw) * ncol(cw)),
        proportion_core = sum(cw > 0) / sum(web > 0))
    } else {
      core_metrics[[yr]] <- list(n_core_lower = length(cl), n_core_higher = length(ch),
                                 core_links = 0, core_connectance = 0, proportion_core = 0)
    }
  }
  list(core_lower = core_lower, core_higher = core_higher, core_metrics = core_metrics)
}

#' Spearman correlation of each metric against ordered time
#'
#' With few time points these tests have almost no power: the smallest
#' two-sided p-value an exact Spearman test can return is 2 / n!
#' (0.33 for n = 3, 0.083 for n = 4, 0.017 for n = 5), reported as
#' `min_possible_p`.
#'
#' @param metrics_df A data frame from [metrics_to_df()].
#' @param time Optional numeric time axis (default: period order, 1..n).
#' @return A named list (one entry per non-constant numeric metric) with
#'   `rho`, `p_value`, `n` and `min_possible_p`.
#' @export
temporal_trend_tests <- function(metrics_df, time = NULL) {
  if (is.null(time)) time <- seq_len(nrow(metrics_df))
  num <- metrics_df[, vapply(metrics_df, is.numeric, logical(1)), drop = FALSE]
  out <- list()
  for (m in names(num)) {
    ok <- !is.na(num[[m]])
    n <- sum(ok)
    if (n >= 3 && length(unique(num[[m]][ok])) > 1) {
      ct <- suppressWarnings(stats::cor.test(time[ok], num[[m]][ok],
                                             method = "spearman"))
      out[[m]] <- list(rho = unname(ct$estimate), p_value = ct$p.value, n = n,
                       min_possible_p = min(1, 2 / factorial(n)))
    }
  }
  if (nrow(metrics_df) < 5)
    message("Only ", nrow(metrics_df), " time points: Spearman trends cannot ",
            "reach p < 0.05; treat them as descriptive.")
  out
}
