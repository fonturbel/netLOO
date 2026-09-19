# =====================================================================
# wrappers.R -- standalone entry points for single analysis components
# =====================================================================
#
# These skip run_network_pipeline() (and therefore skip null models and
# every other block) while still accepting a bd_pipeline_config, so
# existing cfg-based scripts keep working unchanged.

#' Connectance-corrected NODF without running the full pipeline
#'
#' @param x A bipartite matrix, a path to a matrix file, a named list of
#'   matrices, or a `bd_pipeline_config` (see [pipeline_config()]).
#' @param label When `x` is a config or a list of networks, the name/index of
#'   the single network to compute NODF_c for. If NULL, NODF_c is computed
#'   for every network and returned as a named list.
#' @param sep Field separator when `x` is a file path (default ";").
#' @param quality Passed to [calc_NODFc()].
#' @return The list from [calc_NODFc()] for a single network, or a named list
#'   of such results when computed over several networks.
#' @export
nodf_c <- function(x, label = NULL, sep = ";", quality = 2) {
  if (inherits(x, "bd_pipeline_config")) {
    nets <- load_networks(x)
    if (!is.null(label)) return(calc_NODFc(nets[[label]], quality = quality))
    return(stats::setNames(lapply(nets, calc_NODFc, quality = quality), names(nets)))
  }
  if (is.list(x) && !is.matrix(x)) {
    if (!is.null(label)) return(calc_NODFc(x[[label]], quality = quality))
    return(stats::setNames(lapply(x, calc_NODFc, quality = quality), names(x)))
  }
  if (is.character(x) && length(x) == 1) x <- read_network_matrix(x, sep)
  calc_NODFc(x, quality = quality)
}

#' Leave-one-out contribution analysis without running the full pipeline
#'
#' Skips null models, rarefaction, beta diversity and per-network metrics --
#' only builds the metaweb (and leave-one-out metawebs) and runs
#' [loo_contribution_analysis()].
#'
#' @param x A named list of matrices, or a `bd_pipeline_config` (networks are
#'   loaded via [load_networks()]).
#' @param metaweb Optional pooled metaweb; built automatically with
#'   [build_metaweb()] when NULL. If `x` is a config whose `cfg$metaweb` file
#'   exists, that file is used instead unless `metaweb` is supplied explicitly.
#' @return A data frame, see [loo_contribution_analysis()].
#' @export
loo_only <- function(x, metaweb = NULL) {
  if (inherits(x, "bd_pipeline_config")) {
    networks <- load_networks(x)
    if (is.null(metaweb)) {
      metaweb <- if (!is.null(x$metaweb) && file.exists(x$metaweb))
        read_network_matrix(x$metaweb, x$sep) else build_metaweb(networks)
    }
  } else {
    networks <- x
    if (is.null(metaweb)) metaweb <- build_metaweb(networks)
  }
  loo_contribution_analysis(networks, metaweb)
}
