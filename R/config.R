# =====================================================================
# config.R -- package setup helpers and pipeline configuration
# =====================================================================

# Core packages always needed; optional ones loaded on demand per module
.bd_require <- function(pkgs, who = "this step") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("Package(s) needed for %s not installed: %s\nInstall with install.packages() (or devtools::install_github for maxnodf).",
                 who, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(
    require(p, character.only = TRUE))))
}

#' Is a web binary?
#'
#' @param web A matrix.
#' @return TRUE if every non-missing cell is 0 or 1.
#' @export
is_binary_web <- function(web) all(web[!is.na(web)] %in% c(0, 1))

#' Build a pipeline configuration object
#'
#' @param manifest  Path to a CSV with columns `label` and `file` giving the
#'                  ordered time points and their matrix files. If NULL, files
#'                  matching `file_pattern` in `data_dir` are auto-discovered.
#' @param data_dir  Folder holding the matrix files (and where files in the
#'                  manifest are resolved relative to).
#' @param file_pattern Regex to auto-discover matrices when manifest is NULL.
#' @param sep       Field separator of the matrix CSVs (default ";").
#' @param metaweb   Optional path to a pooled/metaweb matrix. If NULL it is
#'                  built automatically as the union of all networks. A
#'                  supplied metaweb is checked against that union (see
#'                  [check_metaweb()]).
#' @param taxonomy  Optional path to a CSV with columns `species` and `group`
#'                  mapping higher-level (column) species to taxonomic groups.
#' @param n_null    Number of null models for significance testing.
#' @param null_method Null model, see [generate_null_webs()]. A character
#'                  string is passed to `vegan::nullmodel()`; the default
#'                  `"quasiswap"` keeps every row/column total and the number
#'                  of links fixed, the appropriate null for binary webs. A
#'                  number (1-5) uses `bipartite::nullmodel()`, which
#'                  generates weighted matrices (quantitative webs only).
#' @param core_threshold Persistence fraction to call a species "core".
#' @param robustness_reps Random-removal replicates for robustness curves.
#' @param n_cores   Cores for parallel null models (NULL = auto).
#' @param do_nodfc  Compute connectance-corrected NODF (needs maxnodf).
#' @param nodfc_quality `maxnodf` optimisation level (0 = greedy, 1 = + hill
#'                  climbing, 2 = + simulated annealing; 2 is slowest and
#'                  has crashed R on some matrices).
#' @param seed      RNG seed so null models, robustness and modularity are
#'                  reproducible. NULL leaves the RNG untouched.
#' @param lower_name,higher_name Friendly labels for the two trophic levels.
#' @param weight_scale Deprecated and ignored. H2' is no longer computed on
#'                  rescaled binary webs (see [calc_H2()]).
#' @return A `bd_pipeline_config` object to pass to [run_network_pipeline()].
#' @export
pipeline_config <- function(manifest        = NULL,
                            data_dir        = "data",
                            file_pattern    = "\\.csv$",
                            sep             = ";",
                            metaweb         = NULL,
                            taxonomy        = NULL,
                            n_null          = 1000,
                            null_method     = "quasiswap",
                            core_threshold  = 0.5,
                            robustness_reps = 100,
                            n_cores         = NULL,
                            do_nodfc        = TRUE,
                            nodfc_quality   = 1,
                            seed            = 1,
                            lower_name      = "Plants",
                            higher_name     = "Pollinators",
                            weight_scale    = NULL) {
  if (!is.null(weight_scale))
    warning("`weight_scale` is deprecated and ignored: H2' is not defined ",
            "for binary webs and rescaling them gives arbitrary values.",
            call. = FALSE)
  structure(list(
    manifest = manifest, data_dir = data_dir, file_pattern = file_pattern,
    sep = sep, metaweb = metaweb, taxonomy = taxonomy,
    n_null = n_null, null_method = null_method,
    core_threshold = core_threshold, robustness_reps = robustness_reps,
    n_cores = n_cores, do_nodfc = do_nodfc, nodfc_quality = nodfc_quality,
    seed = seed, lower_name = lower_name, higher_name = higher_name
  ), class = "bd_pipeline_config")
}
