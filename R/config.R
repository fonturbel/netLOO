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
#'                  built automatically as the union of all networks.
#' @param taxonomy  Optional path to a CSV with columns `species` and `group`
#'                  mapping higher-level (column) species to taxonomic groups.
#' @param weight_scale Multiplier applied before H2' (bipartite needs integers
#'                  for binary webs; 5 mirrors the original analysis).
#' @param n_null    Number of null models for significance testing.
#' @param null_method bipartite::nullmodel method id (1 = non-sequential).
#' @param core_threshold Persistence fraction to call a species "core".
#' @param robustness_reps Random-removal replicates for robustness curves.
#' @param n_cores   Cores for parallel null models (NULL = auto).
#' @param do_nodfc  Compute connectance-corrected NODF (needs maxnodf).
#' @param lower_name,higher_name Friendly labels for the two trophic levels.
#' @return A `bd_pipeline_config` object to pass to [run_network_pipeline()].
#' @export
pipeline_config <- function(manifest        = NULL,
                            data_dir        = "data",
                            file_pattern    = "\\.csv$",
                            sep             = ";",
                            metaweb         = NULL,
                            taxonomy        = NULL,
                            weight_scale    = 5,
                            n_null          = 1000,
                            null_method     = 1,
                            core_threshold  = 0.5,
                            robustness_reps = 100,
                            n_cores         = NULL,
                            do_nodfc        = TRUE,
                            lower_name      = "Plants",
                            higher_name     = "Pollinators") {
  structure(list(
    manifest = manifest, data_dir = data_dir, file_pattern = file_pattern,
    sep = sep, metaweb = metaweb, taxonomy = taxonomy,
    weight_scale = weight_scale, n_null = n_null, null_method = null_method,
    core_threshold = core_threshold, robustness_reps = robustness_reps,
    n_cores = n_cores, do_nodfc = do_nodfc,
    lower_name = lower_name, higher_name = higher_name
  ), class = "bd_pipeline_config")
}
