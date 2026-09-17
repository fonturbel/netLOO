# =====================================================================
# pipeline.R -- master driver
# =====================================================================

#' Run the full multi-temporal bipartite network pipeline
#'
#' Loads networks, computes per-network structural metrics, species
#' persistence/turnover/core, pairwise beta diversity, leave-one-out
#' contribution, null models and rarefaction, and returns everything in one
#' `bd_pipeline_results` object.
#'
#' @param cfg A `bd_pipeline_config` (see [pipeline_config()]).
#' @return A `bd_pipeline_results` list; safe to `save()` and feed to plotting
#'   helpers such as [plot_richness()] and [plot_loo_contribution()].
#' @export
run_network_pipeline <- function(cfg) {
  stopifnot(inherits(cfg, "bd_pipeline_config"))
  cat("=== Multi-temporal bipartite network pipeline ===\n")

  networks <- load_networks(cfg)
  taxonomy <- load_taxonomy(cfg)

  metaweb <- if (!is.null(cfg$metaweb) && file.exists(cfg$metaweb))
    read_network_matrix(cfg$metaweb, cfg$sep) else build_metaweb(networks)

  cat("\n-- per-network metrics --\n")
  temporal_metrics <- stats::setNames(
    lapply(names(networks), function(yr) {
      cat("  ", yr, "\n")
      calculate_network_metrics(networks[[yr]], yr, taxonomy,
                                cfg$do_nodfc, cfg$robustness_reps)
    }), names(networks))
  metrics_df <- metrics_to_df(temporal_metrics)

  cat("-- persistence / turnover / core --\n")
  persistence <- analyze_species_persistence(networks)
  core <- identify_core_species(networks, cfg$core_threshold, persistence)
  trends <- temporal_trend_tests(metrics_df)

  cat("-- beta diversity --\n")
  beta <- tryCatch(pairwise_betadiversity(networks),
                   error = function(e) { warning(e$message); NULL })

  cat("-- leave-one-out contribution --\n")
  loo <- loo_contribution_analysis(networks, metaweb)

  cat("-- null models --\n")
  null_models <- stats::setNames(
    lapply(names(networks), function(yr)
      tryCatch(enhanced_null_models(networks[[yr]], cfg$n_null,
                                    cfg$n_cores, cfg$null_method),
               error = function(e) list(error = e$message))),
    names(networks))

  cat("-- rarefaction --\n")
  rarefaction <- tryCatch(rarefaction_analysis(networks, metaweb),
                          error = function(e) { warning(e$message); NULL })

  cat("=== done ===\n")
  structure(list(
    config = cfg, networks = networks, metaweb = metaweb, taxonomy = taxonomy,
    temporal_metrics = temporal_metrics, metrics_df = metrics_df,
    persistence = persistence, core = core, trends = trends,
    beta = beta, loo = loo, null_models = null_models,
    rarefaction = rarefaction), class = "bd_pipeline_results")
}

#' Print a bd_pipeline_results object
#'
#' @param x A `bd_pipeline_results` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.bd_pipeline_results <- function(x, ...) {
  cat("<bd_pipeline_results>\n")
  cat("  networks :", paste(names(x$networks), collapse = ", "), "\n")
  cat("  metrics  :", nrow(x$metrics_df), "periods x",
      ncol(x$metrics_df), "metrics\n")
  cat("  beta     :", if (is.null(x$beta)) "NA" else nrow(x$beta), "pairwise\n")
  cat("  blocks   : temporal_metrics, persistence, core, trends, beta, loo,",
      "null_models, rarefaction\n")
  invisible(x)
}
