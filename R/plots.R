# =====================================================================
# plots.R -- plotting helpers for bd_pipeline_results
# =====================================================================

#' Temporal trajectories of chosen metrics
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @param metrics Character vector of column names in `results$metrics_df`.
#' @return A list of ggplot objects (one per metric present and non-NA),
#'   or NULL entries for metrics that don't exist.
#' @export
plot_temporal_metrics <- function(results,
                                  metrics = c("connectance", "nestedness",
                                              "modularity", "robustness")) {
  .bd_require("ggplot2", "plots")
  df <- results$metrics_df
  df$order <- seq_len(nrow(df))
  lapply(metrics, function(m) {
    if (!m %in% names(df) || all(is.na(df[[m]]))) return(NULL)
    ggplot2::ggplot(df, ggplot2::aes(x = order, y = .data[[m]])) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = df$order, labels = df$label) +
      ggplot2::labs(title = paste(tools::toTitleCase(m), "over time"),
                    x = NULL, y = m) +
      ggplot2::theme_minimal()
  })
}

#' Species richness of both trophic levels over time
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @return A ggplot object.
#' @export
plot_richness <- function(results) {
  .bd_require(c("ggplot2", "tidyr"), "plots")
  df <- results$metrics_df; df$order <- seq_len(nrow(df))
  cfg <- results$config
  long <- tidyr::pivot_longer(df[, c("label", "order", "n_lower", "n_higher")],
                              c("n_lower", "n_higher"),
                              names_to = "level", values_to = "n")
  long$level <- ifelse(long$level == "n_lower", cfg$lower_name, cfg$higher_name)
  ggplot2::ggplot(long, ggplot2::aes(order, n, color = level)) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = df$order, labels = df$label) +
    ggplot2::labs(title = "Species richness over time", x = NULL,
                  y = "Number of species", color = NULL) +
    ggplot2::theme_minimal()
}

#' Turnover between consecutive periods: Jaccard + gains/losses
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @return A list of two ggplot objects (Jaccard similarity; gains/losses).
#' @export
plot_turnover <- function(results) {
  .bd_require(c("ggplot2", "tidyr"), "plots")
  cfg <- results$config; td <- results$persistence$turnover
  td$transition <- paste(td$from, "->", td$to)

  jac <- tidyr::pivot_longer(td, c("lower_jaccard", "higher_jaccard"),
                             names_to = "level", values_to = "jaccard")
  jac$level <- ifelse(jac$level == "lower_jaccard", cfg$lower_name, cfg$higher_name)
  p1 <- ggplot2::ggplot(jac, ggplot2::aes(transition, jaccard, fill = level)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(title = "Jaccard similarity between consecutive periods",
                  x = NULL, y = "Jaccard", fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  gl <- tidyr::pivot_longer(
    td[, c("transition", "lower_gains", "lower_losses",
           "higher_gains", "higher_losses")],
    -transition, names_to = "type", values_to = "count")
  gl$level <- ifelse(grepl("lower", gl$type), cfg$lower_name, cfg$higher_name)
  gl$change <- ifelse(grepl("gain", gl$type), "Gains", "Losses")
  p2 <- ggplot2::ggplot(gl, ggplot2::aes(transition, count, fill = change)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~level) +
    ggplot2::scale_fill_manual(values = c(Gains = "darkgreen", Losses = "darkred")) +
    ggplot2::labs(title = "Species gains and losses", x = NULL, y = "Species",
                  fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  list(p1, p2)
}

#' Leave-one-out contribution plots: expected vs observed, efficiency
#'
#' Base-graphics 2x2 panel: (a) species contribution, (b) interaction
#' contribution, (c) contribution efficiency, (d) unique species by level.
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @return NULL, invisibly; called for its plotting side effect.
#' @export
plot_loo_contribution <- function(results) {
  loo <- results$loo
  op <- graphics::par(mfrow = c(2, 2)); on.exit(graphics::par(op))
  bar2 <- function(top, bottom, main, ylab, leg, ymax = NULL) {
    m <- rbind(top, bottom)
    graphics::barplot(m, beside = TRUE, names.arg = loo$network,
                      col = c("gray70", "gray25"), legend.text = leg,
                      args.legend = list(x = "topright", cex = 0.8),
                      main = main, ylab = ylab,
                      ylim = c(0, if (is.null(ymax)) max(m) * 1.3 else ymax))
  }
  bar2(loo$exp_species_pct, loo$obs_unique_species_pct,
       "(a) Species contribution (% metaweb)", "%", c("Expected", "Observed unique"))
  bar2(loo$exp_interactions_pct, loo$obs_unique_interactions_pct,
       "(b) Interaction contribution (% metaweb)", "%", c("Expected", "Observed unique"))
  bar2(loo$species_efficiency, loo$interaction_efficiency,
       "(c) Contribution efficiency", "Obs/Exp", c("Species", "Interactions"), ymax = 1)
  graphics::abline(h = 0.5, lty = 2, col = "firebrick2")
  bar2(loo$unique_lower, loo$unique_higher,
       "(d) Unique species by level", "n",
       c(results$config$lower_name, results$config$higher_name))
  invisible(NULL)
}

#' Taxonomic composition for one period (needs a taxonomy lookup)
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @param label Period label (default: the first period).
#' @return A ggplot object, or NULL (with a message) if no taxonomy was
#'   supplied.
#' @export
plot_taxonomic_composition <- function(results, label = NULL) {
  .bd_require("ggplot2", "plots")
  if (is.null(results$taxonomy)) { message("No taxonomy supplied."); return(NULL) }
  if (is.null(label)) label <- names(results$networks)[1]
  tr <- results$temporal_metrics[[label]]$taxonomic_richness
  ggplot2::ggplot(tr, ggplot2::aes(stats::reorder(Group, Species_Count),
                                   Species_Count, fill = Group)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title = paste("Taxonomic richness (", label, ")"),
                  x = NULL, y = "Species count") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
}


# ---------------------------------------------------------------------
# NETWORK VISUALIZATIONS
# ---------------------------------------------------------------------

#' Bipartite plotweb (two-level interaction diagram)
#'
#' Pass a matrix directly, or a `bd_pipeline_results` + period `label`.
#' Uses the bipartite >= 2.24 `plotweb()` API when available, falling back
#' to the pre-2.24 argument names so the original orange/green styling
#' still works on older installs.
#'
#' @param x A matrix, or a `bd_pipeline_results` object.
#' @param label Period label (used when `x` is a `bd_pipeline_results`);
#'   `"metaweb"` plots the pooled metaweb.
#' @param sorting Passed to `bipartite::plotweb()`.
#' @param col_low,col_high,col_interaction Fill colors for lower level,
#'   higher level, and interaction links.
#' @param labels Show species labels.
#' @param ... Passed to `bipartite::plotweb()`.
#' @return NULL, invisibly; called for its plotting side effect.
#' @export
plot_network_web <- function(x, label = NULL, sorting = "ca",
                             col_low = "chartreuse3", col_high = "darkorange1",
                             col_interaction = "azure3", labels = TRUE, ...) {
  .bd_require("bipartite", "plotweb")
  web <- .bd_pick_web(x, label)
  new_api <- "lower_color" %in% names(formals(bipartite::plotweb))
  if (new_api) {
    bipartite::plotweb(web, sorting = sorting, srt = 90,
                       lower_labels = labels, higher_labels = labels,
                       link_color = col_interaction, link_border = "azure4",
                       lower_color = col_low, lower_border = "chartreuse4",
                       higher_color = col_high, higher_border = "darkorange4", ...)
  } else {
    ll <- if (labels) NULL else 0
    bipartite::plotweb(web, method = sorting, text.rot = 90,
                       high.lablength = ll, low.lablength = ll,
                       col.interaction = col_interaction,
                       col.high = col_high, col.low = col_low,
                       bor.col.interaction = "azure4",
                       bor.col.high = "darkorange4", bor.col.low = "chartreuse4", ...)
  }
  invisible(NULL)
}

#' Visweb (nested presence/absence matrix)
#'
#' @param x A matrix, or a `bd_pipeline_results` object.
#' @param label Period label (used when `x` is a `bd_pipeline_results`);
#'   `"metaweb"` plots the pooled metaweb.
#' @param type Passed to `bipartite::visweb()`.
#' @param ... Passed to `bipartite::visweb()`.
#' @return NULL, invisibly; called for its plotting side effect.
#' @export
plot_visweb <- function(x, label = NULL, type = "nested", ...) {
  .bd_require("bipartite", "visweb")
  web <- .bd_pick_web(x, label)
  bipartite::visweb(web, type = type, square = "interaction",
                    text = "none", frame = NULL, ...)
  invisible(NULL)
}

# internal: accept a matrix directly, or (results, label)
.bd_pick_web <- function(x, label) {
  if (is.matrix(x)) return(x)
  if (inherits(x, "bd_pipeline_results")) {
    if (is.null(label)) label <- names(x$networks)[1]
    if (identical(label, "__metaweb__") || identical(label, "metaweb"))
      return(x$metaweb)
    return(x$networks[[label]])
  }
  stop("plot needs a matrix or a bd_pipeline_results (+ label).")
}

#' betalink two-network overlay for a pair of periods
#'
#' Draws shared vs exclusive species/links between two periods. Pass
#' `b = "metaweb"` to compare a single period against the pooled metaweb.
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @param a,b Period labels in `results` (or `"metaweb"` for the pooled web).
#' @param na,nb,ns Colors for network `a`, network `b`, and shared nodes.
#' @param show_metrics Return the `betalink::betalink()` metrics (default
#'   TRUE).
#' @return If `show_metrics`, the `betalink::betalink()` result, invisibly;
#'   otherwise NULL.
#' @export
plot_beta_pair <- function(results, a, b,
                           na = "dodgerblue2", nb = "goldenrod2",
                           ns = "grey", show_metrics = TRUE) {
  .bd_require("betalink", "betalink plots")
  wa <- .bd_pick_web(results, a)
  wb <- .bd_pick_web(results, b)
  pl <- betalink::prepare_networks(stats::setNames(list(wa, wb), c(a, b)),
                                   directed = FALSE)
  betalink::network_betaplot(pl[[1]], pl[[2]], na = na, nb = nb, ns = ns,
                             vertex.label = NA)
  graphics::title(main = paste(a, "vs", b))
  if (show_metrics) {
    m <- betalink::betalink(pl[[1]], pl[[2]], bf = betalink::B01)
    return(invisible(m))
  }
  invisible(NULL)
}

#' betalink overlays for every consecutive pair of periods
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @param palette Optional color vector, one per period (default:
#'   `grDevices::hcl.colors()`).
#' @return A named list of [plot_beta_pair()] results, invisibly.
#' @export
plot_beta_consecutive <- function(results, palette = NULL) {
  labs <- names(results$networks)
  if (is.null(palette))
    palette <- grDevices::hcl.colors(length(labs), "Dynamic")
  out <- list()
  for (i in seq_len(length(labs) - 1)) {
    out[[paste(labs[i], labs[i + 1])]] <-
      plot_beta_pair(results, labs[i], labs[i + 1],
                     na = palette[i], nb = palette[i + 1])
  }
  invisible(out)
}

#' betalink overlay of each period against the pooled metaweb
#'
#' @param results A `bd_pipeline_results` from [run_network_pipeline()].
#' @param palette Optional color vector, one per period (default:
#'   `grDevices::hcl.colors()`).
#' @return A named list of [plot_beta_pair()] results, invisibly.
#' @export
plot_beta_vs_metaweb <- function(results, palette = NULL) {
  labs <- names(results$networks)
  if (is.null(palette))
    palette <- grDevices::hcl.colors(length(labs), "Dynamic")
  out <- list()
  for (i in seq_along(labs)) {
    out[[labs[i]]] <- plot_beta_pair(results, labs[i], "metaweb",
                                     na = palette[i], nb = "grey40")
  }
  invisible(out)
}
