# =====================================================================
# loo.R -- leave-one-out contribution analysis
# =====================================================================

#' Metaweb built from all networks except one focal network
#'
#' @param network_list A named list of matrices, ordered by time point.
#' @param exclude_index Index of the network to leave out.
#' @return A matrix: the union metaweb of every network except
#'   `network_list[[exclude_index]]`.
#' @export
build_loo_metaweb <- function(network_list, exclude_index) {
  build_metaweb(network_list[-exclude_index])
}

#' Species and interaction overlap of a focal network with a reference metaweb
#'
#' @param focal_net A bipartite matrix.
#' @param reference_meta A reference metaweb matrix, e.g. from
#'   [build_metaweb()] or [build_loo_metaweb()].
#' @return A list with species/interaction counts and their overlap with
#'   `reference_meta`.
#' @export
calc_contribution <- function(focal_net, reference_meta) {
  fl <- c(rownames(focal_net), colnames(focal_net))
  rl <- c(rownames(reference_meta), colnames(reference_meta))
  f_links <- which(focal_net > 0, arr.ind = TRUE)
  f_int <- paste(rownames(focal_net)[f_links[, 1]],
                 colnames(focal_net)[f_links[, 2]], sep = "_")
  r_links <- which(reference_meta > 0, arr.ind = TRUE)
  r_int <- paste(rownames(reference_meta)[r_links[, 1]],
                 colnames(reference_meta)[r_links[, 2]], sep = "_")
  list(n_species = length(unique(fl)), n_interactions = length(f_int),
       shared_species = length(intersect(fl, rl)),
       unique_species = length(setdiff(fl, rl)),
       unique_lower  = length(setdiff(rownames(focal_net), rownames(reference_meta))),
       unique_higher = length(setdiff(colnames(focal_net), colnames(reference_meta))),
       shared_interactions = length(intersect(f_int, r_int)),
       unique_interactions = length(setdiff(f_int, r_int)))
}

#' Leave-one-out contribution analysis
#'
#' For each network, compares its expected contribution to the pooled
#' metaweb against its observed-unique contribution relative to the
#' leave-one-out metaweb (i.e. what it adds that no other period covers).
#' An efficiency index (observed unique / expected) summarizes how much of
#' each network's apparent contribution is actually unique to it.
#'
#' @param networks A named list of matrices, ordered by time point.
#' @param metaweb The pooled metaweb, e.g. from [build_metaweb()].
#' @return A data frame, one row per network, with expected/observed-unique
#'   species and interaction counts (raw and as a % of the metaweb) and
#'   `species_efficiency` / `interaction_efficiency`.
#' @export
loo_contribution_analysis <- function(networks, metaweb) {
  labs <- names(networks)
  meta_species <- length(rownames(metaweb)) + length(colnames(metaweb))
  ml <- which(metaweb > 0, arr.ind = TRUE)
  meta_int <- length(paste(rownames(metaweb)[ml[, 1]], colnames(metaweb)[ml[, 2]]))
  out <- data.frame()
  for (i in seq_along(networks)) {
    exp_c <- calc_contribution(networks[[i]], metaweb)
    obs_c <- calc_contribution(networks[[i]], build_loo_metaweb(networks, i))
    out <- rbind(out, data.frame(
      network = labs[i],
      exp_species = exp_c$n_species,
      exp_species_pct = round(100 * exp_c$n_species / meta_species, 1),
      exp_interactions = exp_c$n_interactions,
      exp_interactions_pct = round(100 * exp_c$n_interactions / meta_int, 1),
      obs_unique_species = obs_c$unique_species,
      obs_unique_species_pct = round(100 * obs_c$unique_species / meta_species, 1),
      obs_unique_interactions = obs_c$unique_interactions,
      obs_unique_interactions_pct = round(100 * obs_c$unique_interactions / meta_int, 1),
      unique_lower = obs_c$unique_lower, unique_higher = obs_c$unique_higher,
      total_lower = nrow(networks[[i]]), total_higher = ncol(networks[[i]]),
      stringsAsFactors = FALSE))
  }
  out$species_efficiency <- round(out$obs_unique_species / out$exp_species, 3)
  out$interaction_efficiency <- round(out$obs_unique_interactions / out$exp_interactions, 3)
  out
}
