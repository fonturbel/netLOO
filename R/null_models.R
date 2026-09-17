# =====================================================================
# null_models.R -- parallel z-score / p-value null models
# =====================================================================

#' Parallel null-model significance test for one network
#'
#' Compares observed connectance, nestedness and modularity against a null
#' distribution generated with `bipartite::nullmodel()`.
#'
#' @param web A bipartite matrix.
#' @param n_null Number of null replicates.
#' @param n_cores Cores for parallel evaluation (NULL = auto, capped at 14).
#' @param null_method `bipartite::nullmodel()` method id (1 = non-sequential).
#' @return A named list (one entry per metric) with `observed`, `null_mean`,
#'   `null_sd`, `z_score`, `p_value`, `n_valid`.
#' @export
enhanced_null_models <- function(web, n_null = 1000, n_cores = NULL,
                                 null_method = 1) {
  .bd_require(c("bipartite", "parallel", "foreach", "doParallel"), "null models")
  if (is.null(n_cores)) n_cores <- min(parallel::detectCores() - 1, 14)
  cat("Null models on", n_cores, "cores...\n")

  obs <- list(
    connectance = tryCatch(bipartite::networklevel(web, index = "connectance"),
                           error = function(e) sum(web > 0) / (nrow(web) * ncol(web))),
    nestedness  = tryCatch(bipartite::networklevel(web, index = "NODF"),
                           error = function(e) NA),
    modularity  = tryCatch(bipartite::computeModules(web)@likelihood,
                           error = function(e) NA))

  nulls <- tryCatch(bipartite::nullmodel(web, N = n_null, method = null_method),
                    error = function(e) NULL)
  if (is.null(nulls)) return(list(error = "null model generation failed"))

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  `%dopar%` <- foreach::`%dopar%`

  res <- foreach::foreach(i = seq_len(n_null), .packages = "bipartite",
                          .combine = rbind, .errorhandling = "pass") %dopar% {
    nw <- nulls[[i]]
    data.frame(
      connectance = sum(nw > 0) / (nrow(nw) * ncol(nw)),
      nestedness  = tryCatch(bipartite::networklevel(nw, index = "NODF"),
                             error = function(e) NA),
      modularity  = tryCatch(bipartite::computeModules(nw)@likelihood,
                             error = function(e) NA))
  }

  out <- list()
  for (metric in names(obs)) {
    o <- as.numeric(obs[[metric]][1]); nv <- res[[metric]]; nv <- nv[!is.na(nv)]
    if (length(nv) > 10 && !is.na(o)) {
      z <- (o - mean(nv)) / stats::sd(nv)
      out[[metric]] <- list(observed = o, null_mean = mean(nv),
                            null_sd = stats::sd(nv), z_score = z,
                            p_value = 2 * stats::pnorm(-abs(z)), n_valid = length(nv))
    } else {
      out[[metric]] <- list(observed = o, error = "insufficient valid nulls")
    }
  }
  out
}
