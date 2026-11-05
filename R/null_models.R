# =====================================================================
# null_models.R -- null webs and observed-vs-null significance tests
# =====================================================================

#' Generate null webs
#'
#' A character `method` is passed to `vegan::nullmodel()`. For binary webs
#' use a fixed-fill method such as `"quasiswap"` (default) or `"curveball"`,
#' which keep every row and column total, and so the number of links, fixed.
#' A numeric `method` (1-5) uses `bipartite::nullmodel()`; those methods
#' generate weighted matrices and are meant for quantitative webs. On a
#' binary web, Patefield (`method = 1`) nulls have fewer links than the
#' observed web, which biases nestedness and modularity z-scores.
#'
#' @param web A bipartite matrix.
#' @param n Number of null webs.
#' @param method See Details.
#' @param burnin,thin Used by sequential vegan methods (e.g. `"curveball"`).
#' @return A list of `n` matrices with the dimnames of `web`.
#' @export
generate_null_webs <- function(web, n, method = "quasiswap",
                               burnin = 1000, thin = 100) {
  if (is.character(method)) {
    .bd_require("vegan", "null models")
    web_in <- if (is_binary_web(web)) web else {
      warning("Binary null model applied to a weighted web; using presence/absence.",
              call. = FALSE)
      (web > 0) * 1
    }
    nm  <- vegan::nullmodel(web_in, method)
    sim <- if (isTRUE(nm$commsim$isSeq))
      stats::simulate(nm, nsim = n, burnin = burnin, thin = thin)
    else stats::simulate(nm, nsim = n)
    lapply(seq_len(dim(sim)[3]), function(i)
      matrix(sim[, , i], nrow(web), ncol(web), dimnames = dimnames(web)))
  } else {
    .bd_require("bipartite", "null models")
    if (is_binary_web(web))
      warning("bipartite::nullmodel(method = ", method, ") produces weighted ",
              "matrices with fewer links than a binary web; prefer a vegan ",
              "binary method such as \"quasiswap\".", call. = FALSE)
    bipartite::nullmodel(web, N = n, method = method)
  }
}

#' Null-model significance test for one network
#'
#' Compares observed connectance, nestedness (NODF) and modularity (Q)
#' against null webs from [generate_null_webs()], optionally in parallel.
#'
#' `p_value` is the two-sided empirical p-value,
#' (1 + 2 * min(#null >= obs, #null <= obs)) / (n + 1), capped at 1.
#' `p_value_z` is the normal approximation from the z-score, kept for
#' comparison with netLOO <= 0.0.3. A metric that the null model holds fixed
#' (connectance under `"quasiswap"`) gets NA z and p and a note.
#'
#' One seed per null web is drawn up front, so results do not depend on the
#' number of cores.
#'
#' @param web A bipartite matrix.
#' @param n_null Number of null replicates.
#' @param n_cores Cores for parallel evaluation (NULL = auto, capped at 14;
#'   1 runs sequentially without a cluster).
#' @param null_method See [generate_null_webs()].
#' @param seed Optional RNG seed.
#' @return A named list (one entry per metric) with `observed`, `null_mean`,
#'   `null_sd`, `z_score`, `p_value`, `p_value_z`, `n_valid`, `note`.
#' @export
enhanced_null_models <- function(web, n_null = 1000, n_cores = NULL,
                                 null_method = "quasiswap", seed = NULL) {
  .bd_require("bipartite", "null models")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(n_cores)) n_cores <- max(1, min(parallel::detectCores() - 1, 14))

  obs <- list(
    connectance = sum(web > 0) / (nrow(web) * ncol(web)),
    nestedness  = tryCatch(unname(bipartite::networklevel(web, index = "NODF")),
                           error = function(e) NA),
    modularity  = tryCatch(bipartite::computeModules(web)@likelihood,
                           error = function(e) NA))

  nulls <- tryCatch(generate_null_webs(web, n_null, null_method),
                    error = function(e) { warning(conditionMessage(e)); NULL })
  if (is.null(nulls)) return(list(error = "null model generation failed"))

  task_seeds <- sample.int(.Machine$integer.max, length(nulls))
  one_null <- function(i) {
    set.seed(task_seeds[i])
    nw <- nulls[[i]]
    data.frame(
      connectance = sum(nw > 0) / (nrow(nw) * ncol(nw)),
      nestedness  = tryCatch(unname(bipartite::networklevel(nw, index = "NODF")),
                             error = function(e) NA),
      modularity  = tryCatch(bipartite::computeModules(nw)@likelihood,
                             error = function(e) NA))
  }

  cat("Null models (", if (is.character(null_method)) null_method else
        paste("bipartite method", null_method), ") on ", n_cores, " core(s)...\n",
      sep = "")
  if (n_cores > 1) {
    .bd_require(c("foreach", "doParallel"), "parallel null models")
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    res <- foreach::foreach(i = seq_along(nulls), .packages = "bipartite",
                            .combine = rbind, .errorhandling = "remove") %dopar%
      one_null(i)
  } else {
    res <- do.call(rbind, lapply(seq_along(nulls), one_null))
  }

  out <- list()
  for (metric in names(obs)) {
    o <- as.numeric(obs[[metric]][1]); nv <- res[[metric]]; nv <- nv[!is.na(nv)]
    if (length(nv) > 10 && !is.na(o)) {
      s <- stats::sd(nv)
      fixed <- s < 1e-12
      z <- if (fixed) NA_real_ else (o - mean(nv)) / s
      p_emp <- if (fixed) NA_real_ else
        min(1, (1 + 2 * min(sum(nv >= o), sum(nv <= o))) / (length(nv) + 1))
      out[[metric]] <- list(observed = o, null_mean = mean(nv), null_sd = s,
                            z_score = z, p_value = p_emp,
                            p_value_z = if (fixed) NA_real_ else 2 * stats::pnorm(-abs(z)),
                            n_valid = length(nv),
                            note = if (fixed) "fixed by null model" else "")
    } else {
      out[[metric]] <- list(observed = o, error = "insufficient valid nulls")
    }
  }
  out
}
