# =====================================================================
# rarefaction.R -- species/interaction accumulation curves
# =====================================================================

#' Species accumulation per network (and pooled metaweb)
#'
#' Wraps `vegan::specpool()` / `vegan::specaccum()` for each network, plus
#' the pooled metaweb when supplied.
#'
#' @param networks A named list of matrices.
#' @param metaweb Optional pooled metaweb matrix; included as `"__metaweb__"`
#'   when supplied.
#' @param method Passed to `vegan::specaccum()` (default `"random"`).
#' @return A named list (one entry per network / metaweb) with `pool` and
#'   `accum`, or NULL for any network where accumulation failed.
#' @export
rarefaction_analysis <- function(networks, metaweb = NULL, method = "random") {
  .bd_require("vegan", "rarefaction")
  webs <- networks
  if (!is.null(metaweb)) webs[["__metaweb__"]] <- metaweb
  lapply(webs, function(w) tryCatch(
    list(pool = vegan::specpool(w),
         accum = vegan::specaccum(w, method = method)),
    error = function(e) NULL))
}
