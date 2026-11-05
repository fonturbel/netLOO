# =====================================================================
# data.R -- loading bipartite matrices and building the pooled metaweb
# =====================================================================

#' Read a single bipartite matrix
#'
#' Rows = lower trophic level, columns = higher trophic level.
#'
#' @param file Path to a delimited matrix file.
#' @param sep  Field separator (default ";").
#' Row and column names are trimmed of surrounding/duplicated whitespace, and
#' a warning is raised if that leaves duplicated names.
#'
#' @return A numeric matrix, or NULL (with a warning) if `file` doesn't exist.
#' @export
read_network_matrix <- function(file, sep = ";") {
  if (!file.exists(file)) {
    warning("File not found: ", file)
    return(NULL)
  }
  m <- as.matrix(utils::read.table(file, header = TRUE, sep = sep,
                            row.names = 1, na.strings = "NA",
                            check.names = FALSE))
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0
  # Stray whitespace silently splits one species into two across files
  rownames(m) <- trimws(gsub("\\s+", " ", rownames(m)))
  colnames(m) <- trimws(gsub("\\s+", " ", colnames(m)))
  for (lv in list(list(rownames(m), "row"), list(colnames(m), "column"))) {
    d <- unique(lv[[1]][duplicated(lv[[1]])])
    if (length(d))
      warning(sprintf("%s: duplicated %s name(s) after trimming: %s",
                      basename(file), lv[[2]], paste(d, collapse = ", ")),
              call. = FALSE)
  }
  m
}

#' Load all networks named by time-point label, in order
#'
#' Returns a named list of matrices. Order is preserved from the manifest
#' (or alphabetical filename order when auto-discovering).
#'
#' @param cfg A `bd_pipeline_config` (see [pipeline_config()]).
#' @return A named list of matrices, ordered by time point.
#' @export
load_networks <- function(cfg) {
  if (!is.null(cfg$manifest) && file.exists(cfg$manifest)) {
    man <- utils::read.csv(cfg$manifest, stringsAsFactors = FALSE)
    if (!all(c("label", "file") %in% names(man)))
      stop("Manifest must have columns 'label' and 'file'.")
    files  <- ifelse(file.exists(man$file), man$file,
                     file.path(cfg$data_dir, man$file))
    labels <- as.character(man$label)
  } else {
    files  <- sort(list.files(cfg$data_dir, pattern = cfg$file_pattern,
                              full.names = TRUE))
    labels <- tools::file_path_sans_ext(basename(files))
  }

  nets <- stats::setNames(
    lapply(seq_along(files), function(i) {
      m <- read_network_matrix(files[i], cfg$sep)
      if (!is.null(m))
        cat(sprintf("Loaded %-12s: %d %s x %d %s\n", labels[i],
                    nrow(m), cfg$lower_name, ncol(m), cfg$higher_name))
      m
    }), labels)

  nets <- nets[!vapply(nets, is.null, logical(1))]
  if (length(nets) < 3)
    warning("Fewer than 3 networks loaded; temporal trends need >= 3 points.")
  nets
}

#' Build the union metaweb from a list of networks
#'
#' @param networks A named list of matrices (e.g. from [load_networks()]).
#' @param binary TRUE (default) returns presence/absence; FALSE returns, for
#'   each pair, the number of networks in which it occurs.
#' @return A matrix: the union metaweb.
#' @export
build_metaweb <- function(networks, binary = TRUE) {
  all_rows <- unique(unlist(lapply(networks, rownames)))
  all_cols <- unique(unlist(lapply(networks, colnames)))
  meta <- matrix(0, length(all_rows), length(all_cols),
                 dimnames = list(all_rows, all_cols))
  for (net in networks) {
    meta[rownames(net), colnames(net)] <-
      meta[rownames(net), colnames(net)] + (net > 0)
  }
  if (binary) meta <- (meta > 0) * 1
  meta
}

#' Compare a supplied metaweb with the union of the networks
#'
#' A metaweb cleaned separately from the per-period matrices often ends up
#' with different species names or links. This lists the differences.
#'
#' @param metaweb A pooled metaweb matrix.
#' @param networks A named list of matrices.
#' @return NULL when they agree; otherwise a list with the species found only
#'   in the networks or only in the metaweb (per level) and the number of
#'   mismatched cells among shared species.
#' @export
check_metaweb <- function(metaweb, networks) {
  u <- build_metaweb(networks)
  diff <- list(
    lower_only_in_networks  = setdiff(rownames(u), rownames(metaweb)),
    lower_only_in_metaweb   = setdiff(rownames(metaweb), rownames(u)),
    higher_only_in_networks = setdiff(colnames(u), colnames(metaweb)),
    higher_only_in_metaweb  = setdiff(colnames(metaweb), colnames(u)))
  r <- intersect(rownames(u), rownames(metaweb))
  k <- intersect(colnames(u), colnames(metaweb))
  diff$mismatched_cells <- sum((u[r, k, drop = FALSE] > 0) !=
                               (metaweb[r, k, drop = FALSE] > 0))
  if (all(lengths(diff[1:4]) == 0) && diff$mismatched_cells == 0) NULL else diff
}
