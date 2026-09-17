# =====================================================================
# taxonomy.R -- optional external species -> group lookup
# =====================================================================

#' Load a taxonomy lookup table (species -> group)
#'
#' @param cfg A `bd_pipeline_config` (see [pipeline_config()]).
#' @return A data frame with columns `species`, `group`, or NULL if no
#'   taxonomy file was configured / found.
#' @export
load_taxonomy <- function(cfg) {
  if (is.null(cfg$taxonomy) || !file.exists(cfg$taxonomy)) {
    if (!is.null(cfg$taxonomy))
      warning("Taxonomy file not found: ", cfg$taxonomy, " (skipping groups).")
    return(NULL)
  }
  tax <- utils::read.csv(cfg$taxonomy, stringsAsFactors = FALSE)
  names(tax) <- tolower(names(tax))
  if (!all(c("species", "group") %in% names(tax)))
    stop("Taxonomy CSV must have columns 'species' and 'group'.")
  tax[, c("species", "group")]
}

#' Classify species names against a taxonomy lookup
#'
#' Matching is exact first; if `partial = TRUE` it also tries a genus/substring
#' match so a lookup keyed by genus still classifies "Genus species" names.
#' Unmatched species become `"Other"`.
#'
#' @param species  Character vector of species names to classify.
#' @param taxonomy A lookup data frame from [load_taxonomy()], or NULL.
#' @param partial  Also try genus/substring matching for names not found
#'   exactly (default TRUE).
#' @return A character vector of group labels, same length as `species`.
#' @export
classify_species <- function(species, taxonomy, partial = TRUE) {
  if (is.null(taxonomy)) return(rep("Unclassified", length(species)))
  grp <- taxonomy$group[match(species, taxonomy$species)]
  if (partial) {
    unmatched <- which(is.na(grp))
    for (i in unmatched) {
      hit <- which(vapply(taxonomy$species, function(s)
        grepl(s, species[i], fixed = TRUE), logical(1)))
      if (length(hit) > 0) grp[i] <- taxonomy$group[hit[1]]
    }
  }
  grp[is.na(grp)] <- "Other"
  grp
}
