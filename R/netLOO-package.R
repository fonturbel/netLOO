#' @keywords internal
"_PACKAGE"

# NSE column names used inside ggplot2::aes()/dplyr-style pipelines below;
# silences R CMD check "no visible binding for global variable" notes.
utils::globalVariables(c(
  "order", "level", "jaccard", "transition", "type", "count", "change",
  "Group", "Species_Count", "n", ".data", "i"
))
