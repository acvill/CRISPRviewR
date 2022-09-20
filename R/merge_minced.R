#' merge read_minced() objects into a single tibble
#'
#' Given multiple objects from read_minced(), this function merges
#' them into a single tibble and creates an ID column from 'names'.
#' If a vector of names is not supplied, the object names are used as labels.
#'
#' @param ... One or more objects created with read_minced().
#' @param names A vector of unique names  to label merged objects. If no names are given, the object names are used.
#' @return A tibble containing the merged data.
#' @importFrom tibble lst
#' @importFrom stats setNames
#' @importFrom dplyr bind_rows
#' @export

merge_minced <- function(..., names = NULL) {

  mlist <- tibble::lst(...)
  if (!is.null(names)) {
    if (length(mlist) == length(names)) {
      merged <- mlist |>
        stats::setNames(nm = names) |>
        dplyr::bind_rows(.id = 'id')
    } else {
      message("Error: names vector and object list have different lengths")
    }
  } else {
    merged <- mlist |>
      dplyr::bind_rows(.id = 'id')
  }
  merged

}
