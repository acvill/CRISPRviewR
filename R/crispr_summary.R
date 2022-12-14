#' summarize CRISPR arrays
#'
#' Given an object created with \code{\link{read_minced}}, this function
#' gives the length and number of spacers in each CRISPR array.
#'
#' @param data An object created with \code{\link{read_minced}}.
#' @return A summary tibble.
#' @import dplyr
#' @export

crispr_summary <- function(data) {

  data |> dplyr::group_by(array) |>
    dplyr::mutate(length = max(end) - min(start)) |>
    dplyr::mutate(spacer_count = dplyr::n()) |>
    dplyr::select(array, contig, length, spacer_count) |>
    unique()

}
