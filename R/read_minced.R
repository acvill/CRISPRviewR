#' Read and transform MinCED data
#'
#' MinCED is a program used to identify CRISPRs in genomes
#' or contigs. By default, MinCED outputs both a human-readable
#' text file and a gff-formatted annotation file. This function
#' parses both the .txt and .gff outputs of MinCED to create a
#' single tibble.
#'
#' @param txt Path to the minced .txt output
#' @param gff Path to the minced .gff output
#' @return A tibble containing the minced data
#' @export

read_minced <- function(txt, gff) {

  seqs <- readr::read_tsv(file = txt,
                          col_names = "mess",
                          comment = "--") |>
    dplyr::filter(!grepl('POSITION', mess)) |>
    dplyr::filter(!grepl('Repeats', mess)) |>
    dplyr::mutate(mess = stringr::str_replace(mess, "CRISPR ", "CRISPR")) |>
    dplyr::mutate(mess = stringr::str_replace(mess, "Range.*", "")) |>
    dplyr::group_by(grp = cumsum(stringr::str_detect(mess, "CRISPR"))) |>
    dplyr::mutate(array = dplyr::first(mess)) |>
    dplyr::filter(!grepl('CRISPR', mess)) |>
    dplyr::ungroup() |>
    tidyr::separate(col = mess, remove = TRUE, sep = "\t", fill = "right",
                    into = c("start","drop1","rep","spacer","drop2")) |>
    dplyr::mutate(array = stringr::str_trim(array)) |>
    dplyr::select(array, start, rep, spacer) |>
    dplyr::mutate(start = as.character(start))

  coords <- readr::read_tsv(file = gff,
                            comment = "##",
                            col_names = c("contig","drop1",
                                          "drop2","start","end","drop3",
                                          "drop4","drop5","parse")) |>
    dplyr::select(!dplyr::contains("drop")) |>
    dplyr::filter(!grepl("rpt_type", parse)) |>
    tidyr::separate(col = parse, sep = ";",
                    into = c("array","id"), remove = TRUE) |>
    dplyr::mutate(array = stringr::str_remove(pattern = "Parent=",
                                              string = array)) |>
    dplyr::select(array, contig, start, end) |>
    dplyr::mutate(start = as.character(start))

  dplyr::left_join(seqs, coords, by = c("array", "start")) |>
    dplyr::mutate(start = as.numeric(start)) |>
    dplyr::select(array, rep, start, end, spacer, contig)

}
