#' group CRISPR arrays from a \code{\link{merge_minced}} object by repeat similarity
#'
#' Given a dataset created with \code{\link{merge_minced}}, this function attempts to
#' identify CRISPR arrays across individual profiles / samples that have identical
#' repeat consensus sequences. As the assembly orientation of a CRISPR array
#' may be different across samples, this function compares reverse complements, too.
#'
#' @param data A merged tibble dataset created with \code{\link{merge_minced}}.
#' @return A tibble of CRISPR arrays grouped by repeat consensus sequence.
#' @import dplyr
#' @importFrom purrr as_vector
#' @importFrom tidyr replace_na
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings consensusString
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @export

group_arrays <- function(data) {

  # helper functions to get consensus and reverse complement
  consensus <- function(tib) {
    purrr::as_vector(tib) |>
      Biostrings::DNAStringSet() |>
      Biostrings::consensusString()
  }

  consensus_revcomp <- function(tib) {
    consensus(tib) |>
      Biostrings::DNAString() |>
      Biostrings::reverseComplement() |>
      as.character()
  }

  # generate consensus for each sample + array group
  # replace NA with characters in DNA_ALPHABET
  # alphabetically compare rep_consensus and rep_consensus_revcomp, keep one
  cdat <- data |>
    dplyr::group_by(id, array) |>
    dplyr::mutate(rep_consensus = consensus(rep)) |>
    dplyr::mutate(rep_consensus_revcomp = consensus_revcomp(rep)) |>
    tidyr::replace_na(list(spacer = ".....")) |>
    dplyr::mutate(flip = dplyr::case_when(rep_consensus > rep_consensus_revcomp ~ "yes",
                                          rep_consensus <= rep_consensus_revcomp ~ "no")) |>
    dplyr::mutate(rep = ifelse(flip == "yes",
                               yes = Biostrings::DNAStringSet(rep) |>
                                 Biostrings::reverseComplement() |>
                                 as.character(),
                               no = rep)) |>
    dplyr::mutate(spacer = ifelse(flip == "yes",
                                  yes = Biostrings::DNAStringSet(spacer) |>
                                    Biostrings::reverseComplement() |>
                                    as.character(),
                                  no = spacer)) |>
    dplyr::mutate(repcon = ifelse(flip == "yes",
                                  yes = rep_consensus_revcomp,
                                  no = rep_consensus)) |>
    dplyr::mutate(grp_order = ifelse(flip == "yes",
                                     yes = as.integer(1 + dplyr::n() - dplyr::row_number()),
                                     no = dplyr::row_number())) |>
    dplyr::arrange(grp_order, .by_group = TRUE) |>
    dplyr::mutate(spacer =
                  if (dplyr::first(spacer) == ".....") {
                    spacer[dplyr::row_number() %% dplyr::n() + 1]
                  }
                  else {
                    spacer
                  }) |>
    dplyr::mutate(spacer = gsub(pattern = "\\.\\.\\.\\.\\.",
                                replacement = NA,
                                x = spacer)) |>
    dplyr::ungroup() |>
    dplyr::select(-c(rep_consensus, rep_consensus_revcomp)) |>
    dplyr::arrange(repcon) |>
    dplyr::group_by(repcon) |>
    dplyr::mutate(repcon_grp = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    dplyr::select(-c(flip, grp_order))

  cdat

}
