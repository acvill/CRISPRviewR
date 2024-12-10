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
#' @param fix_repeats Attempts to resolve portions of repeats that were
#' mistakenly assigned to spacers. Default = FALSE.
#' (see https://github.com/ctSkennerton/minced/issues/36)
#' @param window Nucleotide context in each direction
#' used for the rolling average when fix_repeats = TRUE. Default = 2.
#' @param cutoff Threshold value for rolling average conservation
#' across repeats used to determine repeat-spacer break
#' when fix_repeats = TRUE. Default = 0.7.
#' @return A tibble containing the minced data
#' @import dplyr
#' @importFrom readr read_tsv
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stringr str_trim
#' @importFrom stringr str_remove
#' @importFrom tidyr separate
#' @importFrom Biostrings consensusMatrix
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr as_vector
#' @importFrom stats na.omit
#' @export

read_minced <- function(txt,
                        gff,
                        fix_repeats = FALSE,
                        window = 2,
                        cutoff = 0.7) {

  ####################################################
  ## helper functions related to fix_repeats option ##
  ####################################################

  # get average in rolling window
  roll <- function(vector, window1) {
    vc <- (Biostrings::consensusMatrix(vector)[1:4,]/(length(vector))) |>
      t() |> apply(1, max)
    vl <- length(vc)
    sapply(seq_along(vc),
           function(v, w, i) {
             lt <- ifelse(i - w < 1, 1, i - w)
             rt <- ifelse(i + w > vl, vl, i + w)
             mean(v[lt:rt])
           },
           v = vc,
           w = window1)
  }

  # reverse all strings in character vector
  strReverse <- function(x) {
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }

  checkSpacers <- function(dat, window2, cutoff1) {
    # get unique array IDs
    arrays <- dat |> dplyr::select(array) |>
      unique() |> purrr::as_vector()
    out <- dat
    # loop through arrays if number of spacers >= 3 (+1 NA)
    for (arr in arrays) {
      if (nrow(dplyr::filter(dat, array == arr)) > 3) {
        # get the rolling average conservation score for fwd orientation
        spc_fwd <- dat |> dplyr::filter(array == arr) |>
          dplyr::select(spacer) |> stats::na.omit() |>
          purrr::as_vector() |> Biostrings::DNAStringSet() |>
          roll(window1 = window2)
        # get the rolling average conservation score for rev orientation
        spc_rev <- dat |> dplyr::filter(array == arr) |>
          dplyr::select(spacer) |> stats::na.omit() |>
          purrr::as_vector() |> strReverse() |>
          Biostrings::DNAStringSet() |>
          roll(window1 = window2)
        # read rolling average vectors from left
        # when value falls below cutoff, record
        ifwd <- 0
        while (spc_fwd[ifwd+1] >= cutoff1) {
          ifwd <- ifwd + 1
        }
        irev <- 0
        while (spc_rev[irev+1] >= cutoff1) {
          irev <- irev + 1
        }
        # if spacers have more conservation than expected (cutoff)
        # then append conserved sequence to repeat sequence
        if (ifwd > 0 & irev == 0) {
          out <- out |>
            dplyr::filter(array == arr) |>
            dplyr::mutate(spacer = ifelse(test = is.na(spacer),
                                          yes = "",
                                          no = spacer)) |>
            dplyr::mutate(rep = paste0(rep,
                                       substr(x = spacer,
                                              start = 1,
                                              stop = ifwd)),
                          spacer = substr(x = spacer,
                                          start = ifwd + 1,
                                          stop = nchar(spacer))) |>
            dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(.x,""))) |>
            dplyr::mutate(end = ifelse(test = is.na(spacer),
                                       yes = end,
                                       no = end + ifwd)) |>
            rbind(out |> dplyr::filter(array != arr))
        } else if (ifwd == 0 & irev > 0) {
          out <- out |>
            dplyr::filter(array == arr) |>
            dplyr::mutate(spacer = dplyr::lag(spacer)) |>
            dplyr::mutate(spacer = ifelse(is.na(spacer), "", spacer)) |>
            dplyr::mutate(rep = paste0(substr(x = spacer,
                                              start = nchar(spacer) - irev,
                                              stop = nchar(spacer)),
                                       rep),
                          spacer = substr(x = spacer,
                                          start = 1,
                                          stop = nchar(spacer) - irev - 1)) |>
            dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(.x,""))) |>
            dplyr::mutate(start = ifelse(test = is.na(spacer),
                                         yes = start ,
                                         no = start - irev)) |>
            dplyr::mutate(spacer = dplyr::lead(spacer)) |>
            rbind(out |> dplyr::filter(array != arr))
        }

      }

    }

    out |> dplyr::arrange(as.numeric(gsub(x = array,
                                          pattern = "CRISPR",
                                          replacement = "")))

  }

  ####################################################

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

  merged <- dplyr::left_join(seqs, coords, by = c("array", "start")) |>
    dplyr::mutate(start = as.numeric(start)) |>
    dplyr::select(array, rep, start, end, spacer, contig)

  if (isTRUE(fix_repeats)) {

    merged <- checkSpacers(dat = merged,
                           window2 = window,
                           cutoff1 = cutoff)

  }

  merged

}
