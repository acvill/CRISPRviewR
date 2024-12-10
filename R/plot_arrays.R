#' Given a set of grouped CRISPR arrays, plot them
#'
#' generates a graphical representation of array sets associated by shared repeat consensus sequences
#'
#' @param group A subset of the object created by \code{\link{group_arrays}} containing only a single consensus repeat.
#' @param palette A vector of hexadecimal color codes used to distinguish unique spacers.
#' Default is \code{grDevices::palette.colors(palette = "Polychrome 36")[-c(1:2)]}.
#' @param cdist A numeric scalar in \code{[0,1]} representing a proportion of the maximum Levenshtein distance
#' for a set of spacers. Default = 0.1.
#' @param plot_logo If FALSE, do not plot the repeat sequence logo. Default = TRUE.
#' @param number_spacers If TRUE, include cluster number label in each spacer. Default = FALSE.
#' @return A plot comparing CRISPR arrays with a shared consensus repeat.
#' @import dplyr
#' @import ggplot2
#' @importFrom S4Vectors isSingleNumber
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @importFrom Biostrings DNA_BASES
#' @importFrom Biostrings DNAStringSet
#' @importFrom BiocGenerics score
#' @importFrom pwalign stringDist
#' @importFrom pwalign pairwiseAlignment
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom grDevices rgb
#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom ggpubr ggarrange
#' @export

plot_arrays <- function(group,
                        palette = NULL,
                        cdist = 0.1,
                        plot_logo = TRUE,
                        number_spacers = FALSE) {

  # helper function to create asymmetric nucleotide substitution matrix
  # see https://github.com/Bioconductor/Biostrings/pull/77
  asymNSM <- function(match = 1, mismatch = 0,
                      baseOnly = FALSE, type = "DNA",
                      symmetric = TRUE) {
    "%safemult%" <- function(x, y) ifelse(is.infinite(x) & y == 0, 0, x * y)
    type <- match.arg(type, c("DNA", "RNA"))
    if (!S4Vectors::isSingleNumber(match) || !S4Vectors::isSingleNumber(mismatch)) {
      stop("'match' and 'mismatch' must be non-missing numbers")
    }
    if (baseOnly) {
      letters <- Biostrings::IUPAC_CODE_MAP[Biostrings::DNA_BASES]
    }
    else {
      letters <- Biostrings::IUPAC_CODE_MAP
    }
    if (type == "RNA") {
      names(letters) <- chartr("T", "U", names(letters))
    }
    nLetters <- length(letters)
    splitLetters <- strsplit(letters, split = "")
    submat <- matrix(0, nrow = nLetters, ncol = nLetters, dimnames = list(names(letters), names(letters)))
    if (symmetric) {
      for (i in 1:nLetters) {
        for (j in i:nLetters) {
          submat[i,j] <- submat[j,i] <- base::mean(outer(splitLetters[[i]], splitLetters[[j]], "=="))
        }
      }
    }
    else {
      for (i in 1:nLetters) {
        for (j in i:nLetters) {
          submat[i,j] <- base::mean(outer(splitLetters[[i]], splitLetters[[j]], "%in%"))
          submat[j,i] <- base::mean(outer(splitLetters[[j]], splitLetters[[i]], "%in%"))
        }
      }
    }
    abs(match) * submat - abs(mismatch) %safemult% (1 - submat)
  }

  # if adjusted pairwise levenshtein distance is less than cutoff ratio,
  # put spacers in same cluster
  uniq_spacers <- unique(na.omit(group$spacer))
  min_sp_len <- min(nchar(uniq_spacers))
  repcon_len <- nchar(unique(group$repcon))
  dmat <- uniq_spacers |>
    Biostrings::DNAStringSet() |>
    pwalign::stringDist(method = "levenshtein") / min_sp_len
  clusters <- data.frame(spacer = uniq_spacers,
                         cluster = stats::hclust(dmat, method = "complete") |>
                           stats::cutree(h = cdist))
  group <- dplyr::left_join(group, clusters)

  # get distance of each repeat to the consensus
  subst <- asymNSM()

  ## can't use asymmetric matrix with stringDist, use pairwiseAlignment instead
  group <- dplyr::mutate(group,
                         rep_dist = pwalign::pairwiseAlignment(pattern = rep,
                                                               subject = repcon,
                                                               substitutionMatrix = subst) |>
                           BiocGenerics::score() |> abs()
  ) |>
    dplyr::mutate(rep2con_sim = (rep_dist / nchar(repcon)) |> round(digits = 2))

  # get grayscale color for repeat, based on similarity to consensus repeat
  ## uses 0.75 similarity as floor, unless lowest similarity < 0.75, then use that
  minsim <- min(group$rep2con_sim)
  if (minsim < 0.75) {
    floor_ <- minsim
  } else {
    floor_ <- 0.75
  }

  group <- dplyr::mutate(group, temp = (1 - rep2con_sim) / (1 - floor_)) |>
    dplyr::mutate(r_color = grDevices::rgb(temp,temp,temp)) |>
    dplyr::select(id, rep, spacer, cluster,
                  r_color, rep2con_sim)

  # get palette for spacers
  ## if no palette given, generate default
  nclusters <- group$cluster |> na.omit() |> max()
  if (is.null(palette)) {
    default_palette <- c("#F6222E","#FE00FA","#16FF32","#3283FE",
                         "#FEAF16","#B00068","#1CFFCE","#90AD1C",
                         "#2ED9FF","#DEA0FD","#AA0DFE","#F8A19F",
                         "#325A9B","#C4451C","#1C8356","#85660D",
                         "#B10DA1","#FBE426","#1CBE4F","#FA0087",
                         "#FC1CBF","#F7E1A0","#C075A6","#782AB6",
                         "#AAF400","#BDCDFF","#822E1C","#B5EFB5",
                         "#7ED7D1","#1C7F93","#D85FF7","#683B79",
                         "#66B0FF","#3B00FB")
    newpal <- default_palette |>
      {\(x) rep_len(x, nclusters)}()
  } else {
    newpal <- palette |>
      {\(x) rep_len(x, nclusters)}()
  }

  # number samples and spacers for plotting order
  group <- dplyr::left_join(group,
                            data.frame(id = group$id |> unique(),
                                       y_order = seq(1:length(group$id |> unique())))) |>
    dplyr::group_by(id) |>
    dplyr::mutate(x_order = seq(1:dplyr::n())) |>
    dplyr::ungroup()

  # generate x and y coords for repeats and spacers
  ymax <- max(group$y_order)
  group <- group |>
    dplyr::mutate(xfactor = (2 * x_order) - 1) |>
    dplyr::mutate(yfactor = (ymax - y_order) + 1) |>
    dplyr::mutate(xrep = purrr::map(base::strsplit(paste(xfactor,
                                                         xfactor + 0.5,
                                                         xfactor,
                                                         xfactor - 0.5,
                                                         sep = ","),
                                                   split = ","),
                                    as.numeric)) |>
    dplyr::mutate(yrep = purrr::map(base::strsplit(paste(yfactor + 0.4,
                                                         yfactor,
                                                         yfactor - 0.4,
                                                         yfactor,
                                                         sep = ","),
                                                   split = ","),
                                    as.numeric)) |>
    dplyr::mutate(xspa = purrr::map(base::strsplit(paste(xfactor + 1.5,
                                                         xfactor + 1.5,
                                                         xfactor + 0.5,
                                                         xfactor + 0.5,
                                                         sep = ","),
                                                   split = ","),
                                    as.numeric)) |>
    dplyr::mutate(yspa = purrr::map(base::strsplit(paste(yfactor + 0.4,
                                                         yfactor - 0.4,
                                                         yfactor - 0.4,
                                                         yfactor + 0.4,
                                                         sep = ","),
                                                   split = ","),
                                    as.numeric)) |>
    tidyr::unnest(c(xrep,yrep,xspa,yspa)) |>
    dplyr::group_by(x_order, y_order) |>
    dplyr::mutate(plot_id = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    dplyr::mutate(xspa = ifelse(is.na(spacer), NA, xspa),
                  yspa = ifelse(is.na(spacer), NA, yspa),
                  xfactor = ifelse(is.na(spacer), NA, xfactor),
                  yfactor = ifelse(is.na(spacer), NA, yfactor))

  pcr <- ggplot2::ggplot(group) +
    ggplot2::geom_polygon(mapping = ggplot2::aes(x = xspa,
                                                 y = yspa,
                                                 group = plot_id,
                                                 fill = as.character(cluster)),
                          color = "black") +
    ggplot2::scale_fill_manual(values = newpal,
                               guide = "none") +
    {if (isTRUE(number_spacers)) {list(
      ggplot2::geom_point(mapping = ggplot2::aes(x = xfactor + 1,
                                                 y = yfactor),
                          size = 7,
                          pch = 19,
                          fill = "black"),
      ggplot2::geom_text(mapping = ggplot2::aes(x = xfactor + 1,
                                                y = yfactor,
                                                label = as.character(cluster)),
                         check_overlap = TRUE,
                         size = 4,
                         color = "white")
    )}} +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_polygon(mapping = ggplot2::aes(x = xrep,
                                                 y = yrep,
                                                 group = plot_id,
                                                 fill = rep2con_sim),
                          color = "black") +
    ggplot2::scale_fill_gradient(low = '#FFFFFF',
                                 high = '#000000',
                                 name = "similarity to \nrepeat consensus",
                                 breaks = c(floor_, (1 + floor_)/2, 1),
                                 labels = c(floor_, (1 + floor_)/2, 1),
                                 limits = c(floor_, 1)) +
    ggplot2::geom_text(mapping = ggplot2::aes(y = yfactor,
                                              x = -0.5,
                                              label = id),
                       check_overlap = TRUE,
                       size = 7,
                       hjust = 1) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 10)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "right",
                                                   title.vjust = 1,
                                                   frame.colour = "black",
                                                   ticks.colour = "black"))

  if (isTRUE(plot_logo)) {

    rlg <- ggseqlogo::ggseqlogo(group$rep, method = "prob") +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 10)) +
      ggplot2::scale_x_continuous(breaks = seq(5, repcon_len, by = 5))

    h_ <- length(unique(group$id)) / 2.5

    cra <- ggpubr::ggarrange(pcr, ggplot2::ggplot() + ggplot2::theme_void(), rlg,
                             heights = c(h_, h_/10, 0.75),
                             nrow = 3)

  } else {

    cra <- pcr

  }

  cra

}

