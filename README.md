

# CRISPRviewR: an R package for visualizing CRISPRs across environmental datasets
![CRISPRviewR example output](https://user-images.githubusercontent.com/22378512/191805950-8c40f366-4481-4191-b03f-e5b3d82331cb.svg)

## Background

CRISPRviewR uses the output from minCED to associate, compare, and visualize CRISPR arrays across environmental samples. To get a sense for the shape of minCED data, check out the [example files](https://github.com/acvill/CRISPRviewR/tree/master/example_data_minced).  

[![ctSkennerton/minced - GitHub](https://gh-card.dev/repos/ctSkennerton/minced.svg)](https://github.com/ctSkennerton/minced)

This package relies on the functions of other packages for data cleaning and plotting, including the following:
- dplyr, ggplot2, and other components of the [tidyverse](https://www.tidyverse.org/)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [ggpubr](https://github.com/kassambara/ggpubr)
- [ggnewscale](https://github.com/eliocamp/ggnewscale)
- [ggseqlogo](https://github.com/omarwagih/ggseqlogo)

## Installation

Future versions will be available on CRAN or Bioconductor. For now, you can install the development version from GitHub:
```
devtools::install_github("acvill/CRISPRviewR")
```

## Example

Please see the [CRISPRviewR vignette](https://albertvill.com/CRISPRviewR-vignette.html) for a suggested workflow.

## *Caveat emptor*

The CRISPRviewR functions make no assumptions about the completeness of the CRISPR arrays annotated by minCED or the structure of the underlying assembly. 
In that regard, users of CRISPRviewR should be aware of the following possibilities.  
 1. Due to their abundance of direct repeats, CRISPR arrays are often misassembled, particularly in the absence of sufficient coverage.
 2. minCED does not predict orientation of arrays, which requires either identification of *cas* genes or the annotation of a leader sequence. Therefore, CRISPRviewR plots may be backwards with respect to the direction of transcription. If this is a problem, use a tool like [CRISPRleader](https://doi.org/10.1093/bioinformatics/btw454) to get strand orientation, then export your plots and invert manually. 
 3. For fragmented assemblies, CRISPR arrays may occur at contig boundaries.
 4. For time-course metagenomic assemblies, differences in CRISPR array structure through time may be attributed to strain variation as opposed to array expansion / recombination.
 5. Depending on your minCED parameters for the expected lengths of repeats and spacers, portions of repeat sequences may be mistakenly associated with spacers. This results in spacers with stretches of sequences that are more similar than expected. For an example, note the similarity of the first 16 nt of each spacer in the `CRISPR1` array from [example data `s1`](https://github.com/acvill/CRISPRviewR/tree/master/example_data_minced):
```r
> read_minced(txt = url("https://raw.githubusercontent.com/acvill/CRISPRviewR/master/example_data_minced/s1.txt"),
>             gff = url("https://raw.githubusercontent.com/acvill/CRISPRviewR/master/example_data_minced/s1.gff")) |>
>   dplyr::filter(array == "CRISPR1") |> dplyr::select(spacer)

# A tibble: 7 × 1
  spacer                                        
  <chr>                                         
1 CACTTGCTAATACAGCTGTGGTTGAGCCAAACAATGAGATGGTAAT
2 TACTTGCTAATACAGCGCACGCGAGACCTTCACGCGACTAGGACGG
3 TACTTGCTAATACAGCCACGAGCCTCATCACGCGAACTCTCATCAC
4 TACTTGCTAATCCAGCCGAATTATTGCAACGCTTATCCTCGCCTCG
5 CACTTGCTAACACAGCATAAAAACGACGACGACACGACCGACAGGT
6 CACTTGCTAATACAGCTCGGAGGAGTGAAGAATAGCCAGCACCTCG
7 NA 
```

## Bugs and notes

- CRISPRviewR has only been tested with the output from minCED v0.4.2
- If you find a bug or want to suggest a new feature, please [open an issue](https://github.com/acvill/CRISPRviewR/issues/new/choose) or [make a pull request](https://github.com/acvill/CRISPRviewR/pulls).
