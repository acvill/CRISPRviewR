# CRISPRviewR

#### An R package for parsing and visualizing minCED output
[![ctSkennerton/minced - GitHub](https://gh-card.dev/repos/ctSkennerton/minced.svg)](https://github.com/ctSkennerton/minced)

## Installation

Future versions will be available on CRAN or Bioconductor. For now, you can install the development version from GitHub:
```
devtools::install_github("acvill/CRISPRviewR")
```

## Example workflow
 

## *Caveat emptor*

The CRISPRviewR functions make no assumptions about the completeness of the CRISPR arrays annotated by minCED or the structure of the underlying assembly. 
In that regard, users of CRISPRviewR should be aware of the following possibilities.  
 1. Due to their abundance of direct repeats, CRISPR arrays are often misassembled, particularly in the absence of sufficient coverage.
 2. minCED does not predict orientation of arrays, which requires either identification of *cas* genes or the annotation of a leader sequence. Therefore, CRISPRviewR plots may be backwards with respect to the direction of transcription. 
 3. For fragmented assemblies, CRISPR arrays may occur at contig boundaries.
 4. For time-course metagenomic assemblies, differnces in CRISPR array structure through time may be attributed to strain variation as opposed to array expansion / recombination.

## Notes

- CRISPRviewR has only been tested with the output from minCED v0.4.2
- future implementations will be submitted to CRAN or Bioconductor after I handle some [annoying notes](https://stackoverflow.com/q/9439256/7976890) from `R CMD check`
