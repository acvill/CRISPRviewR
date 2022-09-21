# CRISPRviewR

## Description

## Workflow

## *Caveat emptor*

The CRISPRviewR functions make no assumptions about the completeness of the CRISPR arrays annotated by minCED or the structure of the underlying assembly. 
In that regard, users of CRISPRviewR should be aware of the following possibilities.  
 1. Due to their abundance of direct repeats, CRISPR arrays are often misassembled, particularly in the absence of sufficient coverage.
 2. minCED does not predict orientation of arrays, which requires either identification of *cas* genes or the annotation of a leader sequence. Therefore, CRISPRviewR plots may be backwards.
 3. For fragmented assemblies, CRISPR arrays may occur at contig boundaries.
 4. For time-course metagenomic assemblies, differnces in CRISPR array structure through time may be attributed to strain variation as opposed to array expansion / recombination.
