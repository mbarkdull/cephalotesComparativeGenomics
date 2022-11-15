Reference-guided genome assembly and comparative genomics workflow
================
Megan Barkdull

## Introduction

The goal of this project is to sequence five new genomes from the ant
genus *Cephalotes*, assemble them to the reference *C. varians* genome,
and then use these data to better understand the genetic basis of caste
polymorphism.

## Citing

This workflow implements many different tools, all of which should be
cited. I have noted appropriate citations throughout the text of this
ReadMe. Please properly credit these researchers’ hard work\!

In addition, if this workflow has been useful to you, please cite the
Github
repository:

[![DOI](https://zenodo.org/badge/520512524.svg)](https://zenodo.org/badge/latestdoi/520512524)

## Reference-guided genome assembly

Reference-guided genome assembly has the following steps:

1.  Trim raw reads with Trimmomatic.
2.  Map the trimmed reads back to the reference genome (I’ll probably
    use NextGenMap for this).
3.  Call SNPs and indels with SAMtools.
4.  Filter the calls by quality.
5.  Assemble a consensus sequence.

These steps are all implemented by the scripts
`./Scripts/referenceGuidedAssembly`, which requires two arguments:
first, a list of samples, and two, the number of threads to use for
multithreaded steps.

Note that this script currently has a lot of hard-coded paths, including
the paths to each program. The location and names of your raw sequencing
files will also be different than those in this script, so update
accordingly.

If you use this script, please cite:

  - Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a
    flexible trimmer for Illumina sequence data, Bioinformatics, Volume
    30, Issue 15, 1 August 2014, Pages 2114–2120,
    <https://doi.org/10.1093/bioinformatics/btu170>
  - Fritz J. Sedlazeck, Philipp Rescheneder, Arndt von Haeseler,
    NextGenMap: fast and accurate read mapping in highly polymorphic
    genomes, Bioinformatics, Volume 29, Issue 21, 1 November 2013, Pages
    2790–2791, <https://doi.org/10.1093/bioinformatics/btt468>
  - Twelve years of SAMtools and BCFtools. Petr Danecek, James K
    Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O
    Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M
    Davies, Heng Li. GigaScience, Volume 10, Issue 2, February 2021,
    giab008, <https://doi.org/10.1093/gigascience/giab008>

## Processing assembled genomes

Once genomes are assembled, they require some processing before genes
are identified.

Repetitive regions of the genome must be soft-masked prior to
annotation. To do this, we must first predict/identify transposable
element sequences in the genomes, which we will do with RepeatModeler2,
implemented in the Bash script `./Scripts/handlingRepeats`.
