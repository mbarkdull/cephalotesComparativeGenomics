Reference-guided genome assembly and comparative genomics workflow
================
Megan Barkdull

## Introduction

The goal of this project is to sequence five new genomes from the ant
genus *Cephalotes*, assemble them to the reference *C. varians* genome,
and then use these data to better understand the genetic basis of caste
polymorphism.

## Reference-guided genome assembly

Reference-guided genome assembly has the following steps:

1.  Trim raw reads with Trimmomatic.
2.  Map the trimmed reads back to the reference genome (I’ll probably
    use NextGenMap for this).
3.  Call SNPs and indels with SAMtools.
4.  Filter the calls by quality.
5.  Assemble a consensus sequence.
