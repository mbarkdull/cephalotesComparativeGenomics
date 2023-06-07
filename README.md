Reference-guided genome assembly and comparative genomics workflow
================
Megan Barkdull

## Introduction

The goal of this project is to sequence five new genomes from the ant
genus *Cephalotes*, assemble them to the reference *C. varians* genome,
and then use these data to better understand the genetic basis of caste
polymorphism.

The general workflow of this project is as
follows:

<img src="README_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

## Citing

This workflow implements many different tools, all of which should be
cited. I have noted appropriate citations throughout the text of this
ReadMe. Please properly credit these researchers’ hard work\!

In addition, if this workflow has been useful to you, please cite the
Github
repository:

[![DOI](https://zenodo.org/badge/520512524.svg)](https://zenodo.org/badge/latestdoi/520512524)

## Reference-guided genome assembly

Reference-guided genome assembly involves mapping raw reads to a
reference genome and generating a consensus sequence. We will do this
iteratively using the tool
[pseudo-it](https://github.com/goodest-goodlab/pseudo-it#citation),
implemented in the bash script `./Scripts/runPseudoit`. This script
requires two arguments:

1.  The full path to your reference genome fasta file
2.  The full path to a directory containing directories for each of the
    newly sequenced genomes you wish to assemble.

The script will pre-index your reference genome with samtools, BWA, and
Picard, then run pseudo-it. You will get two important outputs, (1) a
final genome assembly and (2) a chain file that can be used to lift over
annotations from the reference genome to your new genome.

If you use this script, please cite pseudo-it and all of the tools it
implements:

  - Sarver BAJ, Keeble S, Cosart T, Tucker PK, Dead MD, Good JM. 2017.
    Phylogenomic insights into mouse evolution using a pseudoreference
    approach. Genome Biology and Evololution.
    <https://doi.org/10.1093/gbe/evx034>.
  - Li H. (2013) Aligning sequence reads, clone sequences and assembly
    contigs with BWA-MEM. arXiv:1303.3997v2 \[q-bio.GN\].
  - Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V.,
    Pollard, M. O., … & Li, H. (2021). Twelve years of SAMtools and
    BCFtools. Gigascience, 10(2), giab008.
  - McKenna, Aaron, et al. “The Genome Analysis Toolkit: a MapReduce
    framework for analyzing next-generation DNA sequencing data.” Genome
    research 20.9 (2010): 1297-1303.
  - Quinlan, Aaron R., and Ira M. Hall. “BEDTools: a flexible suite of
    utilities for comparing genomic features.” Bioinformatics 26.6
    (2010): 841-842.
  - Li, Heng. “Tabix: fast retrieval of sequence features from generic
    TAB-delimited files.” Bioinformatics 27.5 (2011): 718-719.

## Processing assembled genomes

Once genomes are assembled to the reference, we want to identify gene
sequences in those new genomes. To do this, we will lift over the
reference annotation onto the new genomes, and extract transcript
sequences. You can do this using the R script `./Scripts/gffread/`,
which uses liftOver to lift over the annotation, and gffread to extract
transcript sequences.

If you use this script, please cite:

  - Hinrichs, Angela S., et al. “The UCSC genome browser database:
    update 2006.” Nucleic acids research 34.suppl\_1 (2006): D590-D598.
  - Pertea, Geo, and Mihaela Pertea. “GFF utilities: GffRead and
    GffCompare.” F1000Research 9 (2020): 304.

## Identifying evolutionary changes associated with phenotypes

Now, we can move on to using our gene sequences to understand how
genomes are evolving in association with a phenotype of interest.

### Identifying orthogroups with OrthoFinder

The first step in doing so is to identify groups of orthologous genes,
eg. genes which are descended from a single gene in the last common
ancestor of your species. We’ll do this with OrthoFinder.

#### Prepping inputs

OrthoFinder requires amino acid sequences as input, whereas we have
nucleotide sequences. We’ll use TransDecoder to translate our nucleotide
sequences to amino acids. To do this, use the script
`./Scripts/translatingCodingSequences`, which requires, as input, the
path to your list of sample IDs.

If you use this script, please cite:

  - Haas, B., and A. Papanicolaou. “TransDecoder.” (2017).

We will also run OrthoFinder with a user-specific species tree, to
improve accuracy. If possible, obtain one for your group of species, by
doing something like trimming a larger, high-quality phylogeny from the
group you’re working with. You can find an example R script for this
process in `./Scripts/trimmingTree.R`. If no phylogeny for your species
exists, you can run OrthoFinder without the option `-s`; this will allow
OrthoFinder to infer a species tree for you.

### Running OrthoFinder

You can run OrthoFinder on your list of samples with the script
`./Scripts/runningOrthofinder`, which takes two inputs:

  - the maximum number of threads to use on computer
  - the full path to species tree

If you use this script, please cite:

  - Emms D.M. & Kelly S. (2019), Genome Biology 20:238

If you use the species tree in your work then please also cite:

  - Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
  - Emms D.M. & Kelly S. (2018), bioRxiv
    <https://doi.org/10.1101/267914>

### Inferring gene trees

Orthofinder will automatically infer gene trees for any orthogroup with
at least four sequences; however, we may want to analyze smaller
orthogroups too. You can do this by running the script
`./Scripts/treesForSmallOrthogroups`. This script does not take any
arguments; however, you’ll need to edit the path to the Orthofinder
results to match your setup.

### Evolutionary analyses in Hyphy

#### Prepping data for Hyphy

The gene trees from Orthofinder have species labels appended to each
gene name, causing them not to match the gene names in the orthogroup
sequences files. To correct this, run the R script
`./Scripts/removingTreePrefixes.R`.

The Hyphy analyses will need labelled gene trees for each orthogroup,
with the tips (and internal nodes, if appropriate) that possess your
focal trait labelled as being the foreground. You can produce these
labelled trees using the R script `./Scripts/labellingPhylogenies.R`,
which runs a built-in Hyphy utility to label trees across all of the
trees in the directory `./allGeneTrees/`.

#### Assessing positive selection with BUSTED-PH

#### Assessing shifts in selection intensity with RELAX

### Exploring evolution in non-coding elements

<https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15982>
