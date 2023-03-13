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
2.  Map the trimmed reads back to the reference genome, with NextGenMap.
3.  Call SNPs and indels with SAMtools.
4.  Filter the calls by quality.
5.  Assemble a consensus sequence.

These steps are all implemented by the scripts
`./Scripts/referenceGuidedAssembly`, which requires two arguments:
first, a list of samples, and two, the number of threads to use for
multithreaded steps. In other words, run this script with the command
`./Scripts/referenceGuidedAssembly [path to file listing samples]
[number of threads to use when possible]`.

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

Repetitive regions of the genome must be identified, then soft-masked,
prior to annotation. To accomplish this task, we will:

  - Identify repeats in each assembled genome, with RepeatModeler2;
  - Soft mask our identified repeats with RepeatMasker;
  - Soft mask known Hexapod repetitive elements with RepeatMasker, using
    the RepeatMasker repeat database for Hexapoda.

All of these steps are implemented in the script
`./Scripts/handlingRepeats`, which can be run with the command
`./Scripts/handlingRepeats [path to file listing samples] [number of
threads to use when possible]`. Again, note that this script uses
hard-coded paths to RepeatModeler and RepeatMasker on Cornell’s BioHPC,
so you’ll likely need to change those to reflect your setup.

Outputs of this step will appear in the directory `./handlingRepeats`.

If you use this script, please cite:

  - Flynn, Jullien M., et al. “RepeatModeler2 for automated genomic
    discovery of transposable element families.” Proceedings of the
    National Academy of Sciences 117.17 (2020): 9451-9457.
  - Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015
    <http://www.repeatmasker.org>.

## Identifying genes in each genome

Since we assembled our new genomes to a reference genome, genes can be
identified quite simply using the genome annotation which accompanies
our reference.

### Aligning scaffolds

Our genomes were assembled into scaffolds, mapped to the scaffolds of
the reference. We now need to align the new genomes to our reference, so
that we can find gene sequences based on the positions given by the
reference genome annotation file. Because we are working with whole
genomes that contain very large scaffolds, we will use the aligner
[Progressivecactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md),
as implemented in the script `./Scripts/alignWithProgressiveCactus`.

As input, Progressivecactus requires a file that contains:

1.  A Newick-formatted tree of the sequences to be aligned. This tree
    can be generated from an existing species tree using the R script
    `prepForProgressiveCactus.R`.
2.  A two-column list of (a) the sequences to be aligned and (b) the
    paths to the files containing those sequences.

`./Scripts/alignWithProgressiveCactus` will construct this input file,
run Progressivecactus on the genomes, and then produce an aligned FASTA
file for each genome.

If you use this script, please cite Progressivecactus and its
dependency, hal:

  - Armstrong, Joel, et al. “Progressive Cactus is a multiple-genome
    aligner for the thousand-genome era.” Nature 587.7833 (2020):
    246-251.
  - Glenn Hickey, Benedict Paten, Dent Earl, Daniel Zerbino, David
    Haussler, HAL: a hierarchical format for storing and analyzing
    multiple genome alignments, Bioinformatics, Volume 29, Issue 10, May
    2013, Pages 1341–1342,
    <https://doi.org/10.1093/bioinformatics/btt128>

If your scaffolds are fairly small, you can alternatively align them
with the R script `./Scripts/aligningScaffolds.R`, which will select out
homologous scaffolds across your reference and new genomes, align them,
then recombine them by genome into individual, aligned genomes. This
process is parallelized using the excellent R package `furrr`. Note that
if your scaffolds are large, this script will run into memory issues and
fail.

### Remapping the annotation file

When we align our genomes, we may be shifting our reference genome
scaffolds relative to the position of gene features in the annotation
file. For instance, this might happen if our multiple sequence
alignments indicate that our reference genome has a deletion relative to
a newly sequenced genome. This means that the coordinates in the
original annotation will no longer correctly point to the features we
care about.

To solve this problem, we need to remap our genome annotation to the
new, aligned version of the genome. We can do this using a tool called
flo, developed by the Wurm lab for their work on the fire ant genome.

To do this:

1.  First, we need to filter the reference annotation to include only
    CDS features. You can do this with the R script,
    `./Scripts/annotationFilteringToCDS.R`. The script will output the
    filtered annotation in a new directory, `./liftoverAnnotation/`.
2.  Next, run the script `./Scripts/runningFlo`. This script does
    several things, including:
      - install flo and its dependencies.
      - create the configuration file, `flo_opts.yaml`.
      - process the filtered annotation file we just produced using
        genometools into a format that flo will accept.
      - run flo\!

Because this script creates the configuration file for you, you’ll have
to provide several pieces of information as command-line options:

1.  The path to your genometools install, something like
    `/programs/genometools-1.5.9/bin/`
2.  The path to your reference genome assembly, something like
    `/local/workdir/mb2337/cephalotesComparativeGenomics/CVAR/CVAR_genome_v1.0.fasta`
3.  The path to your new, aligned genome assembly, something like
    `/local/workdir/mb2337/cephalotesComparativeGenomics/alignedGenomes/CVAR_alignedScaffolds.fasta`
4.  The number of CPU cores used to parallelise flo; should be near or
    equal to the number of cores on your machine.
5.  The minimum identity parameter used by BLAT; I used `95`

If you use this script, please cite both flo itself, and GNU Parallel:

  - The fire ant social chromosome supergene variant Sb shows low
    diversity but high divergence from SB. 2017. R Pracana, A Priyam, I
    Levantis, Y Wurm. Molecular Ecology, doi: 10.1111/mec.14054. +Tange
    O. GNU Parallel - the command-line power tool. login. 2011;36:42–7.

### Extracting gene sequences

Now that scaffolds have been aligned, we can extract out the sequence of
each gene that was identified in the reference genome from our new
genomes. To do this, use the script `./Scripts/extractCDSwithGFF`, which
requires two inputs:

1.  the path to your list of sample IDs
2.  the path to the GFF-formatted reference annotation.

This script will export a fasta file of gene sequences for each sample
you sequenced, in the output directory `./codingSequences/`.

If you use this script, please cite:

  - Quinlan, Aaron R., and Ira M. Hall. “BEDTools: a flexible suite of
    utilities for comparing genomic features.” Bioinformatics 26.6
    (2010): 841-842.

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

### Assessing positive selection with BUSTED-PH

### Assessing shifts in selection intensity with RELAX

### Exploring evolution in non-coding elements

<https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15982>
