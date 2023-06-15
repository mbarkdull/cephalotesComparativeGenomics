#!/usr/bin/env Rscript

# Make a directory for outputs:
dir.create("./filteredCDS")

# Load packages:
library(phylotools)
library(plyr)
library(tidyverse)

# Read in the input data:
# Get the list of species abbreviations by reading the fasta files in the Orthofinder input directory:
abbreviations <- list.files(path = "./orthoFinder/fasta/") %>%
  str_split_i(pattern = "_", 
              i = 2) %>%
  str_split_i(pattern = "\\.",
              i = 1) %>%
  na.omit()

# I've written this function that will produce a subsetted nucleotide sequence file:
# I should also add error messages into this function. 
# The cdsFile argument is the path to the coding sequences output by Transdecoder, something like: ./4_2_TransdecoderCodingSequences/cds_acol_filteredTranscripts.fasta
# The mdsFile argument is the path to the recombined MSA file output by DataMSA.R, something like: ./6_1_SpeciesMSA/acol_msa.fasta
# The msaPrefix argument is 
cdsFiltering <- function(cdsFile, msaFile, msaPrefix, filteredOutput, msaOutput) {
  # Load packages:
  library(phylotools)
  library(plyr)
  library(tidyverse)
  # Read in coding sequences and split the sequence name column by the space:
  CodingSequences <- read.fasta(file = cdsFile)
  CodingSequences <- separate(data = CodingSequences, col = seq.name, into = c("seq.name", "extraInfo"), sep = " ")
  CodingSequences <- select(CodingSequences, -c("extraInfo"))
  # Read in protein sequences and remove the string "mpha_transcripts_" from the sequence names:
  MSAProteinSequences <- read.fasta(file = msaFile)
  MSAProteinSequences$seq.name <- gsub(msaPrefix, '', MSAProteinSequences$seq.name)
  MSAProteinSequences$seq.name <- str_split_i(MSAProteinSequences$seq.name,
                                              pattern = "\\.",
                                              i = 1)
  # Create a column that checks if the coding sequence found in the protein sequence file?
  CodingSequences$inMSA <- CodingSequences$seq.name %in% MSAProteinSequences$seq.name
  # Subset the coding sequences based on the value in that column:
  FilteredCodingSequences <- subset(CodingSequences, inMSA == "TRUE", select = c("seq.name", "seq.text"))
  # Save the output
  dat2fasta(FilteredCodingSequences, outfile = filteredOutput)
  dat2fasta(MSAProteinSequences, outfile = msaOutput)
}

# This for loop iterates the subsetting function on all species:
for (i in abbreviations)
{
  print(i)
  cdsFile <- (paste("./", i, "_genesAlignedToReference.fasta", sep = ""))
  msaFile <- (paste("./PAL2NALInput/", i, "_msa.fasta", sep = ""))
  msaPrefix <- (paste("translated_", i, "_", sep = ""))
  filteredOutput <- (paste("./filteredCDS/filtered_", i, "_cds.fasta", sep = ""))
  print(filteredOutput)
  msaOutput <- (paste("./PAL2NALInput/proteins_", i, ".fasta", sep = ""))
  print(msaOutput)
  print("__________________________________________________")
  
  cdsFiltering(cdsFile = cdsFile, msaFile = msaFile, msaPrefix = msaPrefix, filteredOutput = filteredOutput, msaOutput = msaOutput)
}