#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Must supply path to Orthofinder MSAs.", call.=FALSE)
} 
# The command to run this script is `Rscript ./scripts/DataMSA.R [inputurls file] [path to Orthofinder MSA files]`, for example: `./scripts/DataMSA.R ./scripts/inputurls_partial /workdir/mb2337/FormicidaeMolecularEvolution/4_OrthoFinder/fasta/OrthoFinder/Results_Oct26/MultipleSequenceAlignments`

# This script will reshuffle the Orthofinder MSA files so that instead of one file per orthogroup with all species present, we have one file per species with all orthogroups present. 

library(phylotools)
library(plyr)
library(tidyverse)

# Get the list of species abbreviations by reading the fasta files in the Orthofinder input directory:
speciesInfo <- list.files(path = "./orthoFinder/fasta/") %>%
  str_split_i(pattern = "_", 
              i = 2) %>%
  str_split_i(pattern = "\\.",
              i = 1) %>%
  na.omit()

# First make a working directory and copy the folder with the MSA files there (THIS WILL NEED TO BE A VALUE SET ON THE COMMAND LINE). 
dir.create("./PAL2NALInput")
file.copy(args[1], "./PAL2NALInput", recursive = TRUE)
# Concatenate all of the MSA files into a single file:
setwd("./PAL2NALInput/MultipleSequenceAlignments")
msaFiles <- list.files(pattern = "*.fa")
allMSAFiles <- bind_rows(lapply(msaFiles, read.fasta))

# Create a function to create a subsetted fasta file for each species:
speciesChecking <- function(abbreviation, outFile){
  # Create a column in the concatenated file that tells us if it matches the species abbreviation:
  allMSAFiles$SpeciesMatch <- ifelse(grepl(abbreviation, allMSAFiles$seq.name), "True", "False")
  # Make a dataframe that contains only the rows where $SpeciesMatch == True. 
  SpeciesMSA <- subset(allMSAFiles, SpeciesMatch == "True", select = c("seq.name", "seq.text"))
  # Save that to a .fasta file: ../test.fasta
  dat2fasta(SpeciesMSA, outfile = outFile)
}

# Write a for loop to run the function on all species:
for (i in speciesInfo)
{
  print(i)
  abbreviation <- i
  print(abbreviation)
  outFile <- (paste("../", i, "_msa.fasta", sep = ""))
  print(outFile)
  print("__________________________________________________")
  
  speciesChecking(abbreviation = abbreviation, outFile = outFile)
}

setwd("../")
setwd("../")