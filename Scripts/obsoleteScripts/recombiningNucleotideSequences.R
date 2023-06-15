#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (path to orthoFinder orthogroup sequences).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "./orthoFinder/fasta/OrthoFinder/Results_Mar28_1/Orthogroup_Sequences/"
}

# The command to run this script is `Rscript ./scripts/DataSubsetCDS.R [path to Orthofinder orthogroup sequence files]`, for example: `Rscript ./scripts/DataSubsetCDS.R ./4_OrthoFinder/fasta/OrthoFinder/Results_Nov12/MultipleSequenceAlignments/`
# Get the list of species abbreviations by reading the fasta files in the Orthofinder input directory:
abbreviations <- list.files(path = "./orthoFinder/fasta/") %>%
  str_split_i(pattern = "_", 
              i = 2) %>%
  str_split_i(pattern = "\\.",
              i = 1) %>%
  na.omit()

# Make a directory for this step:
dir.create("./cdsOrthogroups")

# Copy the nucleotide sequences from TransDecoder to the directory:
for (i in abbreviations) {
  alignedNucleotidesFile <- paste("./translatedData/",
                                  i,
                                  "_genesAlignedToReference.fasta.transdecoder.cds",
                                  sep = "")
  file.copy(alignedNucleotidesFile, "./cdsOrthogroups", recursive = TRUE)
}

# Read in all of the aligned nucleotide sequences
# Concatenate all of the cds files into a single file:
cdsFiles <- list.files(path = "./cdsOrthogroups",
                       pattern = "*.cds",
                       full.names = TRUE)
allCDSFiles <- bind_rows(lapply(cdsFiles, read.fasta))
allCDSFiles$seq.name <- str_split_i(allCDSFiles$seq.name,
                                    pattern = " ",
                                    i = 1)

# Copy the orthogroup amino acid sequences to the directory:
file.copy(args[1], "./cdsOrthogroups", recursive = TRUE)
#file.copy("./orthoFinder/fasta/OrthoFinder/Results_Mar28_1/Orthogroup_Sequences/", "./cdsOrthogroups", recursive = TRUE)

# List all of the amino acid sequence files:
allAAfiles <- list.files("./cdsOrthogroups/Orthogroup_Sequences",
                         full.names = TRUE)

# Write a function to select the nucleotide sequences for a single orthogroup:
extractSingleOrthogroupCDS <- function(i) {
  # Read in an orthogroup:
  orthogroup <- read.fasta(file = i)
  # Get the orthogroup number from the filename:
  orthogroupNumber <- str_split_i(i, 
                                  pattern = "/",
                                  i = 4) %>%
    str_split_i(pattern = "\\.",
                i = 1)
  # Extract out the nucleotide sequences for the orthogroup 
  allCDSFiles$orthogroupMatch <- allCDSFiles$seq.name %in% orthogroup$seq.name
  orthogroupCDSSequences <- filter(allCDSFiles, 
                                   orthogroupMatch == "TRUE") %>%
    select(-c("orthogroupMatch"))
  # Export as a fasta file
  phylotools::dat2fasta(orthogroupCDSSequences,
                        outfile = paste("./cdsOrthogroups/",
                                        orthogroupNumber,
                                        "_cdsSequences.fasta",
                                        sep = ""))
}

# Do this for all the orthogroups. 
purrr::map(allAAfiles,
           extractSingleOrthogroupCDS)
