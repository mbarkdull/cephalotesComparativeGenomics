library(tidyverse)
library(phylotools)
library(msa)
library(stringi)
library(spgs)

# List all of the scaffolds
allScaffolds <- list.files(path = "./individualCactusAlignments",
                           pattern = "*.fasta",
                           full.names = TRUE)

# Subset to all of the fastas for one species
species <- str_split_i(allScaffolds,
                       pattern = "_scaf",
                       i = 1) %>%
  unique()

# Read in and concatenate all the fastas for one species
dir.create("./individualCactusAlignments/wholeGenomes")
readingAndConcateningEachSpecies <- function(speciesCode) {
  speciesFastas <- grep(speciesCode, 
                        allScaffolds,
                        value = TRUE)
  readOneFasta <- function(file) {
    singleFasta <- phylotools::read.fasta(file = file)
  }
  allFasta <- purrr::map(speciesFastas,
                         readOneFasta)
  allFasta <- as.data.frame(do.call(rbind, allFasta))
  exportPrefix <- str_split_i(speciesCode,
                              pattern = "/",
                              i = 3)
  exportPrefix <- str_split_i(exportPrefix, 
                              pattern = "-",
                              i = 1)
  phylotools::dat2fasta(allFasta,
                        outfile = paste("./individualCactusAlignments/wholeGenomes/",
                                        exportPrefix,
                                        "_cactusAlignment.fasta",
                                        sep = ""))
  
}

# Export a concatenated fastas for each species
purrr::map(species,
           readingAndConcateningEachSpecies)