#### A script to prepare data for input to ProgressiveCactus ####
library(tidyverse)
library(phylotools)
library(msa)
library(seqinr)

#### Read in and process raw data ####
# Read in the reference: 
print("Reading in reference.")
referenceScaffolds <- phylotools::read.fasta(file = "./CVAR/CVAR_genome_v1.0.fasta")
referenceScaffolds$codedSequenceName <- paste("CVAR",
                                              referenceScaffolds$seq.name,
                                              sep = "-")

# Read in the other scaffold files:
print("Reading in consensus sequences")
scaffolds <- grep(list.files(path = "./consensusSequences",
                             full.names = TRUE), 
                  pattern = '*.fai', 
                  invert = TRUE, 
                  value = TRUE)

# Write a function to read in one consensus sequence and add a species code:
readSingleConsensus <- function(i) {
  speciesName <- str_split(i, "/")[[1]][[3]]
  speciesName <- str_split(speciesName, ".fa")[[1]][[1]]
  speciesName <- str_remove(speciesName, 
                            pattern = "Consensus.fa")
  test <- phylotools::read.fasta(file = i)
  
  # Add a species code to the sequence names:
  test$codedSequenceName <- paste(speciesName,
                                  test$seq.name,
                                  sep = "-")
  return(test)
}

# Map it over the list of consensus sequences:
allConsensus <- purrr::map(scaffolds,
                           readSingleConsensus)
allConsensus <- as.data.frame(do.call(rbind, allConsensus)) 

# Combine all the sequences:
allScaffolds <- plyr::rbind.fill(referenceScaffolds, allConsensus) %>%
  arrange(str_length(seq.text), seq.text)
rm(referenceScaffolds, allConsensus)
gc()

#### Export each individual scaffold's sequence as a fasta file: ####
dir.create("./individualScaffolds")

exportingIndividualScaffolds <- function(i) {
  group <- str_split(i, 
                     pattern = "-")
  group <-  group[[1]][2]
  dir.create(paste("./individualScaffolds/",
                   group,
                   sep = ""))
  scaffold <- filter(allScaffolds, codedSequenceName == i) %>%
    select(-c("seq.name"))
  colnames(scaffold) <- c("seq.text",
                          "seq.name")
  dat2fasta(scaffold, 
            outfile = paste("./individualScaffolds/",
                            group,
                            "/",
                            i,
                            ".fasta",
                            sep = ""))
}
allScaffoldNames <- unique(allScaffolds$codedSequenceName)

purrr::map(allScaffoldNames,
           exportingIndividualScaffolds)

#### Export a tree for each set of scaffolds: ####
# Create a species tree for each scaffold:
masterTree <- read.tree(file = "speciesTree/speciesCodesTree.txt")
# Get rid of the "translated_" prefix:
masterTree[["tip.label"]] <- str_split_i(masterTree[["tip.label"]],
                                         pattern = "_",
                                         i = 2)
# Rescale the tree- ProgressiveCactus takes a max branch length of 25
# Per their FAQ, "even inexact branch lengths should be fine"
#https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#frequently-asked-questions
maximumBranch <- max(masterTree[["edge.length"]])
masterTree[["edge.length"]] <- masterTree[["edge.length"]]/maximumBranch

# Make a tree for the whole-genome, consensus files:
consensusTree <- masterTree
consensusTree[["tip.label"]] <- paste(consensusTree[["tip.label"]],
                                      "Consensus.fa",
                                      sep = "")
ape::write.tree(consensusTree, 
                file = "./consensusTree.txt")

# Create a list of scaffold groups:
group <- str_split_i(allScaffoldNames, 
                   pattern = "-",
                   i = 2) %>%
  base::unique()

generatingScaffoldTrees <- function(group) {
  desiredTipLabels <- list.files(paste("./individualScaffolds/",
                                       group,
                                       sep = ""))
  # Make a dictionary to match the filenames to the master tip labels:
  dictionary <- as.data.frame(desiredTipLabels)
  dictionary$species <- str_split_i(dictionary$desiredTipLabels,
                                    pattern = "-",
                                    i = 1)
  dictionary$species <- gsub("Consensus",
                             "",
                             as.character(dictionary$species))
  # Replace the tip labels using that dictionary:
  test <- masterTree
  test[["tip.label"]] <- dictionary[match(test[["tip.label"]], dictionary$species), 1]
  # Export the tree as a newick file:
  ape::write.tree(test, 
                  file = paste("./individualScaffolds/",
                               group,
                               "/",
                               group,
                               "_tree.txt",
                               sep = ""))
}
purrr::map(group,
           generatingScaffoldTrees)
