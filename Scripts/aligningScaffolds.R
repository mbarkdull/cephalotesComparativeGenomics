#### A script to align scaffolds of my newly sequenced genomes to the reference genome.####
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

### Align homologous scaffolds ####
# Make a directory for the intermediate fasta files produced below. This will let the function restart without losing progress, if something gets interrupted.
dir.create("./intermediateAlignments/")
# Write a function to extract one set of homologous scaffolds and align it:
alignHomologousScaffolds <- function(i) {
  # Extract out sets of matching scaffolds:
  print(paste("Aligning",
              i,
              sep = " "))
  singleSetOfScaffolds <- subset(allScaffolds, 
                                 seq.name == i)
  
  # Convert to a named character vector, the input format for DNAStringSet:
  singleSetOfScaffoldsVector <- as.character(singleSetOfScaffolds$seq.text)
  names(singleSetOfScaffoldsVector) <- singleSetOfScaffolds$codedSequenceName
  
  # Convert them to an DNAStringSet, the input format for msa:
  dnaStringSet <- Biostrings::DNAStringSet(singleSetOfScaffoldsVector,
                                           use.names = TRUE)
  # Align with msa:
  alignedGenes <- msa(dnaStringSet)
  
  # Calculate a conservation score:
  alignedGenesSequinr <- msaConvert(alignedGenes, type = "seqinr::alignment")
  distanceBetweenSequences <- seqinr::dist.alignment(alignedGenesSequinr, "identity")
  output <- c(i, distanceBetweenSequences[1])
  print(output)
  
  # Get the alignment as a dataframe:
  alignment <- data.frame(alignedGenesSequinr[["nam"]],
                          alignedGenesSequinr[["seq"]])
  colnames(alignment) <- c("new.seq.name",
                           "seq.text")
  alignment <- separate(data = alignment,
                        col = new.seq.name,
                        into = c("species", 
                    "seq.name"),
                    sep = "-",
                    remove = FALSE)
  # Get only the sequence names with the species codes, and the nucleotide sequences.
  alignment <- select(alignment,
                      new.seq.name,
                      seq.text)
  # Rename the columns so that we can export this with dat2fasta, which is picky about column names. 
  colnames(alignment) <- c("seq.name",
                           "seq.text")
  # Save the alignment as a .fasta file
  dat2fasta(alignment, 
            outfile = paste("./intermediateAlignments/",
                            i,
                            ".fasta",
                            sep = ""))
}

# Map that function over all of the scaffolds that haven't already been analyzed:
possiblyAlignHomologousScaffolds <- possibly(alignHomologousScaffolds, 
                                             otherwise = "Error.")
# List all the scaffolds:
scaffolds <- unique(allScaffolds$seq.name) 
# Remove any that have already been aligned:
alreadyAligned <- list.files("./intermediateAlignments/") 
alreadyAligned <- str_split_i(alreadyAligned,
                              pattern = ".fasta",
                              i = 1)
scaffolds <- setdiff(scaffolds, alreadyAligned)
# Run the alignment in parallel, using furrr:
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)
allAlignments <- future_map(scaffolds,
                            possiblyAlignHomologousScaffolds,
                            .progress = TRUE)

# Read in all of the fasta files that have been generated, and convert to a dataframe:
allAlignedFastaFiles <- list.files(path = "./intermediateAlignments",
                                   full.names = TRUE)
allAlignedFastas <- purrr::map(allAlignedFastaFiles,
                               phylotools::read.fasta)
allAlignedFastas <- as.data.frame(do.call(rbind, 
                                          allAlignedFastas))
# Create a column to give us the species code again:
allAlignedFastas$species <- str_split_i(allAlignedFastas$seq.name,
                                      pattern = "-",
                                      i = 1)
# Extract sequences from the same species:
individualSpeciesSequencesList <- split(allAlignedFastas, 
                                        f = allAlignedFastas$species)
# Export those sequences to a subdirectory:
dir.create("./alignedGenomes")
exportingFastaFromList <- function(i) {
  fastaDataFrame <- individualSpeciesSequencesList[[i]]
  outputFasta <- paste("./alignedGenomes/",
                       i,
                       "_alignedScaffolds.fasta",
                       sep = "")
  print(paste("Saving results to", outputFasta, sep = " "))
  dat2fasta(fastaDataFrame, 
            outfile = outputFasta)
}
species <- unique(allAlignedFastas$species)
purrr::map(species,
           exportingFastaFromList)

