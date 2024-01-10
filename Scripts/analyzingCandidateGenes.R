library(rBLAST)
library(tidyverse)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/programs/ncbi-blast-2.13.0/bin/", sep= .Platform$path.sep))

# Read in all the genomes, combine them, export:
if (file.exists("allGenomes.fasta")) {
  print("The file exists!")
} else {
  print("The file does not exist.")
  allGenomePaths <- list.files("./02_annotationsAndExons",
                               full.names = TRUE)
  allGenomes <- list.files("./02_annotationsAndExons",
                           full.names = FALSE)
  allGenomes <- paste(allGenomePaths,
                      "/",
                      allGenomes,
                      "_cds.fasta",
                      sep = "")
  genomes <- purrr::map(allGenomes,
                        phylotools::read.fasta)
  genomes <- as.data.frame(do.call(rbind, genomes)) 
  phylotools::dat2fasta(genomes,
                        outfile = "allGenomes.fasta")
  rm(genomes)
}

# Construct the database:
makeblastdb(file = './allGenomes.fasta', dbtype = "nucl")

# Read in the orthogroups file:
orthogroupMembers <- read_delim(file = "./03_OrthoFinder/fasta/OrthoFinder/Results_Jul07/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", 
                                col_names = TRUE,
                                delim = "\t") %>%
  dplyr::select(-c("OG",
                   "Gene Tree Parent Clade",
                   "acol_proteins")) 
orthogroupMembers <- orthogroupMembers %>%
  tidyr::pivot_longer(cols = -c(HOG),
                      names_to = "junk") %>%
  dplyr::select(-junk) %>%
  dplyr::distinct()
orthogroupMembers$HOG <- str_split_i(orthogroupMembers$HOG,
                                     pattern = "\\.",
                                     i = 2)

# Write a function to blast candidate genes against my genomes and identify the matching orthogroup:
fetchingOrthogroups <- function(i) {
  # Read in my query sequence:
  querySequence <- readDNAStringSet(i, format = 'fasta')
  
  # Prep the search:
  blastSearch <- blast(db = './allGenomes.fasta', type = 'blastn')
  # Run the search:
  searchResults <- predict(blastSearch, querySequence, BLAST_args = "-max_target_seqs 1")
  # Filter the search to include only matches with 100% identity:
  searchResults <- dplyr::filter(searchResults, Perc.Ident == max(Perc.Ident))
  searchResults <- head(searchResults, n = 1)
  # Extract the percent identity of the match:
  percentID <- max(searchResults$Perc.Ident)
  
  # Get out the matching sequence name:
  matchingSequence <- searchResults$SubjectID
  
  # Get the orthogroup that contains the sequence from our BLAST:
  orthogroup <- filter(orthogroupMembers, matchingSequence == value)
  if (length(orthogroup$HOG) == 1) {
    orthogroup <- orthogroup$HOG
  } else {
    orthogroup <- "noMatch"
  }
  
  
  #Construct an object that contains that orthogroup, and the file it matches to, so you know which gene of interest goes with each orthogroup:
  results <- c(orthogroup, i, percentID)
  return(results)
}
possiblyFetchingOrthogroups <- possibly(fetchingOrthogroups,
                                        otherwise = "Error")


allRelaxResults <- read_csv(file = "allRelaxResults.csv") 
allBUSTEDPHResults <- read_csv(file = "allBustedPHResults.csv") 


allGenes <- list.files("./genesOfInterest",
                       full.names = TRUE)
allGenesMatches <- purrr::map(allGenes,
                              possiblyFetchingOrthogroups)
allGenesMatches <- as.data.frame(do.call(rbind, allGenesMatches)) 

selectionIntensityCandidateGenes <- left_join(allGenesMatches,
                                              allRelaxResults,
                                              by = c("V1" = "orthogroup")) %>%
  filter(pValueFDR <= 0.05)
selectionIntensityCandidateGenes <- selectionIntensityCandidateGenes%>%
  mutate(selectionType = case_when(kValue < 1 ~ "Relaxed",
                                   kValue > 1 ~ "Intensified"))

positiveSelectionForegroundCandidateGenes <- left_join(allGenesMatches,
                                                       allBUSTEDPHResults,
                                                       by = c("V1" = "orthogroup")) %>%
  filter(testPvalueFDR <= 0.05 &
           backgroundPvalueFDR >= 0.05 &
           differencePvalueFDR <= 0.05)

positiveSelectionBackgroundCandidateGenes <- left_join(allGenesMatches,
                                                       allBUSTEDPHResults,
                                                       by = c("V1" = "orthogroup")) %>%
  filter(testPvalueFDR >= 0.05 &
           backgroundPvalueFDR <= 0.05 &
           differencePvalueFDR <= 0.05)

