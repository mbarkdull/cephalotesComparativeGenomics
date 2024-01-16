library(rBLAST)
library(tidyverse)
library(googlesheets4)
library(rentrez)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/programs/ncbi-blast-2.13.0/bin/", sep= .Platform$path.sep))

#### Download the candidate genes from NCBI based on my google sheet: ####
# Read in the information from my public google sheet, getting IDs without the "LOC" prefix:
gs4_deauth()
candidates <- read_sheet("https://docs.google.com/spreadsheets/d/1JFIVkdyz9xX_N1ywrEN9uD_KNDuLpaHm8dx_HDHrKRo/edit?usp=sharing")
candidates$`gene symbol` <- as.character(candidates$`gene symbol`)
ncbiIDs <- candidates$`gene symbol`
ncbiIDs <- gsub(pattern = "LOC",
                replacement = "",
                ncbiIDs)

# Get the list of already-downloaded proteins:
alreadyDownloaded <- list.files(path = "./18_CandidateGeneAnalyses/",
                                pattern = "*fasta") %>% 
  str_split_i(pattern = "\\.",
              i = 1)
alreadyDownloaded <- gsub(pattern = "LOC",
                          replacement = "",
                          alreadyDownloaded)

# Get a list of only proteins that still need to be downloaded:
ncbiIDs <- base::setdiff(ncbiIDs,
                         alreadyDownloaded)

##### Get protein sequences from NCBI ####
# Create an output directory:
dir.create("./18_CandidateGeneAnalyses/")
gettingProteinSequence <- function(ID) {
  # Link the gene ID to a protein ID:
  geneToProtein <- entrez_link(dbfrom = 'gene', 
                               id = ID, 
                               db = 'protein')
  proteinID <- geneToProtein[["links"]][["gene_protein_refseq"]][[1]]
  # Get that protein sequence as a fasta:
  entrez_fetch(db = "protein", 
               id = proteinID, 
               rettype = "fasta") %>%
    write(file = paste("./18_CandidateGeneAnalyses/LOC",
                       ID,
                       ".fasta",
                       sep = ""))
}

# Make a safe version with possibly:
possiblyGettingProteinSequence <- purrr::possibly(gettingProteinSequence,
                                                  otherwise = "Error")
# Download the list of protein sequences:
purrr::map(ncbiIDs,
           possiblyGettingProteinSequence)

#### Combine all the proteomes against which we'll BLAST ####
# Use an if-else statement so this step doesn't have to be repeated:
if (file.exists("allProteomes.fasta")) {
  print("The file exists!")
} else {
  print("The file does not exist.")
  # Construct a list of the proteomes:
  allGenomePaths <- list.files("./02_annotationsAndExons",
                               full.names = TRUE)
  allGenomes <- list.files("./02_annotationsAndExons",
                           full.names = FALSE)
  allGenomes <- paste(allGenomePaths,
                      "/",
                      allGenomes,
                      "_proteins.fasta",
                      sep = "")
  # Read in all the proteomes and combine them:
  genomes <- purrr::map(allGenomes,
                        phylotools::read.fasta)
  genomes <- as.data.frame(do.call(rbind, genomes)) 
  phylotools::dat2fasta(genomes,
                        outfile = "allProteomes.fasta")
  rm(genomes)
}

#### BLAST the candidates against my proteomes ####
# Construct the database:
makeblastdb(file = './allProteomes.fasta', dbtype = "prot")

# Process the orthogroups so they can be identified:
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
  querySequence <- readAAStringSet(i, format = 'fasta')
  
  # Prep the search:
  blastSearch <- blast(db = './allProteomes.fasta', type = 'blastp')
  
  # Run the search:
  searchResults <- predict(blastSearch, querySequence, BLAST_args = "-max_target_seqs 1")
  
  # Make a column that is the product of length and identity for selecting the best match:
  searchResults$bestMatch <- searchResults$Perc.Ident*searchResults$Alignment.Length
  
  # Filter the search to include only matches with 100% identity:
  searchResults <- dplyr::filter(searchResults, bestMatch == max(bestMatch))
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

# Make a safe version with possibly:
possiblyFetchingOrthogroups <- possibly(fetchingOrthogroups,
                                        otherwise = "Error")

# List all of the genes to potentially BLAST for:
allGenes <- list.files("./18_CandidateGeneAnalyses/",
                       full.names = TRUE)

# Read in any past results, so as not to duplicate analyses:
if(file.exists("./18_CandidateGeneAnalyses/blastSummary.csv") == TRUE){
  pastResults <- read_delim(file = "./18_CandidateGeneAnalyses/blastSummary.csv")
  pastResults <- pastResults$candidateFasta
  
  # Get a list of genes to search:
  genesToSearch <- base::setdiff(allGenes, pastResults)
  
} else {
  genesToSearch <- allGenes
}

# Read in the HYPHY results so they can be matched to the candidate genes:
allRelaxResults <- read_csv(file = "allRelaxResults.csv") 
allBUSTEDPHResults <- read_csv(file = "allBustedPHResults.csv") 

# BLAST for the new candidates:
allGenesMatches <- purrr::map(genesToSearch,
                              possiblyFetchingOrthogroups)
allGenesMatches <- as.data.frame(do.call(rbind, allGenesMatches)) 

# Do the rest of the analysis only IF there are newly-identified candidates:
if(length(colnames(allGenesMatches)) > 1) {
  colnames(allGenesMatches) <- c("orthogroup",
                                 "candidateFasta",
                                 "percentIdentity")
  allGenesMatches$`gene symbol` <- str_split_i(allGenesMatches$candidateFasta,
                                               pattern = "/",
                                               i = 4) %>%
    str_split_i(pattern = "\\.",
                i = 1)
  
  # Join the BLAST results to the info about the candidate genes:
  allGenesMatches <- full_join(allGenesMatches, 
                               candidates, 
                               by = c("gene symbol" = "gene symbol"))
  
  # Export a summary of the BLAST results to obviate analyses being re-run:
  write_delim(allGenesMatches,
              file = "./18_CandidateGeneAnalyses/blastSummary.csv",
              delim = ",")
  
  # Join the BLAST/candidate results to the relax results:
  selectionIntensityCandidateGenes <- left_join(allGenesMatches,
                                                allRelaxResults,
                                                by = c("orthogroup" = "orthogroup")) %>%
    filter(pValueFDR <= 0.05)
  selectionIntensityCandidateGenes <- selectionIntensityCandidateGenes%>%
    mutate(selectionType = case_when(kValue < 1 ~ "Relaxed",
                                     kValue > 1 ~ "Intensified"))
  
  # Join the BLAST/candidate results to the BUSTED-PH results:
  positiveSelectionForegroundCandidateGenes <- left_join(allGenesMatches,
                                                         allBUSTEDPHResults,
                                                         by = c("orthogroup" = "orthogroup")) %>%
    filter(testPvalueFDR <= 0.05 &
             backgroundPvalueFDR >= 0.05 &
             differencePvalueFDR <= 0.05)
  
  if(length(positiveSelectionForegroundCandidateGenes$orthogroup) > 0){
    positiveSelectionForegroundCandidateGenes$selectionType <- "foregroundPositive"
  } else {
    print("No results")
  }
  
  positiveSelectionBackgroundCandidateGenes <- left_join(allGenesMatches,
                                                         allBUSTEDPHResults,
                                                         by = c("orthogroup" = "orthogroup")) %>%
    filter(testPvalueFDR >= 0.05 &
             backgroundPvalueFDR <= 0.05 &
             differencePvalueFDR <= 0.05)
  
  if(length(positiveSelectionBackgroundCandidateGenes$orthogroup) > 0){
    positiveSelectionBackgroundCandidateGenes$selectionType <- "backgroundPositive"
  } else {
    print("No results")
  }
  
  # Combine new results with old results, if any exist:
  if(file.exists("./18_CandidateGeneAnalyses/significantResults.csv") == TRUE){
    pastResults <- read_delim(file = "./18_CandidateGeneAnalyses/significantResults.csv") %>%
      distinct()
    pastResults$percentIdentity <- as.character(pastResults$percentIdentity)
    
    allResults <- dplyr::bind_rows(positiveSelectionBackgroundCandidateGenes, 
                                   positiveSelectionForegroundCandidateGenes, 
                                   selectionIntensityCandidateGenes,
                                   pastResults) %>%
      distinct()
    
    # Export all results:
    write_delim(allResults,
                file = "./18_CandidateGeneAnalyses/significantResults.csv",
                delim = ",")
  } else {
    allResults <- dplyr::bind_rows(positiveSelectionBackgroundCandidateGenes, 
                                   positiveSelectionForegroundCandidateGenes, 
                                   selectionIntensityCandidateGenes) %>%
      distinct()
    write_delim(allResults,
                file = "./18_CandidateGeneAnalyses/significantResults.csv",
                delim = ",")
  }
} else {
  print("No new matches")
}