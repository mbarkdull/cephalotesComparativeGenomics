dir.create("./11_RERconverge")

# Load packages:
library("furrr")
library("future")
library("RERconverge")
library("Matrix")
library("castor")
library("FSA")
library("dplyr")
library("ggplot2")
library("phangorn")
library("weights")
library("Hmisc")
library("RcppArmadillo")
library("Rcpp")
library("knitr")
library("geiger")
library("phytools")
library("maps")
library("ape")
library("gplots")
library("RColorBrewer")
library("stats")
library("graphics")
library("grDevices")
library("utils")
library("datasets")
library("methods")
library("base")

# Set up future for paralellization:
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)

#### 1. Get gene trees for input ####
## Get versions of the species tree with gene names as tip labels ----
dir.create("./11_RERconverge/01_inputTrees")
# Write a function to produce versions of the species tree with gene names as tip labels
renamingSpeciesTreeTips <- function(alignmentFile){
  # Read in the alignment:
  alignment <- phylotools::read.fasta(alignmentFile)
  # Read in the species tree:
  speciesTree <- read.tree("./speciesTree/speciesCodesTree.txt")
  speciesTree[["tip.label"]] <- stringr::str_split_i(speciesTree[["tip.label"]],
                                                     pattern = "_",
                                                     i = 2)
  # Check each species tree tip to see...
  for (i in speciesTree[["tip.label"]]) {
    # Is it in a gene in the alignment?
    presenceInAlignment <- stringr::str_detect(alignment$seq.name,
                                               paste("^",
                                                     i,
                                                     sep = ""),
                                               negate = FALSE)
    # How many times?
    numberOfOccurences <- length(which(presenceInAlignment == TRUE))
    # If appropriate, rename the species tree tip
    if (numberOfOccurences == 0) {
      print(paste("Species absent from alignment, removing",
                  i,
                  "from tree"))
      speciesTree <- ape::drop.tip(speciesTree,
                                   tip = i)
    } else if (numberOfOccurences == 1) {
      print(paste(i,
                  "present a single time in alignment, renaming tree tip"))
      # Get the index of the species tree tip that matches
      speciesTip <- stringr::str_detect(speciesTree[["tip.label"]],
                                        paste("^",
                                              i,
                                              sep = ""),
                                        negate = FALSE)
      
      speciesTree[["tip.label"]][which(speciesTip == TRUE)] <- stringr::str_subset(alignment$seq.name,
                                                                                   paste("^",
                                                                                         i,
                                                                                         sep = ""),
                                                                                   negate = FALSE)
      speciesTree[["tip.label"]]
      
    } else {
      print(paste(i,
                  "present multiple times in alignment, randomly selecting one gene for renaming"))
      # Get the index of the species tree tip that matches
      speciesTip <- stringr::str_detect(speciesTree[["tip.label"]],
                                        paste("^",
                                              i,
                                              sep = ""),
                                        negate = FALSE)
      # Get the indices of the alignment genes that match:
      alignmentGenes <- stringr::str_detect(alignment$seq.name,
                                            paste("^",
                                                  i,
                                                  sep = ""),
                                            negate = FALSE)
      sample(which(alignmentGenes == TRUE), 1)
      # Replace in species tree:
      speciesTree[["tip.label"]][which(speciesTip == TRUE)] <- alignment$seq.name[sample(which(alignmentGenes == TRUE), 1)]
    }
  }
  
  # Save the tree to a file:
  orthogroupNumber <- stringr::str_split_i(string = alignmentFile,
                                           pattern = "/",
                                           i = 4) %>%
    stringr::str_split_i(pattern = "\\.",
                         i = 1)
  ape::write.tree(speciesTree, 
                  file = paste("./11_RERconverge/01_inputTrees/",
                               orthogroupNumber,
                               ".tree",
                               sep = ""))
}

# Make a safe version with possibly:
possiblyrenamingSpeciesTreeTips <- purrr::possibly(renamingSpeciesTreeTips,
                                                   otherwise = "Error")

# List all of the orthogroup alignments:
alignmentsToFix <- list.files(path = "./06_AlignedNucleotideOrthogroups/alignedOrthogroupNucleotideSequences",
                              full.names = TRUE)

# Map over all alignments:
furrr::future_map(alignmentsToFix,
                  possiblyrenamingSpeciesTreeTips)

## Write a function to drop genes from the alignments that aren't in the trees, i.e. in cases where a species has multiple genes in an orthogroup ----
dir.create("./11_RERconverge/02_filteredAlignments")

filteringAlignments <- function(inputTree) {
  # Read in one tree
  orthogroupTree <- read.tree(inputTree)
  # Get the orthogroup number 
  orthogroupNumber <- stringr::str_split_i(string = inputTree,
                                           pattern = "/",
                                           i = 4) %>%
    stringr::str_split_i(pattern = "\\.",
                         i = 1)
  # Read in the alignment:
  orthogroupAlignment <- phylotools::read.fasta(paste("./06_AlignedNucleotideOrthogroups/alignedOrthogroupNucleotideSequences/",
                                                      orthogroupNumber,
                                                      ".fa",
                                                      sep = ""))
  # Drop alignment rows that aren't in the orthogroup tree tip labels:
  filteredAlignment <- filter(orthogroupAlignment, 
                              seq.name %in% orthogroupTree[["tip.label"]])
  # Export the alignment:
  phylotools::dat2fasta(filteredAlignment,
                        outfile = paste("./11_RERconverge/02_filteredAlignments/filtered_",
                                        orthogroupNumber,
                                        ".fasta",
                                        sep = ""))
}
possiblyfilteringAlignments <- purrr::possibly(filteringAlignments,
                                               otherwise = "Error")
orthogroupTrees <- list.files(path = "./11_RERconverge/01_inputTrees",
                              full.names = TRUE)
furrr::future_map(orthogroupTrees,
                  possiblyfilteringAlignments)

## Using those orthogroup trees and filtered alignments, infer trees with branch lengths for RERconverge input ----
dir.create("./11_RERconverge/03_RERconvergeInputTrees")

# Write a function to read in a single orthogroup reference tree and an alignment, and produce a tree for RERconverge from those:
producingInputGeneTrees <- function(alignmentFile) {
  # Read in a single alignment:
  alignment <- read.phyDat(file = alignmentFile, 
                           type = "DNA", 
                           format = "fasta")
  
  # Read in the corresponding gene tree:
  orthogroupNumber <- stringr::str_split_i(string = alignmentFile,
                                           pattern = "/",
                                           i = 4) %>%
    stringr::str_split_i(pattern = "_",
                         i = 2) %>%
    stringr::str_split_i(pattern = "\\.",
                         i = 1)
  geneTree <- read.tree(paste("./11_RERconverge/01_inputTrees/",
                              orthogroupNumber,
                              ".tree",
                              sep = ""))
  
  if (length(geneTree[["tip.label"]]) > 2) {
    # Unroot the tree
    geneTree <- unroot(geneTree)
    # Just in case, set all branches to 1 first (as per RERconverge function)
    geneTree$edge.length <- c(rep(1,
                                  nrow(geneTree$edge)))
    
    #Run distance estimation using the generalized time reversible submodel:
    # Generate an initial pml tree, using default settings from RERconverge:
    initialTree <- phangorn::pml(tree = geneTree, 
                                 data = alignment, 
                                 model = "GTR", 
                                 k = 4, 
                                 rearrangement = "none")
    # Generate a tree
    finalTree <- optim.pml(initialTree,
                           optInv = T,
                           optGamma = T,
                           optEdge = T,
                           rearrangement = "none",
                           model = "GTR")
    # change the tip labels for the tree to just the species codes:
    tree <- finalTree$tree
    tree[["tip.label"]] <- stringr::str_split_i(tree[["tip.label"]], 
                                                pattern = "_", 
                                                i = 1)
    ape::write.tree(tree, 
                    file = paste("./11_RERconverge/03_RERconvergeInputTrees/",
                                 orthogroupNumber,
                                 ".tree",
                                 sep = ""))
  } else {
    print("Too few tips in tree; skipping this orthogroup")
  }
  
}

# Make a safe version with possibly:
possiblyproducingInputGeneTrees <- purrr::possibly(producingInputGeneTrees,
                                                   otherwise = "Error")

# List all of the alignments:
alignmentFiles <- list.files(path = "./11_RERconverge/02_filteredAlignments",
                             full.names = TRUE) %>%
  stringr::str_split_i(pattern = "/",
                       i = 4) %>%
  stringr::str_split_i(pattern = "\\.",
                       i = 1)  %>%
  stringr::str_split_i(pattern = "_",
                       i = 2)

# List any extant tree files:
extantTreeFiles <- list.files(path = "./11_RERconverge/03_RERconvergeInputTrees",
                              full.names = TRUE) %>%
  stringr::str_split_i(pattern = "/",
                       i = 4) %>%
  stringr::str_split_i(pattern = "\\.",
                       i = 1)

# Remove any alignment files that correspond to extant trees:
stillToInferTree <- setdiff(alignmentFiles,
                            extantTreeFiles)
alignmentFiles <- paste("./11_RERconverge/02_filteredAlignments/filtered_",
                        stillToInferTree,
                        ".fasta",
                        sep = "")

# Map over the alignments using future_map:
furrr::future_map(alignmentFiles,
                  possiblyproducingInputGeneTrees)

# Now combine trees into a single tab delimited file as required by RERconverge:
# See https://clark.genetics.utah.edu/wp-content/uploads/2020/11/FullWalkthroughUTD.html#reading-in-gene-trees-with-readtrees
treesToCombine <- list.files("./11_RERconverge/03_RERconvergeInputTrees",
                             full.names = TRUE)
createMasterTreeFile <- function(input) {
  tree <- readr::read_file(input)
  tree <- gsub("\n",
               "",
               tree)
  orthogroupNumber <- stringr::str_split_i(input, 
                                           pattern = "/",
                                           i = 4) %>%
    stringr::str_split_i(pattern = "\\.",
                         i = 1)
  singleResult <- c(orthogroupNumber, tree)
  return(singleResult)
}
possiblycreateMasterTreeFile <- purrr::possibly(createMasterTreeFile,
                                                otherwise = "Error")
masterTreeFile <- purrr::map(treesToCombine,
                             possiblycreateMasterTreeFile)
masterTreeFile <- as.data.frame(do.call(rbind, masterTreeFile)) 
write.table(masterTreeFile, 
            file = "./11_RERconverge/masterTreeFile.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## Run RERconverge analyses ----
# Read in all trees:
rerConvergeTrees <- readTrees("./11_RERconverge/masterTreeFile.txt")

# Estimate relative evolutionary rates:
relativeEvolutionaryRates <- getAllResiduals(rerConvergeTrees,
                                             min.sp = 2,
                                             transform = "sqrt", 
                                             weighted = T, 
                                             scale = T)
# Do some plotting
#make average and gene tree plots
outgroup <- c("CVAR",
              "POW0461")
par(mfrow = c(1, 2))
averageTree <- plotTreeHighlightBranches(rerConvergeTrees$masterTree, 
                                         outgroup = outgroup,
                                         hlspecies = c("CVAR",
                                                       "CSM3441",
                                                       "POW0123"), 
                                         hlcols = c("blue",
                                                    "red"), 
                                         main = "Average tree") 
specificGeneTree <- plotTreeHighlightBranches(rerConvergeTrees$trees$HOG0000032, 
                                              outgroup = outgroup, 
                                              hlspecies=c("CVAR",
                                                          "CSM3441",
                                                          "POW0123"), 
                                              hlcols = c("blue",
                                                         "red"), 
                                              main="HOG0000032 tree") 
#plot RERs
par(mfrow = c(1,1))
# Set the species of interest as the foreground:
foregroundSpecies <- foreground2Paths(c("CVAR",
                                        "CSM3441",
                                        "POW0123"),
                                      rerConvergeTrees,
                                      clade = "terminal")
plotRers(relativeEvolutionaryRates,
         "HOG0000032",
         phenv = foregroundSpecies)
newbend3rers <- treePlotRers(treesObj = rerConvergeTrees, 
                             rermat = relativeEvolutionaryRates, 
                             index = "HOG0000032", 
                             type = "c", 
                             nlevels = 3, 
                             figwid = 10)
# Do binary trait analysis:
# Select foreground species:
foregroundSpecies <- c("CVAR",
                       "CSM3441",
                       "POW0123")
# Set terminal branches leading to those species as the foreground:
foregroundBranches <- foreground2Tree(foregroundSpecies, 
                                      rerConvergeTrees, 
                                      clade = "terminal",
                                      useSpecies=c("CVAR", "CSM3685", "CSM3441", "POW0123", "POW0461", "CSM3677"))
# Generate paths describing trait evolution for gene trees that don't have all six species:
pathsForSubset <- foreground2Paths(foregroundSpecies, 
                                   rerConvergeTrees, 
                                   clade = "all")
# Correlate the trait with polymorphism:
correlationForPolymorphism <- correlateWithBinaryPhenotype(relativeEvolutionaryRates, 
                                                           pathsForSubset, 
                                                           min.sp = 5, 
                                                           min.pos = 1,
                                                           weighted = "auto")
head(correlationForPolymorphism[order(correlationForPolymorphism$P),])
plotRers(relativeEvolutionaryRates,
         "HOG0000032",
         phenv = pathsForSubset) 
hist(correlationForPolymorphism$p.adj, 
     breaks = 15, 
     xlab = "Adjusted Kendall P-value", 
     main = "Adjusted p-values for correlation between genes and polymorphism")

