library(ape)
library(purrr)
library(tidyverse)
library(ggtree)

# List the species codes for our focal species:
focalSpecies <- c("CVAR", 
                  "POW0123", 
                  "CSM3441")

# Create an output directory:
dir.create("./labelledPhylogenies/")

# List all of the tree files:
orthofinderTreeFiles <- list.files(path = "./allGeneTrees",
                                   full.names = TRUE)
# Write a function to label a single tree with Hyphy:
labellingTrees <- function(i) {
  orthogroupNumber <- str_split_i(i,
                                  pattern = "/",
                                  i = 3)
  # Read in a single tree and get the tip labels of that tree:
  tree <- ape::read.tree(file = i)
  tips <- as.data.frame(tree[["tip.label"]])
  
  # Get the list of tip labels with matches to the species of interest:
  tipsToLabel <- filter(tips,
                        str_split_i(tips$`tree[["tip.label"]]`, pattern = "_", i = 2) %in% focalSpecies)
  tipsToLabel <- tipsToLabel$`tree[["tip.label"]]`
  tipsToLabel
  write(tipsToLabel, "tipsToLabel.txt")
  
  # Run the Hyphy labelling script:
  hyphyCommand <- paste("/programs/hyphy-20210923/bin/hyphy hyphy-analyses/LabelTrees/label-tree.bf --tree ",
                        i,
                        " --list tipsToLabel.txt --output ./labelledPhylogenies/labelled_",
                        orthogroupNumber,
                        " --internal-nodes \"Parsimony\"",
                        sep = "")
  cat(hyphyCommand)
  system(hyphyCommand)
}

# Apply that function to all trees:
possiblyLabellingTrees <- possibly(labellingTrees, 
                                   otherwise = "Error")
purrr::map(orthofinderTreeFiles,
           possiblyLabellingTrees)