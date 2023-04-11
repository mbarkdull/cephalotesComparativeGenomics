library(ape)
library(purrr)
library(tidyverse)
library(ggtree)
library(phytools)

# List all of the unlabelled tree files:
treeFiles <- list.files(path = "./allGeneTrees",
                                full.names = TRUE)

cleanTreePrefixes <- function(i) {
  # Read in one tree
  tree <- ape::read.tree(file = i)
  tree[["tip.label"]]
  
  # Get ride of everything before the second underscore (this looks like `translated_CSM3441_`)
  tree[["tip.label"]] <- sub("[^_]*_[^_]*_", 
                             "", 
                             tree[["tip.label"]])
  # Extract the orthogroup number:
  orthogroup <- str_split_i(string = i,
                            pattern = "/",
                            i = 3) %>%
    str_split_i(pattern = "_",
                i = 1)
  #Export the tree
  write.tree(tree, 
             file = paste("./allGeneTrees/cleaned_",
                          orthogroup,
                          "_tree.txt",
                          sep = ""))
}

possiblyCleanTreePrefixes <- possibly(cleanTreePrefixes,
                                      otherwise = "error")


purrr::map(treeFiles,
           possiblyCleanTreePrefixes)
