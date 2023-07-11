library(ape)
library(purrr)
library(tidyverse)
library(ggtree)

# Create an output directory:
dir.create("./07_labelledPhylogenies/")

# List all of the tree files:
orthofinderTreeFiles <- list.files(path = "./05_HOGTrees",
                                   full.names = TRUE)

# Write a function to label just the tips on a single tree, for a single species, with Hyphy:
# The reason we're doing one species at a time is so that internal nodes of species-specific gene families can get labelled, without labelling internal nodes with descendants from multiple foreground species. 
labellingTrees <- function(i, focalSpecies) {
  filename <- str_split_i(i,
                          pattern = "/",
                          i = 3)
  orthogroup <- str_split_i(filename,
                            pattern = "_",
                            i = 1)
  # Copy the unlabelled tree to the output folder; we'll be editing this tree three times. 
  file.copy(from = i,
            to = paste("./07_labelledPhylogenies/labelled_",
                       orthogroup,
                       "_tree.txt",
                       sep = ""))
  # Read in a single tree and get the tip labels of that tree:
  text <- readr::read_file(paste("./07_labelledPhylogenies/labelled_",
                                 orthogroup,
                                 "_tree.txt",
                                 sep = ""))
  # In case hyphy has already run on this, add back in the semicolon that ape requires:
  text <- gsub("\n",
               "",
               text,
               perl = TRUE)
  text <- gsub("$(?<!;)",
               ";",
               text,
               perl = TRUE)
  tree <- ape::read.tree(text = text)
  #plot(tree)
  tips <- as.data.frame(tree[["tip.label"]])
  tips 
  
  # Get the list of tip labels with matches to the species of interest:
  tipsToLabel <- filter(tips,
                        str_split_i(tips$`tree[["tip.label"]]`, pattern = "_", i = 1) %in% focalSpecies)
  tipsToLabel <- tipsToLabel$`tree[["tip.label"]]`
  tipsToLabel
  write(tipsToLabel, 
        paste("tipsToLabel_",
              focalSpecies,
              "_",
              filename,
              sep = ""))
  
  # Run the Hyphy labelling script:
  hyphyCommand <- paste("/programs/hyphy-20210923/bin/hyphy hyphy-analyses/LabelTrees/label-tree.bf --tree ./07_labelledPhylogenies/labelled_",
                        orthogroup,
                        "_tree.txt --list tipsToLabel_",
                        focalSpecies,
                        "_",
                        orthogroup,
                        "_tree.txt --output ./07_labelledPhylogenies/labelled_",
                        orthogroup,
                        "_tree.txt --internal-nodes \"All descendants\"",
                        sep = "")
  cat(hyphyCommand)
  system(hyphyCommand)
  file.remove(paste("tipsToLabel_",
                    focalSpecies,
                    "_",
                    filename,
                    sep = ""))
}

possiblylabellingTrees <- purrr::possibly(labellingTrees,
                                          otherwise = "Error")

# Set up future for paralellization:
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)
# Run the function separately for each focal species:
future_map2(orthofinderTreeFiles,
            "CSM3441",
            possiblylabellingTrees)
future_map2(orthofinderTreeFiles,
            "POW0123",
            possiblylabellingTrees)
future_map2(orthofinderTreeFiles,
            "CVAR",
            possiblylabellingTrees)

### Don't need to do this based on ASR; keeping code for reference. 
# Now that the tips are appropriately labelled, we can label the one internal node that gets reconstructed as being dimorphic:
#q <- "./07_labelledPhylogenies/labelled_HOG0002868_tree.txt"
#treeText <- readr::read_file(q)
#treeText <- paste(treeText, 
# ";",
#  sep = "")

#tree <- ape::read.tree(text = treeText)
#test <- ggtree(tree)
#plot(test)

# Find the two tips, CSM3441 and CSM3685, who's MRCA should be labelled:
#selected_tips <- grep("CSM3441|CSM3685", 
#    tree$tip.label)

## Finding the direct ancestor for each of these tips
#ancestors <- tree$edge[tree$edge[, 2] %in% selected_tips, 1]

## Finding the direct ancestor for each of these tips
#commonAncestors <- tree$edge[tree$edge[, 2] %in% ancestors, 1]
#length(unique(commonAncestors))




# For each node, find all of it's descendants:
# List all unique parent nodes:
#parentNodes <- unique(test[["data"]]$parent)
#relabellingNodes <- function(species) {
#  for (i in parentNodes) {
# Get the descendants of a particular parent, both their number and their label:
#   descendants <- filter(test[["data"]],
#                         parent == i) %>% 
#     select(node, label)
# Get the subset of descendants that match a focal species
#  selectedTips <- grep(species, 
#                       descendants$label)
# Check if all the descendant match a focal species:
#   length(selectedTips) == length(descendants$label)
# If so, label that particular parent node as foreground:
#   if(length(selectedTips) == length(descendants$label)){
#     print("Relabelling")
#     test[["data"]] <- test[["data"]] %>%
#      mutate(label = case_when(node == i ~ paste(species,
#"_", 
#                     node, 
#                     sep = ""),
#   TRUE ~ as.character(label)))
# } else {
#   print("Not all descendants are from the same species")
#   }
# }
#}


#relabellingNodes("POW0123")
#relabellingNodes("CSM3441")
#plot(test) + 
#  geom_tiplab(size = 3.5, 
#            offset = 0, 
#            hjust = 0, 
#           fontface = 'italic') + 
# geom_text(aes(label = node))

