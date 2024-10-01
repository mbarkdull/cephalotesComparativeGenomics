# Script to make a pruned phylogeny of my focal species for Ch. 2:
# Load packages:
library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(diversitree)


# Info from: https://yulab-smu.top/treedata-book/chapter7.html

# Read in the data from Price, Etienne and Powell:
nodeTree <- read.tree(file = "./speciesTree/Price2022/Price_etal_TurtleAntLandBridgeUCEs_DryadUpload/TreeFiles/RAxMLtree.tre")

# Fix tips:
nodeTree[["tip.label"]] <- str_split_i(string = nodeTree[["tip.label"]], 
                                       pattern = "_",
                                       i = 2)

nodeTree[["tip.label"]] <- gsub(pattern = "^Ce|^C",
                                replacement = "Cephalotes_",
                                x = nodeTree[["tip.label"]])
nodeTree[["tip.label"]] <- gsub(pattern = "^Pr|^P",
                                replacement = "Procryptocerus_",
                                x = nodeTree[["tip.label"]])

# Read in soldier presence/absence data for all Cephalotes species:
allSpeciesPhenotypes <- read.csv("./speciesTree/AllCephalotesPhenotypesWithOutgroup.csv") %>%
  filter(Soldier != "?")
allSpeciesPhenotypes$Species <- gsub(pattern = "^C_",
                                     replacement = "Cephalotes_",
                                     x = allSpeciesPhenotypes$Species)

# Combine the large tree with all of the phenotypes:
allSpecies <- ggtree(nodeTree, 
                     layout = "circular") %<+% 
  allSpeciesPhenotypes 

# Define a color scheme:
colorScheme <- c(`Yes` = "#F1A208EE", `No` = "#06A77D")

# Plot the tree for all of the species:
PhenotypePhylogeny <- allSpecies + 
  geom_tiplab(size = 2, 
              offset = 0, 
              hjust = -0.15, 
              fontface = 'italic') +
  geom_tippoint(aes(color = Soldier), 
                size = 3) +
  scale_color_manual(values = colorScheme) +
  theme(legend.position = "left", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA),
        legend.background = element_rect(fill = "transparent",
                                         colour = NA),
        legend.key = element_rect(fill = "transparent",
                                  colour = NA)) 

plot(PhenotypePhylogeny)
#ggsave(filename = "talkToyWholeCephalotesPhylogeny.png", path = "./Plots/Phylogenies", bg = "transparent", width = 16, height = 9)

plainPhylogeny <- allSpecies + 
  geom_tiplab(size = 2, 
              offset = 0, 
              hjust = -0.15, 
              fontface = 'italic') +
  scale_color_manual(values = colorScheme) +
  theme(legend.position = "left", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA),
        legend.background = element_rect(fill = "transparent",
                                         colour = NA),
        legend.key = element_rect(fill = "transparent",
                                  colour = NA)) 
plot(plainPhylogeny)
#ggsave(filename = "plainCephalotesPhylogenies.png", path = "./Plots/Phylogenies/", bg = "transparent", width = 5, height = 8)


# Create the pruned tree for just my focal species:
# Get list of tips to keep:
tipsToKeep <- c("Cephalotes_varians", "Cephalotes_placidus", "Cephalotes_opacus", "Cephalotes_atratus", "Cephalotes_alfaroi", "Cephalotes_unimaculatus")
rawPrunedTree <- keep.tip(nodeTree, tip = tipsToKeep)
prunedTree <- rawPrunedTree

newTipLabels <- c("Cephalotes varians",
                  "Cephalotes unimaculatus",
                  "Cephalotes atratus",
                  "Cephalotes alfaroi",
                  "Cephalotes opacus",
                  "Cephalotes placidus")

# Replace the original tip labels with my new vector of tip labels:
prunedTree$tip.label <- newTipLabels

# Read in the phenotypes:
focalPhenotypes <- allSpeciesPhenotypes %>%
  filter(Species %in% tipsToKeep)
focalPhenotypes$Species <- gsub(pattern = "_",
                                replacement = " ",
                                x = focalPhenotypes$Species)

# Define a color scheme:
colorScheme <- c(`Yes` = "#F1A208EE", 
                 `No` = "#06A77D")

# Combine the tree and phenotype data with the special %<+% operator:
prunedBasePlot <- ggtree(prunedTree) %<+% focalPhenotypes 

# Reorder the levels of the Polymorphism column so that they won't be in alphabetical order in the legend:
prunedBasePlot[["data"]][["Soldier"]] <- factor(prunedBasePlot[["data"]][["Soldier"]], 
                                                     levels = c("Yes", "No"))

# Transparent background from: http://www.nagraj.net/notes/transparent-background-with-ggplot2/
prunedPhenotypePhylogeny <- prunedBasePlot + 
  geom_tiplab(offset = 0, 
              hjust = -0.15, 
              fontface = 'italic') +
  geom_tippoint(aes(fill = Soldier), 
                color = "black", 
                shape = 21,
                size = 6) + 
  scale_fill_manual(values = colorScheme) +
  theme(legend.position = "left", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA),
        legend.background = element_rect(fill = "transparent",
                                         colour = NA),
        legend.key = element_rect(fill = "transparent",
                                  colour = NA)) + 
  labs(title = "Phylogenetic relationships among focal species", 
       caption = "Phylogeny from Price et al. 2022.") +
  xlim(c(0, 0.04))

plot(prunedPhenotypePhylogeny)
ggsave(filename = "binaryPrunedCephalotesPhylogeny.png", 
       path = "./Plots/Phylogenies", 
       bg = "transparent", 
       width = 8, 
       height = 4.5)
