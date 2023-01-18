# This script will generate a species tree for the species in this analysis, from a published, larger phylogeny. 
library(tidyverse)
library(ape)
library(phylotools)

# Read in the full tree:
fullTree <- read.tree(file = "./speciesTree/Cephalotes_node_calibration_PriceEtiennePowell2016.tre")

# Define the tips we want to keep:
tipsToKeep <- c("C_atratus",
                "C_alfaroi",
                "C_opacus",
                "C_placidus",
                "C_unimaculatus",
                "C_varians")

# Drop all other tips:
subsetTree <- keep.tip(fullTree,
                       tipsToKeep)

# Create a dataframe that links species names to their sequencing codes:
codes <- c("translated_CSM3677.fasta", 
           "translated_POW0123.fasta",
           "translated_CSM3441.fasta",
           "translated_CSM3685.fasta",
           "translated_POW0461.fasta",
           "translated_CVAR.fasta")
namesToCodes <- data.frame(tipsToKeep, 
                           codes)

# Replace the species names with the sequencing codes:
subsetTreeCodes <- phylotools::sub.taxa.label(subsetTree, 
                                 namesToCodes)

# Export that tree so it can be used with OrthoFinder:
ape::write.tree(subsetTreeCodes, 
                file = "./speciesTree/speciesCodesTree.txt")
