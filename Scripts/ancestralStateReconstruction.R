# Doing ancestral state reconstruction with an outgroup
# Load packages:
library(ape)
library(phytools)
library(maps)
library(tidyverse)
library(ggtree)
library(diversitree)
library(corHMM)
library(plotrix)

#### Analysis with just Cephalotes and Procryptocerus: ####
#### Combine the Powell 2020 tree with the Nelsen 2018 tree that has our outgroup: ####
# Read in the tree used in Powell, Price and Kronauer 2020:
targetTree <- read.tree(file = "./speciesTree/Cephalotes_node_calibration.tre")
ggtree(targetTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') +
  xlim(-15, 50)

# Read in the tree from Nelsen et al. 2018. This tree has many species, including many Cephalotes species AND Procryptocerus, the sister genus. 
nelsenTree <- read.tree(file = "./nelsen2018Tree/Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
ggtree(nelsenTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') 


nelsenTree[["tip.label"]]

# Get a list of species names from the Powell 2020 paper, formatted to match the tips from Nelsen:
targetSpecies <- str_replace(targetTree[["tip.label"]], 
                             pattern = "C_", 
                             replacement = "Cephalotes_")

# Check that there is some overlap of species in the two trees:
targetSpecies %in% nelsenTree[["tip.label"]]

# Get Cephalotes and Procryptocerus species that are in the Nelsen tree, then prune the tree to just those tips:
filterTips <- nelsenTree[["tip.label"]] %>% 
  str_subset(pattern = "Cephalotes|Procryptocerus")
trimmedNelsenTree <- keep.tip(nelsenTree,
                              tip = filterTips)
# Plot that tree with node numbers:
ggtree(trimmedNelsenTree) + 
  geom_tiplab(size = 3.5, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  geom_text(aes(label = node)) +
  xlim(-15, 200)

# So I want to take everything from node 67 down, the Cephalotes MRCA, and replace it with the Powell tree.
# What is the length of the Nelsen tree from node 67 (the Cephalotes MRCA) to the tips (e.g. node 62, a tip)?
lengthTrimmedNelsenTree <- castor::get_pairwise_distances(trimmedNelsenTree, 
                                                          c(67), 
                                                          c(62), 
                                                          as_edge_counts = FALSE, 
                                                          check_input = TRUE)
# Rescale the Powell tree to be of that total length: (phytools::rescale)
library(geiger)
rescaledTargetTree <- phytools::rescale(targetTree, 
                                        model = c("depth"), 
                                        depth = lengthTrimmedNelsenTree)
# Check to make sure that it's been rescaled:
castor::get_pairwise_distances(rescaledTargetTree, 
                               c(116), 
                               c(115), 
                               as_edge_counts = FALSE, 
                               check_input = TRUE)

# Now combine the Powell tree with the trimmed Nelsen tree. This tree will have some duplicated Cephalotes species that need to be dropped. 
# Combine the trees:
duplicatesCombinedTree <- ape::bind.tree(trimmedNelsenTree, 
                                         rescaledTargetTree, 
                                         where = 67, 
                                         position = 0, 
                                         interactive = FALSE)
# Plot the combined tree:
ggtree(duplicatesCombinedTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)

# Drop Cephalotes tips that came from the Nelsen tree (e.g. those with the genus spelled out), to remove any duplicates.
tipsToDrop <- duplicatesCombinedTree[["tip.label"]] %>% 
  str_subset(pattern = "Cephalotes")
combinedTree <- drop.tip(duplicatesCombinedTree,
                         tipsToDrop)
ggtree(combinedTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)
# Fix the "C_" in the tree to be "Cephalotes_":
combinedTree[["tip.label"]] <- str_replace(combinedTree[["tip.label"]], 
                                           pattern = "C_", 
                                           replacement = "Cephalotes_")
ggtree(combinedTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)

#### Prep phenotype data for ancestral state reconstructions: ####
## Read in the phenotype data:
casteData <- read.csv("./speciesTree/AllCephalotesPhenotypesWithOutgroup.csv",
                      row.names = 1,
                      stringsAsFactors = TRUE) %>%
  filter(Soldier != "?")
rownames(casteData) <- str_replace(rownames(casteData), 
                                   pattern = "C_", 
                                   replacement = "Cephalotes_")
casteData
casteData$test <- rownames(casteData)

## Filter the tree to only species with phenotype data:
combinedTreeData <- keep.tip(combinedTree, 
                             tip = casteData$test)

## Extract soldier presence/absence as a vector:
casteSystem <- setNames(casteData$Soldier,
                        rownames(casteData))
casteSystem <- factor(casteSystem, 
                      levels = c("No",
                                 "Yes"))

## Set colors for plotting
plotColors <- setNames(c("#06A77D",
                         "#F1A208EE"),
                       levels(casteSystem))
## Plot the tree & data
plotTree.datamatrix(combinedTreeData,
                    as.data.frame(casteSystem),
                    colors = list(plotColors),
                    header = FALSE,
                    fsize = 0.45)
## Add a legend
legend("topright",
       legend = levels(casteSystem),
       pch = 22,
       pt.cex = 1.5,
       pt.bg = plotColors,
       bty = "n", 
       cex = 0.8)

#### Assess models for ancestral state reconstructions: ####
# Now we need to assess which model we want to use.
## Fit an all equal-rates (ER) model:
equalRatesModel <- fitMk(tree = combinedTreeData,
                         x = casteSystem,
                         model = "ER")
## Fit an all-rates-different model:
allRatesDifferentModel <- fitMk(tree = combinedTreeData,
                                x = casteSystem,
                                model = "ARD")
## Fit a model that allows an irreversible transition from no soldier to yes soldier:
fit01 <- fitMk(tree = combinedTreeData,
               x = casteSystem,
               model = matrix(c(0,
                                1,
                                0,
                                0),
                              2,
                              2,
                              byrow = TRUE))
## Fit a model that allows an irreversible transition from yes soldier to no soldier:
fit10 <- fitMk(tree = combinedTreeData,
               x = casteSystem,
               model = matrix(c(0,
                                0,
                                1,
                                0),
                              2,
                              2,
                              byrow = TRUE))
## extract AIC values for each model
aicValues <- c(AIC(equalRatesModel),
               AIC(allRatesDifferentModel),
               AIC(fit01),
               AIC(fit10))
## Summarize the AIC values for each model in a data frame (looks like fit01 is best model):
data.frame(model = c("equalRatesModel",
                     "allRatesDifferentModel",
                     "fit01",
                     "fit10"),
           logL = c(logLik(equalRatesModel),
                    logLik(allRatesDifferentModel),
                    logLik(fit01),
                    logLik(fit10)),
           AIC = aicValues,
           delta.AIC = aicValues-min(aicValues))

## Create new data frame of polymorphism phenotypes for corHMM:
casteData <- data.frame(species = names(casteSystem),
                        casteSystem = as.numeric(casteSystem) - 1)
head(casteData)

## Estimate marginal ancestral states under the 01 model
marginalEstimationIrreversible <- corHMM(phy = combinedTreeData,
                                         data = casteData,
                                         node.states = "marginal",
                                         rate.cat = 1,
                                         rate.mat = matrix(c(0,
                                                             1,
                                                             0,
                                                             0),
                                                           2,
                                                           2))

## Estimate marginal ancestral states under the ARD model:
marginalEstimation <- corHMM(phy = combinedTreeData,
                             data = casteData,
                             node.states = "marginal",
                             rate.cat = 1,
                             model = "ARD")

marginalEstimation                    
head(marginalEstimation$states)      

# Plot the tree & data
# Get out the inferred ancestral states:
inferredStates <- as.data.frame(marginalEstimation$states)
inferredStates$node <- rownames(inferredStates)

# Make pie charts for node inferences:
pies <- nodepie(inferredStates, 
                cols = 1:2, 
                color = c("darkorange1",
                          "blue"), 
                alpha = 0.8)

# Try to plot the ACR with ggtree (so far can't figure out nodes):
test <- ggtree(combinedTreeData) + 
  geom_tiplab(size = 2.5, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)
ggtree::inset(tree_view = test,
                     insets = pies,
                     width = 0.1, 
                     height = 0.1, 
                     x = "node")

# Plot the ACR with phytools plotting:
plotTree.datamatrix(combinedTreeData,
                    as.data.frame(casteSystem),
                    colors = list(plotColors),
                    header = FALSE,
                    fsize = 0.45)
# Add a legend:
legend("topright",
       legend = levels(casteSystem),
       pch = 22,
       pt.cex = 1.5,
       pt.bg = plotColors,
       bty = "n",
       cex = 0.8)

# Add the node labels showing marginal ancestral states from corrHM:
nodelabels(pie = marginalEstimation$states,
           piecol = plotColors,
           cex = 0.25)

#### Do ancestral state reconstruction using stochastic character mapping: ####
# Run 1,000 stochastic character maps:
# Equal rates model:
stochasticCharacterMaps <- make.simmap(tree = combinedTreeData,
                                       x = casteSystem,
                                       model = "ER",
                                       nsim = 1000,
                                       Q = "mcmc",
                                       vQ = 0.01,
                                       prior = list(use.empirical = TRUE),
                                       samplefreq = 10)

# All rates different model:
ARD <- make.simmap(combinedTreeData,
                   casteSystem,
                   model = "ARD",
                   nsim = 1000,
                   Q = "mcmc")

# Summarize all of the runs under the all-rates-different model:
ardSummary <- summary(ARD)
ardSummary

## create a plot showing tip states:
png(filename = "./Plots/tipStatesOnly.png",  
    width = 6, 
    height = 6, 
    units = "in", 
    res = 600)
par(fg = "black")
plot(ardSummary,
     colors = plotColors,
     fsize = 0.5,
     ftype = "i",
     lwd = 1,
     offset = 4,
     #ylim = c(-1,Ntip(combinedTreeData)),
     cex = c(0, 
             0.3),
     outline = FALSE,
     type="fan")
dev.off()

## create a plot showing posterior probabilities at all nodes:
png(filename = "./Plots/AncestralStateReconstruction.png",  
    width = 6, 
    height = 6, 
    units = "in", 
    res = 600)
par(fg = "black")
plot(ardSummary,
     colors = plotColors,
     fsize = 0.5,
     ftype = "i",
     lwd = 1,
     offset = 5,
     #ylim = c(-1,Ntip(combinedTreeData)),
     cex = c(0.3, 
             0.3),
     type="fan")

## add a legend
legend("bottomleft",
       legend = levels(casteSystem),
       pch = 22,
       pt.cex = 1.5,
       pt.bg = plotColors,
       bty = "n",
       cex = 0.8)
dev.off()

# Plot not in a circular way:
plot(ardSummary,
     colors = plotColors,
     fsize = 0.5,
     ftype = "i",
     lwd = 1,
     offset = 5,
     #ylim = c(-1,Ntip(combinedTreeData)),
     cex = c(0.2, 
             0.2))




#### Analysis with Cephalotes + Procryptocerus and Pheidole: ####
#### Combine the Powell 2020 tree with the Nelsen 2018 tree that has our outgroup: ####
# Read in the tree used in Powell, Price and Kronauer 2020:
targetTree <- read.tree(file = "./speciesTree/Cephalotes_node_calibration.tre")
ggtree(targetTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') +
  xlim(-15, 50)

# Read in the tree from Nelsen et al. 2018. This tree has many species, including many Cephalotes species, Procryptocerus, the sister genus, and Pheidole. 
nelsenTree <- read.tree(file = "./nelsen2018Tree/Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
ggtree(nelsenTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') 

nelsenTree[["tip.label"]]

# Get a list of species names from the Powell 2020 paper, formatted to match the tips from Nelsen:
targetSpecies <- str_replace(targetTree[["tip.label"]], 
                             pattern = "C_", 
                             replacement = "Cephalotes_")

# Check that there is some overlap of species in the two trees:
targetSpecies %in% nelsenTree[["tip.label"]]

# Get Cephalotes, Procryptocerus, and Pheidole species that are in the Nelsen tree, then prune the tree to just those tips:
filterTips <- nelsenTree[["tip.label"]] %>% 
  str_subset(pattern = "Cephalotes|Procryptocerus|Pheidole")

trimmedNelsenTree <- keep.tip(nelsenTree,
                              tip = filterTips)

# Plot that tree with node numbers:
ggtree(trimmedNelsenTree) + 
  geom_tiplab(size = 2.5, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  geom_text(aes(label = node),
            size = 2.5) +
  xlim(-15, 200)

# So I want to take everything from node 162 down, the Cephalotes MRCA, and replace it with the Powell tree.
# What is the length of the Nelsen tree from node 162 (the Cephalotes MRCA) to the tips (e.g. node 48, a tip)?
lengthTrimmedNelsenTree <- castor::get_pairwise_distances(trimmedNelsenTree, 
                                                          c(162), 
                                                          c(48), 
                                                          as_edge_counts = FALSE, 
                                                          check_input = TRUE)
# Rescale the Powell tree to be of that total length:
rescaledTargetTree <- phytools::rescale(targetTree, 
                                        model = c("depth"), 
                                        depth = lengthTrimmedNelsenTree)
# Check to make sure that it's been rescaled:
castor::get_pairwise_distances(rescaledTargetTree, 
                               c(116), 
                               c(115), 
                               as_edge_counts = FALSE, 
                               check_input = TRUE)

# Now combine the Powell tree with the trimmed Nelsen tree. This tree will have some duplicated Cephalotes species that need to be dropped. 
# Combine the trees:
duplicatesCombinedTree <- ape::bind.tree(trimmedNelsenTree, 
                                         rescaledTargetTree, 
                                         where = 162, 
                                         position = 0, 
                                         interactive = FALSE)
# Plot the combined tree:
ggtree(duplicatesCombinedTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)

# Drop Cephalotes tips that came from the Nelsen tree (e.g. those with the genus spelled out), to remove any duplicates.
tipsToDrop <- duplicatesCombinedTree[["tip.label"]] %>% 
  str_subset(pattern = "Cephalotes")
combinedTree <- drop.tip(duplicatesCombinedTree,
                         tipsToDrop)
ggtree(combinedTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)
# Fix the "C_" in the tree to be "Cephalotes_":
combinedTree[["tip.label"]] <- str_replace(combinedTree[["tip.label"]], 
                                           pattern = "C_", 
                                           replacement = "Cephalotes_")
ggtree(combinedTree) + 
  geom_tiplab(size = 2, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)

#### Prep phenotype data for ancestral state reconstructions: ####
## Read in the phenotype data:
casteData <- read.csv("./speciesTree/AllCephalotesPhenotypesWithPheidole.csv",
                      row.names = 1,
                      stringsAsFactors = TRUE) %>%
  filter(Soldier != "?")
rownames(casteData) <- str_replace(rownames(casteData), 
                                   pattern = "C_", 
                                   replacement = "Cephalotes_")
casteData
casteData$test <- rownames(casteData)

## Filter the tree to only species with phenotype data:
combinedTreeData <- keep.tip(combinedTree, 
                             tip = casteData$test)

## Extract soldier presence/absence as a vector:
casteSystem <- setNames(casteData$Soldier,
                        rownames(casteData))
casteSystem <- factor(casteSystem, 
                      levels = c("No",
                                 "Yes"))

## Plot the tree & data
plotTree.datamatrix(combinedTreeData,
                    as.data.frame(casteSystem),
                    colors = list(plotColors),
                    header = FALSE,
                    fsize = 0.45)
## Add a legend
legend("topright",
       legend = levels(casteSystem),
       pch = 22,
       pt.cex = 1.5,
       pt.bg = plotColors,
       bty = "n", 
       cex = 0.8)

#### Assess models for ancestral state reconstructions: ####
# Now we need to assess which model we want to use.
## Fit an all equal-rates (ER) model:
equalRatesModel <- fitMk(tree = combinedTreeData,
                         x = casteSystem,
                         model = "ER")
## Fit an all-rates-different model:
allRatesDifferentModel <- fitMk(tree = combinedTreeData,
                                x = casteSystem,
                                model = "ARD")
## Fit a model that allows an irreversible transition from no soldier to yes soldier:
fit01 <- fitMk(tree = combinedTreeData,
               x = casteSystem,
               model = matrix(c(0,
                                1,
                                0,
                                0),
                              2,
                              2,
                              byrow = TRUE))
## Fit a model that allows an irreversible transition from yes soldier to no soldier:
fit10 <- fitMk(tree = combinedTreeData,
               x = casteSystem,
               model = matrix(c(0,
                                0,
                                1,
                                0),
                              2,
                              2,
                              byrow = TRUE))
## extract AIC values for each model
aicValues <- c(AIC(equalRatesModel),
               AIC(allRatesDifferentModel),
               AIC(fit01),
               AIC(fit10))
## Summarize the AIC values for each model in a data frame (looks like fit01 is best model):
data.frame(model = c("equalRatesModel",
                     "allRatesDifferentModel",
                     "fit01",
                     "fit10"),
           logL = c(logLik(equalRatesModel),
                    logLik(allRatesDifferentModel),
                    logLik(fit01),
                    logLik(fit10)),
           AIC = aicValues,
           delta.AIC = aicValues-min(aicValues))

## Create new data frame of polymorphism phenotypes for corHMM:
casteData <- data.frame(species = names(casteSystem),
                        casteSystem = as.numeric(casteSystem) - 1)
head(casteData)

## Estimate marginal ancestral states under the 01 model
marginalEstimationIrreversible <- corHMM(phy = combinedTreeData,
                                         data = casteData,
                                         node.states = "marginal",
                                         rate.cat = 1,
                                         rate.mat = matrix(c(0,
                                                             1,
                                                             0,
                                                             0),
                                                           2,
                                                           2))

## Estimate marginal ancestral states under the ARD model:
marginalEstimation <- corHMM(phy = combinedTreeData,
                             data = casteData,
                             node.states = "marginal",
                             rate.cat = 1,
                             model = "ARD")

marginalEstimation                    
head(marginalEstimation$states)      

# Plot the tree & data
# Get out the inferred ancestral states:
inferredStates <- as.data.frame(marginalEstimation$states)
inferredStates$node <- rownames(inferredStates)

# Make pie charts for node inferences:
pies <- nodepie(inferredStates, 
                cols = 1:2, 
                color = c("darkorange1",
                          "blue"), 
                alpha = 0.8)

# Try to plot the ACR with ggtree (so far can't figure out nodes):
test <- ggtree(combinedTreeData) + 
  geom_tiplab(size = 2.5, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  xlim(-15, 200)
ggtree::inset(tree_view = test,
              insets = pies,
              width = 0.1, 
              height = 0.1, 
              x = "node")

# Plot the ACR with phytools plotting:
plotTree.datamatrix(combinedTreeData,
                    as.data.frame(casteSystem),
                    colors = list(plotColors),
                    header = FALSE,
                    fsize = 0.45)
# Add a legend:
legend("topright",
       legend = levels(casteSystem),
       pch = 22,
       pt.cex = 1.5,
       pt.bg = plotColors,
       bty = "n",
       cex = 0.8)

# Add the node labels showing marginal ancestral states from corrHM:
nodelabels(pie = marginalEstimation$states,
           piecol = plotColors,
           cex = 0.25)

#### Do ancestral state reconstruction using stochastic character mapping: ####
# Run 1,000 stochastic character maps:
# Equal rates model:
stochasticCharacterMaps <- make.simmap(tree = combinedTreeData,
                                       x = casteSystem,
                                       model = "ER",
                                       nsim = 1000,
                                       Q = "mcmc",
                                       vQ = 0.01,
                                       prior = list(use.empirical = TRUE),
                                       samplefreq = 10)

# All rates different model:
ARD <- make.simmap(combinedTreeData,
                   casteSystem,
                   model = "ARD",
                   nsim = 1000,
                   Q = "mcmc")

# Summarize all of the runs under the all-rates-different model:
ardSummary <- summary(ARD)
ardSummary

# Plot the tree with node numbers so that we can rotate basal clades to be together (rotate around node 126):
# Plot that tree with node numbers:
ggtree(combinedTree) + 
  geom_tiplab(size = 3.5, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  geom_text(aes(label = node)) +
  xlim(-15, 200)

#rotatedPD <- purrr::map(ardSummary[["tree"]],
# ~rotateNodes(.x, nodes = 126))
#rotatedPDObject <- ardSummary
#rotatedPDObject[["tree"]] <- rotatedPD

## create a plot showing tip states:
png(filename = "./Plots/tipStatesOnlyPheidole.png",  
    width = 6, 
    height = 6, 
    units = "in", 
    res = 600)
par(fg = "black")
plot(ardSummary,
     colors = plotColors,
     fsize = 0.5,
     ftype = "i",
     lwd = 1,
     offset = 4,
     #ylim = c(-1,Ntip(combinedTreeData)),
     cex = c(0, 
             0.3),
     outline = FALSE,
     type="fan")
dev.off()

## create a plot showing posterior probabilities at all nodes:
png(filename = "./Plots/AncestralStateReconstructionWithPheidole.png",  
    width = 6, 
    height = 6, 
    units = "in", 
    res = 600)
par(fg = "black")
plot(ardSummary,
     colors = plotColors,
     fsize = 0.5,
     ftype = "i",
     lwd = 1,
     offset = 5,
     #ylim = c(-1,Ntip(combinedTreeData)),
     cex = c(0.3, 
             0.3),
     type="fan")
## add a legend
legend("bottomleft",
       legend = levels(casteSystem),
       title = "Soldier present?",
       pch = 22,
       pt.cex = 1.5,
       pt.bg = plotColors,
       bty = "n",
       cex = 0.8)
dev.off()

# Plot not in a circular way:
plot(ardSummary,
     colors = plotColors,
     fsize = 0.5,
     ftype = "i",
     lwd = 1,
     offset = 5,
     #ylim = c(-1,Ntip(combinedTreeData)),
     cex = c(0.2, 
             0.2))
