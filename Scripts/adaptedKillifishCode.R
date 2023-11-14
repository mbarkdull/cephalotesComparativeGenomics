library(tidyverse)
library(rphast)

#### Create output directories ####
dir.create("./16_phastConsAnalyses/")
dir.create("./16_phastConsAnalyses/singleScaffoldAlignments/")
dir.create("./16_phastConsAnalyses/acceleratedElements/")

#### Get alignments for individual scaffolds: ####
# Read in the alignment, produced by ./Scripts/alignmentsForPhast, as a dataframe:
alignment <- read.table("./12_genomeAlignments/orderedCactusOutput.maf", 
                        header = FALSE, 
                        col.names = paste0("V",
                                           seq_len(7)), 
                        fill = TRUE)
# Get rid of the excess species codes that are prepended to the sequence names:
alignment$V2 <- gsub(pattern = "fasta(.*)_",
                     replacement = "fasta.CVAR_scaf_",
                     alignment$V2)
# Create a column with just the scaffold names:
alignment$scaffold <- str_split_i(alignment$V2,
                                  pattern = "\\.",
                                  i = 3)
# Get every single scaffold number:
listOfScaffolds <- alignment$scaffold %>%
  unique()

filteringAlignments <- function(scaffoldNumber) {
  exportFile <- paste("./16_phastConsAnalyses/singleScaffoldAlignments/",
                      scaffoldNumber,
                      ".maf",
                      sep = "")
  if (file.exists(exportFile) == TRUE) {
    print("already filtered")
  } else {
    # Filter to a single scaffold and reformat:
    singleAlignment <- alignment %>%
      tidyr::fill(scaffold,
                  .direction = "up") %>%
      filter(scaffold == scaffoldNumber) %>% 
      dplyr::select(-c(scaffold))
    
    # Add the appropriate header row:
    headerRow <- c("##maf version=1",
                   "",
                   "",
                   "",
                   "",
                   "",
                   "")
    export <- rbind(headerRow, 
                    singleAlignment)
    
    # Export the .maf for a single scaffold:
    write_delim(export, 
                file = paste("./16_phastConsAnalyses/singleScaffoldAlignments/",
                             scaffoldNumber,
                             ".maf",
                             sep = ""),
                na = "",
                delim = " ",
                quote = "none",
                col_names = FALSE)
  }
}

possiblyFilteringAlignments <- possibly(filteringAlignments,
                                        otherwise = "Error.")
purrr::map(listOfScaffolds,
           possiblyFilteringAlignments)

#### Identify accelerated elements: ####
# Get a list of all scaffold numbers:
length(listOfScaffolds)
listOfScaffoldNumbers <- str_split_i(listOfScaffolds,
                                     pattern = "_",
                                     i = 3)

gettingAcceleratedRegions <- function(scaffoldNumber) {
  # Set the chromosome name
  scaffoldName <- paste("CVAR_CVAR_scaf_", 
                        scaffoldNumber,
                        sep = "")
  
  # Set the reference species
  referenceSpecies <- 'CVAR_genome'
  
  # Define the species tree
  tree <- "((POW0461_genome,CVAR_genome),((CSM3685_genome,CSM3441_genome),(POW0123_genome,CSM3677_genome)))"
  
  # Define the tree with only dimorphic species
  dimorphicTree <- "(CVAR_genome,(CSM3441_genome,CSM3677_genome))"
  dimorphicSpecies <- str_split(dimorphicTree, pattern = "\\,|\\(|\\)")
  dimorphicSpecies <- dimorphicSpecies[[1]] %>%
    stringi::stri_remove_empty()
  
  # Define the tree with only monomorphic species
  monomorphicTree <- "(POW0461_genome,(CSM3685_genome,POW0123_genome))"
  monomorphicSpecies <- str_split(monomorphicTree, pattern = "\\,|\\(|\\)")
  monomorphicSpecies <- monomorphicSpecies[[1]] %>%
    stringi::stri_remove_empty()
  
  # Read in the annotation for one chromosome
  features <- read.feat("./CVAR.gp") %>%
    filter(seqname == scaffoldName)
  # Set the sequence names
  features$seqname <- referenceSpecies;
  table(features$feature)
  
  #### Find regions conserved in the foreground and accelerated in the background: ####
  findingConservedAcceleratedRegions <- function(foregroundSpecies, backgroundSpecies, groupAcceleratedIn) {
    # Read in an alignment with only the foreground species
    if ("CVAR_genome" %in% foregroundSpecies == TRUE) {
      foregroundAlignmentWithReference <- read.msa(filename = paste("./16_phastConsAnalyses/singleScaffoldAlignments/CVAR_scaf_",
                                                                    scaffoldNumber,
                                                                    ".maf",
                                                                    sep = ""), 
                                                   format = "MAF", 
                                                   ordered = T, 
                                                   seqnames = c(foregroundSpecies))
    } else {
      foregroundAlignmentWithReference <- read.msa(filename = paste("./16_phastConsAnalyses/singleScaffoldAlignments/CVAR_scaf_",
                                                                    scaffoldNumber,
                                                                    ".maf",
                                                                    sep = ""), 
                                                   format = "MAF", 
                                                   ordered = T, 
                                                   seqnames = c(foregroundSpecies,
                                                                "CVAR_genome"))
    }
    # Get the fourfold degenerate sites
    fourfoldSitesForeground <- get4d.msa(foregroundAlignmentWithReference, 
                                         features)
    fourfoldSitesForeground <- fourfoldSitesForeground[foregroundSpecies, ]
    # Check the names of the fourfold sites
    print(paste("Sequences in the foreground fourfold sites are ",
                (names(fourfoldSitesForeground)),
                sep = ""))
    
    # Estimate a neutral model from fourfold degenerate sites
    foregroundNeutralModel <- phyloFit(fourfoldSitesForeground, 
                                       tree = tree, 
                                       subst.mod = "REV")
    
    # PhastCons in the R rphast package (Hubisz et al., 2011) 
    #was used to identify conserved genomic regions in the annual species, omitting the nonannual.
    # Run phastcons on only the foreground species, estimating all parameters
    foregroundAlignment <- foregroundAlignmentWithReference[foregroundSpecies, ]
    print(paste("Running phastcons on only the ",
                groupAcceleratedIn,
                " species, estimating all parameters",
                sep = ""))
    foregroundPhastCons <- phastCons(msa = foregroundAlignment, 
                                     mod = foregroundNeutralModel, 
                                     viterbi = TRUE, 
                                     estimate.transitions = TRUE, 
                                     estimate.expected.length = TRUE) 
    # estimate.rho = TRUE)    # estimating rho fails for me
    # Reset the sequence names to be the CVAR genome, since that's the reference and comes up in downstream steps
    foregroundPhastCons$most.conserved$seqname <- "CVAR_genome"
    # Read in an alignment with all of the species
    allSpeciesAlignment <- read.msa(filename = paste("./16_phastConsAnalyses/singleScaffoldAlignments/CVAR_scaf_",
                                                     scaffoldNumber,
                                                     ".maf",
                                                     sep = ""), 
                                    format = "MAF", 
                                    ordered = T)
    
    # Get fourfold degenerate sites for all species:
    fourfoldSitesAll <- get4d.msa(allSpeciesAlignment, 
                                  features)
    
    # And estimate a neutral model:
    allSpeciesNeutralModel <- phyloFit(fourfoldSitesAll, 
                                       tree = tree, 
                                       subst.mod = "REV")
    
    # Check which sites in the all-species alignment have representation of all three background species:
    hasBackgroundSpecies <- informative.regions.msa(allSpeciesAlignment,
                                                    3,
                                                    backgroundSpecies)
    
    # Check which sites in the all-species alignment have representation of 5/6 species:
    hasAtLeast5 <- informative.regions.msa(allSpeciesAlignment, 
                                           5)
    
    # Get the intersection of the regions that are most conserved in foreground, have the background, and has 5/6 species:
    # These are the elements to use in downstream bootstrapping tests
    informativeElements <- coverage.feat(foregroundPhastCons$most.conserved, 
                                         hasBackgroundSpecies, 
                                         hasAtLeast5, 
                                         get.feats = TRUE)
    
    # Filter out informative elements that overlap coding elements:
    codingElements <- filter(features,
                             feature == "CDS")
    noncodingInformativeElements <- coverage.feat(informativeElements,
                                                  codingElements,
                                                  not = c(FALSE,
                                                          TRUE), 
                                                  get.feats = TRUE)
    
    # See what proportion of the most conserved foreground regions are also those informative elements:
    coverage.feat(noncodingInformativeElements)/coverage.feat(foregroundPhastCons$most.conserved)
    
    # Define the length into which conserved elements should be split:
    splitLength <- 50
    
    # Split the elements:
    splitElements <- split.feat(noncodingInformativeElements, 
                                f = splitLength, 
                                drop = TRUE)
    
    # Calculate conservation/acceleration p-values on an alignment and evolutionary model.
    # See how the 50bp conserved elements are evolving in foreground species:
    print(paste("Seeing how the 50bp conserved elements are evolving in ",
                groupAcceleratedIn,
                sep = ""))
    observedPhyloPForeground <- phyloP(allSpeciesNeutralModel, 
                                       msa = allSpeciesAlignment, 
                                       mode = "ACC", 
                                       features = splitElements, 
                                       branches = foregroundSpecies)
    # See how 50bp conserved elements are evolving in background species:
    print("Seeing how 50bp conserved elements are evolving in background species:")
    observedPhyloPBackground <- phyloP(allSpeciesNeutralModel, 
                                       msa = allSpeciesAlignment, 
                                       mode ="ACC", 
                                       features = splitElements, 
                                       branches = backgroundSpecies)
    
    #do non parametric bootstrap:
    # Get the sequences of the informative elements out of the all-species alignment
    informativeElementSequences <- extract.feature.msa(copy.msa(allSpeciesAlignment), 
                                                       noncodingInformativeElements, 
                                                       pointer.only = TRUE)
    # Set the number of repetitions to do:
    nrep <- 100000
    
    # Create a simulated alignment by sampling from the informative element sequences, with replacement:
    simulatedMSA <- sample.msa(x = informativeElementSequences, 
                               size = nrep*splitLength, 
                               replace = TRUE)
    
    # produce features allowing long alignment to be interpreted as concatenation of shorter alignments
    # Create a vector starting at one and incrementing by units of splitLength (50) up to nrep (1000)
    startIdx <- seq(from = 1, 
                    by = splitLength, 
                    length.out = nrep)
    
    # Create simulated features on the simulated alignment:
    simulatedFeatures <- feat(seqname = names.msa(simulatedMSA)[1], 
                              src = "sim", 
                              feat = ".", 
                              start = startIdx, 
                              end = startIdx+splitLength-1)
    
    # On those simulated features, see how 50bp elements are evolving in the fore- and background:
    print("Seeing how 50bp elements are evolving in the fore- and background from simulated alignments")
    simulatedPhyloPForeground <- phyloP(allSpeciesNeutralModel, 
                                        msa = simulatedMSA, 
                                        mode = "ACC", 
                                        features = simulatedFeatures, 
                                        branches = foregroundSpecies)
    simulatedPhyloPBackground <- phyloP(allSpeciesNeutralModel, 
                                        msa = simulatedMSA, 
                                        mode = "ACC", 
                                        features = simulatedFeatures, 
                                        branches = backgroundSpecies)
    
    # Compare the simulated likelihood ratios to the observed likelihood ratios:
    # Plotting likelihood ratios 
    plotObservedToSimulated <- function(simulatedPhyloP, observedPhyloP) {
      layout(matrix(c(1, 2), 
                    nrow = 2), 
             heights =c (0.7, 0.3))
      par(mar = c(4.5, 4, 4, 2), 
          mgp = c(2.5 ,1, 0), 
          cex.axis = 1.5, 
          cex.lab = 1.5)
      qqplot(simulatedPhyloP$lnlratio,
             observedPhyloP$lnlratio,
             xlim = c(0, 15),
             ylim = c(0, 15), 
             xlab = "Simulated likelihood ratio",
             ylab = "Observed likelihood ratio")
      abline(0, 1, lty = 2)
      par(mar = c(4, 4, 1, 2))
      plot(density(observedPhyloP$lnlratio,
                   adjust = 3), 
           lty = 1,
           xlim = c(0, 15),
           xlab = "Likelihood Ratio",
           ylab = "Density",
           main = "", 
           col = "red")
      lines(density(simulatedPhyloP$lnlratio,
                    adjust = 3), 
            lty = 1,
            col = "black",
            xlim = c(0, 15))
    }
    
    # x is every element of the observed phyloP likelihood ratios
    # dist is the entire vector of simulated likelihood ratios.
    estimateEmpiricalPValues <- function(x, dist) {
      sum(x <= dist)/length(dist)
    }
    
    #plotObservedToSimulated(simulatedPhyloP = simulatedPhyloPForeground, 
    #observedPhyloP = observedPhyloPForeground)
    #plotObservedToSimulated(simulatedPhyloP = simulatedPhyloPBackground, 
    #observedPhyloP = observedPhyloPBackground )
    
    estimedPValuesForeground <- sapply(observedPhyloPForeground$lnlratio, 
                                       estimateEmpiricalPValues, 
                                       simulatedPhyloPForeground$lnlratio)
    fdrForeground <- p.adjust(estimedPValuesForeground, 
                              method = "BH")
    
    estimedPValuesBackground <- sapply(observedPhyloPBackground$lnlratio, 
                                       estimateEmpiricalPValues, 
                                       simulatedPhyloPBackground$lnlratio)
    fdrBackground <- p.adjust(estimedPValuesBackground, 
                              method = "BH")
    
    # Get the conserved elements that are significantly accelerated in the background:
    backgroundAcceleratedElements <- splitElements[fdrBackground< 0.05,]
    if (length(backgroundAcceleratedElements$seqname) == 0) {
      print("no accelerated features.")
    } else {
      print("exporting accelerated features.")
      nrow(backgroundAcceleratedElements)
      backgroundAcceleratedElements$feature <- groupAcceleratedIn
      backgroundAcceleratedElements$score <- observedPhyloPBackground$lnlratio[fdrBackground < 0.05]
      backgroundAcceleratedElements$seqname <- scaffoldName
      write.feat(backgroundAcceleratedElements, 
                 paste("./16_phastConsAnalyses/acceleratedElements/",
                       groupAcceleratedIn,
                       "AcceleratedElements_",
                       scaffoldName,
                       ".gff",
                       sep = ""))
    }
  }
  findingConservedAcceleratedRegions(foregroundSpecies = dimorphicSpecies, 
                                     backgroundSpecies = monomorphicSpecies, 
                                     groupAcceleratedIn = "monomorphicSpecies")
  findingConservedAcceleratedRegions(foregroundSpecies = monomorphicSpecies, 
                                     backgroundSpecies = dimorphicSpecies, 
                                     groupAcceleratedIn = "dimorphicSpecies")
}

possiblyGettingAcceleratedRegions <- possibly(gettingAcceleratedRegions,
                                              otherwise = "Error.")
purrr::map(listOfScaffoldNumbers,
           possiblyGettingAcceleratedRegions)

### Read in accelerated regions to compare foreground vs. background ####
acceleratedRegions <- list.files("./16_phastConsAnalyses/acceleratedElements/",
                                 full.names = TRUE)
acceleratedRegions <- purrr::map(acceleratedRegions,
                                 ape::read.gff)
acceleratedRegions <- as.data.frame(do.call(rbind, acceleratedRegions)) 

dimorphicSpeciesCount <- length(which(acceleratedRegions$type == "dimorphicSpecies"))
monomorphicSpeciesCount <- length(which(acceleratedRegions$type == "monomorphicSpecies"))
