library(tidyverse)
library(rphast)
library(splitstackshape)

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
  if (length(features$seqname) == 0) {
    print("no annotations on scaffold")
  } else {
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
      print(paste("Running phastcons on these",
                  cat(foregroundAlignment$names),
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
      noncodingInformativeElementsExport <- noncodingInformativeElements
      print(paste("Exporting features conserved in",
                  cat(foregroundAlignment$names),
                  sep = ""))
      nrow(noncodingInformativeElementsExport)
      noncodingInformativeElementsExport$feature <- paste("conservedExcluding",
                                                          groupAcceleratedIn,
                                                          sep = "")
      noncodingInformativeElementsExport$seqname <- scaffoldName
      write.feat(noncodingInformativeElementsExport, 
                 paste("./16_phastConsAnalyses/acceleratedElements/",
                       "conservedElementsExcluding",
                       groupAcceleratedIn,
                       "_",
                       scaffoldName,
                       ".gff",
                       sep = ""))
      
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
      print("Seeing how 50bp conserved elements are evolving in foreground species:")
      
      observedPhyloPForeground <- phyloP(allSpeciesNeutralModel, 
                                         msa = allSpeciesAlignment, 
                                         mode = "ACC", 
                                         features = splitElements, 
                                         branches = foregroundSpecies)
      # See how 50bp conserved elements are evolving in background species:
      print(paste("Seeing how the 50bp conserved elements are evolving in ",
                  groupAcceleratedIn,
                  sep = ""))
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
}

possiblyGettingAcceleratedRegions <- possibly(gettingAcceleratedRegions,
                                              otherwise = "Error.")
purrr::map(listOfScaffoldNumbers,
           possiblyGettingAcceleratedRegions)

### Read in accelerated regions to compare foreground vs. background ####
acceleratedRegions <- list.files("./16_phastConsAnalyses/acceleratedElements",
                                 full.names = TRUE,
                                 pattern = "*AcceleratedElements*")
acceleratedRegions <- purrr::map(acceleratedRegions,
                                 ape::read.gff)
acceleratedRegions <- as.data.frame(do.call(rbind, acceleratedRegions)) 

dimorphicSpeciesAccelerated <- length(which(acceleratedRegions$type == "dimorphicSpecies"))
monomorphicSpeciesAccelerated <- length(which(acceleratedRegions$type == "monomorphicSpecies"))

referenceRegions <- list.files("./16_phastConsAnalyses/acceleratedElements",
                                 full.names = TRUE,
                                 pattern = "*conservedElements*")
referenceRegions <- purrr::map(referenceRegions,
                                 ape::read.gff)
referenceRegions <- as.data.frame(do.call(rbind, referenceRegions)) 
dimorphicSpeciesConserved <- length(which(referenceRegions$type == "conservedExcludingmonomorphicSpecies"))
monomorphicSpeciesConserved <- length(which(referenceRegions$type == "conservedExcludingdimorphicSpecies"))

percentElementsAcceleratedInMonomorphic <- monomorphicSpeciesAccelerated/dimorphicSpeciesConserved
percentElementsAcceleratedInDimorphic <- dimorphicSpeciesAccelerated/monomorphicSpeciesConserved

allResults <- data.frame(category = c("Dimorphic species",
                                      "Monomorphic species",
                                      "Dimorphic species",
                                      "Monomorphic species"),
                         analysis = c("Relaxed selection",
                                      "Relaxed selection",
                                      "Conserved",
                                      "Conserved"),
                         count = c(dimorphicSpeciesAccelerated,
                                   monomorphicSpeciesAccelerated,
                                   dimorphicSpeciesConserved,
                                   monomorphicSpeciesConserved))

colors <- c("Relaxed selection" = "#f1a208",
            "Conserved" = "#dcdcdc")

ggplot(data = allResults) +
  geom_col(mapping = aes(x = category,
                         y = count,
                         fill = analysis))  +
  scale_fill_manual(values = colors) +
  theme_bw() + 
  labs(x = "",
       y = "Number of regulatory elements",
       fill = "Selective regime") +
  scale_y_continuous(expand = c(0, 0))

ggsave(filename = "./Plots/allNoncodingElements.png",
       width = 6,
       height = 6,
       units = "in",
       dpi = 600)


#### Where are the regulatory elements in the genome? ####
# Read in the whole C. varians annotation:
cvarAnnotation <- ape::read.gff("./CVAR/CVAR_OGS_v1.0_longest_isoform.gff3")

# Make sure that the scaffold names for the noncoding elements will match:
acceleratedRegions$seqid <- gsub(pattern = "^CVAR_",
                                 replacement = "",
                                 acceleratedRegions$seqid)

referenceRegions$seqid <- gsub(pattern = "^CVAR_",
                               replacement = "",
                               referenceRegions$seqid)

# Add attribute info describing the nature of the element:
acceleratedRegions <- acceleratedRegions %>%
  mutate(type = case_when(type == "monomorphicSpecies" ~ "monomorphicAccelerated",
                                type == "dimorphicSpecies" ~ "dimorphicAccelerated"))

referenceRegions <- referenceRegions %>%
  mutate(type = case_when(type == "conservedExcludingdimorphicSpecies" ~ "monomorphicConserved",
                                type == "conservedExcludingmonomorphicSpecies" ~ "dimorphicConserved"))

# Combine the accelerated regions with the whole annotation:
allFeatures <- rbind(acceleratedRegions,
                     referenceRegions,
                     cvarAnnotation)

# Group by scaffold
allFeatures <- allFeatures %>%
  group_by(seqid) 

# Arrange by start and then scaffold:
allFeatures <- allFeatures %>%
  arrange(start) %>%
  arrange(seqid)

# To figure out which noncoding elements are inside genes, first, add intergenic features:
# Read in the annotation as a GRanges object:
annotation <- rtracklayer::import.gff("./CVAR/CVAR_OGS_v1.0_longest_isoform.gff3",
                                      feature.type = "gene")
# Get rid of strand info:
BiocGenerics::strand(annotation) <- '*'

# Get the intergene regions (e.g. gaps between genes):
intergenes <- GenomicRanges::gaps(annotation)

#Turn it into a dataframe and clean it up:
intergenes <- Repitools::annoGR2DF(intergenes)
colnames(intergenes) <- c("seqid",
                          "start",
                          "end",
                          "width")
intergenes <- intergenes %>%
  dplyr::select(-c(width))
intergenes$source <- "intergene"
intergenes$type <- "intergene"
intergenes$score <- NA
intergenes$phase <- NA
intergenes$attributes <- "intergene"
intergenes$strand <- NA
  
# Add the intergenes into the overall annotation with all the feature types:
allFeatures <- rbind(allFeatures,
                     intergenes) %>%
  arrange(start) %>%
  arrange(seqid)

# Test which noncoding elements overlap with intergenes:
library(valr)
intergenesToOverlap <- intergenes %>%
  dplyr::select(seqid,
         start,
         end)
colnames(intergenesToOverlap) <- c("chrom", 
                                           "start",
                                           "end")

allNoncodingRegions <- allFeatures %>%
  filter(type == "monomorphicConserved" |
           type == "dimorphicConserved" |
           type == "monomorphicAccelerated" |
           type == "dimorphicAccelerated")

noncodingRegionsToOverlap <- allNoncodingRegions %>%
  dplyr::select(seqid,
         start,
         end)

colnames(noncodingRegionsToOverlap) <- c("chrom", 
                                           "start",
                                           "end")


overlapInfo <- valr::bed_intersect(intergenesToOverlap, 
                                   noncodingRegionsToOverlap, 
                                   suffix = c("_intergene", 
                                              "_AR")) %>%
  distinct()
overlapsWithIntergenic <- right_join(overlapInfo,
                                     allNoncodingRegions,
                                     by = c("start_AR" = "start",
                                            "end_AR" = "end",
                                            "chrom" = "seqid")) %>%
  dplyr::select(-c(attributes)) %>%
  distinct()

overlapsWithIntergenic$.overlap <- replace_na(overlapsWithIntergenic$.overlap,
                                              replace = 0)
overlapsWithIntergenic$lengthNoncoding <- overlapsWithIntergenic$end_AR - overlapsWithIntergenic$start_AR
overlapsWithIntergenic$percentOverlap <- overlapsWithIntergenic$.overlap / overlapsWithIntergenic$lengthNoncoding
overlapsWithIntergenic <- overlapsWithIntergenic %>%
  mutate(isIntergenic = case_when(percentOverlap == 0 ~ "no",
                                  .overlap == 0 ~ "no",
                                  percentOverlap > 0 & percentOverlap < 1 ~ "partial",
                                  percentOverlap == 1 ~ "yes"))

ggplot(data = overlapsWithIntergenic) +
  geom_bar(mapping = aes(x = type,
                         fill = isIntergenic)) 

# If an element is within a gene, which genes?
insideGene <- overlapsWithIntergenic %>%
  filter(isIntergenic != "yes" &
           lengthNoncoding > 0)


genesToOverlap <- allFeatures %>%
  filter(type == "gene") %>%
  dplyr::select(seqid,
         start,
         end)
colnames(genesToOverlap) <- c("chrom", 
                                   "start",
                                   "end")

insideGeneToOverlap <- insideGene %>%
  dplyr::select(chrom,
         start_AR,
         end_AR)
colnames(insideGeneToOverlap) <- c("chrom", 
                              "start",
                              "end")

overlapWithGenes <- valr::bed_intersect(genesToOverlap, 
                                   insideGeneToOverlap, 
                                   suffix = c("_gene", 
                                              "_noncoding")) %>%
  distinct()

overlapWithGenes <- right_join(overlapWithGenes,
                               filter(allFeatures,
                                      type == "gene"),
                                     by = c("start_gene" = "start",
                                            "end_gene" = "end",
                                            "chrom" = "seqid")) %>%
  distinct() 

insideGeneInfo <- right_join(overlapWithGenes,
             insideGene,
             by = c("start_noncoding" = "start_AR",
                    "end_noncoding" = "end_AR",
                    "chrom" = "chrom")) %>%
  dplyr::select(c("chrom",
           "start_gene",
           "end_gene",
           "start_noncoding",
           "end_noncoding",
           "attributes",
           "type.y",
           "lengthNoncoding",
           "percentOverlap",
           "isIntergenic" )) %>%
  distinct()

insideGeneInfo$geneName <- str_split_i(insideGeneInfo$attributes,
                             pattern = "=",
                             i = 3)
insideGeneInfo$geneName <- gsub(pattern = ";",
                      replacement = "",
                      insideGeneInfo$geneName)







# If an element is outside any genes, which genes is it near?
outsideGene <- overlapsWithIntergenic %>%
  filter(isIntergenic == "yes" &
           lengthNoncoding > 0) %>%
  dplyr::select(c(chrom,
           start_AR,
           end_AR,
           type,
           isIntergenic))

colnames(outsideGene) <- c("chrom",
                           "start",
                           "end",
                           "type",
                           "isIntergenic")

# I can use the function nearest from GRanges to find the closest genes; 
# but to do that I need to convert the dataframe to a GRanges object
outsideGeneGRanges <- GenomicRanges::makeGRangesFromDataFrame(df = outsideGene,
                                        keep.extra.columns = TRUE,
                                        seqnames.field = "chrom",
                                        start.field = "start",
                                        end.field = "end")

nearestGenes <- GenomicRanges::nearest(outsideGeneGRanges, 
                                       subject = annotation, 
                                       select = "all", 
                                       ignore.strand=FALSE) %>%
  as.data.frame()
colnames(nearestGenes) <- c("indexOfNoncoding",
                            "indexOfGene")

infoAboutNoncoding <- cbind(as.data.frame(outsideGeneGRanges@seqnames),
                            as.data.frame(outsideGeneGRanges@elementMetadata),
                            as.data.frame(outsideGeneGRanges@ranges))
infoAboutNoncoding$noncodingIndex <- as.numeric(row.names(infoAboutNoncoding))

infoAboutGenes <- cbind(as.data.frame(annotation@seqnames),
                        as.data.frame(annotation@elementMetadata),
                        as.data.frame(annotation@ranges))
infoAboutGenes$geneIndex <- as.numeric(row.names(infoAboutGenes))

# Connect the pieces:
combinedInfo <-right_join(infoAboutNoncoding,
                          nearestGenes, 
                          by = c("noncodingIndex" = "indexOfNoncoding"))
combinedInfo <-left_join(combinedInfo,
                          infoAboutGenes, 
                          by = c("indexOfGene" = "geneIndex")) %>%
  dplyr::select(c(value.x,
           type.x,
           start.x,
           end.x,
           Name))

genesNearAcceleratedElements <- combinedInfo %>%
  filter(grepl("Accelerated", type.x))

#### Are accelerated elements disproportionately near differentially expressed genes? ####
# Read in differential expression results:
differentialExpression <- read_csv("../cephalotesDifferentialExpression/finalResults/allDifferentialExpressionResultsAndFunctions.csv") %>%
  dplyr::select(-c(seq.name)) %>%
  distinct() %>%
  mutate(DE = case_when(padj <= 0.05 ~ "differentiallyExpressed",
                        TRUE ~ "No"))

# Join then to the accelerated element info:
acceleratedElementExpression <- full_join(genesNearAcceleratedElements,
                                           differentialExpression,
                                           by = c("Name" = "gene_name")) 
acceleratedElementExpression$type.x <- acceleratedElementExpression$type.x %>% 
  replace_na("No accelerated element")

# Create a column indicating which morph/stage a gene is upregulated in:
acceleratedElementExpression$left <- str_split_i(acceleratedElementExpression$contrast,
                                                 pattern = "vs.",
                                                 i = 1)
acceleratedElementExpression$right <- str_split_i(acceleratedElementExpression$contrast,
                                                 pattern = "vs.",
                                                 i = 2)
acceleratedElementExpression <- acceleratedElementExpression %>%
  mutate(upregulated = case_when(log2FoldChange < 0 ~ paste("Upregulated in", 
                                                            left, 
                                                            "relative to", 
                                                            right),
                                 log2FoldChange > 0 ~ paste("Upregulated in", 
                                                            right, 
                                                            "relative to", 
                                                            left)))

# Do some exploratory data visualization:
# By accelerated element type:
ggplot(data = filter(acceleratedElementExpression,
                     !is.na(contrast))) +
  geom_bar(mapping = aes(x = type.x,
                         fill = DE),
           position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1)) +
  facet_wrap(~upregulated)

# By differential expression pattern:
ggplot(data = filter(acceleratedElementExpression,
                     padj <= 0.05)) +
  geom_bar(mapping = aes(x = upregulated,
                         fill = type.x),
           position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1)) +
  xlab(label = "Differentially expressed genes across contrasts") +
  ylab(label = "Accelerated noncoding elements") 

# Statistical tests for under/overrepresentation:
chiSquaredAcrossContrasts <- function(specificContrast) {
  chiTestData <- filter(acceleratedElementExpression, 
                        padj <= 0.05 &
                          contrast == specificContrast)
  chiSquareTable <- table(chiTestData$type.x, 
                          chiTestData$contrast)
  test <- chisq.test(chiSquareTable)
  result <- c(test$p.value, test$residuals, specificContrast)
  return(result)
}

chiSquaredResults <- purrr::map(unique(na.omit(acceleratedElementExpression$contrast)),
                                chiSquaredAcrossContrasts)

# Generate a mosaic plot to see which groups are under-represented:
ggstatsplot::ggbarstats(data = filter(acceleratedElementExpression,
                         padj <= 0.05),
           x = type.x,
           y = contrast) +
  labs(caption = NULL)

vcd::mosaic(~ upregulated + type.x,
            direction = c("v"),
            data = filter(acceleratedElementExpression,
                          padj <= 0.05),
            shade = TRUE,
            rot_labels =c (20,
                           90,
                           0,
                           70))

#### What GO terms are enriched among genes near accelerated elements? ####
# Read in the annotations from eggnog:
eggnogAnnotations <- read_delim("./keggAnnotations.emapper.annotations",
                                delim = "\t",
                                skip = 4) %>%
  filter(!is.na(seed_ortholog))

# Subset so we just have gene name and GO domain IDs:
longAnnotations <- dplyr::select(eggnogAnnotations,
                                 c("#query", "GOs"))

# Reshape into a long dataframe:
longAnnotations <- cSplit(longAnnotations, 
                          splitCols = "GOs",
                          sep = ",",
                          direction = "long")

# Take out genes without GO terms and get distinct rows:
longAnnotations <- filter(longAnnotations,
                          GOs != "-") %>%
  distinct()

# Rename the column #query to geneName
longAnnotations <- longAnnotations %>%
  dplyr::rename(geneName = `#query`)

# Convert transcripts to genes:
longAnnotations$geneName <- gsub(pattern = "-R.*", 
                                 replacement = "", 
                                 as.character(longAnnotations$geneName))

# Create list with element for each gene, containing vectors with all terms for each gene
wideListAnnotations <- tapply(longAnnotations$GOs, longAnnotations$geneName, function(x)x)


runGOenrichment <- function(elementType) {
  # Select out genes associated with accelerated elements:
  significanceInfo <- unique(longAnnotations$geneName)
  
  focalGenes <- genesNearAcceleratedElements %>%
    filter(type.x == elementType) %>%
    dplyr::select(Name) %>%
    distinct()
  
  # Create a formatted list of significance information:
  geneList <- ifelse(significanceInfo %in% focalGenes$Name, 
                     1, 
                     0)
  # Give geneList names:
  names(geneList) <- significanceInfo
  
  # Analyze Biological Process terms:
  library(topGO)
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, 
                  gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultsFisherBP <- runTest(GOdataBP, 
                             algorithm = "elim", 
                             statistic = "fisher")
  resultsFisherBP
  resultsFisherBPTable <- GenTable(GOdataBP, 
                                   raw.p.value = resultsFisherBP, 
                                   topNodes = length(resultsFisherBP@score),
                                   numChar = 120)

  # Analyze Molecular Function terms:
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, 
                  gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherMF <- runTest(GOdataMF, 
                            algorithm = "elim", 
                            statistic = "fisher")
  resultFisherMF
  resultsFisherMFTable <- GenTable(GOdataMF, 
                                   raw.p.value = resultFisherMF, 
                                   topNodes = length(resultFisherMF@score),
                                   numChar = 120)

  # Analyze Cellular Component terms:
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, 
                  gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherCC <- runTest(GOdataCC, 
                            algorithm = "elim", 
                            statistic = "fisher")
  resultFisherCC
  resultsFisherCCTable <- GenTable(GOdataCC, 
                                   raw.p.value = resultFisherCC, 
                                   topNodes = length(resultFisherCC@score),
                                   numChar = 120)

  # Combine all of the results:
  enrichedGOTerms <- rbind(resultsFisherBPTable, 
                           resultsFisherMFTable, 
                           resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  return(enrichedGOTerms)
}


monomorphicGO <- runGOenrichment("monomorphicAccelerated")
dimorphicGO <- runGOenrichment("dimorphicAccelerated")

sharedTerms <- intersect(dimorphicGO$Term,
                         monomorphicGO$Term)













geneFunctions <- read_csv(file = "../cephalotesDifferentialExpression/finalResults/allDifferentialExpressionResultsAndFunctions.csv") %>%
  select(c(gene_name,
           seq.name))

genesNearAcceleratedElementsWithFunctions <- left_join(genesNearAcceleratedElements,
                                          geneFunctions,
                                          by = c("Name" = "gene_name")) %>%
  distinct()

#### Are the accelerated elements near genes under selection? ####
# Read in HYPHY results:
allBustedphResults <- read_delim(file = "./allBustedPHResults.csv",
                                 delim = ",") %>%
  dplyr::select(-c("inputFile"))
allRelaxResults <- read_delim(file = "./allRelaxResults.csv",
                              delim = ",") %>%
  dplyr::select(-c("inputFile", 
                   "HOGmembers"))
hyphyResults <- full_join(allBustedphResults,
                          allRelaxResults)

# Split the HOG members column and filter to get a column with CVAR transcript IDs for each result:
hyphyResults <- splitstackshape::cSplit(hyphyResults, 
                       splitCols = "HOGmembers",
                       sep = "|",
                       direction = "long") %>% 
  filter(grepl("^CVAR", 
               HOGmembers)) %>%
  distinct()
# Turn all hyphens to underscores:
hyphyResults$HOGmembers <- gsub(pattern = "\\-",
                                replacement = "_",
                                hyphyResults$HOGmembers)
hyphyResults$HOGmembers <- gsub(pattern = "_R(.*)",
                                replacement = "",
                                hyphyResults$HOGmembers)


acceleratedElementsSelection <- full_join(genesNearAcceleratedElements,
                                          hyphyResults,
                                          by = c("Name" = "HOGmembers")) %>%
  distinct()


acceleratedElementsSelection <- acceleratedElementsSelection %>%
  mutate(positiveSelection = case_when(testPvalueFDR <= 0.05 &
                                         backgroundPvalueFDR >= 0.05 &
                                         differencePvalueFDR <= 0.05 ~ "positiveSelectionOnlyInDimorphic",
                                       testPvalueFDR <= 0.05 &
                                         backgroundPvalueFDR >= 0.05 &
                                         differencePvalueFDR >= 0.05 ~ "nonsignificantDimorphic",
                                       testPvalueFDR >= 0.05 &
                                         backgroundPvalueFDR <= 0.05 &
                                         differencePvalueFDR <= 0.05 ~ "positiveSelectionOnlyInMonomorphic",
                                       testPvalueFDR >= 0.05 &
                                         backgroundPvalueFDR <= 0.05 &
                                         differencePvalueFDR >= 0.05 ~ "nonsignificantMonomorphic",
                                       testPvalueFDR <= 0.05 &
                                         backgroundPvalueFDR <= 0.05 &
                                         differencePvalueFDR <= 0.05 ~ "bothWithDifferentRegimes",
                                       testPvalueFDR <= 0.05 &
                                         backgroundPvalueFDR <= 0.05 &
                                         differencePvalueFDR >= 0.05 ~ "broadPositiveSelection",
                                       testPvalueFDR >= 0.05 &
                                         backgroundPvalueFDR >= 0.05 ~ "noPositiveSelection",
                                       TRUE ~ "Other")) %>%
  
  
  mutate(intensity = case_when(pValueFDR <= 0.05 &
                                 kValue < 1 ~ "Relaxed",
                               pValueFDR <= 0.05 &
                                 kValue > 1 ~ "Intensified",
                               TRUE ~ "No shift"))

ggplot(data = acceleratedElementsSelection) +
  geom_bar(mapping = aes(x = type.x,
                         fill = intensity),
           position = "fill")

ggplot(data = filter(acceleratedElementsSelection,
                     positiveSelection != "Other")) +
  geom_bar(mapping = aes(x = type.x,
                         fill = positiveSelection),
           position = "fill")


