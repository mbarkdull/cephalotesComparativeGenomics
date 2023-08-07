library(plyr)
library(tidyverse)
library(RJSONIO)
library(purrr)
library(topGO)

#### Processing RELAX results and running GO enrichment ####
# Read in RELAX results #
# List all of the RELAX results files with size greater than zero:
relaxResults <- list.files(path = "./08_relaxResults",
                           full.names = TRUE)
relaxResults <- relaxResults[sapply(relaxResults, file.size) > 0]

# Write a function to parse a single result JSON file:
parsingRelax <- function(i) {
  singleRelaxResult <- RJSONIO::fromJSON(content = i)
  singleRelaxResultSimple <- c(singleRelaxResult[["input"]][["file name"]],
                               singleRelaxResult[["test results"]][["relaxation or intensification parameter"]],
                               singleRelaxResult[["test results"]][["p-value"]])
  singleRelaxResultSimple <- as.data.frame(t(singleRelaxResultSimple))
}

# Apply that to all the results files:
possiblyparsingRelax <- purrr::possibly(parsingRelax,
                                        otherwise = "Error")
allRelaxResults <- purrr::map(relaxResults,
                              possiblyparsingRelax)
allRelaxResults <- as.data.frame(do.call(rbind, allRelaxResults))   

# Give the results meaningful column names:
colnames(allRelaxResults) <- c("inputFile",
                               "kValue",
                               "pValue")

# Get a column with the orthogroup ID:
allRelaxResults$orthogroup <- stringr::str_split_i(allRelaxResults$inputFile, 
                                                   pattern = "/",
                                                   i = 9)
allRelaxResults$orthogroup <- stringr::str_split_i(allRelaxResults$orthogroup,
                                                   pattern = "_",
                                                   i = 2) %>%
  stringr::str_split_i(pattern = "\\.",
                       i = 1)

# Make sure the numeric columns are really numeric:
allRelaxResults$pValue <- as.numeric(as.character(allRelaxResults$pValue)) 
allRelaxResults$kValue <- as.numeric(as.character(allRelaxResults$kValue)) 

# Do FDR correction on the p-values:
allRelaxResults$pValueFDR <- p.adjust(allRelaxResults$pValue, method='BH') 

# Export the allRelaxResults to a csv:
write_delim(allRelaxResults,
            file = "./allRelaxResults.csv",
            delim = ",",
            quote = "none")

# Get dataframes of the relaxed and intensified genes
significantlyRelaxedGenes <- filter(allRelaxResults, 
                                    kValue < 1 &
                                      pValueFDR <= 0.05)

significantlyIntensifiedGenes <- filter(allRelaxResults, 
                                        pValueFDR <= 0.05 &
                                          kValue > 1)

# Prep annotations for GO enrichment:
GOannotations <- read_delim("./10_InterProScan/kinfin/kinfin_results/cluster_domain_annotation.GO.txt", 
                            delim = "\t")
# Subset so we just have orthogroup name and GO domain IDs:
longAnnotations <- dplyr::select(GOannotations, 
                                 `#cluster_id`,
                                 domain_id)
# Take out genes without GO terms
longAnnotations <- longAnnotations[which(longAnnotations$domain_id != ""),] 
# Rename the column #cluster_id to orthogroup
longAnnotations <- longAnnotations %>%
  dplyr::rename(orthogroup = `#cluster_id`)

# Create list with element for each gene, containing vectors with all terms for each gene
wideListAnnotations <- tapply(longAnnotations$domain_id, longAnnotations$orthogroup, function(x)x)

# Run GO enrichment for the RELAX results:
runningGOEnrichmentRELAX <- function(geneSet) {
  significanceInfo <- dplyr::select(allRelaxResults, 
                                    orthogroup, 
                                    pValueFDR, 
                                    kValue) 
  
  if (geneSet == "relaxed") {
    print("Analyzing GO enrichment for relaxed orthogroups")
    # Set each gene to 1 if signficantly relaxed, otherwise set to 0
    geneList <- ifelse(significanceInfo$pValueFDR <= 0.05 & 
                         significanceInfo$kValue < 1, 
                       1, 
                       0)
  } else if (geneSet == "intensified") {
    print("Analyzing GO enrichment for intensified orthogroups")
    # Set each gene to 1 if signficantly under positive selection in the foreground, otherwise set to 0
    geneList <- ifelse(significanceInfo$pValueFDR <= 0.05 & 
                         significanceInfo$kValue > 1, 
                       1, 
                       0)
  } else {
    print("Must specify geneSet as either `relaxed` or `intensified`")
  }
  
  # Give geneList names:
  names(geneList) <- significanceInfo$orthogroup
  
  # Create the GOdata object:
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultsFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
  resultsFisherBP
  resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultsFisherBP, topNodes = length(resultsFisherBP@score),
                                   numChar = 120)
  
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherMF <- runTest(GOdataMF, algorithm = "elim", statistic = "fisher")
  resultFisherMF
  resultsFisherMFTable <- GenTable(GOdataMF, raw.p.value = resultFisherMF, topNodes = length(resultFisherMF@score),
                                   numChar = 120)
  
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherCC <- runTest(GOdataCC, algorithm = "elim", statistic = "fisher")
  resultFisherCC
  resultsFisherCCTable <- GenTable(GOdataCC, raw.p.value = resultFisherCC, topNodes = length(resultFisherCC@score),
                                   numChar = 120)
  
  # Combine all of the results for GO terms enriched in the foreground:
  goTermsIntensity <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  
  if (geneSet == "relaxed") {
    goTermsIntensity$selectionCategory <- "Relaxed selection associated with the trait."
    goTermsIntensity
  } else if (geneSet == "intensified") {
    goTermsIntensity$selectionCategory <- "Intensified selection associated with the trait."
    goTermsIntensity
  } else {
    print("Must specify geneSet as either `relaxed` or `intensified`")
  }
}

relaxGoResultsRelaxed <- runningGOEnrichmentRELAX("relaxed")
relaxGoResultsIntensified <- runningGOEnrichmentRELAX("intensified")

relaxGoResults <- rbind(relaxGoResultsRelaxed,
                        relaxGoResultsIntensified)
relaxGoResults

# Export the results to a csv:
write_delim(relaxGoResults,
            file = "./allrelaxGoResults.csv",
            delim = ",",
            quote = "none")

#### Read in BUSTED-PH results and do GO enrichment ####
# List all results file with size greater than zero:
bustedphResults <- list.files(path = "./09_bustedPHResults",
                              pattern = "*.json",
                              full.names = TRUE)
bustedphResults <- bustedphResults[sapply(bustedphResults, file.size) > 0]

# Write a function to parse a single results JSON file:
parsingBustedph <- function(i) {
  singleBustedphResult <- RJSONIO::fromJSON(content = i)
  singleBustedphResultSimple <- c(singleBustedphResult[["input"]][["file name"]],
                                  singleBustedphResult[["test results"]][["p-value"]],
                                  singleBustedphResult[["test results background"]][["p-value"]],
                                  singleBustedphResult[["test results shared distributions"]][["p-value"]])
  singleBustedphResultSimple <- as.data.frame(t(singleBustedphResultSimple))
}

# Apply that to all the results:
possiblyparsingBustedph <- purrr::possibly(parsingBustedph,
                                           otherwise = "Error")
allBustedphResults <- purrr::map(bustedphResults,
                                 possiblyparsingBustedph)
allBustedphResults <- as.data.frame(do.call(rbind, allBustedphResults))   

# Give columns meaningful names:
colnames(allBustedphResults) <- c("inputFile",
                                  "testPvalue",
                                  "backgroundPvalue",
                                  "differencePvalue")

# Get a column with just the orthogroup ID:
allBustedphResults$orthogroup <- stringr::str_split_i(allBustedphResults$inputFile, 
                                                      pattern = "/",
                                                      i = 7) %>%
  stringr::str_split_i(pattern = "_", 
                       i = 1)

# Make sure columns are numeric:
allBustedphResults$testPvalue <- as.numeric(as.character(allBustedphResults$testPvalue)) 
allBustedphResults$backgroundPvalue <- as.numeric(as.character(allBustedphResults$backgroundPvalue)) 
allBustedphResults$differencePvalue <- as.numeric(as.character(allBustedphResults$differencePvalue)) 

# Do FDR corrections on p-values:
allBustedphResults$testPvalueFDR <- p.adjust(allBustedphResults$testPvalue, method='BH') %>%
  as.numeric(as.character()) 
allBustedphResults$backgroundPvalueFDR <- p.adjust(allBustedphResults$backgroundPvalue, method='BH') %>%
  as.numeric(as.character()) 
allBustedphResults$differencePvalueFDR <- p.adjust(allBustedphResults$differencePvalue, method='BH') %>%
  as.numeric(as.character()) 

# Export the results to a csv:
write_delim(allBustedphResults,
            file = "./allBustedPHResults.csv",
            delim = ",",
            quote = "none")

# Get dataframes of orthogroups positively selected in the fore- and background:
foregroundPositiveSelection <- filter(allBustedphResults, 
                                      testPvalueFDR <= 0.05,
                                      backgroundPvalueFDR > 0.05,
                                      differencePvalueFDR <= 0.05)

backgroundPositiveSelection <- filter(allBustedphResults, 
                                      testPvalueFDR > 0.05,
                                      backgroundPvalueFDR <= 0.05,
                                      differencePvalueFDR <= 0.05)

# Prep annotations for GO enrichment:
GOannotations <- read_delim("./10_InterProScan/kinfin/kinfin_results/cluster_domain_annotation.GO.txt", 
                            delim = "\t")
# Subset so we just have orthogroup name and GO domain IDs:
longAnnotations <- dplyr::select(GOannotations, 
                                 `#cluster_id`,
                                 domain_id)
# Take out genes without GO terms
longAnnotations <- longAnnotations[which(longAnnotations$domain_id != ""),] 
# Rename the column #cluster_id to orthogroup
longAnnotations <- longAnnotations %>%
  dplyr::rename(orthogroup = `#cluster_id`)

# Create list with element for each gene, containing vectors with all terms for each gene
wideListAnnotations <- tapply(longAnnotations$domain_id, longAnnotations$orthogroup, function(x)x)

# Run GO enrichment for the BUSTED-PH results:
runningGOEnrichmentBUSTEDPH <- function(geneSet) {
  significanceInfo <- dplyr::select(allBustedphResults, 
                                    orthogroup, 
                                    testPvalueFDR, 
                                    backgroundPvalueFDR,
                                    differencePvalueFDR) 
  
  if (geneSet == "foreground") {
    print("Analyzing GO enrichment for orthogroups positively selected in the foreground")
    # Set each gene to 1 if signficantly under positive selection in the foreground, otherwise set to 0
    geneList <- ifelse(significanceInfo$testPvalueFDR <= 0.05 & 
                         significanceInfo$backgroundPvalueFDR > 0.05 & 
                         significanceInfo$differencePvalueFDR <= 0.05, 
                       1, 
                       0)
  } else if (geneSet == "background") {
    print("Analyzing GO enrichment for orthogroups positively selected in the background")
    # Set each gene to 1 if signficantly under positive selection in the foreground, otherwise set to 0
    geneList <- ifelse(significanceInfo$testPvalueFDR > 0.05 & 
                         significanceInfo$backgroundPvalueFDR <= 0.05 & 
                         significanceInfo$differencePvalueFDR <= 0.05, 
                       1, 
                       0)
  } else {
    print("Must specify geneSet as either `foreground` or `background`")
  }
  
  # Give geneList names:
  names(geneList) <- significanceInfo$orthogroup
  
  # Create the GOdata object:
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultsFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
  resultsFisherBP
  resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultsFisherBP, topNodes = length(resultsFisherBP@score),
                                   numChar = 120)
  
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherMF <- runTest(GOdataMF, algorithm = "elim", statistic = "fisher")
  resultFisherMF
  resultsFisherMFTable <- GenTable(GOdataMF, raw.p.value = resultFisherMF, topNodes = length(resultFisherMF@score),
                                   numChar = 120)
  
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = annFUN.gene2GO, gene2GO = wideListAnnotations)
  # Run Fisher's exact test to check for enrichment:
  resultFisherCC <- runTest(GOdataCC, algorithm = "elim", statistic = "fisher")
  resultFisherCC
  resultsFisherCCTable <- GenTable(GOdataCC, raw.p.value = resultFisherCC, topNodes = length(resultFisherCC@score),
                                   numChar = 120)
  
  # Combine all of the results for GO terms enriched in the foreground:
  goTermsPositive <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
    dplyr::filter(raw.p.value <= 0.01)
  
  if (geneSet == "foreground") {
    goTermsPositive$selectionCategory <- "Positive selection associated with the trait."
    goTermsPositive
  } else if (geneSet == "background") {
    goTermsPositive$selectionCategory <- "Positive selection associated with absence of the trait."
    goTermsPositive
  } else {
    print("Must specify geneSet as either `foreground` or `background`")
  }
}

bustedGoResultsForeground <- runningGOEnrichmentBUSTEDPH("foreground")
bustedGoResultsBackground <- runningGOEnrichmentBUSTEDPH("background")
bustedGoResults <- rbind(bustedGoResultsForeground,
                         bustedGoResultsBackground)

# Export the results to a csv:
write_delim(bustedGoResults,
            file = "./allBustedPHGOResults.csv",
            delim = ",",
            quote = "none")
