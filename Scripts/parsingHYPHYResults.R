library(plyr)
library(tidyverse)
library(RJSONIO)
library(purrr)
library(topGO)

#### RELAX Results ####
relaxResults <- list.files(path = "./08_relaxResults",
                           full.names = TRUE)
relaxResults <- relaxResults[sapply(relaxResults, file.size) > 0]

parsingRelax <- function(i) {
  singleRelaxResult <- RJSONIO::fromJSON(content = i)
  singleRelaxResultSimple <- c(singleRelaxResult[["input"]][["file name"]],
                               singleRelaxResult[["test results"]][["relaxation or intensification parameter"]],
                               singleRelaxResult[["test results"]][["p-value"]])
  singleRelaxResultSimple <- as.data.frame(t(singleRelaxResultSimple))
}

possiblyparsingRelax <- purrr::possibly(parsingRelax,
                                        otherwise = "Error")
allRelaxResults <- map(relaxResults,
                       possiblyparsingRelax)
allRelaxResults <- as.data.frame(do.call(rbind, allRelaxResults))   
colnames(allRelaxResults) <- c("inputFile",
                               "kValue",
                               "pValue")
allRelaxResults$orthogroup <- str_split_i(allRelaxResults$inputFile, 
                                          pattern = "/",
                                          i = 9)
allRelaxResults$orthogroup <- str_split_i(allRelaxResults$orthogroup,
                                          pattern = "_",
                                          i = 2) %>%
  str_split_i(pattern = "\\.",
              i = 1)
allRelaxResults$pValue <- as.numeric(as.character(allRelaxResults$pValue)) 
allRelaxResults$kValue <- as.numeric(as.character(allRelaxResults$kValue)) 
allRelaxResults$pValueFDR <- p.adjust(allRelaxResults$pValue, method='BH') 

# Read in the GO term annotations for each orthogroup, from KinFin:
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

# Make a wide version:
annotations <- tapply(longAnnotations$domain_id, longAnnotations$orthogroup, function(x)x)
annotations <- plyr::ldply(annotations, rbind)

# Combine GO annotations with RELAX results:
allRelaxResults <- left_join(allRelaxResults,
                             annotations, 
                             by = c("orthogroup" = ".id"))

# Get dataframes of the relaxed and intensified genes
significantlyRelaxedGenes <- filter(allRelaxResults, 
                                    kValue < 1 &
                                      pValueFDR <= 0.05)

significantlyIntensifiedGenes <- filter(allRelaxResults, 
                                        pValueFDR <= 0.05 &
                                          kValue > 1)

#### BUSTED-PH Results ####
bustedphResults <- list.files(path = "./09_bustedPHResults",
                              pattern = "*.json",
                              full.names = TRUE)
bustedphResults <- bustedphResults[sapply(bustedphResults, file.size) > 0]

parsingBustedph <- function(i) {
  singleBustedphResult <- RJSONIO::fromJSON(content = i)
  singleBustedphResultSimple <- c(singleBustedphResult[["input"]][["file name"]],
                                  singleBustedphResult[["test results"]][["p-value"]],
                                  singleBustedphResult[["test results background"]][["p-value"]],
                                  singleBustedphResult[["test results shared distributions"]][["p-value"]])
  singleBustedphResultSimple <- as.data.frame(t(singleBustedphResultSimple))
}

possiblyparsingBustedph <- purrr::possibly(parsingBustedph,
                                           otherwise = "Error")
allBustedphResults <- map(bustedphResults,
                          possiblyparsingBustedph)
allBustedphResults <- as.data.frame(do.call(rbind, allBustedphResults))   
colnames(allBustedphResults) <- c("inputFile",
                                  "testPvalue",
                                  "backgroundPvalue",
                                  "differencePvalue")
allBustedphResults$orthogroup <- str_split_i(allBustedphResults$inputFile, 
                                             pattern = "/",
                                             i = 7)
allBustedphResults$testPvalue <- as.numeric(as.character(allBustedphResults$testPvalue)) 
allBustedphResults$backgroundPvalue <- as.numeric(as.character(allBustedphResults$backgroundPvalue)) 
allBustedphResults$differencePvalue <- as.numeric(as.character(allBustedphResults$differencePvalue)) 

allBustedphResults$testPvalueFDR <- p.adjust(allBustedphResults$testPvalue, method='BH') %>%
  as.numeric(as.character()) 
allBustedphResults$backgroundPvalueFDR <- p.adjust(allBustedphResults$backgroundPvalue, method='BH') %>%
  as.numeric(as.character()) 
allBustedphResults$differencePvalueFDR <- p.adjust(allBustedphResults$differencePvalue, method='BH') %>%
  as.numeric(as.character()) 

foregroundPositiveSelection <- filter(allBustedphResults, 
                                      testPvalueFDR <= 0.05,
                                      backgroundPvalueFDR > 0.05,
                                      differencePvalueFDR <= 0.05)

backgroundPositiveSelection <- filter(allBustedphResults, 
                                      testPvalueFDR > 0.05,
                                      backgroundPvalueFDR <= 0.05,
                                      differencePvalueFDR <= 0.05)

### Adding KinFin annotations ###
#### Do GO enrichment: #####
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

# Run GO enrichment:
significanceInfo <- dplyr::select(allRelaxResults, 
                                  orthogroup, 
                                  pValueFDR, 
                                  kValue) 
significanceInfo$orthogroup <- str_split_i(significanceInfo$orthogroup,
                                           pattern = "_",
                                           i = 2) %>%
  str_split_i(pattern = "\\.",
              i = 1)
# Set each gene to 1 if adjP < cutoff amd kValue is < 1, otherwise set to 0
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$pValueFDR < pcutoff & significanceInfo$kValue < 1, 1, 0)
geneList <- tmp
rm(tmp)

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
goTermsRelaxed <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
  dplyr::filter(raw.p.value <= 0.01)
goTermsRelaxed$selectionCategory <- "Relaxed selection associated with the trait"

goTermsRelaxed








# Set each gene to 1 if adjP < cutoff amd kValue is < 1, otherwise set to 0
pcutoff <- 0.05 
tmp <- ifelse(significanceInfo$pValueFDR < pcutoff & significanceInfo$kValue > 1, 1, 0)
geneList <- tmp
rm(tmp)

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
goTermsIntensified <- rbind(resultsFisherBPTable, resultsFisherMFTable, resultsFisherCCTable) %>%
  dplyr::filter(raw.p.value <= 0.01)
goTermsIntensified$selectionCategory <- "Intensified selection associated with the trait"

goTermsRelaxed
goTermsIntensified
