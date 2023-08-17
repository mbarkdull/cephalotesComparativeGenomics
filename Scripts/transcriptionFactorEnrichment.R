library(rtracklayer)
library(Gviz)
library(ape)
library(tidyverse)

dir.create("./14_SEAenrichment/")

#### Read in all annotations and use them to create promoter regions: ####
annotation <- read.gff("./CVAR/CVAR_OGS_v1.0_longest_isoform.gff3", 
                       na.strings = c(".",
                                      "?"),
                       GFF3 = TRUE)
# Get an ID column:
annotation$ID <- stringr::str_extract(annotation$attributes,
                                      regex("ID[^;]+", 
                                            ignore_case = T)) %>%
  stringr::str_split_i(pattern = "=",
                       i = 2) %>%
  stringr::str_split_i(pattern = ":",
                       i = 1)
annotation$ID <- gsub(pattern = "_",
                      replacement = "-",
                      annotation$ID)
# Filter out gene annotations:
annotation <- filter(annotation, 
                     type != "gene")

# For every gene, select out the region 2kb upstream as a promoter:
# "We proceeded to examine whether ant genes that harbor conserved TFBSs in their promoters (0–2 kb upstream of ORFs) exhibit evolutionary changes in TF regulation among insects", Simola et al. 2013
identifyingPromoterRegions <- function(gene) {
  relevantAnnotations <- filter(annotation, 
                                ID == gene &
                                  type == "CDS")
  
  if (as.character(unique(relevantAnnotations$strand)) == "+") {
    # Find where the gene starts and calculate 2kb upstream:
    geneStart <- as.numeric(min(relevantAnnotations$start))
    promoter2kbUpstreamOfORF <- geneStart-2000
    # If that start site is not negative, add it as an annotation:
    if (promoter2kbUpstreamOfORF >= 1) {
      # Add this information as an annotation for the promoter:
      relevantAnnotationsWithPromoter <- relevantAnnotations %>% 
        add_row(seqid = as.character(unique(relevantAnnotations$seqid)),
                source = "rCode",
                type = "promoter",
                start = as.numeric(promoter2kbUpstreamOfORF),
                end = as.numeric(geneStart),
                score = NA,
                strand = "+",
                phase = "0",
                attributes = as.character(unique(relevantAnnotations$attributes)),
                ID = as.character(unique(relevantAnnotations$ID))) %>%
        distinct()
    }
    
  } else {
    # Find where the gene starts and calculate 2kb upstream; note this is all backwards because we're on the - strand:
    geneStart <- as.numeric(max(relevantAnnotations$end))
    promoter2kbUpstreamOfORF <- geneStart+2000
    # Add this information as an annotation for the promoter:
    relevantAnnotationsWithPromoter <- relevantAnnotations %>% 
      add_row(seqid = as.character(unique(relevantAnnotations$seqid)),
              source = "rCode",
              type = "promoter",
              start = as.numeric(geneStart),
              end = as.numeric(promoter2kbUpstreamOfORF),
              score = NA,
              strand = "-",
              phase = "0",
              attributes = as.character(unique(relevantAnnotations$attributes)),
              ID = as.character(unique(relevantAnnotations$ID))) %>%
      distinct()
  }
  relevantAnnotationsWithPromoter$start <- as.numeric(relevantAnnotationsWithPromoter$start)
  relevantAnnotationsWithPromoter <- relevantAnnotationsWithPromoter %>% 
    mutate(start = if_else(start < 0, 
                           0, 
                           start))
  readr::write_delim(subset(relevantAnnotationsWithPromoter,
                            select = -c(ID)),
                     file = "./14_SEAenrichment/genesWithPromoters.gff",
                     delim = "\t",
                     quote = "none",
                     col_names = FALSE, 
                     na = ".",
                     append = TRUE)
  
  
  justPromoters <- filter(relevantAnnotationsWithPromoter, 
                          type == "promoter")
  readr::write_delim(subset(justPromoters,
                            select = -c(ID)),
                     file = "./14_SEAenrichment/justPromoters.gff",
                     delim = "\t",
                     quote = "none",
                     col_names = FALSE, 
                     na = ".",
                     append = TRUE)
}

# Make a safe version and apply to all genes:
possiblyIdentifyingPromoterRegions <- purrr::possibly(identifyingPromoterRegions, 
                                                      otherwise = "Error")
genes <- unique(annotation$ID)
purrr::map(genes, 
           possiblyIdentifyingPromoterRegions)

# Read in the results for further use:
genesWithPromoters <- read_delim(file = "./14_SEAenrichment/genesWithPromoters.gff",
                                 delim = "\t",
                                 col_names = c("seqid",
                                               "source",
                                               "type",
                                               "start",
                                               "end",
                                               "score",
                                               "strand",
                                               "phase",
                                               "attributes"))
genesWithPromoters$ID <- stringr::str_extract(genesWithPromoters$attributes,
                                              regex("ID[^;]+", 
                                                    ignore_case = T)) %>%
  stringr::str_split_i(pattern = "=",
                       i = 2) %>%
  stringr::str_split_i(pattern = ":",
                       i = 1)
genesWithPromoters$ID <- gsub(pattern = "_",
                              replacement = "-",
                              genesWithPromoters$ID)

justPromoters <- filter(genesWithPromoters, 
                        type == "promoter")

#### Run this over genes under relaxed, intensified, and positive selection in both directions: ####
shiftedSelectionIntensity <- readr::read_delim("./allRelaxResults.csv") %>%
  filter(pValueFDR <= 0.05) %>%
  mutate(selectionCategory = case_when(kValue < 1 ~ "relaxedSelection",
                                       kValue > 1 ~ "intensifiedSelection",
                                       TRUE ~ "noSelection")) %>%
  dplyr::select(c("orthogroup",
                  "HOGmembers",
                  "selectionCategory"))


positiveSelection <- readr::read_delim("./allBustedPHResults.csv") %>%
  mutate(selectionCategory = case_when(testPvalueFDR <= 0.05 &
                                         backgroundPvalueFDR > 0.05 &
                                         differencePvalueFDR <= 0.05 ~ "foregroundOnly",
                                       testPvalueFDR > 0.05 &
                                         backgroundPvalueFDR <= 0.05 &
                                         differencePvalueFDR <= 0.05 ~ "backgroundOnly",
                                       TRUE ~ "noSelection")) %>%
  dplyr::select(c("orthogroup",
                  "HOGmembers",
                  "selectionCategory"))

genesUnderSelection <- rbind(shiftedSelectionIntensity,
                             positiveSelection) %>%
  filter(selectionCategory != "noSelection")

# Find the gene name for the CVAR member, if possible 
genesUnderSelection$CVARgene <- str_extract(genesUnderSelection$HOGmembers,
                                            pattern = "(\\||^)CVAR(.*?)(\\||$)") 
genesUnderSelection$CVARgene <- gsub(pattern = "\\|",
                                     replacement = "",
                                     genesUnderSelection$CVARgene)

# Now run enrichment analyses for all for selection types:
selectionRegimes <- unique(genesUnderSelection$selectionCategory)

transcriptionFactorEnrichment <- function(type) {
  selectionSubset <- filter(genesUnderSelection,
                            selectionCategory == type)
  
  # Now subset all of the annotations, including promoters, to just those corresponding to genes under selection:
  annotationsUnderSelection <- filter(genesWithPromoters, 
                                      ID %in% selectionSubset$CVARgene) %>%
    subset(select = -c(ID))
  
  promotersUnderSelection <- filter(annotationsUnderSelection, 
                                    type == "promoter")
  
  # Subset all of the annotations, including promoters, to those corresponding to genes NOT under selection:
  annotationsNotUnderSelection <- filter(genesWithPromoters, 
                                         !(ID %in% genesUnderSelection$CVARgene)) %>%
    subset(select = -c(ID))
  
  promotersNotUnderSelection <- filter(annotationsNotUnderSelection, 
                                       type == "promoter")
  
  # Export those data:
  readr::write_delim(promotersUnderSelection,
                     file = paste("./14_SEAenrichment/",
                                  type,
                                  "promotersUnderSelection.gff",
                                  sep = ""),
                     delim = "\t",
                     quote = "none",
                     col_names = FALSE,
                     na = ".")
  readr::write_delim(promotersNotUnderSelection,
                     file = paste("./14_SEAenrichment/",
                                  type,
                                  "promotersNotUnderSelection.gff",
                                  sep = ""),
                     delim = "\t",
                     quote = "none",
                     col_names = FALSE,
                     na = ".")
  # Convert the promoter gff files to bed files:
  system(paste("/programs/bin/bedops/convert2bed --input=gff --do-not-sort < ./14_SEAenrichment/",
               type,
               "promotersUnderSelection.gff > ./14_SEAenrichment/",
               type,
               "promotersUnderSelection.bed",
               sep = ""))
  system(paste("/programs/bin/bedops/convert2bed --input=gff --do-not-sort < ./14_SEAenrichment/",
               type,
               "promotersNotUnderSelection.gff > ./14_SEAenrichment/",
               type,
               "promotersNotUnderSelection.bed",
               sep = ""))
  
  # Extract out .fasta sequences for each type of promoter:
  system(paste("/programs/bin/bedtools/bin/bedtools getfasta -s -fi ./CVAR/CVAR_genome_v1.0.fasta -bed ./14_SEAenrichment/",
               type,
               "promotersUnderSelection.bed -fo ./14_SEAenrichment/",
               type,
               "promotersUnderSelectionSequences.fasta",
               sep = ""))
  
  system(paste("/programs/bin/bedtools/bin/bedtools getfasta -s -fi ./CVAR/CVAR_genome_v1.0.fasta -bed ./14_SEAenrichment/",
               type,
               "promotersNotUnderSelection.bed -fo ./14_SEAenrichment/",
               type,
               "promotersNotUnderSelectionSequences.fasta",
               sep = ""))
  
  # Run enrichment analysis with MEME-SEA:
  system(paste("/programs/meme-5.5.2/bin/sea --p ./14_SEAenrichment/",
               type,
               "promotersUnderSelectionSequences.fasta --m JASPAR2022_CORE_insects_non-redundant_pfms_meme.txt --n ./14_SEAenrichment/",
               type,
               "promotersNotUnderSelectionSequences.fasta --oc ./14_SEAenrichment/",
               type,
               "/",
               sep = ""))
}

possiblytranscriptionFactorEnrichment <- purrr::possibly(transcriptionFactorEnrichment,
                                                         otherwise = "Error.")

library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)

furrr::future_map(selectionRegimes,
                  possiblytranscriptionFactorEnrichment)
