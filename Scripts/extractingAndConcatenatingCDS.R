library(tidyverse)
library(phylotools)
library(msa)
library(stringi)
library(spgs)

# The goal here is to test how well sequences extracted from the genome 
# assembly based on the GFF3 annotation file match to the CVAR/CVAR_OGS_v1.0_longest_isoform.cds.fasta sequences.

# First, we want to extract and concatenate the CDS sequences identified by the Cephalotes varians genome annotation. 

# Read in the actual nucleotide sequences of each feature in the annotation file. 
# This file was generated using bedtools getfasta on the file CVAR/CVAR_genome_v1.0.fasta
#bedtoolsOutput <- phylotools::read.fasta("./codingSequences/CSM3441_allFeatures.fasta")
bedtoolsOutput <- phylotools::read.fasta("testCDSExtraction.fasta")

# Read in the annotation file itself:
annotationGFF3 <- read.gff("./CVAR/CVAR_OGS_v1.0_longest_isoform.gff3", 
                       na.strings = c(".",
                                      "?"),
                       GFF3 = TRUE)
annotationGFF3$Parent <- stringr::str_extract(annotationGFF3$attributes,
                                          regex("Parent[^;]+", 
                                                ignore_case = T))
annotation <- read.gff("floAttempt/run/annotationSolelyCDSTidiedTranscripts/lifted_cleaned.gff", 
                       na.strings = c(".",
                                      "?"),
                       GFF3 = FALSE)
colnames(annotation) <- c("seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes")

# Split the attributes column to get individual columns for ID, Parent, and Name:
annotation$ID <- stringr::str_extract(annotation$attributes,
                                      regex("ID[^;]+", 
                                            ignore_case = T))
annotation$Parent <- stringr::str_extract(annotation$attributes,
                                          regex("Parent[^;]+", 
                                                ignore_case = T))
annotation$Name <- stringr::str_extract(annotation$attributes,
                                        regex("Name[^;]+", 
                                              ignore_case = T))
# For some reason, the sequence start value from the annotation file is 1 greater than the values in the CVAR all features sequences file
# So fix that by subtracting 1 from the start value given in the annotation:
annotation$startValueInAllFeatures <- annotation$start-1
# Create a seq.name that will match to the seq.name in the CVAR all features file produced by bedtools getfasta:
annotation$seq.name <- paste(annotation$seqid, 
                             ":", 
                             annotation$startValueInAllFeatures, 
                             "-", 
                             annotation$end, 
                             sep = "")

# Do a full join on the two dataframes, so that we have the attributes of each feature, together with the sequence of each feature:
mergedData <- full_join(bedtoolsOutput, 
                        annotation, 
                        by = c("seq.name" = "seq.name"))

# Filter to just the CDS features:
justCDS <- subset(mergedData, type == "CDS")

#### Handle the positive strand genes ####
# Get a dataframe with just the positive strand genes (sequences and annotation information):
justPositive <- subset(justCDS, strand == "+")

# Group by their parent, so that all of the exons for each gene are grouped together, then arrange by start site:
justPositiveGrouped <- group_by(justPositive, Parent) %>% 
  arrange(start) %>%
  distinct()

# And paste together the nucleotide sequences for each group, so that the individual exon sequences are combined into a full gene sequence:
positivecdsSequences <- justPositiveGrouped %>%
  summarise(seq.text = paste(seq.text, collapse = ""))
# Give the sequences a name that will allow them to match the name in the CVAR longest isoform CDS file; that way we can align them later and check if our method worked:
positivecdsSequences$seq.name <- stringr::str_remove(positivecdsSequences$Parent, "Parent=") 
# And give them a unique name so we can tell where they came from:
positivecdsSequences$sequenceName <- paste("generated_", positivecdsSequences$seq.name, sep = "")

#### Handle the negative strand genes ####
# Now subset out the negative sense coding sequences. These are a bit trickier to extract and combine. 
justNegative <- subset(justCDS, strand == "-")

# Write a function to reverse complement and then combine the coding sequences for a single gene:
concatenatingMinusSenseGenes <- function(i) {
  # Get one gene as a test case:
  gene00002 <- subset(justNegative, Parent == i)
  
  # Sort by the start value:
  gene00002 <- gene00002 %>% 
    arrange(start) %>%
    distinct()
  
  # Complement each exon:
  for (row in 1:nrow(gene00002)) {
    strand <- gene00002[row, "strand"]
    sequenceRow  <- gene00002[row, "seq.text"]
    # sequence <- as.character(sequenceRow[1,1])
    
    if (strand == "-") {
      gene00002[row, "seq.text"] <- spgs::complement(sequenceRow,
                                                     case = "as is" )
    } else {
      print("Plus")
    }
  }
  
  # Group by their parent and then arrange in descending order within the group:
  gene00002 <- group_by(gene00002, Parent) %>% 
    arrange(start, .by_group = TRUE)
  
  # And paste together the nucleotide sequences for each group:
  gene00002 <- gene00002 %>%
    summarise(seq.text = paste(seq.text, collapse = ""))
  
  # Then reverse:
  gene00002$seq.text <- stringi::stri_reverse(gene00002$seq.text)
  
  # Get the results:
  results <- c(gene00002$Parent, gene00002$seq.text)
  return(results)
}

# Make a safe version with possibly:
possiblyConcatenatingMinusSenseGenes <- possibly(concatenatingMinusSenseGenes, 
                                                 otherwise = "error")

# Map it over all the minus sense genes:
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)
justNegativeSequences <- furrr::future_map(justNegative$Parent,
                                    possiblyConcatenatingMinusSenseGenes,
                                    .progress = TRUE)
justNegativeSequences <- as.data.frame(do.call(rbind, justNegativeSequences)) %>%
  distinct()
colnames(justNegativeSequences) <- c("Parent", "seq.text")
# Give the sequences a name that will allow them to match the name in the CVAR longest isoform CDS file:
justNegativeSequences$seq.name <- stringr::str_remove(justNegativeSequences$Parent, "Parent=") 
# And give them a unique name so we can tell where they came from:
justNegativeSequences$sequenceName <- paste("generated_", justNegativeSequences$seq.name, sep = "")

#### Combine both sets ###
allConcatenatedSequences <- plyr::rbind.fill(justNegativeSequences, positivecdsSequences) %>%
  select(sequenceName, seq.name, seq.text) 

#### Check the new sequences align to the concatenated ones: ####
# Read in the longest isoforms cds file:
desiredOutput <- phylotools::read.fasta("./CVAR/CVAR_OGS_v1.0_longest_isoform.trans.fasta")
desiredOutput$sequenceName <- paste("real_",
                                    desiredOutput$seq.name,
                                    sep = "")
# Merge the two data frames top-to-bottom:
bedtoolsOutputAndCorrectOutput <- plyr::rbind.fill(allConcatenatedSequences, desiredOutput) %>%
  select(sequenceName, seq.name, seq.text) 

# Get a list of unique genes:
uniqueGenesToAlign <- unique(allConcatenatedSequences$seq.name)

# Write a function that will select a set of homologous sequences from the bedtools output and the desired output, then align them and check how similar they are:
checkingSimilarity <- function(i) {
  # Select out a set of homologous sequences:
  homologousSequences <- filter(bedtoolsOutputAndCorrectOutput, seq.name == i) %>%
    select(sequenceName, seq.text) 
  
  # Replace any Xs with Ns:
  homologousSequences$seq.text <- str_replace_all(homologousSequences$seq.text,'X','N')
  homologousSequences$seq.text <- str_replace_all(homologousSequences$seq.text,'x','n')
  
  # Convert to a named character vector, the input format for DNAStringSet:
  homologousSequencesVector <- as.character(homologousSequences$seq.text)
  names(homologousSequencesVector) <- homologousSequences$sequenceName
  
  # Convert them to an DNAStringSet, the input format for msa:
  dnaStringSet <- Biostrings::DNAStringSet(homologousSequencesVector,
                                           use.names = TRUE)
  # Align with msa:
  alignedGenes <- msa(dnaStringSet)
  # Calculate a conservation score:
  alignedGenesSequinr <- msaConvert(alignedGenes, type="seqinr::alignment")
  distanceBetweenSequences <- seqinr::dist.alignment(alignedGenesSequinr, "identity")
  output <- c(i, distanceBetweenSequences[1])
  return(output)
}
possiblyCheckingSimilarity <- possibly(checkingSimilarity, otherwise = "Error")

# Map that function over all of the genes:
subsetOfGenes <- sample(uniqueGenesToAlign, 1000)
similarityValues <- furrr::future_map(subsetOfGenes, possiblyCheckingSimilarity)

similarityValues <- as.data.frame(do.call(rbind, similarityValues))
similarityValues$Name <- paste("Name=",
                                 similarityValues$V1,
                                 sep = "")
similarityValues <- separate(similarityValues, 
                 Name, 
                 into = c("Name"),
                 sep = "-")

similarityValues <- left_join(similarityValues, 
                  filter(annotation,
                         type == "gene"),
                  by = c("Name" = "Name")) %>%
  select(c("V1",
           "V2",
           "strand",
           "phase"))

# Check which genes failed to align altogether, and fix them:
failedGenes <- setdiff(subsetOfGenes, similarityValues$V1)
failedGenes <- furrr::future_map(failedGenes, possiblyCheckingSimilarity)
failedGenes <- as.data.frame(do.call(rbind, failedGenes))
failedGenes$Name <- paste("Name=",
                          failedGenes$V1,
                               sep = "")
failedGenes <- separate(failedGenes, 
                             Name, 
                             into = c("Name"),
                             sep = "-")

failedGenes <- left_join(failedGenes, 
                              filter(annotation,
                                     type == "gene"),
                              by = c("Name" = "Name")) %>%
  select(c("V1",
           "V2",
           "strand",
           "phase"))

# Combine all of the alignment results:
similarityValues <- rbind(similarityValues,
                          failedGenes) %>%
  dplyr::filter(V1 != "Error") %>%
  distinct()

# Plot the distribution of similarity values:
ggplot(data = similarityValues) + 
  geom_histogram(mapping = aes(x = as.numeric(V2)))

# Check what's going on with any that align poorly:
poorAlignments <- filter(similarityValues, 
                         V2 != 0)
poorAlignments <- poorAlignments$V1
poorAlignments <- paste("Parent=",
                        poorAlignments,
                        sep = "")
# See which scaffold they are on in the original annotation:
originalScaffoldOfPoorAlignments <- filter(annotationGFF3,
                                           Parent %in% poorAlignments)
originalScaffoldOfPoorAlignments <- originalScaffoldOfPoorAlignments %>%
  select(seqid) %>%
  distinct()
# See if that scaffold is in the aligned genome yet:
alignedGenome <- phylotools::read.fasta("./alignedGenomes/CVAR_alignedScaffolds.fasta")
alignedScaffolds <- alignedGenome$seq.name
rm(alignedGenome)
intersect(originalScaffoldOfPoorAlignments$seqid,
          alignedScaffolds)
