library(tidyverse)
library(phylotools)
library(msa)
library(stringi)
library(spgs)

# pseudo-it produced a reference-guided genome assembly and a chain file so that you can lift over the reference annotations. 
# Here, we will do that with UCSC's liftOver tool.

# Create an output directory:
dir.create("./annotationsAndExons/")
dir.create(paste("./annotationsAndExons/",
                 i,
                 sep = ""))

# Run liftOver to create the new annotation file:
liftOverCommand <- paste("~/miniconda3/bin/liftOver -gff ./CVAR/CVAR_OGS_v1.0.gff3 ./pseudo-it/aligned",
                         i,
                         "/iter-04/fa/iter-04-softmask-final.chain ./annotationsAndExons/",
                         i,
                         "/",
                         i,
                         "annotation.gff3 ./annotationsAndExons/",
                         i,
                         "/",
                         i,
                         "unmapped -minMatch=0.4",
                        sep = "")
cat(liftOverCommand)
system(liftOverCommand)

# Now use that new annotation to extract transcript sequences:
gffreadCommand <- paste("/programs/gffread-0.9.12/gffread/gffread -w ./annotationsAndExons/",
                        i,
                        "/", 
                        i,
                        "_transcripts.fasta -g ./pseudo-it/aligned",
                        i, 
                        "/iter-04/fa/iter-04-softmask-final.fa ./annotationsAndExons/",
                        i,
                        "/",
                        i,
                        "annotation.gff3",
                        sep = "")
cat(gffreadCommand)
system(gffreadCommand)

# Now we can read those transcripts in and see how they align to the CVAR genome:
newlyGeneratedTranscriptSequences <- phylotools::read.fasta(paste("./annotationsAndExons/",
                                                                  i,
                                                                  "/",
                                                                  i,
                                                                  "_transcripts.fasta",
                                                                  sep = ""))
newlyGeneratedTranscriptSequences$seq.name <- str_split_i(string = newlyGeneratedTranscriptSequences$seq.name, 
                                                          pattern = " ",
                                                          i = 1)
newlyGeneratedTranscriptSequences$sequenceName <- paste(i, 
                                                        "_",
                                                        newlyGeneratedTranscriptSequences$seq.name, 
                                                        sep = "")

cvarSequences <- phylotools::read.fasta("CVAR/CVAR_OGS_v1.0_longest_isoform.trans.fasta")
cvarSequences$sequenceName <- paste("CVAR_",
                                    cvarSequences$seq.name, 
                                    sep = "")
allSequences <- plyr::rbind.fill(newlyGeneratedTranscriptSequences,
                                 cvarSequences)

# Write a function that will select a set of homologous sequences from the bedtools output and the desired output, then align them and check how similar they are:
checkingSimilarity <- function(sequence) {
  # Select out a set of homologous sequences:
  homologousSequences <- filter(allSequences, seq.name == sequence) %>%
    dplyr::select(sequenceName, seq.text) 
  
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
  output <- c(sequence, distanceBetweenSequences[1])
  return(output)
}
possiblyCheckingSimilarity <- possibly(checkingSimilarity, otherwise = "Error")

# Map that function over all of the genes:
subsetOfGenes <- sample(unique(allSequences$seq.name), 500)
similarityValues <- furrr::future_map(subsetOfGenes, possiblyCheckingSimilarity)
similarityValues <- as.data.frame(do.call(rbind, similarityValues))

# Plot the distribution of similarity values:
ggplot(data = similarityValues) + 
  geom_histogram(mapping = aes(x = as.numeric(V2)))






