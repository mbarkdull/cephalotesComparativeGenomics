### This script will go from aligned orthogroup files and cds files and run pal2nal on each gene, individualy, circumventing the errors you get when you run pal2nal on a batch input. ###
# First make a directory for this step:
dir.create("./testAligning")

# Set up future for paralellization:
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)

# Now save every aligned amino acid sequence as a separate fasta:
# Concatenate all of the MSA files into a single file:
aminoAcidAlignments <- list.files("./04_alignedHOGs",
                                  full.names = TRUE)

possiblyRead.fasta <- possibly(phylotools::read.fasta,
                               otherwise = "Error")
allaminoAcidAlignments <- future_map(aminoAcidAlignments,
                                     possiblyRead.fasta)
allaminoAcidAlignments <- as.data.frame(do.call(rbind, allaminoAcidAlignments))

# Get the species abbreviations:
abbreviations <- list.files(path = "./03_OrthoFinder/fasta/") %>%
  str_split_i(pattern = "_",
              i = 1)

# Create a function to create a subsetted fasta file for each species:
speciesChecking <- function(abbreviation){
  outFile <- (paste("./testAligning/", 
                    abbreviation, 
                    "_alignedProteins.fasta", 
                    sep = ""))
  # Create a column with just the species abbreviation for the gene:
  allaminoAcidAlignments$species <- str_split_i(allaminoAcidAlignments$seq.name,
                                                pattern = "_",
                                                i = 1)
  # Create a column in the concatenated file that tells us if it matches the species abbreviation:
  allaminoAcidAlignments$speciesMatch <- ifelse(grepl(abbreviation, 
                                                      allaminoAcidAlignments$species), 
                                                "True", 
                                                "False")
  # Make a dataframe that contains only the rows where $speciesMatch == True. 
  SpeciesMSA <- subset(allaminoAcidAlignments, 
                       speciesMatch == "True", 
                       select = c("seq.name", 
                                  "seq.text"))
  # Replace non-terminal stop codons with Ns:
  SpeciesMSA$seq.text <- gsub("\\*(?!-+$)", # This regex matches an asterisk, UNLESS it is followed by any number of hyphens and then an end-of-string.
                              "N", 
                              SpeciesMSA$seq.text, 
                              perl = TRUE)
  # Replace terminal Ns with stop codons:
  SpeciesMSA$seq.text <- gsub("N$", # This regex matches any N followed by an end-of-string.
                              "\\*", 
                              SpeciesMSA$seq.text, 
                              perl = TRUE)
  # Replace Ns followed by only gaps with stop codons:
  SpeciesMSA$seq.text <- gsub("(?!N)N(?=-+$)", # This regex matches any N followed by any number of hyphens and an end-of-string, UNLESS it is preceded by another N. 
                              "\\*", 
                              SpeciesMSA$seq.text, 
                              perl = TRUE)
  # Save that to a .fasta file: ../test.fasta
  dat2fasta(SpeciesMSA, 
            outfile = outFile)
}

# Apply the function to all species with future_map:
possiblySpeciesChecking <- possibly(speciesChecking,
                                    otherwise = "Error")
future_map(abbreviations,
           possiblySpeciesChecking)












aminoAcidFiles <- paste("./testAligning/",
                        abbreviations,
                        "_alignedProteins.fasta",
                        sep = "")

splittingAminoAcidAlignments <- function(i) {
  sequence_name <- get.fasta.name(i)
  sequence_group <- gsub("^CVAR_",
                         "CVAR_CVAR_",
                         sequence_name)
  sequence_group <- paste("./testAligning/protein_",
                          sequence_group,
                          sep = "")
  group <- data.frame(sequence_name, sequence_group)
  fasta <- read.fasta(i)
  split_dat(fasta, group)
}

possiblysplittingAminoAcidAlignments <- possibly(splittingAminoAcidAlignments,
                                                 otherwise = "Error")
future_map(aminoAcidFiles,
           possiblysplittingAminoAcidAlignments)

# Get every coding sequence as a separate fasta:
abbreviations <- list.files(path = "./03_OrthoFinder/fasta/") %>%
  str_split_i(pattern = "_",
              i = 1)
splittingCDS <- function(abbreviation) {
  cdsFastaFile <- paste("./02_annotationsAndExons/",
                        abbreviation,
                        "/",
                        abbreviation,
                        "_cds.fasta",
                        sep = "")
  sequence_name <- get.fasta.name(cdsFastaFile) 
  
  sequence_name_correct <- str_split_i(sequence_name,
                                       pattern = " ",
                                       i = 1) 
  
  sequence_name_correct <- paste(abbreviation,
                                 "_",
                                 sequence_name_correct,
                                 sep = "")
  sequence_group <- paste("./testAligning/cds_",
                          sequence_name_correct,
                          sep = "")
  group <- data.frame(sequence_name, sequence_group)
  fasta <- read.fasta(cdsFastaFile)
  split_dat(fasta, group)
}
possiblysplittingCDS <- possibly(splittingCDS,
                                 otherwise = "Error")
future_map(abbreviations,
           possiblysplittingCDS)

# Run pal2nal on individual pairs of files:
allGenes <- list.files(path = "./testAligning/",
                       pattern = "cds_*") %>%
  str_split_i(pattern = "\\.",
              i = 1) %>%
  str_split_i(pattern = "cds_",
              i = 2)
dir.create("./testAligning/alignedNucleotideSequences")

# Run PAL2NAL:
runPAL2NAL <- function(gene) {
  alignedProteinsFile <- paste("./testAligning/protein_", 
                               gene, 
                               ".fasta", 
                               sep = "")
  cdsSequencesFile <- paste("./testAligning/cds_", 
                            gene, 
                            ".fasta", 
                            sep = "")
  outputFileName <- paste("./testAligning/alignedNucleotideSequences/", 
                          gene, 
                          "_nucleotideAlignments.fasta", 
                          sep = "")
  # Construct a system command to run pal2nal
  pal2nalCommand <- paste("/programs/pal2nal/pal2nal.pl  ",
                          alignedProteinsFile,
                          "  ",
                          cdsSequencesFile,
                          "  -output fasta > ",
                          outputFileName,
                          sep = "")
  cat(pal2nalCommand)
  system(pal2nalCommand)
}
possiblyrunPAL2NAL <- possibly(runPAL2NAL,
                               otherwise = "Error")
# Run the alignment in parallel, using furrr:
future_map(allGenes,
           possiblyrunPAL2NAL)

