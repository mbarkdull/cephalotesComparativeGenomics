# Save every aligned amino acid sequence as a separate fasta:
# First make a working directory and copy the folder with the MSA files there (THIS WILL NEED TO BE A VALUE SET ON THE COMMAND LINE). 
dir.create("./testAligning")
aminoAcidAlignments <- list.files("./04_alignedHOGs",
                                  full.names = TRUE)

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
purrr::map(aminoAcidAlignments,
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
purrr::map(abbreviations,
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
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)
future_map(allGenes,
           possiblyrunPAL2NAL)

