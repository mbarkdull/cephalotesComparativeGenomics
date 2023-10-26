library(tidyverse)

# Read in the annotation from flo that has orphan CDS features:
failedAnnotation <- ape::read.gff("/workdir/mb2337/cephalotesComparativeGenomics/12_genomeAlignments/liftoverAnnotation/run/annotationSolelyCDS/lifted.gff3")
unique(failedAnnotation$type)

# Make an ID column and a Parent column. Then, every entry in the parent column needs to be present in the ID column:
failedAnnotation$ID <- stringr::str_extract(failedAnnotation$attributes,
                                            regex("ID[^;]+", 
                                                  ignore_case = T)) %>%
  str_split_i(pattern = "=",
              i = 2) 
failedAnnotation$Parent <- stringr::str_extract(failedAnnotation$attributes,
                                                regex("Parent[^;]+", 
                                                      ignore_case = T)) %>%
  str_split_i(pattern = "=",
              i = 2)
failedAnnotation$Name <- stringr::str_extract(failedAnnotation$attributes,
                                              regex("Name[^;]+", 
                                                    ignore_case = T)) %>%
  str_split_i(pattern = "=",
              i = 2)

# Now, for every entry in the parent column, check if that parent is present in the ID column:
failedAnnotation$isParentPresent <- failedAnnotation$Parent %in% failedAnnotation$ID

# Now, for anything that is a CDS type feature and lacks a parent, remove it:
passingAnnotations <- filter(failedAnnotation,
                             type != "CDS" | (type == "CDS" & isParentPresent != "FALSE"))

# Get rid of excess columns:
passingAnnotations <- select(passingAnnotations,
                             -c("ID",
                                "Parent",
                                "Name",
                                "isParentPresent"))

# Now export that as a .gff3 file and it should work with gff3ToGenePred:
readr::write_delim(passingAnnotations,
                   file = "passingAnnotations.gff",
                   delim = "\t",
                   quote = "none",
                   col_names = FALSE,
                   na = ".")
