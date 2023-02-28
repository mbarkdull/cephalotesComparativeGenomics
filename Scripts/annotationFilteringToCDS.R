### I need to filter the CVAR annotation file to just CDS features, to work with flo for annotation remapping. 
library(tidyverse)
library(phylotools)
library(msa)
library(seqinr)

# Read in the annotation:
annotation <- read.gff("./CVAR/CVAR_OGS_v1.0_longest_isoform.gff3", 
                       na.strings = c(".",
                                      "?"),
                       GFF3 = TRUE)

# Filter to just CDS features:
annotationSolelyCDS <- filter(annotation,
                              type == "CDS" |
                                type == "mRNA" |
                                type == "gene")

# Replace NAs in the score column with `.` so that gt doesn't get mad:
annotationSolelyCDS$score <- as.character(annotationSolelyCDS$score)
annotationSolelyCDS$score <- "."
# Likewise, replace NAs in the phase column:
annotationSolelyCDS$phase <- as.character(annotationSolelyCDS$phase)
annotationSolelyCDS$phase <- annotationSolelyCDS$phase %>% replace_na("unknown")
annotationSolelyCDS <- annotationSolelyCDS %>% 
  mutate(phase = str_replace(phase, "unknown", "."))

# Export the filtered annotation:
dir.create("./liftoverAnnotation/")
write.table(annotationSolelyCDS, 
            file = "./liftoverAnnotation/annotationSolelyCDS.gff3", 
            append = FALSE, 
            sep = "\t", 
            dec = ".",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)






