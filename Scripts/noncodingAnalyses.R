library("rphast")

#### Convert the whole genome alignment to an alignment of just one scaffold, in the correct format: ####
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
# Filter to a single scaffold and reformat:
alignment <- alignment %>%
  tidyr::fill(scaffold,
              .direction = "up") %>%
  filter(scaffold == "CVAR_scaf_0") %>% 
  select(-c(scaffold))

# Add the appropriate header row:
headerRow <- c("##maf version=1",
               "",
               "",
               "",
               "",
               "",
               "")
export <- rbind(headerRow, 
                alignment)

# Export the .maf for a single scaffold:
write_delim(export, 
            file = "test.maf",
            na = "",
            delim = " ",
            quote = "none",
            col_names = FALSE)

#### Run PhasCons ####
# read myAlignment
myAlign <- read.msa("test.maf")

# read gene annotations from a UCSC "genepred" file, created by the script ./Scripts/alignmentsForPhast:
myFeats <- read.feat("./12_genomeAlignments/CVAR.gp") %>%
  filter(seqname == "CVAR_scaf_0")
# Replace the scaffold seqname with the name of the whole genome, since that's how things are named in the alignment msa:
myFeats$seqname <- "CVAR_genome"
# Add features corresponding to introns:
myFeats <- add.introns.feat(myFeats)
# Filter out exons:
myFeats <- myFeats[myFeats$feature != "exon",]
# Get a summary of the features on this scaffold:
table(myFeats$feature)

# Define your species tree in Newick format:
myTree <- "((POW0461_genome,CVAR_genome),((CSM3685_genome,CSM3441_genome),(POW0123_genome,CSM3677_genome)));"

# Create a feature corresponding to the whole scaffold:
myWholeChrom <- feat(seq = "CVAR_genome", # This is the sequence, or scaffold, name, which we changed to CVAR_genome
                     src=".", 
                     feature="all",
                     start = 1,
                     end = 1 + ncol.msa(myAlign, 
                                        "CVAR_genome")) # This is the genome in the alignment to use as the reference

# Add annotations for intergenic regions:
myIntergenicFeats <- inverse.feat(myFeats, region.bounds = myWholeChrom)
myIntergenicFeats$feature <- "intergenic"
myFeats <- rbind.feat(myFeats, myIntergenicFeats)

# Extract fourfold degenerate sites from the scaffold alignment:
myAlign4d <- get4d.msa(myAlign, 
                       features = myFeats)
# Using those sites, estimate a neutral model:
myNeutralMod <- phyloFit(myAlign4d, 
                         tree = myTree, 
                         subst.mod = "REV")

# Get conservation scores for each base and predict conserved elements:
myPc <- phastCons(msa = myAlign, 
                  mod = myNeutralMod, 
                  expected.length = 20,
                  target.coverage = 0.125, 
                  viterbi = TRUE)
names(myPc)

# Extract the most conserved elements:
myConsElements <- myPc$most.conserved

# This shows how many bases, and what percentage of bases, are predicted to be conserved
coverage.feat(myConsElements)
coverage.feat(myConsElements)/coverage.feat(myWholeChrom)

# The posterior probabilities for every base are here:
names(myPc$post.prob.wig)
dim(myPc$post.prob.wig)

# And the overall likelihood is here:
myPc$likelihood

# For comparison, we will produce an alternative set of conservation scores using phyloP.
myPp <- phyloP(myNeutralMod, 
               myAlign, 
               method = "LRT", 
               mode = "CONACC")
# the returned object is a data frame giving statistics for every base in the alignment
names(myPp)
dim(myPp)

# Now we can plot genes, conserved elements, and conservation scores for a segment of the scaffold:
# Get the coding sequences and make them a track:
myCodingFeats <- myFeats[myFeats$feature=="CDS",]
myGeneTrack <- as.track.feat(myCodingFeats, 
                             "genes", 
                             is.gene = TRUE)
# Get the conserved elements as a track:
myConsElTrack <- as.track.feat(myConsElements, 
                               "phastCons most conserved", 
                               col = "red") 
# Plot the conservation scores as a track:
myPhastConsScoreTrack <- as.track.wig(wig = myPc$post.prob.wig,
                                      name = "phastCons post prob", 
                                      col = "red", 
                                      ylim = c(0, 1)) 
myPhyloPTrack <- as.track.wig(coord = myPp$coord, 
                              score = myPp$score, 
                              name="phyloP score",
                              col = "blue", 
                              smooth = TRUE, 
                              horiz.line = 0) 
# Plot all of the tracks on a segment of the scaffold from position 0 to 20000:
plot.track(list(myGeneTrack, 
                myConsElTrack, 
                myPhastConsScoreTrack, 
                myPhyloPTrack),
           xlim = c(0, 200000), 
           cex.labels = 1.25,
           cex.axis = 1.25, 
           cex.lab = 1.5)


# Examine the distribution of lengths of the conserved elements:
# Get just the conserved elements:
myCe <- myPc$most.conserved
plot(density.feat(myCe), 
     ylim = c(0, 
            0.018),
     main = "Element Length by Type", 
     xlab = "Length",
     mgp = c(1.5,
             0.5,
             0),
     mar = c(2,
             2,
             2,
             2))
# Get the elements that overlap genes by at least 50 percent:
myCodingConsEle <- overlap.feat(myCe, 
                                myCodingFeats, 
                                min.percent = 0.5)
# obtain elements that overlap by less than 50 percent:
myNoncodingConsEle <- overlap.feat(myCe, 
                                   myCodingFeats, 
                                   min.percent = 0.5,
                                   overlapping = FALSE)
# Plot those two sub-types:
lines(density.feat(myCodingConsEle), 
      col = "red")
lines(density.feat(myNoncodingConsEle), 
      col = "blue")
legend(c("All", 
         "Coding", 
         "Noncoding"), 
       x = "topright", 
       inset = 0.05,
       lty = 1, 
       col = c("black", 
               "red", 
               "blue"))



# Plot and compare the enrichment of conserved elements in different region types:
par(mfrow = c(2, 
              2), 
    cex.main = 1.5, 
    cex.lab = 1.5, 
    cex.axis = 1.5, 
    mar = c(5,
            5,
            4,
            2))
# look at fold-enrichment of each annotation type by conserved element
myEnrich <- enrichment.feat(x = myCe, 
                            annotations = myFeats, 
                            region.bounds = myWholeChrom)
col <- rainbow(nrow(myEnrich))
barplot(myEnrich$enrichment, 
        col = col,
        main = "Enrichment of\nConserved Elements",
        ylab = "Fold Enrichment")
plot.new()
legend(x = "center", 
       legend = myEnrich$type, 
       fill = col,
       cex = 1.5)
# look at the composition of the conserved elements
myComp <- composition.feat(myCe, 
                           myFeats)
pie(myComp$composition, 
    col = rainbow(nrow(myComp)), 
    radius = 1.0,
    main = "Composition of\nConserved Elements", 
    labels = NA)
# compare with background composition
myComp <- composition.feat(myWholeChrom, 
                           myFeats)
pie(myComp$composition, 
    col = rainbow(nrow(myComp)), 
    radius = 1.0,
    main = "Background\nComposition", 
    labels = NA)



# Identify conserved elements again, this time without fixed transition probabilities. 
# instead, those probabilities are estimated by maximum likelihood:
myPcEM <- phastCons(myAlign, 
                    myNeutralMod, 
                    viterbi = TRUE, 
                    estimate.transitions = TRUE)
names(myPcEM)
myPcEM$transition.rates
myPcEM$likelihood

# Now we can compare the coverage of conserved elements between the two estimation methods:
coverage.feat(myPcEM$most.conserved)
coverage.feat(myPcEM$most.conserved, 
              myPc$most.conserved)
coverage.feat(myPcEM$most.conserved, 
              myPc$most.conserved, 
              or = TRUE)
coverage.feat(myPcEM$most.conserved, 
              myPc$most.conserved,
              not = c(FALSE, 
                      TRUE))
coverage.feat(myPcEM$most.conserved, 
              myPc$most.conserved,
              not = c(TRUE, 
                      FALSE))
plot.track(list(as.track.feat(myPc$most.conserved, 
                              name="No estimation"),
                as.track.feat(myPcEM$most.conserved, 
                              name="With estimation")))


plot(density.feat(myPc$most.conserved),
     main = "Distribution of Element Lengths", 
     xlab = "Length", 
     xlim = c(0,
              1000))
lines(density.feat(myPcEM$most.conserved), 
      col = "red")
legend(x = "topright",
       inset = 0.05, 
       c("without estimation", 
         "with estimation"),
       lty = 1, 
       col = c("black", 
               "red"))
