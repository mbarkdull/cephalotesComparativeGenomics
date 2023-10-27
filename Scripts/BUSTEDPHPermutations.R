# List all of the species:
species <- list.files("02_annotationsAndExons/")
allPossibleCombinationsOfSpecies <- combn(species, 3)
allPossibleCombinationsOfSpecies <- apply(allPossibleCombinationsOfSpecies, 
                                          2, 
                                          as.list)
numberOfCombinations <- 1:length(allPossibleCombinationsOfSpecies)

# Label trees for all of the runs:
labelForAllRuns <- function(run) {
  # Create an output directory:
  dir.create(paste("./15_BUSTEDPHPermutations/",
                   run,
                   "_run/labelled/",
                   sep = ""),
             recursive = TRUE)
  # Get the species for that combination:
  speciesSetForThisCombination <- unlist(allPossibleCombinationsOfSpecies[[run]])
  # List all of the tree files:
  orthofinderTreeFiles <- list.files(path = "./05_HOGTrees",
                                     full.names = TRUE)    
  # Create the function that labels the trees for a single run:
  labelAllTreesForOneRun <- function(i) {
    # Get the filename alone:
    filename <- str_split_i(i,
                            pattern = "/",
                            i = 3)
    # Get the orthogroup number:
    orthogroup <- str_split_i(filename,
                              pattern = "_",
                              i = 1)
    # Copy the unlabelled tree to the output folder; we'll be editing this tree three times. 
    file.copy(from = i,
              to = paste("./15_BUSTEDPHPermutations/",
                         run,
                         "_run/labelled/",
                         orthogroup,
                         "_tree.txt",
                         sep = ""))
    # Read in a single tree and get the tip labels of that tree:
    text <- readr::read_file(paste("./15_BUSTEDPHPermutations/",
                                   run,
                                   "_run/labelled/",
                                   orthogroup,
                                   "_tree.txt",
                                   sep = ""))
    # In case hyphy has already run on this, add back in the semicolon that ape requires:
    text <- gsub("\n",
                 "",
                 text,
                 perl = TRUE)
    text <- gsub("$(?<!;)",
                 ";",
                 text,
                 perl = TRUE)
    tree <- ape::read.tree(text = text)
    tips <- as.data.frame(tree[["tip.label"]])
    
    labelSingleSpecies <- function(focalSpecies) {
      # Get the list of tip labels with matches to the species of interest:
      tipsToLabel <- filter(tips,
                            str_split_i(tips$`tree[["tip.label"]]`, pattern = "_", i = 1) %in% focalSpecies)
      tipsToLabel <- tipsToLabel$`tree[["tip.label"]]`
      base::write(x = tipsToLabel, 
                  file = paste(run,
                               "tipsToLabel_",
                               focalSpecies,
                               "_",
                               filename,
                               sep = ""))
      # Run the Hyphy labelling script:
      hyphyCommand <- paste("/programs/hyphy-20210923/bin/hyphy hyphy-analyses/LabelTrees/label-tree.bf --tree ./15_BUSTEDPHPermutations/",
                            run,
                            "_run/labelled/",
                            orthogroup,
                            "_tree.txt --list ",
                            run,
                            "tipsToLabel_",
                            focalSpecies,
                            "_",
                            orthogroup,
                            "_tree.txt --output ./15_BUSTEDPHPermutations/",
                            run,
                            "_run/labelled/",
                            orthogroup,
                            "_tree.txt --internal-nodes \"All descendants\"",
                            sep = "")
      #cat(hyphyCommand)
      system(hyphyCommand)
      file.remove(paste(run,
                        "tipsToLabel_",
                        focalSpecies,
                        "_",
                        filename,
                        sep = ""))
    }
    possiblyLabelSingleSpecies <- possibly(labelSingleSpecies, otherwise = "Error")
    purrr::map(speciesSetForThisCombination,
               possiblyLabelSingleSpecies)
  }      
  possiblyLabelAllTreesForOneRun <- possibly(labelAllTreesForOneRun,
                                             otherwise = "Error")
  # Run that function:
  purrr::map(orthofinderTreeFiles,
             possiblyLabelAllTreesForOneRun)
}
possiblyLabelForAllRuns <- possibly(labelForAllRuns,
                                    otherwise = "Error")
# Set up future for paralellization:
library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)

furrr::future_map(numberOfCombinations,
                  possiblyLabelForAllRuns)

