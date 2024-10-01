# Load required packages:
library(RJSONIO)
library(tidyverse)
library(ggthemes)
library(gt)
library(beepr) 

# Get a list of results files, sort by decreasing size, and remove any that are empty:
#files <- list.files(path = args[1], pattern = "*.json", full.names = TRUE)
files <- list.files(path = "./bustedPHResults/parsimony", pattern = "*.json", full.names = TRUE)
files <- sort(files, decreasing = TRUE)
files <- files[sapply(files, file.size) > 0]

readingBustedPHResults <- function(content) {
  rawResult <- RJSONIO::fromJSON(content = content)
  result <- c(rawResult[["input"]][["file name"]], 
              rawResult[["test results"]][["p-value"]], 
              rawResult[["test results background"]][["p-value"]], 
              rawResult[["test results shared distributions"]][["p-value"]], 
              rawResult[["test results"]][["LRT"]], 
              rawResult[["test results background"]][["LRT"]], 
              rawResult[["test results shared distributions"]][["LRT"]]) %>%
    as.list()
  print(rawResult[["input"]][["file name"]])
  result <- as.data.frame(do.call(cbind, result))   
  colnames(result) <- c("file", 
                        "test results p-value", 
                        "test results background p-value", 
                        "test results shared distributions p-value", 
                        "test results LRT", 
                        "test results background LRT", 
                        "test results shared distributions LRT")
  result$orthogroup <- str_split_i(result$file,
                                   pattern = "/",
                                   i = -1) %>%
    str_split_i(pattern = "_",
                i = 2)
  
  result <- result %>% mutate(selectionOn =
                                case_when(as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "ForegroundOnly",
                                          
                                          as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "SelectionOnBothButDifferent",
                                          
                                          as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                                          
                                          as.numeric(as.character(`test results p-value`)) <= 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                                          
                                          as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "BackgroundOnly",
                                          
                                          as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) <= 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                          
                                          as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) <= 0.05 ~ "NoEvidenceOfSelection",
                                          
                                          as.numeric(as.character(`test results p-value`)) > 0.05 & 
                                            as.numeric(as.character(`test results background p-value`)) > 0.05 &
                                            as.numeric(as.character(`test results shared distributions p-value`)) > 0.05 ~ "NoEvidenceOfSelection"))
  
}

library(furrr)
future::plan(multisession)
options(future.globals.maxSize= +Inf)

allResults <- furrr::future_map(files,
                                readingBustedPHResults)
beep(sound = 5)

allResults <- as.data.frame(do.call(rbind, allResults)) 
allResults <- allResults %>%  
  mutate(foregroundpValueFDR = p.adjust(`test results p-value`, 
                                        method='BH')) %>% 
  mutate(backgroundpValueFDR = p.adjust(`test results background p-value`, 
                                        method='BH')) %>% 
  mutate(sharedpValueFDR = p.adjust(`test results shared distributions p-value`, 
                                        method='BH')) 

allResults <- allResults %>% mutate(fdrSelectionOn =
                                      case_when(as.numeric(as.character(foregroundpValueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) > 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) <= 0.05 ~ "ForegroundOnly",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) <= 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) <= 0.05 ~ "SelectionOnBothButDifferent",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) <= 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) > 0.05 ~ "SelectionOnBothButNoSignificantDifference",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) <= 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) > 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithTraitButNS",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) <= 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) <= 0.05 ~ "BackgroundOnly",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) <= 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) > 0.05 ~ "EvidenceOfSelectionAssociatedWithLackOfTraitButNS",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) > 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) <= 0.05 ~ "NoEvidenceOfSelection",
                                                
                                                as.numeric(as.character(foregroundpValueFDR)) > 0.05 & 
                                                  as.numeric(as.character(backgroundpValueFDR)) > 0.05 &
                                                  as.numeric(as.character(sharedpValueFDR)) > 0.05 ~ "NoEvidenceOfSelection"))


positiveSelectionOnForeground <- filter(allResults, 
                                        fdrSelectionOn == "ForegroundOnly")

positiveSelectionOnBackground <- filter(allResults, 
                                        fdrSelectionOn == "BackgroundOnly")
