# EEOB 546 R assignment Joshua Kemp

setwd(.)
## Create environment with the correct packages
install.packages("tidyverse", "dplyr", "sjmisc")
install.packages("sjmisc")
yes
yes
install.packages("gitr")
library(tidyverse)
library(dplyr)
library(rmarkdown)
library(sjmisc)

## Download files and transpose fang
snppositions <- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/raw/main/assignments/UNIX_Assignment/snp_position.txt", col_names=TRUE, show_col_types = FALSE) %>% select(SNP_ID, Chromosome, Position)
fang.original<- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/blob/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt?raw=true", col_names=FALSE, show_col_types = FALSE) 
#Col.names = FALSE makes sure that the header will be a row so that we have the SNP column when we transpose.  I had added the column using mutate, but this is cleaner.

#seperate by species and traspose
maize.fang <- filter(fang.original, X3 == "ZMMIL" | X3 == "ZMMLR" | X3 =="ZMMMR" | X1 == "Sample_ID") %>% rotate_df()
teosinte.fang <- filter(fang.original, X3 == "ZMPBA" | X3 == "ZMPIL" | X3 =="ZMPJA" | X1 == "Sample_ID") %>% rotate_df()


#Test renaming columns to the sample_ID
fang_labeled <- originalfang
colnames(fang_labeled) <- filter(fang_labeled, V1 == "Sample_ID")


merge.df <- merge.data.frame(snppositions, fang_labeled, by.x = "SNP_ID", by.y = "Sample_ID", all.y = TRUE, sort = FALSE)
labeled.merge.df <- merge.data.frame(snppositions, fang.transposed, by.x = "SNP_ID", by.y = "V1", sort = FALSE)


## workflow 1. add sample names to col names because it is nice 2. merge genofiles with snpostions 3. remove unknown missing multiple and "" positions and write to files 4. loop or lapply through chr #s filtering and sorting
#create functions and lists
  #lists
  species.files.names <- c("maize.fang", "teosinte.fang")
  genotypes.by.species <- list(maize.fang, teosinte.fang)
  chromosomes <- c(1:10)
  problem.locations <- c("unknown", "multiple", NA, "missing")
  unk_miss_mult.locations <- c
  #fuctions
  Merge.with.snppositions <- function(f) merge.data.frame(snppositions, f, by.x = "SNP_ID", by.y = "V1", all.y = TRUE, sort = FALSE)
  
  Firstrow.to.colnames <- function(dfx) colnames(dfx)<-c(filter(dfx, SNP_ID == "Sample_ID"))
  Remove.problem.postions <- function(dfx) postions.problems <- filter(dfx, Chromosome %in% problem.locations)
  
#merged.genotypes.by.species <- merge.data.frame(snppositions, genotypes.by.species, by.x = "SNP_ID", by.y = "V1", all.y = TRUE, sort = FALSE)


merged.genotypes.by.species <-lapply(genotypes.by.species, Merge.with.snppositions)
names(merged.genotypes.by.species) <- c(paste("merged_", species.files.names, sep = ""))
lapply(merged.genotypes.by.species, Firstrow.to.colnames)
problem.postions <- lapply(merged.genotypes.by.species, Remove.problem.postions)
names(problem.postions) <- c(paste(problem.postions, "miss_unk_mult_", sep = ""))

#assign column names

unique(snppositions$Chromosome)


