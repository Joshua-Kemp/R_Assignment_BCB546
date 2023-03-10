library(tidyverse)
library(dplyr)
library(rmarkdown)
library(sjmisc)
library(magrittr)


## Download files and transpose fang
snppositions <- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/raw/main/assignments/UNIX_Assignment/snp_position.txt", col_names=TRUE, show_col_types = FALSE) %>% select(SNP_ID, Chromosome, Position)
fang.original<- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/blob/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt?raw=true", col_names=FALSE, show_col_types = FALSE) 
#Col.names = FALSE makes sure that the header will be a row so that we have the SNP column when we transpose.  I had added the column using mutate, but this is cleaner.

#seperate by species and traspose
maize.fang <- filter(fang.original, X3 == "ZMMIL" | X3 == "ZMMLR" | X3 =="ZMMMR" | X1 == "Sample_ID") %>% rotate_df() 
teosinte.fang <- filter(fang.original, X3 == "ZMPBA" | X3 == "ZMPIL" | X3 =="ZMPJA" | X1 == "Sample_ID") %>% rotate_df() 
colnames(filter(teosinte.fang, V1 == "Sample_ID"), prefix = "teosinte_") -> colnames(teosinte.fang)
colnames(filter(maize.fang, V1 == "Sample_ID"), prefix ="maize_") -> colnames(maize.fang) 
#c(paste("maize_", (maize.fang[1, ]))) -> colnames(maize.fang)
#c(paste("teosinte_", (teosinte.fang[1, ]))) -> colnames(teosinte.fang)

## workflow 1. add sample names to col names because it is nice 2. merge genofiles with snpostions 3. remove unknown missing multiple and "" positions and write to files 4. loop or lapply through chr #s filtering and sorting
#create functions and lists
#lists
species.files.names <- c("maize.fang", "teosinte.fang")
genotypes.by.species <- list(maize.fang, teosinte.fang)
chromosomes <<- c(1:10)
problem.locations <<- c("unknown", "multiple", NA, "missing")

#functions
Merge.with.snppositions <- function(f) merge.data.frame(snppositions, f, by.x = "SNP_ID", by.y = [ ,1], all.y = TRUE, sort = FALSE)
Firstrow.to.colnames <- function(dfx) colnames(dfx)<-c(filter(dfx, SNP_ID == "Sample_ID"))
Remove.problem.postions <- function(dfx) postions.problems <- filter(dfx, Chromosome %in% problem.locations)
Filter.by.chromosome <- function(dfx) {
lapply(chromosomes, function(y) filter(dfx, Chromosome == y )) %>% lapply(function(y) filter(y, !(Position %in% problem.locations)) %>% arrange(as.numeric(Position)))
lapply(chromosomes, function(y) filter(dfx, Chromosome == y )) %>% lapply(function(y) filter(y, !(Position %in% problem.locations)) %>% arrange(desc(as.numeric(Position))))
}
merged.genotypes.by.species <-lapply(genotypes.by.species, Merge.with.snppositions)
names(merged.genotypes.by.species) <- c(paste("merged_", species.files.names, sep = ""))
genotypes.by.species <- lapply(merged.genotypes.by.species, Firstrow.to.colnames)
problem.postions <- lapply(merged.genotypes.by.species, Remove.problem.postions)
names(problem.postions) <- c(paste(problem.postions, "miss_unk_mult_", sep = ""))
split.genotype.files <- lapply(merged.genotypes.by.species, Filter.by.chromosome)
str(split.genotype.files)

#mutate(merged.genotypes.by.species[[1]], Species = "Maize", .before = Chromosome) -> maize.genotypes
#mutate(merged.genotypes.by.species[[2]], Species = "Teosinte", .before = Chromosome) -> teosinte.genotypes
#rbind(maize.genotypes,teosinte.genotypes) -> full.dataframe

