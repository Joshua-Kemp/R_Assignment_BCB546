# EEOB 546 R assignment Joshua Kemp

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
originalfang <- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/blob/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt?raw=true", col_names=FALSE, show_col_types = FALSE)
#Col.names = FALSE makes sure that the header will be a row so that we have the SNP column when we transpose.  I had added the column using mutate, but this is cleaner.
colnames(fang_labeled) <- filter(fang.transposed, V1 == "sampleID")
maize.fang <- filter(fang.labled, "Group" = )

merge.df <- merge.data.frame(snppositions, fang.transposed, by.x = "SNP_ID", by.y = "Group")
labeled.merge.df <- merge.data.frame(snppositions, fang.transposed, by.x = "SNP_ID", by.y = "V1", all.y = TRUE)




