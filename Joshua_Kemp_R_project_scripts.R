#install.packages("tidyverse", "dplyr", "sjmisc")
#install.packages("sjmisc")

library(tidyverse)
library(dplyr)
library(rmarkdown)
library(sjmisc)
library(stats)


## Download files and transpose fang
snppositions <- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/raw/main/assignments/UNIX_Assignment/snp_position.txt", col_names=TRUE, show_col_types = FALSE) %>% select(SNP_ID, Chromosome, Position)
fang.original<- read_tsv("https://github.com/EEOB-BioData/BCB546_Spring2023/blob/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt?raw=true", col_names=FALSE, show_col_types = FALSE) 
#Col.names = FALSE makes sure that the header will be a row so that we have the SNP column when we transpose.  I had added the column using mutate, but this is cleaner.

## Data Inspection
### There are genotypes in the fang file that do not belong to teosinte or maize groups we filter for looks like they come from the related group
976 + 1574
unique(duplicated(fang.original$X1))
unique(fang.original$X3)
unique(duplicated(snppositions$SNP_ID))
unique(snppositions$Chromosome)


### Seperate by species and traspose
maize.fang <- filter(fang.original, X3 == "ZMMIL" | X3 == "ZMMLR" | X3 =="ZMMMR" | X1 == "Sample_ID") %>% rotate_df() 
teosinte.fang <- filter(fang.original, X3 == "ZMPBA" | X3 == "ZMPIL" | X3 =="ZMPJA" | X1 == "Sample_ID") %>% rotate_df() 


## Functions and Lists
#### Lists
species.files.names <- c("maize.fang", "teosinte.fang")
chromosomes <- c(1:10)
problem.locations <- c("unknown", "multiple", NA, "missing")

#### Functions
Firstrow.to.colnames <- function(dfx) {
  newcolnames <- filter(dfx, V1 == "Sample_ID")
  colnames(dfx) <- c(newcolnames)
  dfx }
Merge.with.snppositions <- function(f) merge.data.frame(snppositions, f, by.x = "SNP_ID", by.y = "Sample_ID", all.y = TRUE, sort = FALSE)
Remove.problem.postions <- function(dfx) postions.problems <- filter(dfx, Position %in% problem.locations)
Replace_quest_with_hyphen <- function(dfx) as.data.frame(lapply(dfx, function(x) (gsub(pattern="?", replacement="-", x, fixed = TRUE))) -> dfx)
Filter.by.chromosome <- function(dfx) {
  append(
  lapply(chromosomes, function(y) filter(dfx, Chromosome == y )) %>% lapply(function(x) filter(x, !(Position %in% problem.locations)) %>% arrange(as.numeric(Position))),
  lapply(chromosomes, function(y) filter(Replace_quest_with_hyphen(dfx), Chromosome == y )) %>% lapply(function(x) filter(x, !(Position %in% problem.locations)) %>% arrange(desc(as.numeric(Position))))
  )
  }
Create.Split.Filenames <- {
  sapply(1:2, function(x) {
  c(
    sapply(chromosomes, function(y) c(paste0(species.files.names[x], "_chr-", y, "_ascending"))),
    sapply(chromosomes, function(y) c(paste0(species.files.names[x], "_chr-", y, "_descending"))))
  })
  }
Write_csv <- function(df) walk2(df, paste0(names(df), ".csv"), write_csv)

## Workflow
#### Add Names and Marker locations (Chromosome and Postion) to our genotype files then name them
genotypes.by.species <- list(maize.fang, teosinte.fang)
genotypes.by.species <- lapply(genotypes.by.species, Firstrow.to.colnames)
merged.genotypes.by.species <-lapply(genotypes.by.species, Merge.with.snppositions)
names(merged.genotypes.by.species) <- c(paste("merged_", species.files.names, sep = ""))
#### Remove problematic marker positions
problem.positions <- lapply(merged.genotypes.by.species, Remove.problem.postions)
names(problem.positions) <- paste0(species.files.names, "_miss_unk_mult", sep = "")
problem.positions_hyphen <- (lapply(problem.positions, function(l) Replace_quest_with_hyphen(l)))
names(problem.positions_hyphen) <- paste0(species.files.names, "_hypen_miss_unk_mult", sep = "")
#### Split dataframes by Chromosome
split.genotype.files <- lapply(merged.genotypes.by.species, Filter.by.chromosome)
Create.Split.Filenames -> newsplitnames
#### Flatten List Structure (remove nesting) and Name the split files
flat.split.genos <- list_flatten(split.genotype.files)
names(flat.split.genos) <- newsplitnames
#### write lists to csv files
Write_csv(flat.split.genos)
Write_csv(problem.positions)
Write_csv(problem.positions_hyphen)


## Visualization

### Reformatting

both_species.dataframe <- full_join(merged.genotypes.by.species[[1]], merged.genotypes.by.species[[2]], by = join_by("SNP_ID" == "SNP_ID", "Chromosome" == "Chromosome", "Position" == "Position"), suffix = c("", ""))
transposed_full.dataframe <- rotate_df(both_species.dataframe)
transposed_full.dataframe <- rotate_df(both_species.dataframe) %>% mutate(Species = ifelse(V986 == "ZMMIL" | V986 == "ZMMLR" | V986 =="ZMMMR", "Maize", ifelse(V986 == "ZMPBA" | V986 == "ZMPIL" | V986 =="ZMPJA", "Teosinte", NA)), .before = V1)
transposed_full.dataframe[1,1] <- "Species"
colnames(transposed_full.dataframe) <- filter(transposed_full.dataframe, Species == "Species")
reformatted_one_genotypecall_per_row.df <- filter(transposed_full.dataframe, !(Sample_ID == "NA" | Sample_ID == "Sample_ID")) %>% pivot_longer(cols = -c("Species", "Group", "JG_OTU", "Sample_ID")) %>% 
  merge.data.frame(snppositions, by.x = "name", by.y = "SNP_ID", all.x = TRUE, sort = FALSE) %>% 
  relocate(c(Chromosome, Position, Species, Group), .before = name) %>% 
  mutate(Allele1 = substr(value, 1, 1), .before = value) %>% 
  mutate(Allele2 = substr(value, 3, 3), .before = value) %>% 
  mutate(Homozygous = ifelse(Allele1 == "?" | Allele2 == "?", NA, eval(Allele1 == Allele2)), .before = Allele1) %>% 
  mutate(Missing_call = Allele1 == "?" | Allele2 == "?", .before = Allele1) %>% 
  filter(!(Position %in% problem.locations)) %>% 
  mutate(bin_position = cut(as.integer(Position), 20, labels = FALSE), .after = Position)

### Plotting
#### SNPs per Chromosome
ggplot(data = reformatted_one_genotypecall_per_row.df ) + geom_bar(mapping = aes(x= as.factor(as.integer(Chromosome)), color = Species )) +
  ylab('Number of tested SNPs') +
  xlab('Chromosome') 

#### Distribution accross Chromosomes
ggplot(data = reformatted_one_genotypecall_per_row.df) +
  aes(as.double(Position)/98507715) +
  geom_histogram(aes(y = ..density..), alpha = 0.4, bins = 20) + 
  geom_density(aes(as.numeric(Position)/98507715),) + 
  facet_wrap(~ as.double(Chromosome)) +
  ggtitle("SNPs accross Chromosomes") +
  xlab("Genome Position") +
  ylab('SNP density')


## this ends up being more than can easily plotted and may crash the computer.  That said, there are a number of lines that were very far from inbred, and may not have correct phenotype to genotype matching from pollen contamination. Overall, lines are less inbred than you would expect for this kind of study.
contaminated_lines <- summarise(reformatted_one_genotypecall_per_row.df, Mean_Homozygosity = mean(Homozygous, na.rm = TRUE), .by = c(Sample_ID, Species))
contaminated_lines <- arrange(contaminated_lines, Mean_Homozygosity)
head(contaminated_lines, n = 20)
ggplot(data = contaminated_lines) + geom_point(mapping = aes(x= Species, y= Mean_Homozygosity, color = Species, alpha=0.05)) + ggtitle("Sample Homozygosity by Species")

ggplot(data = contaminated_lines) + geom_point(mapping = aes(x= Sample_ID, y= Mean_Homozygosity, color = Species))


