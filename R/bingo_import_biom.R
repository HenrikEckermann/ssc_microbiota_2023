# options(repr.plot.width=15, repr.plot.height=15)

library(tidyverse)
library(mia)
library(glue)
library(tidySummarizedExperiment)
library(miaViz)
library(scater)
library(readxl)



### BINGO 1

# biom_file <- here::here("data/raw_data/bingo_infancy_SILVA_138.biom")
# file.exists(biom_file)
# # because of some issue in NGT2 we need to replace NaN in the biom file with 0
# biom_lines <- read_lines(biom_file)
# sum(str_count(biom_lines, "NaN"))
# biom_lines <- str_replace_all(biom_lines, "NaN", "0.0")
# sum(str_count(biom_lines, "NaN"))
# write_lines(
#   biom_lines,
#   here::here("data/raw_data/bingo_infancy_SILVA_138_rpl.biom")
# )

# now we should be able to use read_biom 
biom_file <- here::here("data/raw_data/bingo_infancy_SILVA_138_rpl.biom")


xl <- read_excel("bingo_infancy_mapping.xlsx", sheet = "baby metadata good samples")
colnames(xl)
sample_data <- select(
  xl, 
  sample_id = `#SampleID`, 
  id = ID_participants, 
  age = Actual_Age_Collection_days)
bingo1 <- loadFromBiom(biom_file)

# use practical names for the ranks instead 
colnames(rowData(bingo1)) <- c(
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
)

# add sample ids as column
colData(bingo1)[["sample_id"]] <- rownames(colData(bingo1))
colData(bingo1) <- colData(bingo1) %>% as.data.frame() %>%
  left_join(sample_data, by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  mutate(sample_id = rownames(colData(bingo1))) %>%
  DataFrame()


tse_bingo1 <- filter(bingo1, !str_detect(sample_id, "MOCK"))
rns <- str_c(rowData(tse_bingo1)$family, rowData(tse_bingo1)$genus, sep = "_")
rns[duplicated(rns)]
rownames(tse_bingo1) <- rns
save(tse_bingo1, file = here::here("data/tse_bingo1.Rds"))



### BINGO 2


# biom_file <- here::here("data/raw_data/bingo_year13_SILVA_138.biom")
# file.exists(biom_file)
# # because of some issue in NGT2 we need to replace NaN in the biom file with 0
# biom_lines <- read_lines(biom_file)
# sum(str_count(biom_lines, "NaN"))
# biom_lines <- str_replace_all(biom_lines, "NaN", "0.0")
# sum(str_count(biom_lines, "NaN"))
# write_lines(
#   biom_lines,
#   here::here("data/raw_data/bingo_year13_SILVA_138_rpl.biom")
# )

# now we should be able to use read_biom 
biom_file <- here::here("data/raw_data/bingo_year13_SILVA_138_rpl.biom")
bingo2 <- loadFromBiom(biom_file)
colnames(rowData(bingo2)) <- c(
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
)


# add sample ids as column
colData(bingo2)[["sample_id"]] <- rownames(colData(bingo2))
colData(bingo2) <- colData(bingo2) %>% as.data.frame() %>%
  left_join(sample_data, by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  mutate(sample_id = rownames(colData(bingo2))) %>%
  DataFrame()


tse_bingo2 <- filter(bingo2, !str_detect(sample_id, "MOCK"))
rns <- str_c(rowData(tse_bingo2)$family, rowData(tse_bingo2)$genus, sep = "_")
rns[duplicated(rns)]
rownames(tse_bingo2) <- rns
save(tse_bingo2, file = here::here("data/tse_bingo2.Rds"))