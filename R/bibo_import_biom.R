# options(repr.plot.width=15, repr.plot.height=15)

library(tidyverse)
library(mia)
library(glue)
library(tidySummarizedExperiment)
library(miaViz)
library(scater)
library(readxl)

# list.files("data/raw_data")
# biom_file <- here::here("data/raw_data/bibo_SILVA_138.biom")
# file.exists(biom_file)
# # because of some issue in NGT2 we need to replace NaN in the biom file with 0
# biom_lines <- read_lines(biom_file)
# sum(str_count(biom_lines, "NaN"))
# biom_lines <- str_replace_all(biom_lines, "NaN", "0.0")
# sum(str_count(biom_lines, "NaN"))
# write_lines(
#   biom_lines,
#   here::here("data/raw_data/bibo_SILVA_138_rpl.biom")
# )

# now we should be able to use read_biom 
biom_file <- here::here("data/raw_data/bibo_SILVA_138_rpl.biom")
bibo <- loadFromBiom(biom_file)
colnames(rowData(bibo)) <- c(
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
)
colData(bibo)[["sample_id"]] <- rownames(colData(bibo))
tse_bibo <- filter(bibo, !str_detect(sample_id, "MOCK"))


rns <- str_c(rowData(tse_bibo)$family, rowData(tse_bibo)$genus, sep = "_")
rns[duplicated(rns)]
rownames(tse_bibo) <- rns
rownames(tse_bibo)


save(tse_bibo, file = here::here("data/tse_bibo.Rds"))

