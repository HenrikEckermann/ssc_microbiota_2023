library(tidyverse)
library(omixerRpm)

# load the KEGG orthologue table 
ko <- read.delim(
    "picrust2/data/KO_metagenome_out/pred_metagenome_unstrat.tsv",
    header = TRUE
    )

# pick most recent database 
listDB()
db <- loadDB(name = listDB()[1])
head(ko)
# calculate GBM abundance and store in df
# threads = 16, java.mem = 16
gbm <- rpm(x = ko, module.db = db, )
gbm <- asDataFrame(gbm, "abundance")
head(gbm)


# swap ids back to the sample ids I use in the other datasets
mapfiles <- list.files(here::here("data/map_files"), full.names = TRUE)
id_swap <- map_dfr(mapfiles, ~readxl::read_excel(.x)) %>%
  select(sample_id = Seq_ID, internal_sample_id) %>%
  mutate(
    sample_id = str_replace(sample_id, "SKIPPY_", ""),
    sample_id = str_replace(sample_id, "_2", "_1"),
    sample_id = str_replace(sample_id, "_52", "_3"),
    sample_id = str_replace(sample_id, "_5", "_2"),
    sample_id = str_replace(sample_id, "245_1_15-1-17", "245_1"),
    sample_id = str_replace(sample_id, "245_1_22-1-17", "245_1_2"),
    sample_id = str_replace(sample_id, "269_1_15_mei_2017", "269_1"),
    sample_id = str_replace(sample_id, "269_1_18_mei_2017", "269_1_2")
    )


modules <- select(gbm, Module, Description)
gbm <- select(gbm, -contains("MOCK"), -Description)
colnames(gbm) <- str_remove(colnames(gbm), "X")
gbm <- column_to_rownames(gbm, "Module") %>% 
    t() %>% 
    as.data.frame() %>%
    mutate_all(function(x) as.integer(x)) %>%
    rownames_to_column("internal_sample_id") %>%
    left_join(id_swap, by = "internal_sample_id") %>%
    select(-internal_sample_id)
head(gbm)
dim(gbm)


# store the gbm in a tse to then perform DAA using similar script as before
load(here::here("data/data_imp.Rds"))
load(here::here("data/data.Rds"))

fvars <- c("siblings", "condition")
# add metadata to tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  select(-siblings) %>%
  left_join(select(d, age, sample_id, siblings), by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  mutate(across(all_of(fvars), function(x) as.factor(x))) %>%
  DataFrame()
colData(tse)$age <- colData(tse)$age + as.numeric(colData(tse)$week) * 7
colData(tse)$age_s <- scale(colData(tse)$age)[, 1]
gbm <- gbm %>%
    filter(sample_id %in% colnames(tse)) %>%
    column_to_rownames("sample_id") %>%
    t()

colnames(gbm)[!colnames(gbm) %in% colnames(tse)]
colnames(tse)[!colnames(tse) %in% colnames(gbm)]
md <- colData(tse) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    arrange(factor(sample_id, levels = colnames(gbm))) %>%
    mutate(sid = sample_id) %>%
    column_to_rownames("sid") %>%
    DataFrame()

tse_pc <- TreeSummarizedExperiment(assays = gbm, colData = md)
assayNames(tse_pc) <- "counts"
save(tse_pc, modules, file = here::here("data/data_pc.Rds"))


gbm %>% dim()
