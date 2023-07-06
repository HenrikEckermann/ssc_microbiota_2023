# tse_bingo, tse_bibo and tse_skippy must go through same steps in terms
# of preprocessing
library(tidyverse)
library(mia)
library(glue)
library(tidySummarizedExperiment)
library(scater)
library(readxl)


# biom_file <- here::here("data/raw_data/skippy_non_concat_SILVA_138.biom")
# file.exists(biom_file)
# # because of some issue in NGT2 we need to replace NaN in the biom file with 0
# biom_lines <- read_lines(biom_file)
# sum(str_count(biom_lines, "NaN"))
# biom_lines <- str_replace_all(biom_lines, "NaN", "0.0")
# sum(str_count(biom_lines, "NaN"))
# write_lines(
#   biom_lines,
#   here::here("data/raw_data/skippy_non_concat_SILVA_138_rpl.biom")
# )

# now we should be able to use read_biom 
biom_file <- here::here("data/raw_data/skippy_non_concat_SILVA_138_rpl.biom")


# extract skippy id and week from map files
raw_path <- here::here("data/map_files")
map_files <- list.files(raw_path, pattern = ".xlsx$")
map_files %>% length()
lib_nums <- str_extract(map_files, "\\d\\d\\d\\d_00\\d\\d")
LibraryNumber_merged <- c(1:6)
map_file_merged <- map2_dfr(map_files[1:6], LibraryNumber_merged, function(filename, lnm) {
  # extract library number
  lib_num <- str_extract(filename, "\\d\\d\\d\\d_00\\d\\d")
  # read provided map file
  xl <- read_excel(glue("{raw_path}/{filename}"))
  map_file <- xl %>% 
    mutate(
      LibraryNumber = lnm,
      ProjectName = ifelse(is.na(ProjectName), "unspecified", ProjectName)
    ) %>% 
    filter(!ProjectName == "Shime-uncooled_ileal_digesta") %>%
    select(
      "sample_id" = internal_sample_id,
      id = Seq_ID,
      ProjectName   
    ) %>%
      filter(!is.na(`sample_id`))
    
    return(map_file)
})
# make sample_ids unique as they are in the tse object later 
map_file_merged <- map_file_merged %>% 
  mutate(
  namealt = glue("{sample_id}_{1:dim(map_file_merged)[1]}"),
  sample_id = ifelse(str_detect(sample_id, "MOCK"), namealt, sample_id))

# extract skippy id and week
# the following provided duplicates and need extra labeling later 
# when I create the new sample names (see further below)
duplicates <- c(
  "245_2_22-1-17", 
  "240_52_11-2-18",
  "269_2_18_mei_2017", 
  "276_5_26052017"
)

sample_data <- str_match(map_file_merged$id, ".*(\\d\\d\\d)(_|-)(\\d+).*") %>%
  as_tibble() %>%
  select(id = V1, skippy_id = V2, week = V4) %>%
  full_join(map_file_merged, by = "id")  %>%
  filter(!is.na(sample_id)) %>%
  mutate(duplicated = ifelse(id %in% duplicates, TRUE, FALSE)) %>%
  select(wur_id = id, skippy_id, sample_id, week, duplicated)
filter(sample_data, duplicated)


skippy <- loadFromBiom(biom_file)
colnames(rowData(skippy)) <- c(
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
)

# add sample ids as column
colData(skippy)[["sample_id"]] <- rownames(colData(skippy))
rownames(colData(skippy))
colData(skippy)[["sample_id"]]



colData(skippy) <- colData(skippy) %>% as.data.frame() %>%
  left_join(sample_data, by = "sample_id") %>%
  mutate(neg_control = as.factor(ifelse(is.na(week), 1, 0))) %>%
  column_to_rownames("sample_id") %>%
  mutate(sample_id = rownames(colData(skippy))) %>%
  DataFrame()

# we will now store MOCK samples in a separate tse 
mock_samples <- filter(skippy, str_detect(sample_id, "MOCK"))

# also check they cluster together closely when using ordination
colData(skippy) %>% dim()
skippy <- mutate(skippy, 
  mock = as.factor(ifelse(str_detect(sample_id, "MOCK"), 1, 0))
)
colData(skippy)[["library_size"]] <- colSums(assay(skippy))
# perform NMDS coordination method

# MOCK samples do cluster together but particularly 1 outlier is there
tse_skippy <- filter(skippy, !str_detect(sample_id, "MOCK"))


# library sizes (we need to first get orignal IDs in order to filter out 
# negative control)
colData(tse_skippy) %>% as.data.frame() %>% filter(!week %in% c(2, 5, 52))

# the low ls are negative controls as we can see in above table
colData(tse_skippy) %>% as.data.frame() %>%
  arrange(library_size)
# exclude negative controls
tse_skippy <- filter(tse_skippy, neg_control == 0) %>%
  mutate(sample_ids_new = glue("{skippy_id}_{ifelse(week == 2, 1, ifelse(week == 5, 2, ifelse(week == 52, 3, 'CHECK')))}{ifelse(duplicated, '_2', '')}"))



# we retain 331 samples at this step 
# lets change row and column names to better descriptions 
colnames(tse_skippy) <- colData(tse_skippy)[, "sample_ids_new"]
colData(tse_skippy)[, "sample_ids_new"]
tse_skippy <- filter(skippy, !str_detect(sample_id, "MOCK"))
rns <- str_c(rowData(tse_skippy)$family, rowData(tse_skippy)$genus, sep = "_")
rns[duplicated(rns)]
rownames(tse_skippy) <- rns
save(tse_skippy, file = here::here("data/tse_skippy.Rds"))

