library(tidyverse)
library(mia)
library(glue)
library(tidySummarizedExperiment)
library(miaViz)
library(scater)
library(readxl)
library(lubridate)

# extract skippy id and week from map files
raw_path <- here::here("data/map_files")
map_files <- list.files(raw_path, pattern = "2021_\\d{4}.*.xlsx$")
map_files %>% length()
lib_nums <- str_extract(map_files, "\\d\\d\\d\\d_00\\d\\d")
library_number_merged <- c(1:6)
map_file_merged <- map2_dfr(
  map_files,
library_number_merged,
function(filename, lnm) {
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


# the following provided duplicates and need extra labeling later
# when I create the new sample names (see further below)
duplicates <- c(
  "245_2_22-1-17",
  "240_52_11-2-18",
  "269_2_18_mei_2017",
  "276_5_26052017"
)

# create a dataframe that has id from wur, our id, an indicator whether we have
# 2 samples from an id
sample_data <- str_match(map_file_merged$id, ".*(\\d\\d\\d)(_|-)(\\d+).*") %>%
  as_tibble() %>%
  select(id = V1, skippy_id = V2, week = V4) %>%
  full_join(map_file_merged, by = "id")  %>%
  filter(!is.na(sample_id)) %>%
  mutate(duplicated = ifelse(id %in% duplicates, TRUE, FALSE)) %>%
  select(wur_id = id, skippy_id, sample_id, week, duplicated)
filter(sample_data, duplicated)

# import the biome file to combine with above metadata

# read_biom does not work with the output of NGT2 in its current version
# we must replace certain NaNs with 0 (not relevant for our analysis)
biom_file <- here::here("data/raw_data/skippy_non_concat_SILVA_138_rpl.biom")
if (!file.exists(biom_file)) {
  biom_file_old <- here::here("data/raw_data/skippy_non_concat_SILVA_138.biom")
  biom_lines <- read_lines(biom_file_old)
  biom_lines <- str_replace_all(biom_lines, "NaN", "0.0")
  write_lines(
    biom_lines,
    here::here("data/raw_data/skippy_non_concat_SILVA_138_rpl.biom")
  )
}

skippy <- loadFromBiom(biom_file)
# add phylogenetic tree made with NGT2 (clustal omega)
tree_file <- here::here(
  "data/raw_data/skippy_non_concat_tree_files/all_otus.tree"
)
tree <- ape::read.tree(tree_file)
# tree can only be added to tse but not se object
skippy <- as(skippy, "TreeSummarizedExperiment")
rowTree(skippy) <- tree
# use practical names for the ranks instead
colnames(rowData(skippy)) <- c(
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
)

# add sample ids as column
colData(skippy)[["sample_id"]] <- rownames(colData(skippy))
# now we can add our sample data we created above to the tse
colData(skippy) <- colData(skippy) %>%
  as.data.frame() %>%
  left_join(sample_data, by = "sample_id") %>%
  mutate(neg_control = as.factor(ifelse(is.na(week), 1, 0))) %>%
  column_to_rownames("sample_id") %>%
  mutate(sample_id = rownames(colData(skippy))) %>%
  DataFrame()

# we will now store MOCK samples in a separate tse
mock_samples <- filter(skippy, str_detect(sample_id, "MOCK"))
mock_samples <- transformSamples(mock_samples, method = "clr", pseudocount = 1)
# check correlations between mock samples
assay(mock_samples, "clr") %>% cor()


# also check they cluster together closely when using ordination
skippy <- mutate(skippy,
  mock = as.factor(ifelse(str_detect(sample_id, "MOCK"), 1, 0))
)
colData(skippy)[["library_size"]] <- colSums(assay(skippy))
# perform NMDS coordination method
skippy <- runNMDS(
  skippy,
  FUN = vegan::vegdist,
  name = "NMDS"
)
# plot results of a 2-component NMDS on tse,
plotReducedDim(skippy, "NMDS", colour_by = "mock")
plotReducedDim(skippy, "NMDS", colour_by = "library_size")
# MOCK samples do cluster together but particularly 1 outlier is there
# now drop the MOCK samples
tse <- filter(skippy, !str_detect(sample_id, "MOCK"))
# fix names of taxa 
rowData(tse) <- rowData(tse) %>%
  as_tibble() %>%
  mutate(across(everything(), function(x) str_remove(x, ".*[kpcofgs]__"))) %>%
  DataFrame()

# which are the negative control samples?
colData(tse) %>% as.data.frame() %>% filter(!week %in% c(2, 5, 52))
# the low ls are negative controls as we can see in above table
colData(tse) %>%
as.data.frame() %>%
  arrange(library_size)
# exclude negative controls and create unique sample ids indicative of week 
tse <- filter(tse, neg_control == 0) %>%
  mutate(
    sample_ids_new = glue(
      "{skippy_id}_{ifelse(week == 2, 1, 
        ifelse(week == 5, 2, ifelse(week == 52, 3, 'CHECK')))}{
          ifelse(duplicated, '_2', '')}"))

# we retain 331 samples at this step
# lets change row and column names to better descriptions
colnames(tse) <- colData(tse)[, "sample_ids_new"]
save(tse, file = here::here("data/tse.Rds"))



# store a file that lists all samples we have
colData(tse) %>%
  as.data.frame() %>%
  select(sample_ids_new, skippy_id, week, wur_id) %>%
  arrange(sample_ids_new) %>%
  write_excel_csv(here::here("data/processed_samples.csv"))

# test if we need to rarefy or not: do the methods in ad_tjazi and mia 
# produce different results?
tse <- estimateDiversity(
  tse,
  assay_name = "counts",
  index = c("shannon", "faith", "inverse_simpson"),
  name = c("shannon", "faith", "inverse_simpson")
)

ad_tjazi <- Tjazi::get_asymptotic_alpha(assay(tse)) %>%
  rownames_to_column("sample_id")

colData(tse) %>%
  as.data.frame() %>%
  select(
    library_size,
    sample_id = sample_ids_new,
    shannon,
    inverse_simpson,
    faith) %>%
  full_join(ad_tjazi, by = "sample_id") %>%
  select_if(is_numeric) %>%
  cor()

  # colData(tse) %>% as.data.frame() %>%
  #   select(library_size, sample_id = sample_ids_new, shannon, inverse_simpson, faith) %>%
  #   full_join(ad_tjazi, by = "sample_id") %>%
  #   select_if(is_numeric) %>%
  #   select(-library_size, -Chao1) %>%
  #   pivot_longer(everything(), names_to = "index", values_to = "diversity") %>%
  #   ggplot(aes(index, diversity)) +
  #     geom_boxplot() +
  #     geom_jitter()

# no correlation with library size + perfect correlation of shannon -->
# we can pick any method 


# here I check if if end up with all the samples as per documented by Kelly:
samples <- tibble(sample_id = colnames(tse)) %>%
  mutate(
    week = str_match(sample_id, "_(\\d+)$")[, 2],
    id = str_match(sample_id, "^(\\d\\d\\d)")[, 2]
  ) %>%
  arrange(id, week)
doc <- read_excel(here::here("data/StoolSamples_20162017_SKIPPY.xlsx"))
doc <- select(doc,
  id = ID,
  w2 = "Stoolsample week 2",
  w5 = "Stoolsample week 5",
  w52 = "Stoolsample week 52")  %>%
  pivot_longer(matches("w\\d"), names_to = "week", values_to = "planned") %>%
  mutate(
    week = str_extract(week, "\\d+"),
    week = ifelse(
      week == 52, 3, ifelse(
        week == 5, 2, ifelse(
          week == 2, 1, NA)))
  )
samples$id <- as.double(samples$id)
samples$week <- as.double(samples$week)
full_join(samples, doc, by = c("id", "week")) %>%
 filter(is.na(sample_id), planned == 1) %>%
 arrange(id, week)
missing_ids <- full_join(samples, doc, by = c("id", "week")) %>%
  filter(is.na(sample_id), planned == 1) %>%
  arrange(id, week) %>%
  .$id
missing_ids %in% samples$id
filter(samples, id %in% missing_ids) %>%
  arrange(id, week)


# there is a discrepancy in the samples that should be there and that
# are there. I will search through all map files once more but this time
# without filtering to make sure that I didnt filter them out before

raw_path <- here::here("data/map_files")
map_files <- list.files(raw_path, pattern = "2021_\\d{4}.*.xlsx$")
map_files %>% length()
lib_nums <- str_extract(map_files, "\\d\\d\\d\\d_00\\d\\d")
library_number_merged <- c(1:6)
map_file_merged <- map2_dfr(
  map_files,
  library_number_merged,
  function(filename, lnm) {
  # extract library number
  lib_num <- str_extract(filename, "\\d\\d\\d\\d_00\\d\\d")
  # read provided map file
  xl <- read_excel(glue("{raw_path}/{filename}"))
  map_file <- xl %>%
    mutate(
      LibraryNumber = lnm,
      ProjectName = ifelse(is.na(ProjectName), "unspecified", ProjectName)
    ) %>%
    select(
      "sample_id" = internal_sample_id,
      id = Seq_ID,
      ProjectName
    ) %>%
      filter(!is.na(`sample_id`))

    return(map_file)}) %>%
  mutate(id2 = str_match(id, ".*(\\d\\d\\d).*")[, 2])
filter(map_file_merged, id2 %in% missing_ids)

# according to documentation all samples except 226 week 52 are correctly
# missing as they have never been sent to WUR

n_total <- length(colnames(tse))
duplicate_ids <- colnames(tse)[map_int(colnames(tse), ~str_length(.x)) > 5]
n_dupl <- length(colnames(tse)[map_int(colnames(tse), ~str_length(.x)) > 5])
n_total - n_dupl
# how many unique samples per time point?
colData(tse) %>% 
  as.data.frame() %>%
  filter(!sample_ids_new %in% duplicate_ids) %>%
  count(week)


# after Kelly's latest prereg email I check now which of the ITT individuals 
# are in the data

library(foreign)
sav <- read.spss(here::here("data/kelly_documents/data_itt_pp_dr.sav"), to.data.frame = TRUE) %>% select(id = ID, itt = ITT) %>%
  filter(itt != "Excl") %>%
  mutate(id = as.character(id))
  
df <- colData(tse) %>%
  as.data.frame() %>%
  select(wur_id, id = skippy_id, week, duplicated) %>%
  full_join(sav, by = "id")
ids <- filter(df, is.na(wur_id)) %>%
  .$id
cat(glue("{ids},"))
filter(df, id %in% ids)

df %>% filter(id %in% sav$id, !duplicated) %>% count(week)


# after all has been checked, I will now double check that only the correct
# ids remain in the tse + that I add any necessary variables for the project
load(file = here::here("data/tse.Rds"))

colData(tse) %>%
  as.data.frame() %>%
  filter(skippy_id %in% sav$id, !duplicated) %>%
  count(week)
tse <- filter(tse, skippy_id %in% sav$id, !duplicated)

# apparently, there were several date files and this file has more complete data
# than the mdata below for the birthdate variable:
xl <- read_excel(here::here("data/Birthdates_SKIPPY.xlsx")) %>%
  mutate(birthdate = dmy(GeboortedatumBaby)) %>%
  select(id = ID, birthdate) %>%
  distinct()
# fix wrong year entry for id 272 (typo)
xl[xl$id == 272, "birthdate"] <- ymd("2017-04-24")


# import metadata
mdata <- read_csv2(here::here("data/skippy_stool_data_cleaned.csv")) %>%
  mutate_all(function(x) ifelse(x == 99999, NA, ifelse(x == 88888, NA, x))) %>%
  filter(id %in% sav$id) %>%
  left_join(xl, by = "id") %>%
  # dates must be fixed in order to calculate age
  mutate(across(contains("dat_"), function(x) {
    x_new <- str_replace(x, "okt", "october")
    x_new <- str_replace(x_new, "mrt", "march")
    x_new <- str_replace(x_new, "mei", "may")
    x_new <- str_replace(x_new, "15-jan-17/22-jan-17", "18/jan/17")
    x_new <- dmy(x_new)
    return(x_new)
  }))

# again a typo for dat_week2. Easier to fix manually in base R
# dat_week2 = ifelse(id == 300, dmy("25-06-2017"), dmy(dat_week2)),
mdata[mdata$id == 300, "dat_week2"] <- ymd("2017-06-25")

mdata <- mutate(mdata,
  datbirth_infant = dmy(datbirth_infant),
  # now that dates are fine we can calculate age
  checkold_week2 = as.numeric(dat_week2 - datbirth_infant),
  checkold_week5 = as.numeric(dat_week5 - datbirth_infant),
  checkold_week52 = as.numeric(dat_1year - datbirth_infant),
  age_week2 = as.numeric(dat_week2 - birthdate),
  age_week5 = as.numeric(dat_week5 - birthdate),
  age_week52 = as.numeric(dat_1year - birthdate),
  id = as.character(id)) %>%
  select(
    skippy_id = id, condition, csection, birthweight,
    siblings, sex, bfexcl, contains("bf"), everything())

filter(
  mdata,
  age_week2 != checkold_week2 |
    age_week5 != checkold_week5 |
      age_week52 != checkold_week52) %>%
  select(skippy_id, birthdate, datbirth_infant)

# now we add all the meta variables to the tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  rownames_to_column("sid") %>%
  left_join(mdata, by = "skippy_id") %>%
  select(-wur_id, -duplicated, -neg_control, -mock, -sample_id) %>%
  select(everything(), sample_id = sample_ids_new) %>%
  column_to_rownames("sid") %>%
  DataFrame()


# this resulting file is our starting point for the analyses
# i will also save it with the mdata and sav as we need that for multiple
# imputation
save(tse, sav, mdata, file = here::here("data/data.Rds"))
select(
  mdata, 
  birthdate, 
  datbirth_infant, 
  skippy_id, 
  age_week2, 
  age_week5, 
  age_week52
)