set.seed(1)
library(mia)
library(LinDA)
library(tidyverse)
library(tidySummarizedExperiment)
library(glue)


###############################################################################
#########################           1. ITT       ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data_pc.Rds"))
load(file = here::here("data/data_imp.Rds"))

# for analyses we apply prevalence fitlering
tse_pc <- subsetByPrevalentTaxa(tse_pc, detection = 0, prevalence = 0.1)

####################### 1.1 Complete Case Analysis ############################


# model includes random intercepts and all samples 
asv_tab <- as.data.frame(assay(tse_pc, "counts"))
vars <- c("age_s", "condition", "siblings", "constipation", "diarrhea", "skippy_id")
meta <- colData(tse_pc) %>% as.data.frame() %>%
  select(all_of(vars))
linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + (1|skippy_id)', alpha = 0.1,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)


filter(linda_obj$output$condition1, reject)
filter(linda_obj$output$siblings1, reject)
 

# model includes random intercepts and excludes 1 year samples
tse_pc_infancy <- filter(tse_pc, week != 52)
asv_tab <- as.data.frame(assay(tse_pc_infancy, "counts"))
meta <- colData(tse_pc_infancy) %>% as.data.frame() %>%
  select(all_of(vars))
dim(na.omit(meta))
linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + (1|skippy_id)', alpha = 0.1,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
# linda.plot(linda_obj, c('condition'), 
#            titles = c('Condition: n v.s. y'), alpha = 0.1, lfc.cut = 1,
#            legend = TRUE, directory = NULL, width = 11, height = 8)

filter(linda_obj$output$condition1, reject)
filter(linda_obj$output$siblings1, reject)



# model includes 1 years samples only 
tse_pc_year1 <- filter(tse_pc, week == 52)
asv_tab <- as.data.frame(assay(tse_pc_year1, "counts"))
meta <- colData(tse_pc_year1) %>% as.data.frame() %>%
  select(all_of(vars))
linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings', alpha = 0.1,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
res <- filter(linda_obj$output$condition1, reject) %>%
  rownames_to_column("taxid")
rd <- rowData(tse_pc_year1) %>% 
  as.data.frame() %>% 
  rownames_to_column("taxid") %>%
  filter(taxid %in% res$taxid)
res <- left_join(res, rd, by = "taxid")
filter(linda_obj$output$condition1, reject)
filter(linda_obj$output$siblings1, reject)



######################## 1. 2Multiple imputation  #############################


load(here::here("data/data_pc.Rds"))
load(file = here::here("data/data_imp.Rds"))
# for analyses we apply prevalence fitlering
tse_pc <- subsetByPrevalentTaxa(tse_pc, detection = 0, prevalence = 0.1)

models_imp <- map2(implist, 1:length(implist), function(dimp, imp) {
  tse_pc_map <- tse_pc 
  fvars <- c("constipation", "siblings", "diarrhea", "condition")
  # add metadata to tse_pc
  colData(tse_pc_map) <- colData(tse_pc_map) %>%
    as.data.frame() %>%
    select(-all_of(fvars), -contains("age")) %>%
    left_join(select(dimp, age, sample_id, constipation, diarrhea, siblings, condition), by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(across(all_of(fvars), function(x) as.factor(x))) %>%
    DataFrame()
  colData(tse_pc_map)$age <- colData(tse_pc_map)$age + as.numeric(colData(tse_pc_map)$week) * 7
  colData(tse_pc_map)$age_s <- scale(colData(tse_pc_map)$age)[, 1]

  # all samples
  asv_tab <- as.data.frame(assay(tse_pc_map, "counts"))
  meta <- colData(tse_pc_map) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars))
    asv_tab <- asv_tab[, rownames(meta)]

  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + (1|skippy_id)', alpha = 0.4,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
  
  res_all <- linda_obj$output$condition1
  res_all2 <- linda_obj$output[["condition1:age_s"]]
  
  # model includes random intercepts and excludes 1 year samples
  tse_pc_infancy <- filter(tse_pc_map, week != 52)
  asv_tab <- as.data.frame(assay(tse_pc_infancy, "counts"))
  meta <- colData(tse_pc_infancy) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars))
  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + (1|skippy_id)', alpha = 0.1,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)


  res_infancy <- linda_obj$output$condition1
  res_infancy2 <- linda_obj$output[["condition1:age_s"]]
  
  # model includes 1 years samples only 
  tse_pc_year1 <- filter(tse_pc_map, week == 52)
  asv_tab <- as.data.frame(assay(tse_pc_year1, "counts"))
  meta <- colData(tse_pc_year1) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars))
  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings', alpha = 0.1,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)

  
  res_year1 <- linda_obj$output$condition1
  res_year12 <- linda_obj$output[["condition1:age_s"]]

  list(res_all, res_infancy, res_year1, res_all2, res_infancy2, res_year12)
})

tables_linda_itt <- map(models_imp, function(x) {
  x[[1]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3)))
})
tables_linda_itt

tables_linda_itt2 <- map(models_imp, function(x) {
  x[[4]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3)))
})
tables_linda_itt2





tables_linda_itt_infancy <- map(models_imp, function(x) {
  x[[2]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3)))
})
 

tables_linda_itt_year1 <- map(models_imp, function(x) {
  x[[3]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3)))
})

tables_linda_itt
tables_linda_itt2
tables_linda_itt_infancy
tables_linda_itt_year1

linda_identified <- map2_dfr(
  list(
    tables_linda_itt,
    tables_linda_itt2,
    tables_linda_itt_infancy,
    tables_linda_itt_year1
    ),
  c(
    "main effect",
    "interaction term",
    "early",
    "late"
  ),
  function(table, term) {
    filter(table[[1]], padj <= 0.2)  %>%
      mutate(term = term)
  }) %>%
  arrange(desc(abs(log2FoldChange)), padj)

# # delete duplicate row
# linda_identified <- linda_identified[-4, ]
linda_identified


save(
  tables_linda_itt, 
  tables_linda_itt2, 
  linda_identified, 
  file = here::here("data/tables_linda_itt_pc.Rds")
)

save(
  tables_linda_itt_infancy, 
  file = here::here("data/tables_linda_itt_infancy_pc.Rds")
)

save(
  tables_linda_itt_year1, 
  file = here::here("data/tables_linda_itt_year1_pc.Rds")
) 








models_imp <- map2(implist, 1:length(implist), function(dimp, imp) {
  tse_pc_map <- tse_pc 
  fvars <- c("constipation", "siblings", "diarrhea", "condition")
  # add metadata to tse_pc
  colData(tse_pc_map) <- colData(tse_pc_map) %>%
    as.data.frame() %>%
    select(-all_of(fvars), -contains("age"), -bfexcl) %>%
    left_join(select(dimp, age, sample_id, constipation, diarrhea, siblings, bfexcl, condition), by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(across(all_of(fvars), function(x) as.factor(x))) %>%
    DataFrame()
  colData(tse_pc_map)$age <- colData(tse_pc_map)$age + as.numeric(colData(tse_pc_map)$week) * 7
  colData(tse_pc_map)$age_s <- scale(colData(tse_pc_map)$age)[, 1]
  
  
  # model includes 1 years samples only 
  tse_pc_year1 <- filter(tse_pc_map, week == 52)
  asv_tab <- as.data.frame(assay(tse_pc_year1, "counts"))
  meta <- colData(tse_pc_year1) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars), bfexcl)
  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + bfexcl', alpha = 0.1,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
  res <- filter(linda_obj$output$condition1, reject) %>%
    rownames_to_column("taxid")
  if (dim(res)[1] >= 1) {
    rd <- rowData(tse_pc_year1) %>% 
     as.data.frame() %>% 
     rownames_to_column("taxid") %>%
     filter(taxid %in% res$taxid)
    res_year1 <- left_join(res, rd, by = "taxid") %>% 
     mutate(time = "year1")
  }
  
  out <- list()
  if(exists("res_year1")) {
    out[[1]] <- res_year1 
  } 
  out
})


models_imp












































###############################################################################
#########################           2. PP        ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data_pc.Rds"))
load(file = here::here("data/data_imp.Rds"))

tse <- subsetByPrevalentTaxa(tse_pc, detection = 0, prevalence = 0.1)


# obtain ids that were selected for PP analyses 
pp_indicator <- foreign::read.spss(here::here("data/raw_data/kelly141022/Data_ITT_PP_ExploratoryDRselections.sav"), to.data.frame = TRUE)
pp_indicator <- select(pp_indicator, skippy_id = ID, pp = PP)
# add pp info to existing data 
if (!"pp" %in% colnames(d)) {
  d <- left_join(d, pp_indicator, by = "skippy_id")
}
# 60 that are in PP and condition 0; 18 that are condition 1 and pp. Fits...
d_pp <- filter(d, pp == 1)
implist_pp <- map(implist, function(dimp) {
  dimp_new <- left_join(dimp, pp_indicator, by = "skippy_id") %>%
                filter(pp == 1)
  dimp_new
})


fvars <- c("constipation", "siblings", "diarrhea", "condition", "csection", "sex")
# add metadata to tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  select(-ges_age, -birthweight, -edlevel, -csection, -all_of(fvars), -age, -age_s) %>%
  left_join(select(d_pp, age, sample_id, constipation, condition, diarrhea, siblings, pp, csection, sex, ges_age, edlevel, birthweight), by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  mutate(
    across(all_of(fvars), function(x) as.factor(x)),
    ges_age_s = scale(ges_age)[, 1],
    birthweight_s = scale(birthweight)[, 1],
    edlevel_s = scale(edlevel)[, 1]
    ) %>%
  DataFrame()
colData(tse)$age <- colData(tse)$age + as.numeric(colData(tse)$week) * 7
colData(tse)$age_s <- scale(colData(tse)$age)[, 1]
tse <- filter(tse, pp == 1)

####################### 2.1 Complete Case Analysis ############################


# model includes random intercepts and all samples 
asv_tab <- as.data.frame(assay(tse, "counts"))
meta <- colData(tse) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars), ges_age_s, birthweight_s, edlevel_s)
asv_tab <- asv_tab[, rownames(meta)]
linda_obj <- linda(
  asv_tab, 
  meta, 
  formula = '~condition * age_s + siblings + csection + sex + ges_age_s + birthweight_s + edlevel_s + (1|skippy_id)', alpha = 0.4,
  prev.cut = 0.1, 
  lib.cut = 1000, 
  winsor.quan = 0.97
)


filter(linda_obj$output$condition1, reject)
filter(linda_obj$output$siblings1, reject)


# model includes random intercepts and excludes 1 year samples
tse_infancy <- filter(tse, week != 52)
asv_tab <- as.data.frame(assay(tse_infancy, "counts"))
meta <- colData(tse_infancy) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars), ges_age_s, edlevel_s, birthweight_s)
dim(na.omit(meta))
asv_tab <- asv_tab[, rownames(meta)]

linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + edlevel_s + ges_age_s + birthweight_s + sex + csection + (1|skippy_id)', alpha = 0.4,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
# linda.plot(linda_obj, c('condition'), 
#            titles = c('Condition: n v.s. y'), alpha = 0.4, lfc.cut = 1,
#            legend = TRUE, directory = NULL, width = 11, height = 8)

filter(linda_obj$output$condition1, reject)
filter(linda_obj$output$siblings1, reject)



# model includes 1 years samples only 
tse_year1 <- filter(tse, week == 52)
asv_tab <- as.data.frame(assay(tse_year1, "counts"))
meta <- colData(tse_year1) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars), ges_age_s, birthweight_s, edlevel_s)
asv_tab <- asv_tab[, rownames(meta)]
linda_obj <- linda(
  asv_tab, 
  meta, 
  formula = '~condition * age_s + siblings + csection + sex + ges_age_s + birthweight_s + edlevel_s', alpha = 0.4,
  prev.cut = 0.1, 
  lib.cut = 1000, 
  winsor.quan = 0.97
)
res <- filter(linda_obj$output$condition1, reject) %>%
  rownames_to_column("taxid")
rd <- rowData(tse_year1) %>% 
  as.data.frame() %>% 
  rownames_to_column("taxid") %>%
  filter(taxid %in% res$taxid)
res <- left_join(res, rd, by = "taxid")
filter(linda_obj$output$siblings1, reject)
filter(linda_obj$output$condition1, reject)





######################## 2.2 Multiple imputation  #############################

load(here::here("data/data_pc.Rds"))
load(file = here::here("data/data_imp.Rds"))
# for analyses we apply prevalence fitlering
tse <- subsetByPrevalentTaxa(tse_pc, detection = 0, prevalence = 0.1)


models_imp <- map2(implist_pp, 1:length(implist), function(dimp, imp) {
  tse_map <- tse 
  fvars <- c("constipation", "siblings", "diarrhea", "condition", "csection", "sex")
  # add metadata to tse
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(-all_of(fvars), -age, -age_s, -ges_age, -birthweight, -edlevel, -csection) %>%
    left_join(select(dimp, age, sample_id, condition, constipation, diarrhea, siblings, pp, csection, sex, ges_age_s, edlevel_s, birthweight_s), by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(across(all_of(fvars), function(x) as.factor(x))) %>%
    DataFrame()
  colData(tse_map)$age <- colData(tse_map)$age + as.numeric(colData(tse_map)$week) * 7
  colData(tse_map)$age_s <- scale(colData(tse_map)$age)[, 1]
  tse_map <- filter(tse_map, pp == 1)
  
  # all samples
  asv_tab <- as.data.frame(assay(tse_map, "counts"))
  meta <- colData(tse_map) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars), edlevel_s, birthweight_s, ges_age_s)
  asv_tab <- asv_tab[, rownames(meta)]
  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + csection + sex + ges_age_s + birthweight_s + edlevel_s + (1|skippy_id)', alpha = 0.4,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
  

  res_all <- linda_obj$output$condition1
  res_all2 <- linda_obj$output[["condition1:age_s"]]
  
  # model includes random intercepts and excludes 1 year samples
  tse_infancy <- filter(tse_map, week != 52)
  asv_tab <- as.data.frame(assay(tse_infancy, "counts"))
  meta <- colData(tse_infancy) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars), edlevel_s, birthweight_s, ges_age_s)
  asv_tab <- asv_tab[, rownames(meta)]
  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + csection + sex + ges_age_s + birthweight_s + edlevel_s + (1|skippy_id)', alpha = 0.4,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
  
  res_infancy <- linda_obj$output$condition1
  res_infancy2 <- linda_obj$output[["condition1:age_s"]]

  
  
  # model includes 1 years samples only 
  tse_year1 <- filter(tse_map, week == 52)
  asv_tab <- as.data.frame(assay(tse_year1, "counts"))
  meta <- colData(tse_year1) %>% as.data.frame() %>%
    select(skippy_id, age_s, all_of(fvars), edlevel_s, birthweight_s, ges_age_s)
  asv_tab <- asv_tab[, rownames(meta)]
  linda_obj <- linda(asv_tab, meta, formula = '~condition * age_s + siblings + csection + sex + ges_age_s + birthweight_s + edlevel_s', alpha = 0.4,
                     prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)

  res_year1 <- linda_obj$output$condition1
  res_year12 <- linda_obj$output[["condition1:age_s"]]

  list(res_all, res_infancy, res_year1, res_all2, res_infancy2, res_year12)
})

tables_linda_pp <- map(models_imp, function(x) {
  x[[1]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3))) 
})
tables_linda_pp

save(tables_linda_pp, file = here::here("data/tables_linda_pp_pc.Rds"))



tables_linda_pp_infancy <- map(models_imp, function(x) {
  x[[2]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3))) 
})
tables_linda_pp_infancy

save(tables_linda_pp_infancy, file = here::here("data/tables_linda_pp_infancy_pc.Rds"))



tables_linda_pp_year1 <- map(models_imp, function(x) {
  x[[3]] %>%
    rownames_to_column("taxon") %>%
    select(taxon, log2FoldChange, lfcSE, pvalue, padj) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(across(where(is.numeric), function(x) round(x, 3))) 
})
tables_linda_pp_year1

save(tables_linda_pp_year1, file = here::here("data/tables_linda_pp_year1_pc.Rds"))

