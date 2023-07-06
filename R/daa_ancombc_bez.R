set.seed(1)
library(mia)
library(ANCOMBC)
library(tidyverse)
library(tidySummarizedExperiment)
library(glue)




###############################################################################
#########################           1. ITT       ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

# for analyses we apply prevalence fitlering and analyze at genus level
tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)
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

####################### 1.1 Complete Case Analysis ############################



# model includes random intercepts and all samples 
if (!file.exists(here::here("data/daa_cc_mlm.Rds"))) {
  output <- ancombc2(
    data = tse, 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s",
    rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl()
    )
    save(output, file = here::here("data/daa_cc_mlm.Rds"))
 } else {
   load(here::here("data/daa_cc_mlm.Rds"))
}
   
# structural zeros 
tab_zero = output$zero_ind
sum(tab_zero[, 2]) + sum(tab_zero[, 3])
# sensitivity scores 
tab_sens = output$pseudo_sens_tab
head(tab_sens)
pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
  filter(value >5)

res_prim <- output$res
colnames(res_prim)
effects <- c("age_s", "condition1", "condition1:age_s")
map(effects, function(effect) {
  select(res_prim, 
    taxon, 
    lfc = glue::glue("lfc_{effect}"),
    se = glue::glue("se_{effect}"),
    indicator = glue::glue("diff_{effect}")) %>%
  filter(indicator) %>%
  mutate(effect = effect)
})

res_prim %>% select(taxon, contains("condition1")) %>%
  arrange(q_condition1, p_condition1) %>%
  head(15)


# model includes random intercepts and excludes 1 year samples
if (!file.exists(here::here("data/daa_cc_mlm_infancy.Rds"))) {
  output <- ancombc2(
    data = filter(tse, week != 52), 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s",
    rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl()
    )
    save(output, file = here::here("data/daa_cc_mlm_infancy.Rds"))
 } else {
   load(here::here("data/daa_cc_mlm_infancy.Rds"))
}
   
# structural zeros 
tab_zero = output$zero_ind
sum(tab_zero[, 2]) + sum(tab_zero[, 3])
# sensitivity scores 
tab_sens = output$pseudo_sens_tab
head(tab_sens)
pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
 filter(value >5)

res_prim <- output$res
colnames(res_prim)
effects <- c("age_s", "condition1", "condition1:age_s")
map(effects, function(effect) {
 select(res_prim, 
   taxon, 
   lfc = glue::glue("lfc_{effect}"),
   se = glue::glue("se_{effect}"),
   indicator = glue::glue("diff_{effect}")) %>%
 filter(indicator) %>%
 mutate(effect = effect)
})

res_prim %>% select(taxon, contains("condition1")) %>%
  arrange(p_condition1) %>%
  head(15)

# model includes 1 years samples only 
if (!file.exists(here::here("data/daa_cc_year1.Rds"))) {
  output <- ancombc2(
    data = filter(tse, week == 52), 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s",
    # rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl()
    )
    save(output, file = here::here("data/daa_cc_year1.Rds"))
 } else {
   load(here::here("data/daa_cc_year1.Rds"))
 }
   
# structural zeros 
tab_zero = output$zero_ind
sum(tab_zero[, 2]) + sum(tab_zero[, 3])
 # sensitivity scores 
tab_sens = output$pseudo_sens_tab
head(tab_sens)
pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
   filter(value >5)
 
res_prim <- output$res
colnames(res_prim)
effects <- c("age_s", "condition1", "condition1:age_s")
map(effects, function(effect) {
   select(res_prim, 
     taxon, 
     lfc = glue::glue("lfc_{effect}"),
     se = glue::glue("se_{effect}"),
     indicator = glue::glue("diff_{effect}")) %>%
   filter(indicator) %>%
   mutate(effect = effect)
})





######################## 1.2 Multiple imputation  #############################



load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
# for analyses we apply prevalence fitlering
tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)

models_imp <- map2(implist, length(implist), function(dimp, imp) {
  tse_map <- tse
  fvars <- c("siblings", "condition")
  # add metadata to tse
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(-siblings) %>%
    left_join(select(dimp, age, sample_id, siblings), by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(across(all_of(fvars), function(x) as.factor(x))) %>%
    DataFrame()
  colData(tse_map)$age <- colData(tse_map)$age + as.numeric(colData(tse_map)$week) * 7
  colData(tse_map)$age_s <- scale(colData(tse_map)$age)[, 1]


  # model includes random intercepts and excludes 1 year samples
  if (!file.exists(here::here(glue("data/daa_cc_mlm_all_imp{imp}2_bez.Rds")))) {
    output <- ancombc2(
      data = filter(tse_map), 
      assay_name = "counts", 
      tax_level = "genus",
      fix_formula = "condition * age_s",
      rand_formula = "(1 | skippy_id)",
      p_adj_method = "holm", 
      pseudo = 0,
      pseudo_sens = TRUE,
      prv_cut = 0.10,
      lib_cut = 1000,
      s0_perc = 0.05,
      group = "condition",
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.2,
      n_cl = 2, 
      verbose = TRUE,
      global = FALSE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
        node = list(2), 
        solver = "ECOS",B = 100)
      )
      save(output, file = here::here(glue("data/daa_cc_mlm_all_imp{imp}2_bez.Rds")))
   } else {
     load(here::here(glue("data/daa_cc_mlm_all_imp{imp}2_bez.Rds")))
  }
     
  # structural zeros 
  tab_zero_all = output$zero_ind
  # sensitivity scores 
  tab_sens_all = output$pseudo_sens_tab
  sens_scores_all <- pivot_longer(tab_sens_all, -taxon, names_to = "contrast", values_to = "value") %>%
     filter(value >5)
  res_prim_all <- output$res
  
  effects_all <- map(effects, function(effect) {
     select(res_prim_all, 
       taxon, 
       lfc = glue::glue("lfc_{effect}"),
       se = glue::glue("se_{effect}"),
       indicator = glue::glue("diff_{effect}")) %>%
     filter(indicator) %>%
     mutate(effect = effect)
   })


  
  # model includes random intercepts and excludes 1 year samples
  if (!file.exists(here::here(glue("data/daa_cc_mlm_infancy_imp{imp}_bez.Rds")))) {
    output <- ancombc2(
      data = filter(tse_map, week != 52), 
      assay_name = "counts", 
      tax_level = "genus",
      fix_formula = "condition * age_s",
      rand_formula = "(1 | skippy_id)",
      p_adj_method = "holm", 
      pseudo = 0, 
      pseudo_sens = TRUE,
      prv_cut = 0.10, 
      lib_cut = 1000, 
      s0_perc = 0.05,
      group = "condition", 
      struc_zero = TRUE, 
      neg_lb = TRUE,
      alpha = 0.2, 
      n_cl = 2, 
      verbose = TRUE,
      global = FALSE, 
      pairwise = FALSE, 
      dunnet = FALSE, 
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
        node = list(2), 
        solver = "ECOS",B = 100)
      )
      save(output, file = here::here(glue("data/daa_cc_mlm_infancy_imp{imp}_bez.Rds")))
   } else {
     load(here::here(glue("data/daa_cc_mlm_infancy_imp{imp}_bez.Rds")))
  }
     
  # structural zeros 
  tab_zero_infancy = output$zero_ind
  # sensitivity scores 
  tab_sens_infancy = output$pseudo_sens_tab
  sens_scores_infancy <- pivot_longer(tab_sens_infancy, -taxon, names_to = "contrast", values_to = "value") %>%
     filter(value >5)
  res_prim_infancy <- output$res
  
  effects_infancy <- map(effects, function(effect) {
     select(res_prim_infancy, 
       taxon, 
       lfc = glue::glue("lfc_{effect}"),
       se = glue::glue("se_{effect}"),
       indicator = glue::glue("diff_{effect}")) %>%
     filter(indicator) %>%
     mutate(effect = effect)
   })
  
  
  # model includes 1 years samples only 
  if (!file.exists(here::here(glue("data/daa_cc_year1_imp{imp}_bez.Rds")))) {
    output <- ancombc2(
      data = filter(tse_map, week == 52), 
      assay_name = "counts", 
      tax_level = "genus",
      fix_formula = "condition * age_s",
      # rand_formula = "(1 | skippy_id)",
      p_adj_method = "holm", 
      pseudo = 0, 
      pseudo_sens = TRUE,
      prv_cut = 0.10, 
      lib_cut = 1000, 
      s0_perc = 0.05,
      group = "condition", 
      struc_zero = TRUE, 
      neg_lb = TRUE,
      alpha = 0.2, 
      n_cl = 2, 
      verbose = TRUE,
      global = FALSE, 
      pairwise = FALSE, 
      dunnet = FALSE, 
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
        node = list(2), 
        solver = "ECOS",B = 100)
      )
      save(output, file = here::here(glue("data/daa_cc_year1_imp{imp}_bez.Rds")))
   } else {
     load(here::here(glue("data/daa_cc_year1_imp{imp}_bez.Rds")))
   }
     
  # structural zeros 
  tab_zero_year1 = output$zero_ind
  # sensitivity scores 
  tab_sens_year1 = output$pseudo_sens_tab
  sens_scores_year1 <- pivot_longer(tab_sens_year1, -taxon, names_to = "contrast", values_to = "value") %>%
     filter(value >5)
  res_prim_year1 <- output$res
  effects <- c("age_s", "condition1")
  effects_year1 <- map(effects, function(effect) {
     select(res_prim_year1, 
       taxon, 
       lfc = glue::glue("lfc_{effect}"),
       se = glue::glue("se_{effect}"),
       indicator = glue::glue("diff_{effect}")) %>%
     filter(indicator) %>%
     mutate(effect = effect)
   })

  list(
      all = list(
        tab_zero_all,
        sens_scores_all,
        res_prim_all,
        effects_all
    ),
    infancy = list(
      tab_zero_infancy,
      sens_scores_infancy,
      res_prim_infancy,
      effects_infancy
    ),
    year1 = list(
      tab_zero_year1,
      sens_scores_year1,
      res_prim_year1,
      effects_year1
    )
  )
})

# you must have effects defined before running the map command
models_imp[[1]]$all$effects_all
models_imp[[1]]$infancy$effects_infancy

ancom_tables_itt <- models_imp[[1]]$all[[3]] %>% 
  select(taxon, lfc_condition1, se_condition1, p_condition1, q_condition1) %>%
  arrange(q_condition1, desc(abs(lfc_condition1))) %>%
  mutate(across(where(is.numeric), round, 3))
colnames(ancom_tables_itt) <- str_remove(colnames(ancom_tables_itt), "_condition1") 

ancom_tables_itt2 <- models_imp[[1]]$all[[3]] %>% 
  select(taxon, "lfc_condition1:age_s", "se_condition1:age_s", "p_condition1:age_s", "q_condition1:age_s") %>%
  arrange(`q_condition1:age_s`, desc(abs(`lfc_condition1:age_s`))) %>%
  mutate(across(where(is.numeric), round, 3))

colnames(ancom_tables_itt2) <- str_remove(colnames(ancom_tables_itt2), "_condition1:age_s") 

save(ancom_tables_itt, file = here::here("data/ancom_tables_itt.Rds"))

ancom_tables_itt
ancom_tables_itt2
# ancom did not identify any as significant but effect sizes are in line







t1 <- select(
  models_imp[[1]]$all[[3]],
  taxon,
  lfc_age_s) %>%
  pivot_longer(
    contains("lfc"), 
    names_to = "variable", 
    values_to = "lfc", 
    names_prefix = "lfc_")
  

t2 <- select(
  models_imp[[1]]$all[[3]],
  taxon,
  p_age_s
  ) %>%
  pivot_longer(
    contains("p_"), 
    names_to = "variable", 
    values_to = "p", 
    names_prefix = "p_")
t2
t3 <- select(
  models_imp[[1]]$all[[3]],
  taxon,
  q_age_s) %>%
  pivot_longer(
    contains("q_"), 
    names_to = "variable", 
    values_to = "q", 
    names_prefix = "q_")
t3

ancombc_remaining <- full_join(t1, t2, by = c("taxon", "variable")) %>%
  full_join(t3, by = c("taxon", "variable")) %>%
  filter(p <= 0.05) %>%
  mutate(
    across(where(is.numeric), round, 3),
    taxon = str_remove(taxon, "genus:")
    ) %>%
  arrange(variable, q) 
save(ancombc_remaining, file = here::here("data/ancombc_remaining.Rds"))




# same for early infancy samples:
ancom_tables_itt_infancy <- models_imp[[1]]$infancy[[3]] %>% 
  select(taxon, lfc_condition1, se_condition1, p_condition1, q_condition1) %>%
  arrange(q_condition1, desc(abs(lfc_condition1))) %>%
  mutate(across(where(is.numeric), round, 3))
colnames(ancom_tables_itt_infancy) <- str_remove(colnames(ancom_tables_itt), "_condition1") 

ancom_tables_itt2_infancy <- models_imp[[1]]$infancy[[3]] %>% 
  select(taxon, "lfc_condition1:age_s", "se_condition1:age_s", "p_condition1:age_s", "q_condition1:age_s") %>%
  arrange(`q_condition1:age_s`, desc(abs(`lfc_condition1:age_s`))) %>%
  mutate(across(where(is.numeric), round, 3))

colnames(ancom_tables_itt2_infancy) <- str_remove(colnames(ancom_tables_itt2), "_condition1:age_s") 

save(ancom_tables_itt_infancy, file = here::here("data/ancom_tables_itt_infancy.Rds"))

ancom_tables_itt_infancy




t1_infancy <- select(
  models_imp[[1]]$infancy[[3]],
  taxon,
  lfc_age_s) %>%
  pivot_longer(
    contains("lfc"), 
    names_to = "variable", 
    values_to = "lfc", 
    names_prefix = "lfc_")


t2_infancy <- select(
  models_imp[[1]]$infancy[[3]],
  taxon,
  p_age_s
) %>%
  pivot_longer(
    contains("p_"), 
    names_to = "variable", 
    values_to = "p", 
    names_prefix = "p_")
t2
t3_infancy <- select(
  models_imp[[1]]$infancy[[3]],
  taxon,
  q_age_s) %>%
  pivot_longer(
    contains("q_"), 
    names_to = "variable", 
    values_to = "q", 
    names_prefix = "q_")
t3

ancombc_remaining_infancy <- full_join(t1_infancy, t2_infancy, by = c("taxon", "variable")) %>%
  full_join(t3_infancy, by = c("taxon", "variable")) %>%
  filter(p <= 0.05) %>%
  mutate(
    across(where(is.numeric), round, 3),
    taxon = str_remove(taxon, "genus:")
  ) %>%
  arrange(variable, q) 
save(ancombc_remaining_infancy, file = here::here("data/ancombc_remaining_infancy.Rds"))



# same for year1 samples:
ancom_tables_itt_year1 <- models_imp[[1]]$year1[[3]] %>% 
  select(taxon, lfc_condition1, se_condition1, p_condition1, q_condition1) %>%
  arrange(q_condition1, desc(abs(lfc_condition1))) %>%
  mutate(across(where(is.numeric), round, 3))
colnames(ancom_tables_itt_year1) <- str_remove(colnames(ancom_tables_itt), "_condition1") 

ancom_tables_itt2_year1 <- models_imp[[1]]$year1[[3]] %>% 
  select(taxon, "lfc_condition1:age_s", "se_condition1:age_s", "p_condition1:age_s", "q_condition1:age_s") %>%
  arrange(`q_condition1:age_s`, desc(abs(`lfc_condition1:age_s`))) %>%
  mutate(across(where(is.numeric), round, 3))

colnames(ancom_tables_itt2_year1) <- str_remove(colnames(ancom_tables_itt2), "_condition1:age_s") 

save(ancom_tables_itt_year1, file = here::here("data/ancom_tables_itt_year1.Rds"))

ancom_tables_itt_year1




t1_year1 <- select(
  models_imp[[1]]$year1[[3]],
  taxon,
  lfc_age_s) %>%
  pivot_longer(
    contains("lfc"), 
    names_to = "variable", 
    values_to = "lfc", 
    names_prefix = "lfc_")


t2_year1 <- select(
  models_imp[[1]]$year1[[3]],
  taxon,
  p_age_s
) %>%
  pivot_longer(
    contains("p_"), 
    names_to = "variable", 
    values_to = "p", 
    names_prefix = "p_")
t2
t3_year1 <- select(
  models_imp[[1]]$year1[[3]],
  taxon,
  q_age_s) %>%
  pivot_longer(
    contains("q_"), 
    names_to = "variable", 
    values_to = "q", 
    names_prefix = "q_")
t3

ancombc_remaining_year1 <- full_join(t1_year1, t2_year1, by = c("taxon", "variable")) %>%
  full_join(t3_year1, by = c("taxon", "variable")) %>%
  filter(p <= 0.05) %>%
  mutate(
    across(where(is.numeric), round, 3),
    taxon = str_remove(taxon, "genus:")
  ) %>%
  arrange(variable, q) 
save(ancombc_remaining_year1, file = here::here("data/ancombc_remaining_year1.Rds"))






















###############################################################################
#########################           2. PP        ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

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


# for analyses we apply prevalence fitlering
tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)
fvars <- c("siblings", "condition")
# add metadata to tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  select(-siblings, -ges_age, -birthweight, -edlevel, -csection, -sex) %>%
  left_join(select(d_pp, age, sample_id, siblings, pp, csection, sex, ges_age, edlevel, birthweight), by = "sample_id") %>%
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
if (!file.exists(here::here("data/daa_pp_mlm_bez.Rds"))) {
  output <- ancombc2(
    data = tse, 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
    rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl(),
    # mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
    # trend_control = list(
    #   contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
    #   node = list(2), 
    #   solver = "ECOS",B = 100)
    # 
    )
    save(output, file = here::here("data/daa_pp_mlm_bez.Rds"))
 } else {
   load(here::here("data/daa_pp_mlm_bez.Rds"))
 }
   
# structural zeros 
tab_zero = output$zero_ind
sum(tab_zero[, 2]) + sum(tab_zero[, 3])
# sensitivity scores 
tab_sens = output$pseudo_sens_tab
head(tab_sens)
pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
  filter(value >5)

res_prim <- output$res
colnames(res_prim)
effects <- c("age_s", "condition1", "condition1:age_s")
map(effects, function(effect) {
  select(res_prim, 
    taxon, 
    lfc = glue::glue("lfc_{effect}"),
    se = glue::glue("se_{effect}"),
    indicator = glue::glue("diff_{effect}")) %>%
  filter(indicator) %>%
  mutate(effect = effect)
})


# model excludes random intercepts and includes all samples
if (!file.exists(here::here("data/daa_pp_bez.Rds"))) {
  output <- ancombc2(
    data = tse, 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
    #rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl()
    )
    save(output, file = here::here("data/daa_pp_bez.Rds"))
 } else {
   load(here::here("data/daa_pp_bez.Rds"))
 }
   
 # structural zeros 
 tab_zero = output$zero_ind
 sum(tab_zero[, 2]) + sum(tab_zero[, 3])
 # sensitivity scores 
 tab_sens = output$pseudo_sens_tab
 head(tab_sens)
 pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
   filter(value >5)
 
 res_prim <- output$res
 colnames(res_prim)
 effects <- c("age_s", "condition1", "condition1:age_s")
 map(effects, function(effect) {
   select(res_prim, 
     taxon, 
     lfc = glue::glue("lfc_{effect}"),
     se = glue::glue("se_{effect}"),
     indicator = glue::glue("diff_{effect}")) %>%
   filter(indicator) %>%
   mutate(effect = effect)
 })

# model includes random intercepts and excludes 1 year samples
if (!file.exists(here::here("data/daa_pp_mlm_infancy_bez.Rds"))) {
  output <- ancombc2(
    data = filter(tse, week != 52), 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
    rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl()
    )
    save(output, file = here::here("data/daa_pp_mlm_infancy_bez.Rds"))
 } else {
   load(here::here("data/daa_pp_mlm_infancy_bez.Rds"))
}
   
# structural zeros 
tab_zero = output$zero_ind
sum(tab_zero[, 2]) + sum(tab_zero[, 3])
# sensitivity scores 
tab_sens = output$pseudo_sens_tab
head(tab_sens)
pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
 filter(value >5)

res_prim <- output$res
colnames(res_prim)
effects <- c("age_s", "condition1", "condition1:age_s")
map(effects, function(effect) {
 select(res_prim, 
   taxon, 
   lfc = glue::glue("lfc_{effect}"),
   se = glue::glue("se_{effect}"),
   indicator = glue::glue("diff_{effect}")) %>%
 filter(indicator) %>%
 mutate(effect = effect)
})


# model includes 1 years samples only 
if (!file.exists(here::here("data/daa_pp_year1_bez.Rds"))) {
  output <- ancombc2(
    data = filter(tse, week == 52), 
    assay_name = "counts", 
    tax_level = "genus",
    fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
    # rand_formula = "(1 | skippy_id)",
    p_adj_method = "holm", 
    pseudo = 0, 
    pseudo_sens = TRUE,
    prv_cut = 0.10, 
    lib_cut = 1000, 
    s0_perc = 0.05,
    group = "condition", 
    struc_zero = TRUE, 
    neg_lb = TRUE,
    alpha = 0.2, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE, 
    pairwise = FALSE, 
    dunnet = FALSE, 
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = lme4::lmerControl()
    )
    save(output, file = here::here("data/daa_pp_year1_bez.Rds"))
 } else {
   load(here::here("data/daa_pp_year1_bez.Rds"))
 }
   
# structural zeros 
tab_zero = output$zero_ind
sum(tab_zero[, 2]) + sum(tab_zero[, 3])
 # sensitivity scores 
tab_sens = output$pseudo_sens_tab
head(tab_sens)
pivot_longer(tab_sens, -taxon, names_to = "contrast", values_to = "value") %>%
   filter(value >5)
 
res_prim <- output$res
colnames(res_prim)
effects <- c("age_s", "condition1", "condition1:age_s")
map(effects, function(effect) {
   select(res_prim, 
     taxon, 
     lfc = glue::glue("lfc_{effect}"),
     se = glue::glue("se_{effect}"),
     indicator = glue::glue("diff_{effect}")) %>%
   filter(indicator) %>%
   mutate(effect = effect)
})





######################## 2.2 Multiple imputation  #############################



load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
# for analyses we apply prevalence fitlering
tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)

models_imp <- map2(implist_pp, 1:length(implist), function(dimp, imp) {
  tse_map <- tse 
  fvars <- c("siblings", "condition", "csection", "sex")
  # add metadata to tse
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(-siblings, -ges_age, -birthweight, -edlevel, -csection, -sex) %>%
    left_join(select(dimp, age, sample_id, siblings, pp, csection, sex, ges_age_s, edlevel_s, birthweight_s), by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(across(all_of(fvars), function(x) as.factor(x))) %>%
    DataFrame()
  colData(tse_map)$age <- colData(tse_map)$age + as.numeric(colData(tse_map)$week) * 7
  colData(tse_map)$age_s <- scale(colData(tse_map)$age)[, 1]
  
  tse_map <- filter(tse_map, pp == 1)
  
  # model includes random intercepts and excludes 1 year samples
  if (!file.exists(here::here(glue("data/daa_pp_mlm_imp{imp}_bez.Rds")))) {
    output <- ancombc2(
      data = tse_map, 
      assay_name = "counts", 
      tax_level = "genus",
      fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
      rand_formula = "(1 | skippy_id)",
      p_adj_method = "holm", 
      pseudo = 0, 
      pseudo_sens = TRUE,
      prv_cut = 0.10, 
      lib_cut = 1000, 
      s0_perc = 0.05,
      group = "condition", 
      struc_zero = TRUE, 
      neg_lb = TRUE,
      alpha = 0.2, 
      n_cl = 2, 
      verbose = TRUE,
      global = FALSE, 
      pairwise = FALSE, 
      dunnet = FALSE, 
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
        node = list(2), 
        solver = "ECOS",B = 100)
    )
    save(output, file = here::here(glue("data/daa_pp_mlm_imp{imp}_bez.Rds")))
  } else {
    load(here::here(glue("data/daa_pp_mlm_imp{imp}_bez.Rds")))
  }
  
  # structural zeros 
  tab_zero_all = output$zero_ind
  # sensitivity scores 
  tab_sens_all = output$pseudo_sens_tab
  sens_scores_all <- pivot_longer(tab_sens_all, -taxon, names_to = "contrast", values_to = "value") %>%
    filter(value >5)
  res_prim_all <- output$res
  effects_all <- map(effects, function(effect) {
    select(res_prim_all, 
           taxon, 
           lfc = glue::glue("lfc_{effect}"),
           se = glue::glue("se_{effect}"),
           indicator = glue::glue("diff_{effect}")) %>%
      filter(indicator) %>%
      mutate(effect = effect)
  })
  
  
  
  # model includes random intercepts and excludes 1 year samples
  if (!file.exists(here::here(glue("data/daa_pp_mlm_infancy_imp{imp}_bez.Rds")))) {
    output <- ancombc2(
      data = filter(tse_map, week != 52), 
      assay_name = "counts", 
      tax_level = "genus",
      fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
      rand_formula = "(1 | skippy_id)",
      p_adj_method = "holm", 
      pseudo = 0, 
      pseudo_sens = TRUE,
      prv_cut = 0.10, 
      lib_cut = 1000, 
      s0_perc = 0.05,
      group = "condition", 
      struc_zero = TRUE, 
      neg_lb = TRUE,
      alpha = 0.2, 
      n_cl = 2, 
      verbose = TRUE,
      global = FALSE, 
      pairwise = FALSE, 
      dunnet = FALSE, 
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
        node = list(2), 
        solver = "ECOS",B = 100)
      )
      save(output, file = here::here(glue("data/daa_pp_mlm_infancy_imp{imp}_bez.Rds")))
   } else {
     load(here::here(glue("data/daa_pp_mlm_infancy_imp{imp}_bez.Rds")))
  }
     
  # structural zeros 
  tab_zero_infancy = output$zero_ind
  # sensitivity scores 
  tab_sens_infancy = output$pseudo_sens_tab
  sens_scores_infancy <- pivot_longer(tab_sens_infancy, -taxon, names_to = "contrast", values_to = "value") %>%
     filter(value >5)
  res_prim_infancy <- output$res
  effects_infancy <- map(effects, function(effect) {
     select(res_prim_infancy, 
       taxon, 
       lfc = glue::glue("lfc_{effect}"),
       se = glue::glue("se_{effect}"),
       indicator = glue::glue("diff_{effect}")) %>%
     filter(indicator) %>%
     mutate(effect = effect)
   })
  
  
  # model includes 1 years samples only 
  if (!file.exists(here::here(glue("data/daa_pp_year1_imp{imp}_bez.Rds")))) {
    output <- ancombc2(
      data = filter(tse_map, week == 52), 
      assay_name = "counts", 
      tax_level = "genus",
      fix_formula = "condition * age_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex",
      # rand_formula = "(1 | skippy_id)",
      p_adj_method = "holm", 
      pseudo = 0, 
      pseudo_sens = TRUE,
      prv_cut = 0.10, 
      lib_cut = 1000, 
      s0_perc = 0.05,
      group = "condition", 
      struc_zero = TRUE, 
      neg_lb = TRUE,
      alpha = 0.2, 
      n_cl = 2, 
      verbose = TRUE,
      global = FALSE, 
      pairwise = FALSE, 
      dunnet = FALSE, 
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)), 
        node = list(2), 
        solver = "ECOS",B = 100)
      )
      save(output, file = here::here(glue("data/daa_pp_year1_imp{imp}_bez.Rds")))
   } else {
     load(here::here(glue("data/daa_pp_year1_imp{imp}_bez.Rds")))
   }
     
  # structural zeros 
  tab_zero_year1 = output$zero_ind
  # sensitivity scores 
  tab_sens_year1 = output$pseudo_sens_tab
  sens_scores_year1 <- pivot_longer(tab_sens_year1, -taxon, names_to = "contrast", values_to = "value") %>%
     filter(value >5)
  res_prim_year1 <- output$res
  effects_year1 <- map(effects, function(effect) {
     select(res_prim_year1, 
       taxon, 
       lfc = glue::glue("lfc_{effect}"),
       se = glue::glue("se_{effect}"),
       indicator = glue::glue("diff_{effect}")) %>%
     filter(indicator) %>%
     mutate(effect = effect)
   })

  list(
    all = list(
      tab_zero_all,
      sens_scores_all,
      res_prim_all,
      effects_all
    ),
    infancy = list(
      tab_zero_infancy,
      sens_scores_infancy,
      res_prim_infancy,
      effects_infancy
    ),
    year1 = list(
      tab_zero_year1,
      sens_scores_year1,
      res_prim_year1,
      effects_year1
    )
  )
  
})





ancom_tables_pp <- models_imp[[1]]$all[[3]] %>% 
  select(taxon, lfc_condition1, se_condition1, p_condition1, q_condition1) %>%
  arrange(q_condition1, desc(abs(lfc_condition1))) %>%
  mutate(across(where(is.numeric), round, 3))
colnames(ancom_tables_pp) <- str_remove(colnames(ancom_tables_pp), "_condition1") 

ancom_tables_pp2 <- models_imp[[1]]$all[[3]] %>% 
  select(taxon, "lfc_condition1:age_s", "se_condition1:age_s", "p_condition1:age_s", "q_condition1:age_s") %>%
  arrange(`q_condition1:age_s`, desc(abs(`lfc_condition1:age_s`))) %>%
  mutate(across(where(is.numeric), round, 3))

colnames(ancom_tables_pp2) <- str_remove(colnames(ancom_tables_pp2), "_condition1:age_s") 

save(ancom_tables_pp, file = here::here("data/ancom_tables_pp.Rds"))

ancom_tables_pp
ancom_tables_pp2








t1 <- select(
  models_imp[[1]]$all[[3]],
  taxon,
  lfc_age_s) %>%
  pivot_longer(
    contains("lfc"), 
    names_to = "variable", 
    values_to = "lfc", 
    names_prefix = "lfc_")


t2 <- select(
  models_imp[[1]]$all[[3]],
  taxon,
  p_age_s
) %>%
  pivot_longer(
    contains("p_"), 
    names_to = "variable", 
    values_to = "p", 
    names_prefix = "p_")
t2
t3 <- select(
  models_imp[[1]]$all[[3]],
  taxon,
  q_age_s) %>%
  pivot_longer(
    contains("q_"), 
    names_to = "variable", 
    values_to = "q", 
    names_prefix = "q_")
t3

ancombc_remaining_pp <- full_join(t1, t2, by = c("taxon", "variable")) %>%
  full_join(t3, by = c("taxon", "variable")) %>%
  filter(p <= 0.05) %>%
  mutate(
    across(where(is.numeric), round, 3),
    taxon = str_remove(taxon, "genus:")
  ) %>%
  arrange(variable, q) 
save(ancombc_remaining_pp, file = here::here("data/ancombc_remaining_pp.Rds"))




# same for early infancy samples:
ancom_tables_pp_infancy <- models_imp[[1]]$infancy[[3]] %>% 
  select(taxon, lfc_condition1, se_condition1, p_condition1, q_condition1) %>%
  arrange(q_condition1, desc(abs(lfc_condition1))) %>%
  mutate(across(where(is.numeric), round, 3))
colnames(ancom_tables_pp_infancy) <- str_remove(colnames(ancom_tables_pp), "_condition1") 

ancom_tables_pp2_infancy <- models_imp[[1]]$infancy[[3]] %>% 
  select(taxon, "lfc_condition1:age_s", "se_condition1:age_s", "p_condition1:age_s", "q_condition1:age_s") %>%
  arrange(`q_condition1:age_s`, desc(abs(`lfc_condition1:age_s`))) %>%
  mutate(across(where(is.numeric), round, 3))

colnames(ancom_tables_pp2_infancy) <- str_remove(colnames(ancom_tables_pp2), "_condition1:age_s") 

save(ancom_tables_pp_infancy, file = here::here("data/ancom_tables_pp_infancy.Rds"))

ancom_tables_pp_infancy




t1_infancy <- select(
  models_imp[[1]]$infancy[[3]],
  taxon,
  lfc_age_s) %>%
  pivot_longer(
    contains("lfc"), 
    names_to = "variable", 
    values_to = "lfc", 
    names_prefix = "lfc_")


t2_infancy <- select(
  models_imp[[1]]$infancy[[3]],
  taxon,
  p_age_s
) %>%
  pivot_longer(
    contains("p_"), 
    names_to = "variable", 
    values_to = "p", 
    names_prefix = "p_")
t2
t3_infancy <- select(
  models_imp[[1]]$infancy[[3]],
  taxon,
  q_age_s) %>%
  pivot_longer(
    contains("q_"), 
    names_to = "variable", 
    values_to = "q", 
    names_prefix = "q_")
t3

ancombc_remaining_infancy_pp <- full_join(t1_infancy, t2_infancy, by = c("taxon", "variable")) %>%
  full_join(t3_infancy, by = c("taxon", "variable")) %>%
  filter(p <= 0.05) %>%
  mutate(
    across(where(is.numeric), round, 3),
    taxon = str_remove(taxon, "genus:")
  ) %>%
  arrange(variable, q) 
save(ancombc_remaining_infancy_pp, file = here::here("data/ancombc_remaining_infancy_pp.Rds"))



# same for year1 samples:
ancom_tables_pp_year1 <- models_imp[[1]]$year1[[3]] %>% 
  select(taxon, lfc_condition1, se_condition1, p_condition1, q_condition1) %>%
  arrange(q_condition1, desc(abs(lfc_condition1))) %>%
  mutate(across(where(is.numeric), round, 3))
colnames(ancom_tables_pp_year1) <- str_remove(colnames(ancom_tables_pp), "_condition1") 

ancom_tables_pp2_year1 <- models_imp[[1]]$year1[[3]] %>% 
  select(taxon, "lfc_condition1:age_s", "se_condition1:age_s", "p_condition1:age_s", "q_condition1:age_s") %>%
  arrange(`q_condition1:age_s`, desc(abs(`lfc_condition1:age_s`))) %>%
  mutate(across(where(is.numeric), round, 3))

colnames(ancom_tables_pp2_year1) <- str_remove(colnames(ancom_tables_pp2), "_condition1:age_s") 

save(ancom_tables_pp_year1, file = here::here("data/ancom_tables_pp_year1.Rds"))

ancom_tables_pp_year1




t1_year1 <- select(
  models_imp[[1]]$year1[[3]],
  taxon,
  lfc_age_s) %>%
  pivot_longer(
    contains("lfc"), 
    names_to = "variable", 
    values_to = "lfc", 
    names_prefix = "lfc_")


t2_year1 <- select(
  models_imp[[1]]$year1[[3]],
  taxon,
  p_age_s
) %>%
  pivot_longer(
    contains("p_"), 
    names_to = "variable", 
    values_to = "p", 
    names_prefix = "p_")
t2
t3_year1 <- select(
  models_imp[[1]]$year1[[3]],
  taxon,
  q_age_s) %>%
  pivot_longer(
    contains("q_"), 
    names_to = "variable", 
    values_to = "q", 
    names_prefix = "q_")
t3

ancombc_remaining_year1_pp <- full_join(t1_year1, t2_year1, by = c("taxon", "variable")) %>%
  full_join(t3_year1, by = c("taxon", "variable")) %>%
  filter(p <= 0.05) %>%
  mutate(
    across(where(is.numeric), round, 3),
    taxon = str_remove(taxon, "genus:")
  ) %>%
  arrange(variable, q) 
save(ancombc_remaining_year1_pp, file = here::here("data/ancombc_remaining_year1_pp.Rds"))













###############################################################################
################################# PLOTS #######################################
###############################################################################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
tse <- transformCounts(tse, method = "relabundance")

asv_tab <- t(assay(tse, "relabundance")) %>% as.data.frame() %>%
  select(all_of(c("genus:Parabacteroides", "genus:Flavonifractor", "genus:Faecalibacterium"))) %>%
  rownames_to_column("sample_id")
colnames(asv_tab) <- str_remove(colnames(asv_tab), "genus:")
meta <- colData(tse) %>% as.data.frame() %>%
  mutate(
    condition_label = ifelse(condition == 1, "SSC", ifelse(condition == 0, "CAU", NA)),
    week_label = glue::glue("Week {week}")
    ) %>%
  select(sample_id, condition_label, week_label)

ds <- meta %>% full_join(asv_tab, by = "sample_id") %>%
  pivot_longer(
    all_of(c("Parabacteroides", "Flavonifractor", "Faecalibacterium")),
    names_to = "Genus", values_to = "Abundance"  
  ) %>%
  group_by(Genus) %>%
  nest()


daa_plots <- map2(ds[[1]], ds[[2]], function(genus, d) {
    p <- d %>%
    ggplot(aes_string("condition_label", "Abundance", fill = "condition_label")) +
      geom_boxplot(outlier.alpha = 0) +
      #ggbeeswarm::geom_beeswarm(size = 3, cex = 1) +
      geom_jitter(width = 0.1, size = 2.5, alpha = 0.5) +
      scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
      theme_bw(base_size = 25) +
      theme(
        legend.position = "none") +
      xlab("") + ylab(str_to_title(genus))
  if (genus == "Parabacteroides") {
    p <- p + facet_wrap(~week_label, strip.position = "bottom") +
      theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank())
  }
  p
})
daa_plots
daa_plots <- map2(ds[[1]], ds[[2]], function(genus, d) {
    if (genus == "Parabacteroides") {
      df <- d 
    } else {
      df <- filter(d, week_label == "Week 52")
    }
    p <- df %>%
    ggplot(aes_string("condition_label", "Abundance", fill = "condition_label")) +
      geom_boxplot(outlier.alpha = 0) +
      #ggbeeswarm::geom_beeswarm(size = 3, cex = 1) +
      geom_jitter(width = 0.1, size = 2.5, alpha = 0.5) +
      scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
      theme_bw(base_size = 25) +
      theme(
        legend.position = "none") +
      xlab("") + ylab(str_to_title(genus))
  if (genus == "Parabacteroides") {
    p <- p + facet_wrap(~week_label, strip.position = "bottom") +
      theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank())
  }
  p
})
daa_plots[[1]]

library(patchwork)
daa_plot <-  daa_plots[[1]] /
  (daa_plots[[2]] + daa_plots[[3]]) 

daa_plot
save(daa_plot, file = here::here("data/daa_plots.Rds"))



