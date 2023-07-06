set.seed(1)
library(mia)
library(Maaslin2)
library(tidyverse)
library(tidySummarizedExperiment)
library(glue)




###############################################################################
#########################           1. ITT       ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))


# for analyses we apply prevalence fitlering
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



asv_tab <- t(assay(tse))
meta <- colData(tse) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars))
asv_tab <- asv_tab[rownames(meta),]

# you can specifiy different GLMs/normalizations/transforms. We used similar
# settings as in Nearing et al. (2021) here:
fit_data <- Maaslin2(
  asv_tab,
  meta,
  output = here::here("data/maaslin/1"),
  transform = "AST",
  fixed_effects = c("condition", "age_s"),
  random_effects = "skippy_id", 
  reference = "0",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
filter(fit_data$results, qval <= 0.5, metadata == "condition")

# model includes random intercepts and excludes 1 year samples
tse_infancy <- filter(tse, week != 52)
asv_tab <- t(assay(tse))
meta <- colData(tse_infancy) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars))
asv_tab <- asv_tab[rownames(meta), ]

fit_data <- Maaslin2(
  asv_tab,
  meta,
  output = here::here("data/maaslin/1"),
  transform = "AST",
  fixed_effects = c("condition", "age_s"),
  random_effects = "skippy_id", 
  reference = "0",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
filter(fit_data$results, qval <= 0.4, metadata == "condition")




# model includes 1 years samples only 
tse_year1 <- filter(tse, week == 52)
asv_tab <- t(assay(tse))
meta <- colData(tse_year1) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars))
asv_tab <- asv_tab[rownames(meta), ]


fit_data <- Maaslin2(
  asv_tab,
  meta,
  output = here::here("data/maaslin/1"),
  transform = "AST",
  fixed_effects = c("condition", "age_s"),
  #random_effects = "skippy_id", 
  reference = "0",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
filter(fit_data$results, qval <= 0.4, metadata == "condition")




######################## 1. 2Multiple imputation  #############################


load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
# for analyses we apply prevalence fitlering
tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)


if (!file.exists(here::here("data/maaslin2_itt_mi_bez.Rds"))) {
  models_imp <- map2(implist, 1:length(implist), function(dimp, imp) {
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

    # all samples
    asv_tab <- t(assay(tse_map))
    meta <- colData(tse_map) %>% as.data.frame() %>%
      select(skippy_id, age_s, all_of(fvars))
    asv_tab <- asv_tab[rownames(meta), ]

    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = c("condition", "age_s"),
      random_effects = "skippy_id", 
      reference = "0",  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )
    
    res_all <- fit_data$results



    
    # model includes random intercepts and excludes 1 year samples
    tse_infancy <- filter(tse_map, week != 52)
    asv_tab <- t(assay(tse_infancy))
    meta <- colData(tse_infancy) %>% as.data.frame() %>%
      select(skippy_id, age_s, all_of(fvars))
    asv_tab <- asv_tab[rownames(meta), ]


    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = c("condition", "age_s"),
      random_effects = "skippy_id", 
      reference = "0",  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )

    res_inf <- fit_data$results


    # model includes 1 years samples only 
    tse_year1 <- filter(tse_map, week == 52)
    asv_tab <- t(assay(tse_year1))
    meta <- colData(tse_year1) %>% as.data.frame() %>%
      select(skippy_id, age_s, all_of(fvars))
    asv_tab <- asv_tab[rownames(meta), ]

    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = c("condition", "age_s"),
      #random_effects = "skippy_id", 
      reference = "0",  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )

    res_year1 <- fit_data$results

    list(res_all, res_inf, res_year1)
  })
  save(models_imp, file = here::here("data/maaslin2_itt_mi_bez.Rds"))
 } else {
  load(here::here("data/maaslin2_itt_mi_bez.Rds"))
 }


# change [[1]] to 2-5 to inspect the other imputations
maaslin2_tables_itt <- map(models_imp[[1]], function(x) {
    x %>% filter(metadata == "condition") %>%
      select(feature, coef, stderr, pval, qval) %>%
      arrange(qval, desc(abs(coef))) %>%
      mutate(
        across(where(is.numeric), round, 3),
        feature = str_replace(feature, "\\.\\.", "."),
        feature = str_replace(feature, "\\.", ":")
        )
})
maaslin2_tables_itt

# since MaAsLin does not support interactions I need to evaluate early and late infancy separately
maaslin2_identified <- map2_dfr(maaslin2_tables_itt, c("all", "early", "late"), function(table, term) {
  filter(table, qval <= 0.2) %>%
  mutate(term = term) %>%
  arrange(coef, qval)
})

maaslin2_identified
save(maaslin2_tables_itt, maaslin2_identified, file = here::here("data/maaslin2_tables_itt.Rds"))




































###############################################################################
#########################           2. PP        ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)

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


fvars <- c("siblings", "condition", "csection", "sex")
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
asv_tab <- t(assay(tse))
meta <- colData(tse) %>% as.data.frame() %>%
  select(skippy_id, age_s, all_of(fvars))
asv_tab <- asv_tab[rownames(meta), ]

fit_data <- Maaslin2(
  asv_tab,
  meta,
  output = here::here("data/maaslin/1"),
  transform = "AST",
  fixed_effects = c(
    "condition", 
    "siblings", 
    "age_s",
    "csection",
    "sex",
    "ges_age_s",
    "birthweight_s",
    "edlevel_s"
    ),
  random_effects = "skippy_id", 
  reference = "0",
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
filter(fit_data$results, qval <= 0.1, metadata == "condition")







######################## 2.2 Multiple imputation  #############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
# for analyses we apply prevalence fitlering
tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentTaxa(tse, detection = 0.001, prevalence = 0.1)

if (!file.exists(here::here("data/maaslin2_pp_mi_bez.Rds"))) {
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
    
    # all samples
    asv_tab <- t(assay(tse_map))
    meta <- colData(tse_map) %>% as.data.frame() %>%
      select(skippy_id, age_s, all_of(fvars), edlevel_s, birthweight_s, ges_age_s)
    asv_tab <- asv_tab[rownames(meta), ]

    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = c(
        "condition", 
        "siblings", 
        "age_s",
        "csection",
        "sex",
        "ges_age_s",
        "birthweight_s",
        "edlevel_s"
        ),
      random_effects = "skippy_id", 
      reference = "0",  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )
      
    res_all <- fit_data$results


    # model includes random intercepts and excludes 1 year samples
    tse_infancy <- filter(tse_map, week != 52)
    asv_tab <- t(assay(tse_infancy))
    meta <- colData(tse_infancy) %>% as.data.frame() %>%
      select(skippy_id, age_s, all_of(fvars), edlevel_s, birthweight_s, ges_age_s)
    asv_tab <- asv_tab[rownames(meta), ]


    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = c(
        "condition", 
        "siblings", 
        "age_s",
        "csection",
        "sex",
        "ges_age_s",
        "birthweight_s",
        "edlevel_s"
        ),
      random_effects = "skippy_id", 
      reference = "0",  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )

    res_inf <- fit_data$results
    
    
      # model includes 1 years samples only 
    tse_year1 <- filter(tse_map, week == 52)
    asv_tab <- t(assay(tse_year1))
    meta <- colData(tse_year1) %>% as.data.frame() %>%
      select(skippy_id, age_s, all_of(fvars), edlevel_s, birthweight_s, ges_age_s)
    asv_tab <- asv_tab[rownames(meta), ]

    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = c(
        "condition", 
        "siblings", 
        "age_s",
        "csection",
        "sex",
        "ges_age_s",
        "birthweight_s",
        "edlevel_s"
        ),
      #random_effects = "skippy_id", 
      reference = "0",  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )

    res_year1 <- fit_data$results
    
    list(res_all, res_inf, res_year1)
  })
  save(models_imp, file = here::here("data/maaslin2_pp_mi_bez.Rds"))
 } else {
  load(here::here("data/maaslin2_pp_mi_bez.Rds"))
 }





# switch the 1 to 2-5 to check other imputations
maaslin2_tables_itt_pp <- map(models_imp[[1]], function(x) {
    x %>% filter(metadata == "condition") %>%
      select(feature, coef, stderr, pval, qval) %>%
      arrange(qval, desc(abs(coef))) %>%
      mutate(
        across(where(is.numeric), round, 3),
        feature = str_replace(feature, "\\.\\.", "."),
        feature = str_replace(feature, "\\.", ":")
        )
})
maaslin2_tables_itt_pp
save(maaslin2_tables_itt_pp, file = here::here("data/maaslin2_tables_itt_pp.Rds"))