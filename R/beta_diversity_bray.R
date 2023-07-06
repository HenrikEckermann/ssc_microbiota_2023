# thanks to Gavin for helping out with this part:
# https://stats.stackexchange.com/questions/590510/repeated-measures-permanova-nowhere-to-find

set.seed(1)
library(mia)
library(tidyverse)
library(tidySummarizedExperiment)
library(vegan)
library(permute)



# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
tse <- transformSamples(x = tse, method = "relabundance", pseudocount = 1, name = "relabundance")
# add metadata to tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  left_join(select(d, age, sample_id), by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  DataFrame()

colData(tse)$age <- colData(tse)$age + as.numeric(as.character(colData(tse)$week)) * 7
colData(tse)$age_s <- scale(colData(tse)$age)[, 1]

###############################################################################
#########################           1. ITT       ##############################
###############################################################################


####################### 1.1 Complete Case Analysis ############################

### First I fit a model to all samples 


# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse, "counts"))
asv <- asv[rownames(asv) %in% meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "bray",
                     # does not work if trend is in data (therefore use 999)
                     permutations = 999 
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "bray",
                        permutations = 999)
permanova2



### Now split models by infancy and 1 year olds 

## first infancy
tse_inf <- filter(tse, week != 52)
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_inf) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_inf, "counts"))
asv <- asv[rownames(asv) %in% meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "bray",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "bray",
                        permutations = 999)
permanova2




## then 1 year olds 
tse_y <- filter(tse, week == 52)
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_y) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings) %>%
                  na.omit()

# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_y, "counts"))
asv <- asv[rownames(asv) %in% meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "bray",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "bray",
                        permutations = 999)
permanova2




######################## 1.2 Multiple imputation  #############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
permanovas <- map2_dfr(implist, 1:length(implist), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  tse_map <- transformSamples(x = tse, method = "relabundance", pseudocount = 1, name = "relabundance")
  # step 1
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(sample_id) %>%
    left_join(
      select(dimp, condition, siblings, age, sample_id, week, skippy_id), 
      by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(
      age = age + as.numeric(as.character(week)),
      age_s = scale(age)[, 1]
    ) %>%
    DataFrame()

    ## first all samples
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_map) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_map, "counts"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~ condition + age_s + condition:age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
                       )

    
    # step 2
    ## first infancy
    tse_inf <- filter(tse_map, week != 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_inf) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_inf, "counts"))
    asv <- asv[rownames(asv) %in% meta$sample_id, ]

    # fit and inspect model
    permanova_inf <- adonis2(asv ~ condition + age_s + condition:age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
                       )


    ## then 1 year olds
    tse_y <- filter(tse_map, week == 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_y) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings) %>%
                      na.omit()

    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_y, "counts"))
    asv <- asv[rownames(asv) %in% meta$sample_id, ]

    # fit and inspect model 
    permanova_y <- adonis2(asv ~ condition + age_s + condition:age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
                       )

    permanova_all$time <- "all"
    permanova_inf$time <- "infancy"
    permanova_y$time <- "year1"
    permanova <- bind_rows(
      as.data.frame(permanova_all) %>% rownames_to_column("parameter"),
      as.data.frame(permanova_inf) %>% rownames_to_column("parameter"), 
      as.data.frame(permanova_y) %>% rownames_to_column("parameter")
      )
    permanova$imp <- imp
    permanova
    #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
})

permanovas




######################## 1.2 WITH BREASTFEEDING  #############################



load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

permanovas_bf <- map2_dfr(implist, 1:length(implist), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  
  # step 1
  # we use Aitchison distance
  tse_map <- agglomerateByRank(tse, rank = "genus")
  tse_map <- transformSamples(x = tse_map, method = "relabundance", pseudocount = 1, name = "relabundance")
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(sample_id) %>%
    left_join(
      select(
        dimp, condition, siblings, age, 
        sample_id, 
        week, skippy_id, bfexcl), 
      by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(
      age = age + as.numeric(as.character(week)),
      age_s = scale(age)[, 1]
    ) %>%
    DataFrame()
    

    # step 2

    ## first all samples
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_map) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings,
                      bfexcl) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_map, "counts"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~ bfexcl + siblings  + condition + age_s + condition:age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999,
                         by = "margin"
                       )


    ## first infancy
    tse_inf <- filter(tse_map, week != 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_inf) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_inf, "counts"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_inf <- adonis2(asv ~ condition + age_s + condition:age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999
                       )

    ## then 1 year olds 
    tse_y <- filter(tse_map, week == 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_y) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings) %>%
                      na.omit()

    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_y, "counts"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_y <- adonis2(asv ~ condition + age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999
                       )
    permanova_all$time <- "all"
    permanova_inf$time <- "infancy"
    permanova_y$time <- "year1"
    permanova <- bind_rows(
      as.data.frame(permanova_all) %>% rownames_to_column("parameter"), 
      as.data.frame(permanova_inf) %>% rownames_to_column("parameter"), 
      as.data.frame(permanova_y) %>% rownames_to_column("parameter")
      )
    permanova$imp <- imp 
     permanova
    # #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
})
permanovas_bf





















###############################################################################
#########################         2. PP          ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
tse <- transformSamples(x = tse, method = "relabundance", pseudocount = 1, name = "relabundance")

# obtain ids that were selected for PP analyses 
pp_indicator <- foreign::read.spss(here::here("data/raw_data/kelly141022/Data_ITT_PP_ExploratoryDRselections.sav"), to.data.frame = TRUE)
count(pp_indicator, PP, SSC)
pp_indicator <- select(pp_indicator, skippy_id = ID, pp = PP)
# add pp info to existing data 
if (!"pp" %in% colnames(d)) {
  d <- left_join(d, pp_indicator, by = "skippy_id")
}
# 60 that are in PP and condition 0; 18 that are condition 1 and pp. Fits...
count(d, condition, pp)
d_pp <- filter(d, pp == 1)
implist_pp <- map(implist, function(dimp) {
  dimp_new <- left_join(dimp, pp_indicator, by = "skippy_id") %>%
                filter(pp == 1)
  dimp_new
})

# add metadata to tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  left_join(select(d_pp, age, sample_id, pp), by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  DataFrame()

colData(tse)$age <- colData(tse)$age + as.numeric(as.character(colData(tse)$week)) * 7
colData(tse)$age_s <- scale(colData(tse)$age)[, 1]
tse_pp <- filter(tse, pp == 1)

####################### 2.1 Complete Case Analysis ############################

### First I fit a model to all samples 



# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_pp) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_pp, "counts"))
asv <- asv[rownames(asv) %in% meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~  condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "bray",
                     # does not work if trend is in data (therefore use 999)
                     permutations = 999 
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~  condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "bray",
                        permutations = 999)
permanova2



### Now split models by infancy and 1 year olds 

## first infancy
tse_inf <- filter(tse_pp, week != 52)
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_inf) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_inf, "counts"))
asv <- asv[rownames(asv) %in% meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~  condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "bray",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~  condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "bray",
                        permutations = 999)
permanova2




## then 1 year olds 
tse_y <- filter(tse_pp, week == 52)
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_y) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings) %>%
                  na.omit()

# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_y, "counts"))
asv <- asv[rownames(asv) %in% meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~  condition + age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "bray",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~  condition + age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "bray",
                        permutations = 999)
permanova2

# no effects in PP either!


######################## 2.2 Multiple imputation  #############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
permanovas_pp <- map2_dfr(implist_pp, 1:length(implist_pp), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  
  # step 1
  # we use Aitchison distance
  tse_map <- transformSamples(x = tse, method = "relabundance", pseudocount = 1, name = "relabundance")
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(sample_id) %>%
    left_join(
      select(dimp, condition, siblings, age, birthweight_s, ges_age_s, edlevel_s, csection, sex, sample_id, week, skippy_id, pp), 
      by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(
      age = age + as.numeric(as.character(week)),
      age_s = scale(age)[, 1]
    ) %>%
    DataFrame()
  
  tse_map <- filter(tse_map, pp == 1)
  
  
  # step 2
  ## first all samples
  # extract relevant meta data and omit na as adonis doesnt accept them.
  meta <- colData(tse_map) %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, skippy_id, condition, age_s, siblings,
           birthweight_s, ges_age_s,
           edlevel_s, csection, sex) %>%
    na.omit()
  # we need to account for non-independence of data in the infancy model 
  ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
    .$skippy_id
  h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_map, "counts"))
  asv <- asv[meta$sample_id, ]
  
  # fit and inspect model 
  permanova_all <- adonis2(asv ~ siblings  + condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex,
                           # by = "margin", # each term analyzed individually
                           data = meta,
                           method = "bray",
                           permutations = 999
  )
  
  ## then infancy
  tse_inf <- filter(tse_map, week != 52)
  # extract relevant meta data and omit na as adonis doesnt accept them.
  meta <- colData(tse_inf) %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, skippy_id, condition, age_s, siblings,
           birthweight_s, ges_age_s, 
           edlevel_s, csection, sex) %>%
    na.omit()
  # we need to account for non-independence of data in the infancy model 
  ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
    .$skippy_id
  h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_inf, "counts"))
  asv <- asv[meta$sample_id, ]
  
  # fit and inspect model 
  permanova_inf <- adonis2(asv ~ siblings  + condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex,
                           # by = "margin", # each term analyzed individually
                           data = meta,
                           method = "bray",
                           permutations = 999
  )
  
  
  ## then 1 year olds 
  tse_y <- filter(tse_map, week == 52)
  # extract relevant meta data and omit na as adonis doesnt accept them.
  meta <- colData(tse_y) %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, skippy_id, condition, age_s, siblings,
           birthweight_s, ges_age_s,
           edlevel_s, csection, sex) %>%
    na.omit()
  
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_y, "counts"))
  asv <- asv[meta$sample_id, ]
  
  # fit and inspect model 
  permanova_y <- adonis2(asv ~ siblings  + condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
  )
  permanova_all$time <- "all"
  permanova_inf$time <- "infancy"
  permanova_y$time <- "year1"
  permanova <- bind_rows(
    as.data.frame(permanova_all) %>% rownames_to_column("parameter"), 
    as.data.frame(permanova_inf) %>% rownames_to_column("parameter"), 
    as.data.frame(permanova_y) %>% rownames_to_column("parameter")
  )
  permanova$imp <- imp 
  permanova
  #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
})

permanovas_pp





###############################################################################
################################ DR analyses ##################################
###############################################################################


load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

dr <- foreign::read.spss(
  here::here("data/kelly_documents/data_itt_pp_dr.sav"),
  to.data.frame = TRUE
  ) %>%
  select(skippy_id = ID, ITT, SSC = TotalSSCwk1wk5) %>%
  mutate(SSC_s = scale(SSC)[, 1])

implist <- map(implist, function(imp) {
  imp_new <- imp %>%
    left_join(
      select(dr, skippy_id, SSC_s),
      by = "skippy_id")
  #mice::complete(mice::mice(imp_new))
})

permanovas <- map2_dfr(implist, 1:length(implist), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  tse_map <- transformSamples(x = tse, method = "relabundance", pseudocount = 1, name = "relabundance")
  # step 1
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(sample_id) %>%
    left_join(
      select(dimp, condition, siblings, age, sample_id, week, skippy_id, SSC_s), 
      by = "sample_id") %>%
    column_to_rownames("sample_id") %>%
    mutate(
      age = age + as.numeric(as.character(week)),
      age_s = scale(age)[, 1]
    ) %>%
    DataFrame()

    ## first all samples
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_map) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings,
                      SSC_s) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_map, "counts"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~  SSC_s + age_s + SSC_s:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
                       )

   
    # step 2
    ## first infancy
    tse_inf <- filter(tse_map, week != 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_inf) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings,
                      SSC_s) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_inf, "counts"))
    asv <- asv[rownames(asv) %in% meta$sample_id, ]

    # fit and inspect model
    permanova_inf <- adonis2(asv ~  SSC_s + age_s + SSC_s:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
                       )


    ## then 1 year olds
    tse_y <- filter(tse_map, week == 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_y) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings,
                      SSC_s) %>%
                      na.omit()

    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_y, "counts"))
    asv <- asv[rownames(asv) %in% meta$sample_id, ]

    # fit and inspect model 
    permanova_y <- adonis2(asv ~  SSC_s + age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "bray",
                         permutations = 999
                       )

    permanova_all$time <- "all"
    permanova_inf$time <- "infancy"
    permanova_y$time <- "year1"
    permanova <- bind_rows(
      as.data.frame(permanova_all) %>% rownames_to_column("parameter"),
      as.data.frame(permanova_inf) %>% rownames_to_column("parameter"), 
      as.data.frame(permanova_y) %>% rownames_to_column("parameter")
      )
    permanova$imp <- imp
    permanova
    #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
})

permanovas
