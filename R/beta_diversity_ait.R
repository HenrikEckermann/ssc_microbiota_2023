# thanks to Gavin for helping out with this part:
# https://stats.stackexchange.com/questions/590510/repeated-measures-permanova-nowhere-to-find

set.seed(1)
library(mia)
library(tidyverse)
library(tidySummarizedExperiment)
library(vegan)
library(permute)
library(glue)


# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
# add metadata to tse
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  left_join(select(d, age, sample_id), by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  DataFrame()

colData(tse)$age <- colData(tse)$age + as.numeric(as.character(colData(tse)$week)) * 7
colData(tse)$age_s <- scale(colData(tse)$age)[, 1]

tse <- agglomerateByRank(tse, rank = "genus")

###############################################################################
#########################           1. ITT       ##############################
###############################################################################


####################### 1.1 Complete Case Analysis ############################

### First I fit a model to all samples 

# we use Aitchison distance
tse <- transformSamples(x = tse, method = "clr", pseudocount = 1, name = "clr")
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
asv <- t(assay(tse, "clr"))
asv <- asv[meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     # h does not work if trend is in data (therefore use 999), see Gavins post
                     permutations = 999
                   )

permanova



# Perform dbRDA
dbrda <- dbrda(asv ~ age_s + age_s*condition, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "euclidean",
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
asv <- t(assay(tse_inf, "clr"))
asv <- asv[meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "euclidean",
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
asv <- t(assay(tse_y, "clr"))
asv <- asv[meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "euclidean",
                        permutations = 999)
permanova2




######################## 1. 2Multiple imputation  #############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))


permanovas <- map2_dfr(implist, 1:length(implist), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  
  # step 1
  # we use Aitchison distance
  tse <- agglomerateByRank(tse, rank = "genus")
  tse_map <- transformSamples(x = tse, method = "clr", pseudocount = 1, name = "clr")
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
    

    # step 2

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
    asv <- t(assay(tse_map, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~ condition + age_s + age_s:condition,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999
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
    asv <- t(assay(tse_inf, "clr"))
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
    asv <- t(assay(tse_y, "clr"))
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
    
    #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
    permanova
})
permanovas




######################## 1. 2 WITH BREASTFEEDING  #############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

permanovas_bf <- map2_dfr(implist, 1:length(implist), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  
  # step 1
  # we use Aitchison distance
  tse <- agglomerateByRank(tse, rank = "genus")
  tse_map <- transformSamples(x = tse, method = "clr", pseudocount = 1, name = "clr")
  colData(tse_map) <- colData(tse_map) %>%
    as.data.frame() %>%
    select(sample_id) %>%
    left_join(
      select(
        dimp, condition, siblings, age, sample_id, 
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
              select(sample_id, skippy_id, condition, age_s, siblings, bfexcl) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_map, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~ bfexcl + condition + age_s + condition:age_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999,
                         # by = "margin"
                       )


    ## first infancy
    tse_inf <- filter(tse_map, week != 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_inf) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings, bfexcl) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_inf, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_inf <- adonis2(asv ~ bfexcl + condition + age_s + condition:age_s,
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
              select(sample_id, skippy_id, condition, age_s, siblings, bfexcl) %>%
                      na.omit()

    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_y, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_y <- adonis2(asv ~ bfexcl + condition + age_s,
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

    #list(all = permanova_all, infancy = permanova_inf, year1 = permanova_y, imp = imp)
    permanova
})
permanovas_bf










###############################################################################
#########################         2. PP          ##############################
###############################################################################

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")

# obtain ids that were selected for PP analyses 
pp_indicator <- foreign::read.spss(here::here("data/raw_data/kelly141022/Data_ITT_PP_ExploratoryDRselections.sav"), to.data.frame = TRUE)
count(pp_indicator, PP, SSC)
pp_indicator <- select(pp_indicator, skippy_id = ID, pp = PP)
# add pp info to existing data 
if (!"pp" %in% colnames(d)) {
  d <- left_join(d, pp_indicator, by = "skippy_id") %>%
    mutate(
      birthweight_s = scale(birthweight)[, 1],
      ges_age_s = scale(ges_age)[, 1],
      edlevel_s = scale(edlevel)[, 1],
    )
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
  select(sample_id) %>%
  left_join(select(d_pp, skippy_id, siblings, week, age, condition, birthweight_s, ges_age_s, edlevel_s, csection, sex, sample_id, pp), by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  DataFrame()
colnames(colData(tse))
colData(tse)$age <- colData(tse)$age + as.numeric(as.character(colData(tse)$week)) * 7
colData(tse)$age_s <- scale(colData(tse)$age)[, 1]
tse_pp <- filter(tse, pp == 1)

colData(tse_pp)$csection
####################### 2.1 Complete Case Analysis ############################

### First I fit a model to all samples 

# we use Aitchison distance
tse_pp <- transformSamples(x = tse_pp, method = "clr", pseudocount = 1, name = "clr")
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_pp) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings, birthweight_s, ges_age_s,
                  edlevel_s, csection, sex) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_pp, "clr"))
asv <- asv[meta$sample_id, ]


# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     # does not work if trend is in data (therefore use 999)
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "euclidean",
                        permutations = 999)
permanova2



### Now split models by infancy and 1 year olds 

## first infancy
tse_inf <- filter(tse_pp, week != 52)
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_inf) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings, birthweight_s, ges_age_s,
                  edlevel_s, csection, sex) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_inf, "clr"))
asv <- asv[meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s * condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "euclidean",
                        permutations = 999)
permanova2




## then 1 year olds 
tse_y <- filter(tse_pp, week == 52)
# extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_y) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings, birthweight_s, ges_age_s,
                  edlevel_s, csection, sex) %>%
                  na.omit()

# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_y, "clr"))
asv <- asv[meta$sample_id, ]

# fit and inspect model 
permanova <- adonis2(asv ~ condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     permutations = 999
                   )

permanova

# Perform dbRDA
dbrda <- dbrda(asv ~ condition + age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings, data = meta)
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term analyzed individually
                        method = "euclidean",
                        permutations = 999)
permanova2

# same as above


######################## 2.2 Multiple imputation  #############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
permanovas_pp <- map2_dfr(implist_pp, 1:length(implist_pp), function(dimp, imp) {
  # Steps are repeated as from the beginning in the script above
  
  # step 1
  # we use Aitchison distance
  tse_map <- transformSamples(x = tse, method = "clr", pseudocount = 1, name = "clr")
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
              select(sample_id, skippy_id, condition, age_s, siblings, birthweight_s, ges_age_s,
                      edlevel_s, csection, sex) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_map, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~ condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999
                       )

    ## then infancy
    tse_inf <- filter(tse_map, week != 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_inf) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings, birthweight_s, ges_age_s, 
                      edlevel_s, csection, sex) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_inf, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_inf <- adonis2(asv ~ condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
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
              select(sample_id, skippy_id, condition, age_s, siblings, birthweight_s, ges_age_s,
                      edlevel_s, csection, sex) %>%
                      na.omit()

    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_y, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_y <- adonis2(asv ~ condition + age_s + condition:age_s + birthweight_s + ges_age_s + edlevel_s + csection + sex + siblings,
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
    #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
})

permanovas_pp

# same as above






##############################################################################
#################################PLOT & TABLE ################################
##############################################################################


# for the tables I use first imputed dataset of ITT analyses
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")

tse_plot <- transformSamples(x = tse, method = "clr", pseudocount = 1, name = "clr")
colData(tse_plot) <- colData(tse_plot) %>%
  as.data.frame() %>%
  select(sample_id) %>%
  left_join(
    select(
      implist[[1]], 
      condition, siblings, age,
      sample_id, week, skippy_id, bfexcl), 
    by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  mutate(
    age = age + as.numeric(as.character(week)),
    age_s = scale(age)[, 1]
  ) %>%
  DataFrame()
  


## first all samples
# extract relevant meta data and omit na as adonis doesnt accept them.
meta_plot <- colData(tse_plot) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings, bfexcl) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta_plot %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
            nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_plot, "clr"))
asv <- asv[meta_plot$sample_id, ]

# fit and inspect model 
permanova_all <- adonis2(asv ~ condition * age_s,
                      # by = "margin", # each term analyzed individually
                      data = meta_plot,
                      method = "euclidean",
                      permutations = 999
                    )
beta_table <- permanova_all %>% 
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(
    `R2 (%)` = glue("{round(R2 * 100, 2)}"),
    P = as.character(`Pr(>F)`),
    P = str_remove(P, "^0"),
    parameter = str_to_title(parameter),
    parameter = ifelse(parameter == "Condition", "SSC", ifelse(
      parameter == "Age_s", "Age", ifelse(
        parameter == "Condition:age_s", "SSC x Age", parameter
      )))
    ) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>%
  select(parameter, Df, SumOfSqs, "R2 (%)", F, P)

colnames(beta_table) <- str_to_title(colnames(beta_table))

# fit and inspect model 
permanova_all <- adonis2(
  asv ~ bfexcl + condition * age_s,
  # by = "margin", # each term analyzed individually
  data = meta_plot,
  method = "euclidean",
  permutations = 999)

beta_table_bf <- permanova_all %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(
    `R2 (%)` = glue("{round(R2 * 100, 2)}"),
    P = as.character(`Pr(>F)`),
    P = str_remove(P, "^0"),
    parameter = str_to_title(parameter),
    parameter = ifelse(parameter == "Condition", "SSC", ifelse(
      parameter == "Age_s", "Age", ifelse(
        parameter == "Condition:age_s", "SSC x Age", ifelse(
          parameter == "Bfexcl", "Breastfeeding", parameter
      ))))
    ) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>%
  select(parameter, Df, SumOfSqs, "R2 (%)", F, P)

colnames(beta_table_bf) <- str_to_title(colnames(beta_table_bf))

beta_table
beta_table_bf








# same but split by infants and year1 

# for the tables I use first imputed dataset of ITT analyses
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
tse_inf <- filter(tse, week != 52)

tse_plot <- transformSamples(x = tse_inf, method = "clr", pseudocount = 1, name = "clr")
colData(tse_plot) <- colData(tse_plot) %>%
  as.data.frame() %>%
  select(sample_id) %>%
  left_join(
    select(
      implist[[1]], 
      condition, siblings, age,
      sample_id, week, skippy_id, bfexcl), 
    by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  mutate(
    age = age + as.numeric(as.character(week)),
    age_s = scale(age)[, 1]
  ) %>%
  DataFrame()
  

# step 2

# extract relevant meta data and omit na as adonis doesnt accept them.
meta_plot <- colData(tse_plot) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings, bfexcl) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta_plot %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
            nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_plot, "clr"))
asv <- asv[meta_plot$sample_id, ]

# fit and inspect model 
permanova_inf <- adonis2(asv ~ condition * age_s,
                      # by = "margin", # each term analyzed individually
                      data = meta_plot,
                      method = "euclidean",
                      permutations = 999
                    )
beta_table_inf <- permanova_inf %>% 
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(
    `R2 (%)` = glue("{round(R2 * 100, 2)}"),
    P = as.character(`Pr(>F)`),
    P = str_remove(P, "^0"),
    parameter = str_to_title(parameter),
    parameter = ifelse(parameter == "Condition", "SSC", ifelse(
      parameter == "Age_s", "Age", ifelse(
        parameter == "Condition:age_s", "SSC x Age", parameter
      )))
    ) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>%
  select(parameter, Df, SumOfSqs, "R2 (%)", F, P)

colnames(beta_table_inf) <- str_to_title(colnames(beta_table_inf))

# fit and inspect model 
permanova_inf <- adonis2(
  asv ~ bfexcl + condition * age_s,
  # by = "margin", # each term analyzed individually
  data = meta_plot,
  method = "euclidean",
  permutations = 999)

beta_table_bf_inf <- permanova_inf %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(
    `R2 (%)` = glue("{round(R2 * 100, 2)}"),
    P = as.character(`Pr(>F)`),
    P = str_remove(P, "^0"),
    parameter = str_to_title(parameter),
    parameter = ifelse(parameter == "Condition", "SSC", ifelse(
      parameter == "Age_s", "Age", ifelse(
        parameter == "Condition:age_s", "SSC x Age", ifelse(
          parameter == "Bfexcl", "Breastfeeding", parameter
      ))))
    ) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>%
  select(parameter, Df, SumOfSqs, "R2 (%)", F, P)

colnames(beta_table_bf_inf) <- str_to_title(colnames(beta_table_bf_inf))

beta_table_inf
beta_table_bf_inf





# now lastly year 1
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
tse <- agglomerateByRank(tse, rank = "genus")
tse_year1 <- filter(tse, week == 52)

tse_plot <- transformSamples(x = tse_year1, method = "clr", pseudocount = 1, name = "clr")
colData(tse_plot) <- colData(tse_plot) %>%
  as.data.frame() %>%
  select(sample_id) %>%
  left_join(
    select(
      implist[[1]], 
      condition, siblings, age,
      sample_id, week, skippy_id, bfexcl), 
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
meta_plot <- colData(tse_plot) %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          select(sample_id, skippy_id, condition, age_s, siblings, bfexcl) %>%
                  na.omit()
# we need to account for non-independence of data in the infancy model 
ids <- meta_plot %>% mutate(skippy_id = as.factor(skippy_id)) %>%
  .$skippy_id
h <- how(plots = Plots(strata = ids, type = "none"),
            nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_plot, "clr"))
asv <- asv[meta_plot$sample_id, ]

# fit and inspect model 
permanova_year1 <- adonis2(asv ~ condition + age_s,
                      # by = "margin", # each term analyzed individually
                      data = meta_plot,
                      method = "euclidean",
                      permutations = 999
                    )
beta_table_year1 <- permanova_year1 %>% 
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(
    `R2 (%)` = glue("{round(R2 * 100, 2)}"),
    P = as.character(`Pr(>F)`),
    P = str_remove(P, "^0"),
    parameter = str_to_title(parameter),
    parameter = ifelse(parameter == "Condition", "SSC", ifelse(
      parameter == "Age_s", "Age", ifelse(
        parameter == "Condition:age_s", "SSC x Age", parameter
      )))
    ) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>%
  select(parameter, Df, SumOfSqs, "R2 (%)", F, P)

colnames(beta_table_year1) <- str_to_title(colnames(beta_table_year1))

# fit and inspect model 
permanova_year1 <- adonis2(
  asv ~ bfexcl + condition + age_s,
  # by = "margin", # each term analyzed individually
  data = meta_plot,
  method = "euclidean",
  permutations = 999)

beta_table_bf_year1 <- permanova_year1 %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(
    `R2 (%)` = glue("{round(R2 * 100, 2)}"),
    P = as.character(`Pr(>F)`),
    P = str_remove(P, "^0"),
    parameter = str_to_title(parameter),
    parameter = ifelse(parameter == "Condition", "SSC", ifelse(
      parameter == "Age_s", "Age", ifelse(
        parameter == "Condition:age_s", "SSC x Age", ifelse(
          parameter == "Bfexcl", "Breastfeeding", parameter
      ))))
    ) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>%
  select(parameter, Df, SumOfSqs, "R2 (%)", F, P)

colnames(beta_table_bf_year1) <- str_to_title(colnames(beta_table_bf_year1))

beta_table_year1
beta_table_bf_year1




load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")


pseq <- makePhyloseqFromTreeSE(tse)
pseq_clr <- microbiome::transform(pseq, transform = "clr")
sample_data(pseq_clr) <- sd_to_df(pseq_clr) %>%
  mutate(
    Infancy = ifelse(week == 52, "Late", "Early"),
    Condition = ifelse(condition == 1, "SSC", ifelse(condition == 0, "CAU", NA))
  ) %>%
  df_to_sd()


# all samples
bp <- biplot(
  pseq_clr, 
  color = "Condition", 
  point_size = 5, 
  otu_alpha = 0,
  colors = c("#909090", "#000000"),
  shape = "Infancy"
)

save(bp, beta_table, beta_table_bf, beta_table_inf, beta_table_bf_inf, beta_table_year1, beta_table_bf_year1, file = here::here("data/beta_plot_table.Rds"))


sids <- colData(tidySummarizedExperiment::filter(tse, week != "52")) %>% rownames()
bp <- biplot(
  pseq_clr, 
  color = "condition", 
  point_size = 5, 
  otu_alpha = 0,
  colors = c("#909090", "#000000"),
  filter_samples = sids
)

bp[[1]]

sids <- colData(tidySummarizedExperiment::filter(tse, week == "52")) %>% rownames()
bp <- biplot(
  pseq_clr, 
  color = "condition", 
  point_size = 5, 
  otu_alpha = 0,
  colors = c("#909090", "#000000"),
  filter_samples = sids
  )

bp[[2]]


sids <- colData(tidySummarizedExperiment::filter(tse, week != "52")) %>% rownames()
sids
filter <- dplyr::filter
bp_series <- biplot(
  pseq_clr,
  color = "condition", 
  point_size = 5,
  otu_alpha = 0,
  connect_series = "week",
  subject_id = "skippy_id",
  colors = c("#909090", "#000000"),
  filter_samples = sids
  )


bp_series[[1]]



# make also a plot with Bray Curtis

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

library(scater)
# Bray-Curtis is usually applied to relative abundances
tse_ra <- transformCounts(tse, method = "relabundance")
tse_ra <- mutate(
    tse_ra, 
    Infancy = ifelse(week == 52, "Late", "Early"),
    Condition = factor(ifelse(condition == 1, "SSC", ifelse(condition == 0, "CAU", NA)), levels = c("SSC", "CAU"))
)
# Perform PCoA
tse_ra <- runMDS(tse_ra, FUN = vegan::vegdist, method = "bray", name = "PCoA_BC", exprs_values = "relabundance")
# Create ggplot object
p <- plotReducedDim(
  tse_ra, 
  "PCoA_BC", 
  colour_by = "Condition",
  shape_by = "Infancy", 
  point_size = 5, 
  point_alpha = 1,
  theme_size = 25)
p
?"scater-plot-args"
plotReducedDim(tse_ra, "PCoA_BC", colour_by = "Condition", other_fields = list(point_size = 3))

# Add explained variance for each axis
e <- attr(reducedDim(tse, "PCoA_BC"), "eig");
rel_eig <- e/sum(e[e>0])          
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))

print(p)





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
  
  # step 1
  # we use Aitchison distance
  tse <- agglomerateByRank(tse, rank = "genus")
  tse_map <- transformSamples(x = tse, method = "clr", pseudocount = 1, name = "clr")
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
    

    # step 2

    ## first all samples
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_map) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings, SSC_s) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_map, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_all <- adonis2(asv ~ SSC_s + age_s + age_s:SSC_s,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         permutations = 999
                       )

    ## first infancy
    tse_inf <- filter(tse_map, week != 52)
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_inf) %>% as.data.frame() %>%
              rownames_to_column("sample_id") %>%
              select(sample_id, skippy_id, condition, age_s, siblings, SSC_s) %>%
                      na.omit()
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(skippy_id = as.factor(skippy_id)) %>%
      .$skippy_id
    h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_inf, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_inf <- adonis2(asv ~ SSC_s + age_s + SSC_s:age_s,
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
              select(sample_id, skippy_id, condition, age_s, siblings, SSC_s) %>%
                      na.omit()

    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_y, "clr"))
    asv <- asv[meta$sample_id, ]

    # fit and inspect model 
    permanova_y <- adonis2(asv ~ SSC_s + age_s,
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
    
    #list(infancy = permanova_inf, year1 = permanova_y, imp = imp)
    permanova
})
permanovas
