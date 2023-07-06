library(tidyverse)
library(mia)
library(glue)
library(vegan) 
library(brms)
library(bayesplot)
library(posterior)
library(tidybayes)

# import of biomfile and meta data can be found in the import script
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

tse <- transformSamples(tse, method = "clr", name = "clr", pseudocount = 1)
asv <- t(assay(tse, "clr"))
ait <- vegdist(asv, method = "euclidean")
ait <- as.matrix(ait)
ids <- sort(unique(d$skippy_id))
vol <- map_dfr(ids, function(id) {
  
  # each subject has maximal 3 samples:
  s1 <- glue("{id}_1")
  s2 <- glue("{id}_2")
  s3 <- glue("{id}_3")
  
  # we can only obtain vol if there is no missing sample in a pair 
  vol1 <- ifelse(s1 %in% rownames(ait) & s2 %in% rownames(ait), ait[s1, s2], NA)
  vol2 <- ifelse(s2 %in% rownames(ait) & s3 %in% rownames(ait), ait[s2, s3], NA)
  
  tibble(
    skippy_id = id,
    time = c("2-5", "5-52"),
    vol = c(vol1, vol2)
  )
})


vol_by_comp <- group_by(vol, time) %>% nest()
# visualize distributions of volatility per time point pair 
voldist <- map(vol_by_comp[[2]], function(df) {
  ggplot(df, aes(vol)) +
    geom_density()
})
voldist[[2]]
vol1 <- vol_by_comp[[2]][[1]]
vol2 <- vol_by_comp[[2]][[2]]
colnames(vol1) <- c("skippy_id", "vol1")
colnames(vol2) <- c("skippy_id", "vol2")




###############################################################################
#########################           1. ITT       ##############################
###############################################################################




####################### 1.1 Complete Case Analysis ############################

d <- dplyr::left_join(d, vol1, by = "skippy_id") %>%
      dplyr::left_join(vol2, by = "skippy_id") %>%
      filter(week == 2) %>%
      mutate(across(contains("vol"), function(x) scale(x)[, 1]))
fvars <- c("siblings", "condition")
d <- mutate(
  d, 
  across(all_of(fvars), function(x) as.factor(x)),
  age = week * 7 + age,
  ges_age_s = scale(ges_age)[, 1],
  edlevel_s = scale(edlevel)[, 1],
  birthweight_s = scale(birthweight)[, 1],
  condition_label = ifelse(condition == 0, "CAU", ifelse(
    condition == 1, "SSC", NA))
)

implist_vol <- map(implist, function(dimp) {
  df <- dplyr::left_join(dimp, vol1, by = "skippy_id") %>%
        dplyr::left_join(vol2, by = "skippy_id") %>%
        filter(week == 2) %>%
        mutate(across(contains("vol"), function(x) scale(x)[, 1]))
  # the vol columns will be imputed
  imp <- mice::mice(df, m = 1)
  complete(imp)
})
dim(implist_vol[[1]])

# make a plot
d %>%
  pivot_longer(contains("vol"), names_to = "Time", values_to = "Volatility") %>%
  mutate(Time = ifelse(Time == "vol1", "2-5 weeks", "5-52 weeks")) %>%
  ggplot(aes(condition_label, Volatility, fill = condition_label)) +
    geom_boxplot(outlier.alpha = 0) +
    # ggbeeswarm::geom_beeswarm(alpha = 0.4) +
    geom_jitter(width = 0.1) +
    #stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red") +
    scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
    #scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
    theme_bw(base_size = 25) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.background = element_blank()) +
    xlab("")


vol_plot <- d %>%
  pivot_longer(contains("vol"), names_to = "Time", values_to = "Volatility") %>%
  mutate(Time = ifelse(Time == "vol1", "2-5 weeks", "5-52 weeks")) %>%
  ggplot(aes(condition_label, Volatility, fill = condition_label)) +
    geom_boxplot(outlier.alpha = 0) +
    # ggbeeswarm::geom_beeswarm(alpha = 0.4) +
    geom_jitter(width = 0.1, size = 2) + 
    #stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red") +
    facet_wrap(~Time, strip.position = "bottom") +
    scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
    #scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
    theme_bw(base_size = 25) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.background = element_blank()) +
    xlab("") 
vol_plot

save(vol_plot, file = here::here("data/vol_out.Rds"))






# determine model structure


# I consider the following variables for selection 
coefs <- c(
  "csection", 
  "birthweight_s", 
  "birthweight_s + ges_age_s",
  "siblings",
  "sex",
  "apgar_5_s",
  "ges_age_s",
  "edlevel_s",
  "age_s"
)


# I will use the following algorithm:
# For each dataset in implist:
  # Calculate LOO for base model
  # Then for each var in coefs:
    # calculate LOO for base model + coef 
    # If LOO indicates the new model is a better fit
    # keep that coef in list 

loo_comp <- map_dfr(coefs, function(coef) {
  map_dfr(1:5, function(i) {
    # base model 
    f1 <- bf(Volatility ~ condition + (1|skippy_id))
    m1 <- brm(
      data = pivot_longer(implist_vol[[i]], contains("vol"), names_to = "Time", values_to = "Volatility"),
      formula = f1,
      file = here::here(glue("data/m1_imp{i}_vol.Rds"))
    )
    loo_m1 <- add_criterion(
      m1,
      "loo",
      file = here::here(glue("data/loo_m1_imp{i}_vol")),
      moment_match = FALSE
    )

    f2 <- bf(glue("Volatility ~ condition + {coef} + (1|skippy_id)"))
    m2 <- brm(
      data = pivot_longer(implist_vol[[i]], contains("vol"), names_to = "Time", values_to = "Volatility"),
      formula = f2,
      file = here::here(glue("data/m2_imp{i}_{coef}_vol.Rds"))
    )

    loo_m2 <- add_criterion(
      m2,
      "loo",
      file = here::here(glue("data/loo_m2_imp{i}_{coef}_vol")),
      moment_match = FALSE
    )
    lcomp <- loo_compare(loo_m2, loo_m1)
    score <- ifelse(rownames(lcomp)[1] == "loo_m1", 0, 1)

  tibble(
    coef = coef,
    imp = i,
    score = score
  )
  }) 
})

group_by(loo_comp, coef) %>%
  summarise(ss = sum(score))



model1 <- brm(
  family = student(),
  formula = vol1 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = d,
  file = here::here("data/vol1_cc.Rds")
)
summary(model1)
hypothesis(model1, "condition1 < 0")




model1_bez <- brm(
  family = student(),
  formula = vol1 ~ condition,
  data = d,
  file = here::here("data/vol1_cc_bez.Rds")
)
summary(model1_bez)
hypothesis(model1_bez, "condition1 < 0")



model2 <- brm(
  family = student(),
  formula = vol2 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = d,
  file = here::here("data/vol2_cc.Rds")
)

summary(model2)
hypothesis(model2, "condition1 < 0")



mlm <- brm(
  family = student(),
  formula = Volatility ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + (1 | skippy_id),
  data = pivot_longer(d, contains("vol"), names_to = "Time", values_to = "Volatility"),
  file = here::here("data/vol1_cc_mlm.Rds")
)
summary(mlm)
hypothesis(mlm, "condition1 < Intercept")



################### 1.2 Multiple Imputation Analysis ##########################



model1 <- brm_multiple(
  family = student(),
  formula = vol1 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = implist_vol,
  file = here::here("data/vol1_imp.Rds")
)
summary(model1)
hypothesis(model1, "condition1 < 0")

hypothesis(model1, "ges_age_s < 0")
post <- posterior_samples(model1)
mean(post$b_ges_age_s<0)


model1_bez <- brm_multiple(
  family = student(),
  formula = vol1 ~ condition,
  data = implist_vol,
  file = here::here("data/vol1_imp_bez.Rds")
)
summary(model1_bez)
hypothesis(model1_bez, "condition1 < 0")


post <- posterior_samples(model1)
mean(post$b_condition1<0)

model1_bf <- brm_multiple(
  family = student(),
  formula = vol1 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + bfexcl,
  data = implist_vol,
  file = here::here("data/vol1_imp_bf.Rds")
)

summary(model1_bf)
hypothesis(model1_bf, "condition1 < 0")

model2 <- brm_multiple(
  family = student(),
  formula = vol2 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = implist_vol,
  file = here::here("data/vol2_imp.Rds")
)

summary(model2)
hypothesis(model2, "condition1 < 0")

model2_bf <- brm_multiple(
  family = student(),
  formula = vol2 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + bfexcl,
  data = implist_vol,
  file = here::here("data/vol2_imp_bf.Rds")
)

summary(model2_bf)
hypothesis(model2_bf, "bfexcl < 0")
post <- posterior_samples(model2_bf)
mean(post$b_bfexcl<0)


mlm <- brm_multiple(
  family = student(),
  formula = Volatility ~ condition + ges_age_s + birthweight_s + csection + edlevel_s + (1 | skippy_id),
  data = map(implist_vol, ~pivot_longer(.x, contains("vol"), names_to = "Time", values_to = "Volatility")),
  file = here::here("data/vol1_imp_mlm.Rds")
)

summary(mlm)
hypothesis(mlm, "condition1 < 0")
# file.remove(here::here("data/vol1_imp_mlm.Rds"))

mlm2 <- brm_multiple(
  family = student(),
  formula = Volatility ~ condition + (1 | skippy_id),
  data = map(implist_vol, ~pivot_longer(.x, contains("vol"), names_to = "Time", values_to = "Volatility")),
  file = here::here("data/vol1_imp_mlm2.Rds")
)

summary(mlm2)
hypothesis(mlm2, "condition1 < 0")



save(model1, model1_bez, model1_bf, model2, model2_bf, file = here::here("data/volmodels_itt.Rds"))






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
  d <- dplyr::left_join(d, pp_indicator, by = "skippy_id")
}



####################### 2.1 Complete Case Analysis ############################

d <- dplyr::left_join(d, vol1, by = "skippy_id") %>%
      dplyr::left_join(vol2, by = "skippy_id") %>%
      filter(week == 2) %>%
      mutate(across(contains("vol"), function(x) scale(x)[, 1]))
fvars <- c("siblings", "condition", "csection", "sex")
d <- mutate(
  d, 
  across(all_of(fvars), function(x) as.factor(x)),
  age = week * 7 + age,
  ges_age_s = scale(ges_age)[, 1],
  edlevel_s = scale(edlevel)[, 1],
  birthweight_s = scale(birthweight)[, 1]) %>% 
  filter(pp == 1)

model1 <- brm(
  family = student(),
  formula = vol1 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = d,
  file = here::here("data/vol1_cc_pp.Rds")
)
summary(model1)
hypothesis(model1, "condition1 < 0")



model2 <- brm(
  family = student(),
  formula = vol2 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = d,
  file = here::here("data/vol2_cc_pp.Rds")
)

summary(model2)
hypothesis(model2, "condition1 < 0")



mlm <- brm(
  family = student(),
  formula = Volatility ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + (1 | skippy_id),
  data = pivot_longer(d, contains("vol"), names_to = "Time", values_to = "Volatility"),
  file = here::here("data/vol1_cc_mlm_pp.Rds")
)
summary(mlm)
hypothesis(mlm, "condition1 < Intercept")


################### 2.2 Multiple Imputation Analysis ##########################

implist_pp <- map(implist, function(dimp) {
  dimp_new <- dplyr::left_join(dimp, pp_indicator, by = "skippy_id") %>%
                filter(pp == 1)
  dimp_new
})

implist_vol <- map(implist_pp, function(dimp) {
  df <- dplyr::left_join(dimp, vol1, by = "skippy_id") %>%
        dplyr::left_join(vol2, by = "skippy_id") %>%
        filter(week == 2, pp == 1) %>%
        mutate(across(contains("vol"), function(x) scale(x)[, 1]))
  # the vol columns will be imputed
  imp <- mice::mice(df, m = 1)
  complete(imp)
})

model1 <- brm_multiple(
  family = student(),
  formula = vol1 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = implist_vol,
  file = here::here("data/vol1_imp_pp.Rds")
)


dim(implist_vol[[1]])
summary(model1)
hypothesis(model1, "condition1 < 0")


model1_bf <- brm_multiple(
  family = student(),
  formula = vol1 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + bfexcl,
  data = implist_vol,
  file = here::here("data/vol1_imp_bf_pp.Rds")
)

summary(model1_bf)
hypothesis(model1_bf, "condition1 < 0")

model2 <- brm_multiple(
  family = student(),
  formula = vol2 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection,
  data = implist_vol,
  file = here::here("data/vol2_imp_pp.Rds")
)

summary(model2)
hypothesis(model2, "condition1 < 0")

model2_bf <- brm_multiple(
  family = student(),
  formula = vol2 ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + bfexcl,
  data = implist_vol,
  file = here::here("data/vol2_imp_bf_pp.Rds")
)

summary(model2_bf)


mlm <- brm_multiple(
  family = student(),
  formula = Volatility ~ condition + ges_age_s + birthweight_s + edlevel_s + csection + (1 | skippy_id),
  data = map(implist_vol, ~pivot_longer(.x, contains("vol"), names_to = "Time", values_to = "Volatility")),
  file = here::here("data/vol1_imp_mlm_pp.Rds")
)

summary(mlm)
hypothesis(mlm, "condition1 < 0")

save(model1, model1_bf, model2, model2_bf, file = here::here("data/volmodels_pp.Rds"))

mlm2 <- brm_multiple(
  family = student(),
  formula = Volatility ~ condition + (1 | skippy_id),
  data = map(implist_vol, ~pivot_longer(.x, contains("vol"), names_to = "Time", values_to = "Volatility")),
  file = here::here("data/vol1_imp_mlm2_pp.Rds")
)

summary(mlm2)
hypothesis(mlm2, "condition1 < 0")








# we can argue that there is evidence for an effect of SSC on MBA
# therefore we will look at the DR analysis as well here:

dr <- foreign::read.spss(
  here::here("data/kelly_documents/data_itt_pp_dr.sav"),
  to.data.frame = TRUE
  ) %>%
  select(skippy_id = ID, ITT, SSC = TotalSSCwk1wk5) %>%
  mutate(SSC_s = scale(SSC)[, 1])


implist_vol <- map(implist, function(dimp) {
  df <- dplyr::left_join(dimp, vol1, by = "skippy_id") %>%
        dplyr::left_join(vol2, by = "skippy_id") %>%
        filter(week == 2) %>%
        mutate(across(contains("vol"), function(x) scale(x)[, 1])) %>%
        dplyr::left_join(select(dr, skippy_id, SSC_s), by = "skippy_id")
  # the vol columns will be imputed
  imp <- mice::mice(df, m = 1)
  complete(imp)
})
colnames(implist_vol[[1]])

model1 <- brm_multiple(
  family = student(),
  formula = vol1 ~ SSC_s + ges_age_s + birthweight_s + edlevel_s + csection,
  data = implist_vol,
  file = here::here("data/vol1_imp_dr.Rds")
)

summary(model1)
hypothesis(model1, "SSC_s < 0")


model2 <- brm_multiple(
  family = student(),
  formula = vol2 ~ SSC_s + ges_age_s + birthweight_s + edlevel_s + csection,
  data = implist_vol,
  file = here::here("data/vol2_imp_dr.Rds")
)


summary(model2)
dim(model2$data)


# only within SSC

formula <- bf(vol1 ~ SSC_s + ges_age_s + birthweight_s + edlevel_s + csection)
model_vol1_within <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_vol, ~filter(.x, condition == 1)),
  file = here::here("data/vol1_dr_within_cov")
)

save(model_vol1_within, file = here::here("data/volmodels_dr.Rds"))

summary(model_vol1_within)

formula <- bf(vol1 ~ SSC_s + ges_age_s + birthweight_s + edlevel_s + csection)
model_vol1_out <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_vol, ~filter(.x, condition == 0)),
  file = here::here("data/vol1_dr_out_cov")
)

summary(model_vol1_out)


implist_vol <- map(implist_vol, ~mutate(.x, upper = ifelse(SSC_s <= -0.5, "no", "yes")))

colnames(implist_vol[[1]])
formula <- formula <- bf(vol1 ~ upper + ges_age_s + birthweight_s + edlevel_s + csection)
model_year1_cat <- brm_multiple(
  family = student(),
  formula = formula,
  data = implist_vol,
  file = here::here("data/vol1_dr_cat_cov")
)

summary(model_year1_cat)
hypothesis(model_year1_cat, "upperyes < 0")
dim(implist_vol[[1]])

