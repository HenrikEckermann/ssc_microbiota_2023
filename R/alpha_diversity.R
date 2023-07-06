# for reproducibility
set.seed(1)
# this script contains step 1 of https://osf.io/s45mu (alpha diversity)
library(mia)
library(glue)
library(tidyverse)
library(rstanarm)
library(mice)
library(brms)
library(here)
library(HDInterval)
source(here("R/helper_functions.R"))




# Data Preparation
################################################################################
if (!file.exists(here::here("data/data_imp.Rds"))) {


  # import of biomfile and meta data can be found in the import script
  load(here::here("data/data.Rds"))

  mdata$siblings <- ifelse(
    mdata$siblings == 0, 0, ifelse(
      mdata$siblings >= 1, 1, NA))
  # before estimating AD, I bring the mdata in long format where appropriate 
  l_mdata <- pivot_longer(
    mdata, 
    matches("age_week\\d+"), 
    names_to = "week",
    names_prefix = "age_week",
    values_to = "age") %>%
    mutate(sample_id = ifelse(week == "2", glue("{skippy_id}_1"), ifelse(week == "5", glue("{skippy_id}_2"), ifelse(week == "52", glue("{skippy_id}_3"), NA))))


  # obtain AD and combine dfs
  tse <- estimateDiversity(
    tse,
    assay_name = "counts",
    index = c("shannon", "faith", "inverse_simpson"),
    name = c("shannon", "faith", "inverse_simpson")
  )
  tse <- estimateRichness(
    tse,
    assay_name = "counts",
    index = "chao1",
    name = "chao1"
  )
  ad <- colData(tse) %>% 
    as.data.frame() %>%
    select(sample_id, chao1, inverse_simpson, shannon, faith)

  df <- ad %>%
    dplyr::full_join(l_mdata, by = "sample_id") %>%
    arrange(skippy_id, week)



  # variables I choose to retain for imputation based on expected predictive 
  # value for microbiota data imputation
  vars <- c(
    "sample_id",
    "skippy_id",
    "week",
    "age",
    "chao1",
    "inverse_simpson",
    "shannon",
    "faith",
    "condition",
    "csection",
    "birthweight",
    "siblings",
    "sex",
    "bfexcl",
    "bfperc_w1",
    "bfperc_w2",
    "bfperc_w3",
    "bfperc_w4",
    "bfperc_w5",
    "apgar_5",
    "ges_age",
    "edlevel",
    "parity",
    "weaning",
    "antibiotic_1year",
    "antibiotic_week2",
    "antibiotic_week5"
  )
  mdata <- dplyr::rename(
    mdata, 
    constipation_week52 = constipation_1year,
    diarrhea_week52 = diarrhea_1year,
    antibiotic_week52 = antibiotic_1year
  )



  # furthermore I include constipation and diarrhea but first i need to put 
  # them into long format
  c_and_d <- select(
    mdata,
    skippy_id, 
    matches("constipation_week\\d+$"), 
    matches("diarrhea_week\\d+$"),
    matches("antibiotic_week\\d+$")
    ) %>%
    pivot_longer(
    cols = contains("week"),
    names_to = c(".value", "week"),
    names_pattern = "(\\w+)_week(\\d+)")



  d <- dplyr::left_join(
  select(df, all_of(vars)),
  c_and_d, 
  by = c("skippy_id", "week")) %>%
  arrange(sample_id) %>%
  select(
    sample_id, skippy_id, week, age, chao1, inverse_simpson, shannon, faith, condition,
    csection, birthweight, siblings, sex, contains("bf"), constipation, antibiotic, 
    diarrhea, everything(), -contains("antibiotic_w"))

  save(d, file = here::here("data/table1.Rds"))
  # now we are ready for imputation, we start with m = 5 to make the script run 
  # but for the final estimates we will increase number of imputations
  # if skippy_id and week are numeric the imputation model function better, other-
  # wise the age var imputation is not good.
  nvars <- c("skippy_id", "week")
  d <- mutate(d, 
    across(all_of(nvars), function(x) as.numeric(x)),
    age = age - week * 7
  )
  imp <- mice(d, m = 50)

  # for the analysis we need to change dtypes for some vars
  # variables to standardize

  svars <- c(
    "age", 
    "age_dev",
    "chao1",
    "inverse_simpson",
    "faith",
    "shannon", 
    "birthweight", 
    "ges_age",
    "apgar_5",
    "edlevel"
  )
  fvars <- c(
    "week", 
    "condition", 
    "csection", 
    "siblings", 
    "sex", 
    "antibiotic"
  )
  implist <- map(1:5, function(x) {
    dtemp <- complete(imp, x)
    dtemp <- mutate(dtemp, 
      across(all_of(fvars), function(x) as.factor(x)),
      # across(all_of(svars), function(x) scale(x)[, 1]),
      age_dev = age,
      age = as.numeric(levels(week)) * 7 + age_dev,
      age_s = scale(age)[, 1],
      age_y = age/365,
      age_dev_s = scale(age_dev)[, 1],
      chao1_s = scale(chao1)[, 1],
      faith_s = scale(faith)[, 1],
      inverse_simpson_s = scale(inverse_simpson)[, 1],
      shannon_s = scale(shannon)[, 1],
      birthweight_s = scale(birthweight)[, 1],
      ges_age_s = scale(ges_age)[, 1],
      apgar_5_s = scale(apgar_5)[, 1],
      edlevel_s = scale(edlevel)[, 1],
      skippy_id = as.integer(skippy_id)
    )
    return(dtemp)
  })

  save(d, implist, file = here::here("data/data_imp.Rds"))
 } else {
  load(file = here::here("data/data_imp.Rds"))
}



# Figure out model structure for the ITT analyses 
################################################################################


# I decide which time variable we will use based on LOO 
f1 <- bf(shannon_s ~ week * condition + (1|skippy_id))
m1 <- brm_multiple(
  data = implist,
  formula = f1,
  file = here::here("data/m1.Rds")
)
loo_m1 <- add_criterion(
  m1,
  "loo",
  file = here::here("data/loo_m1"),
  moment_match = FALSE
)
loo_m1

f2 <- bf(shannon_s ~ age_y * condition + (1|skippy_id))
m2 <- brm_multiple(
  data = implist,
  formula = f2,
  file = here::here("data/m2.Rds")
)

loo_m2 <- add_criterion(
  m2,
  "loo",
  file = here::here("data/loo_m2"),
  moment_match = FALSE
)
f3 <- bf(shannon_s ~ week * condition + age_dev + (1|skippy_id))
m3 <- brm_multiple(
  data = implist,
  formula = f3,
  file = here::here("data/m3.Rds")
)
loo_m3 <- add_criterion(
  m3,
  "loo",
  file = here::here("data/loo_m3"),
  moment_match = FALSE
)


loo_comp <- loo_compare(loo_m1, loo_m2, loo_m3)
loo_comp

# the model with age has best fit. Therefore, we will use these as our main
# models. We will double check if results differ as compared to using week 
# But especially for the frequentist models this setup makes it easier as well.


# first we must decide which covariates we want to include. For ITT we can 
# decide totally based on model fit as long as it does not interfere with 
# causal chain of our DAG. E.g. we must exclude BF although it prob help to 
# predict AD as SSC --> BF and we want ACE of SSC. However, variables such as 
# age, csection, siblings etc. can be included based on whether they improve 
# model fit. We can use loo for this to get a quick answer: 


# I consider the following variables for selection 
coefs <- c(
  "csection", 
  "birthweight_s", 
  "siblings",
  "sex",
  "apgar_5_s",
  "ges_age_s",
  "edlevel_s"
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
    map_dfr(list(c(2, 5), 52), function(n_week) {
      if (52 %in% n_week) {
        # base model 
        f1 <- bf(shannon_s ~ age * condition)
        m1 <- brm(
          data = filter(implist[[i]], week == 52) %>%
                    mutate(age = scale(age)[, 1]),
          formula = f1,
          file = here::here(glue("data/m1_age_imp{i}_{n_week[1]}.Rds"))
        )
        loo_m1 <- add_criterion(
          m1,
          "loo",
          file = here::here(glue("data/loo_m1_age_imp{i}_{n_week[1]}")),
          moment_match = FALSE
        )

        f2 <- bf(glue("shannon_s ~ age * condition + {coef}"))
        m2 <- brm(
          data = filter(implist[[i]], week == 52) %>%
                    mutate(age = scale(age)[, 1]),
          formula = f2,
          file = here::here(glue("data/m2_age_imp{i}_{coef}_{n_week[1]}.Rds"))
        )

        loo_m2 <- add_criterion(
          m2,
          "loo",
          file = here::here(glue("data/loo_m2_age_imp{i}_{coef}_{n_week[1]}")),
          moment_match = FALSE
        )
        lcomp <- loo_compare(loo_m2, loo_m1)
        score <- ifelse(rownames(lcomp)[1] == "loo_m1", 0, 1)
      } else {
        # base model 
        f1 <- bf(shannon_s ~ age * condition + (1|skippy_id))
        m1 <- brm(
          data = filter(implist[[i]], week != 52) %>%
                    mutate(age = scale(age)[, 1]),
          formula = f1,
          file = here::here(glue("data/m1_age_imp{i}_{n_week[1]}.Rds"))
        )
        loo_m1 <- add_criterion(
          m1,
          "loo",
          file = here::here(glue("data/loo_m1_age_imp{i}_{n_week[1]}")),
          moment_match = FALSE
        )

        f2 <- bf(glue("shannon_s ~ age * condition + {coef} + (1|skippy_id)"))
        m2 <- brm(
          data = filter(implist[[i]], week != 52) %>%
                    mutate(age = scale(age)[, 1]),
          formula = f2,
          file = here::here(glue("data/m2_age_imp{i}_{coef}_{n_week[1]}.Rds"))
        )

        loo_m2 <- add_criterion(
          m2,
          "loo",
          file = here::here(glue("data/loo_m2_age_imp{i}_{coef}_{n_week[1]}")),
          moment_match = FALSE
        )
        lcomp <- loo_compare(loo_m2, loo_m1)
        score <- ifelse(rownames(lcomp)[1] == "loo_m1", 0, 1)
      }
      tibble(
        model = ifelse(n_week == 52, "1year", "2 and 5 weeks"),
        coef = coef,
        imp = i,
        score = score
      )      
    })
  })
})

# siblings was the only covariate that improved out of sample 
# predictions. Therefore, we will only keep this variable
group_by(loo_comp, coef, model) %>%
  summarise(ss = sum(score))

# lastly: should we include random effects for age?
loo_comp <- map_dfr(1:5, function(i) {
  # base model 
  f1 <- bf(shannon_s ~ age * condition + (1|skippy_id))
  m1 <- brm(
    data = implist[[i]],
    formula = f1,
    file = here::here(glue("data/m1_age_imp{i}.Rds"))
  )
  loo_m1 <- add_criterion(
    m1,
    "loo",
    file = here::here(glue("data/loo_m1_age_imp{i}")),
    moment_match = FALSE
  )

  f2 <- bf(glue("shannon_s ~ age * condition + (1 + age|skippy_id)"))
  m2 <- brm(
    data = implist[[i]],
    formula = f2,
    file = here::here(glue("data/m2_age_imp{i}_rftime.Rds"))
  )

  loo_m2 <- add_criterion(
    m2,
    "loo",
    file = here::here(glue("data/loo_m2_age_imp{i}_rftime")),
    moment_match = FALSE
  )
  lcomp <- loo_compare(loo_m2, loo_m1)
  score <- ifelse(rownames(lcomp)[1] == "loo_m1", 0, 1)
  tibble(
    imp = i,
    score = score
  )
})
# results indicate that varying time effects arent necessary
loo_comp





# ITT analyses with MI 
#############################################################################


# start with posterior predictive checks and model diagnostics before we inter 
# pret the results. In the end I will use the final models for interpretation.
# because some loo estimates indicate outliers I will use robust regression.
# Shannon is least prone to non-normality and yet we got some warnings while 
# there is no disadvantage of using robust regression:

# "inverse_simpson_s" cannot be modelled with student()
indeces <- c("shannon_s", "chao1_s", "faith_s")
ad_models <- map(indeces, function(index) {
  f <- bf(glue("{index} ~ age_s * condition + siblings + (1|skippy_id)"))
  m <- brm_multiple(
    family = student(),
    data = implist,
    formula = f,
    file = here::here(glue("data/m_age_{index}_bez.Rds"))
  )
})
ad_models
save(ad_models, file = here::here("data/ad_models.Rds"))
# now for each model we perform posterior predictive checks and have a look 
# at residuals 
pp_checks <- map(ad_models, function(m) {
  p <- pp_check(m)
})
pp_checks
# load helper function to diagnose lms 
source(here::here("R/ml_helper.R"))
names(ad_models) <- indeces
lm_diags <- map(indeces, function(index) {
  lm_diag(ad_models[[index]], ad_models[[index]]$data, index, id = "skippy_id")
})
lm_diags

# check also models split by time points (this only makes really a difference for siblings effect)
ad_models2 <- map(c("2and5", "52"), function(weekchr) {
  if (weekchr == "2and5") {
    ad_models <- map(indeces, function(index) {
      f <- bf(glue("{index} ~ age * condition + siblings + (1|skippy_id)"))
      m <- brm_multiple(
        family = student(),
        data = map(implist, ~filter(.x, week != 52)),
        formula = f,
        file = here::here(glue("data/m_age_{index}_{weekchr}_bez.Rds"))
      )
    })
  } else {
    ad_models <- map(indeces, function(index) {
      f <- bf(glue("{index} ~ condition + age + siblings"))
      m <- brm_multiple(
        family = student(),
        data = map(implist, ~filter(.x, week == 52)),
        formula = f,
        file = here::here(glue("data/m_{index}_{weekchr}_bez.Rds"))
      )
    })
  }
  ad_models
})
ad_models2
post <- posterior_samples(ad_models2[[1]][[1]])
mean(post$b_siblings1 < 0)
length(unlist(ad_models2, recursive = FALSE))
# now for each model we perform posterior predictive checks and have a look 
# at residuals 
ad_models2 <- unlist(ad_models2, recursive = FALSE)

length(ad_models2)
pp_checks2 <- map(ad_models2, function(m) {
  p <- pp_check(m)
})
pp_checks2


lm_diags2 <- map(c("2and5", "52"), function(weekchr) {
  if (weekchr == "2and5") {
    lm_diags <- map(indeces, function(index) {
      lm_diag(
        ad_models2[1:4][[index]], 
        filter(implist[[1]], week != 52), 
        index, id = "skippy_id")
    })
  } else {
    lm_diags <- map(indeces, function(index) {
      lm_diag(
        ad_models2[5:8][[index]], 
        filter(implist[[1]], week== 52), 
        index, id = "skippy_id")
    })
  }
  lm_diags
})
lm_diags2[[1]][[1]]
lm_diags2[[2]][[1]]




# create a plot for AD 
if (!is.factor(d$week)) {
  d <-mutate(d, 
        week = as.factor(week), condition = as.factor(condition),
        group = ifelse(week == 2 & condition == 0, "W2 CAU", ifelse(
          week == 2 & condition == 1, "W2 SSC", ifelse(
            week == 5 & condition == 0, "W5 CAU", ifelse(
              week == 5 & condition == 1, "W5 SSC", ifelse(
                week == 52 & condition == 0, "W52 CAU", ifelse(
                  week == 52 & condition == 1, "W52 SSC", NA)))))),
        group = as.factor(group),
        condition_label = ifelse(condition == 0, "CAU", ifelse(
          condition == 1, "SSC", NA)),
        week_label = glue::glue("Week {week}")
        )
}
  d <-mutate(d, 
        week = as.factor(week), condition = as.factor(condition),
        group = ifelse(week == 2 & condition == 0, "W2 CAU", ifelse(
          week == 2 & condition == 1, "W2 SSC", ifelse(
            week == 5 & condition == 0, "W5 CAU", ifelse(
              week == 5 & condition == 1, "W5 SSC", ifelse(
                week == 52 & condition == 0, "W52 CAU", ifelse(
                  week == 52 & condition == 1, "W52 SSC", NA)))))),
        group = as.factor(group),
        condition_label = ifelse(condition == 0, "CAU", ifelse(
          condition == 1, "SSC", NA)),
        week_label = glue::glue("Week {week}")
        )
# alternative 1
adplots <- map(str_remove(indeces, "_s$"), function(index) {
  d %>% 
    ggplot(aes_string("group", index, fill = "condition")) +
      geom_boxplot() +
      #ggbeeswarm::geom_beeswarm(size = 3, cex = 1) +
      geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
      # scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
      scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
      theme_bw(base_size = 25) +
      theme(legend.position = "none") +
      xlab("") + ylab(str_to_title(index))
})
# alternative 2
adplots <- map(str_remove(indeces, "_s$"), function(index) {
  d %>%
    ggplot(aes_string("condition_label", index, fill = "condition_label")) +
      geom_boxplot() +
      #ggbeeswarm::geom_beeswarm(size = 3, cex = 1) +
      geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
      # scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
      scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
      facet_wrap(~week_label, strip.position = "bottom") +
      #scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
      theme_bw(base_size = 25) +
      theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank()) +
      xlab("") + ylab(str_to_title(index))
})



save(adplots, file = here::here("data/adplots.Rds"))




# ITT analyses complete case 
#############################################################################
if (!"shannon_s" %in% colnames(d)) {
  d <- mutate(d, 
   age_s = scale(age)[, 1],
   chao1_s = scale(chao1)[, 1],
   faith_s = scale(faith)[, 1],
   inverse_simpson_s = scale(inverse_simpson)[, 1],
   shannon_s = scale(shannon)[, 1],
   skippy_id = as.integer(skippy_id),
   siblings = as.factor(siblings)
  )
}


ad_models_cc <- map(c("2and5", "52"), function(weekchr) {
  if (weekchr == "2and5") {
    ad_models <- map(indeces, function(index) {
      f <- bf(glue("{index} ~ age * condition + siblings + (1|skippy_id)"))
      m <- brm(
        family = student(),
        data = filter(d, week != 52) %>% 
                              mutate(age = scale(age)[, 1]),
        formula = f,
        file = here::here(glue("data/m_age_{index}_{weekchr}_cc_bez.Rds"))
      )
    })
  } else {
    ad_models <- map(indeces, function(index) {
      f <- bf(glue("{index} ~ condition + age + siblings"))
      m <- brm(
        family = studetnt(),
        data = filter(d, week == 52) %>% 
                              mutate(age = scale(age)[, 1]),
        formula = f,
        file = here::here(glue("data/m_age_{index}_{weekchr}_cc_bez.Rds"))
      )
    })
  }
  ad_models
})
parameters <- c(
  "b_age", 
  "b_condition1",
  "b_siblings1",
  "b_age:condition1"
)



tbs_cc <- map2_dfr(c("2and5", "52"), ad_models_cc, function(weekchr, models) {
  if (weekchr == "2and5") {
    tb <- map2_dfr(indeces, models, function(index, model) {
      summarise_posterior(model, parameters, 2) %>%
      mutate(model = weekchr, index = index, indicator = (lower <=0 & upper <=0) | (lower >=0 & upper >=0))
    })
  } else {
    tb <- map2_dfr(indeces, models, function(index, model) {
      summarise_posterior(model, parameters[-length(parameters)], 2) %>%
      mutate(model = weekchr, index = index, indicator = (lower <=0 & upper <=0) | (lower >=0 & upper >=0))
    })
  }
})
filter(tbs_cc, indicator)

filter(tbs_cc, indicator)



# results do not differ meaningully between complete case analyses and mi
# for inverse simpson the models are not a good fit in either case. I ran 
# wilcoxon tests here:
filter(d, group %in% c("W52 CAU", "W52 SSC")) %>%
  mutate(group = fct_drop(group)) %>%
  wilcox.test(inverse_simpson ~ group, data = .)




# PP analyses with MI 
#############################################################################

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


covariates <- c(
  "siblings", 
  "birthweight_s", 
  "ges_age_s", 
  "edlevel_s", 
  "csection", 
  "sex"
)

model_str <- "age_s * condition"
for (coef in covariates) {
  model_str <- glue("{model_str} + {coef}")
}

# "inverse_simpson_s" cannot be modelled with student()
indeces <- c("shannon_s", "chao1_s", "faith_s")
ad_models_pp <- map(indeces, function(index) {
  f <- bf(glue("{index} ~ {model_str} + (1|skippy_id)"))
  m <- brm_multiple(
    family = student(),
    data = implist_pp,
    formula = f,
    file = here::here(glue("data/m_age_{index}_pp_bez.Rds"))
  )
})
ad_models_pp

indeces <- c("shannon_s", "inverse_simpson_s", "chao1_s", "faith_s")
ad_models_pp <- map(c("2and5", "52"), function(weekchr) {
  if (weekchr == "2and5") {
    ad_models <- map(indeces, function(index) {
      f <- bf(glue("{index} ~ {model_str} + (1|skippy_id)"))
      m <- brm_multiple(
        family = student(),
        data = map(implist_pp, ~filter(.x, week != 52) %>% 
                              mutate(age = scale(age)[, 1])),
        formula = f,
        file = here::here(glue("data/m_age_{index}_{weekchr}_sn_pp_bez.Rds"))
      )
    })
  } else {
    ad_models <- map(indeces, function(index) {
      f <- bf(glue("{index} ~ {model_str}"))
      if (index %in% c("chao1_s", "inverse_simpon_s")) {
        family <- skew_normal()
      } else {
        family <- student()
      }
      m <- brm_multiple(
        family = family,
        data = map(implist_pp, ~filter(.x, week == 52) %>% 
                              mutate(age = scale(age)[, 1])),
        formula = f,
        file = here::here(glue("data/m_age_{index}_{weekchr}_sn_pp_bez.Rds"))
      )
    })
  }
  ad_models
})



pp_checks_pp <- map(ad_models_pp, function(m) {
  p <- pp_check(m)
})
pp_checks_pp

tbs_pp <- map2_dfr(c("2and5", "52"), ad_models_pp, function(weekchr, models) {
  if (weekchr == "2and5") {
    tb <- map2_dfr(indeces, models, function(index, model) {
      summarise_posterior(model, parameters, 2) %>%
      mutate(model = weekchr, index = index, indicator = (lower <=0 & upper <=0) | (lower >=0 & upper >=0))
    })
  } else {
    tb <- map2_dfr(indeces, models, function(index, model) {
      summarise_posterior(model, parameters[-length(parameters)], 2) %>%
      mutate(model = weekchr, index = index, indicator = (lower <=0 & upper <=0) | (lower >=0 & upper >=0))
    })
  }
})
filter(tbs_pp, indicator)


# results are similar to the ITT analyses. No effects on alpha diversity 


# create a supplementary table with all beta coefficients

ad_stbl <- map2_dfr(c("2and5", "52"), ad_models_pp, function(weekchr, models) {
  if (weekchr == "2and5") {
    tb <- map2_dfr(indeces, models, function(index, model) {
      summarise_posterior(model, contains("b_"), 3) %>%
      mutate(model = weekchr, index = index, indicator = (lower <=0 & upper <=0) | (lower >=0 & upper >=0))
    })
  } else {
    tb <- map2_dfr(indeces, models, function(index, model) {
      summarise_posterior(model, contains("b_"), 3) %>%
      mutate(model = weekchr, index = index, indicator = (lower <=0 & upper <=0) | (lower >=0 & upper >=0))
    })
  }
})
ad_stbl <- group_by(ad_stbl, model, index) %>% nest()

ad_stbl_pp <- pmap(list(ad_stbl[[1]], ad_stbl[[2]], ad_stbl[[3]]), function(time, index, tbl){
  time_t <- ifelse(time == "2and5", "week 2 and 5", "1 year")
  index_t <- str_to_title(str_remove(index, "_s$"))
  caption <- glue("Alpha diversity model using samples obtained at {time_t} and {index_t}")
  knitr::kable(tbl, caption = caption)
})

ad_stbl_pp










