library(glue)
library(mia)
library(readxl)
library(lubridate)
library(here)
library(ranger)
library(tidyverse)
library(brms)
library(ggbeeswarm)
library(ComplexHeatmap)
library(circlize)

full_join <- dplyr::full_join

# double check by using two approaches.
# APPROACH 1: combine to one TSE, then extract data:
for (file in list.files(here("data/"), pattern = "tse_\\w+.Rds")) {
  load(here(glue("data/{file}")))
}

lnsk <- length(rownames(tse_skippy))
lnbib <- length(rownames(tse_bibo))
lnbin1 <- length(rownames(tse_bingo1))
lnbin2 <- length(rownames(tse_bingo2))
rnskippy <- glue("x{1:lnsk}")
rnbibo <- glue("x{(lnsk + 1): (lnsk +lnbib)}")
rnbingo1 <- glue("x{(lnsk + lnbib + 1): (lnsk +lnbib + lnbin1)}")
rnbingo2 <- glue("x{(lnsk + lnbib + lnbin1 + 1): (lnsk +lnbib + lnbin1 + 
  lnbin2)}")

tses <- map2(
  list(rnskippy, rnbibo, rnbingo1, rnbingo2),
  list(tse_skippy, tse_bibo, tse_bingo1, tse_bingo2),
  function(rn, tse) {
    rownames(tse) <- rn
    tse
  }
)

tses[[1]] %>% colnames()
skbib <- mergeSEs(tses[[1]], tses[[2]], missing_values = 0)
bin <- mergeSEs(tses[[3]], tses[[4]], missing_values = 0)
assay(bin)
tse <- mergeSEs(skbib, bin, missing_values = 0)
rowData(tse) %>% dim()
atest <- assay(tse)

plot(density(colSums(atest)))
# now we need to replace "g__", "f__" etc. by NA in order to aggregate correctly
rowData(tse) <- rowData(tse) %>%
  as_tibble() %>%
  mutate(across(everything(), function(x) str_remove(x, ".*[kpcofgs]__"))) %>%
  DataFrame()

test <- agglomerateByRank(tse, rank = "Genus")
plot(density(colSums(assay(tse))))
plot(density(assay(test)[1, ]))

# APPROACH 2: extract data from each tse and then combine
tses <- map(
  list(tse_skippy, tse_bibo, tse_bingo1, tse_bingo2),
  function(tse) {
    new_tse <- tse
    rowData(new_tse) <- rowData(new_tse) %>%
      as_tibble() %>%
      mutate(
        across(everything(),
        function(x) str_remove(x, ".*[kpcofgs]__"))) %>%
      DataFrame()
    new_tse <- agglomerateByRank(new_tse, rank = "genus")
    new_tse <- transformSamples(new_tse, method = "relabundance")
    as.data.frame(assay(new_tse, "relabundance")) %>% rownames_to_column("taxon")
  }
)

# we can see that there are different numbers of genera in each df
map(tses, ~dim(.x))
df <- full_join(tses[[1]], tses[[2]], by = "taxon") %>%
  full_join(tses[[3]], by = "taxon") %>%
  full_join(tses[[4]], by = "taxon") %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x))


# Plots are identical (checked them with the counts but now we need relab)
# the next step is to get the age at sample collection for each study

# SKIPPY
load(file = here::here("data/data_imp.Rds"))
skage <- select(d, skippy_id, sample_id, week, age) %>%
  mutate(age = age + as.numeric(as.character(week)) * 7) %>%
  select(sample_id, age)
mapfiles <- list.files(here::here("data/map_files"), full.names = TRUE)
id_swap <- map_dfr(mapfiles, ~read_excel(.x)) %>%
  select(sample_id = Seq_ID, internal_sample_id) %>%
  mutate(
    sample_id = str_replace(sample_id, "SKIPPY_", ""),
    sample_id = str_replace(sample_id, "_2", "_1"),
    sample_id = str_replace(sample_id, "_52", "_3"),
    sample_id = str_replace(sample_id, "_5", "_2"),
    sample_id = str_replace(sample_id, "245_1_15-1-17", "245_5"),
    sample_id = str_replace(sample_id, "245_1_22-1-17", "245_5_2"),
    sample_id = str_replace(sample_id, "269_1_15_mei_2017", "269_5"),
    sample_id = str_replace(sample_id, "269_1_18_mei_2017", "269_5_2")
    )

id_swap %>% filter(str_detect(sample_id, "269_"))
skage <- full_join(skage, id_swap, by = "sample_id") %>%
  select(sample_id = internal_sample_id, age) %>%
  na.omit()
mean(skage$sample_id %in% id_swap$sample_id)


# BINGO
load(here::here("data/tse_bingo1.Rds"))
bin_age0 <- colData(tse_bingo1) %>%
  as.data.frame() %>%
  mutate(
    # id = str_extract(id, "\\d\\d\\d"),
    time = ifelse(
      str_detect(sample_id, "k1"), "week2", ifelse(
        str_detect(sample_id, "k2"), "week6", ifelse(
          str_detect(sample_id, "k3"), "week12", NA
        )
      )
    )
  ) %>%
  select(sample_id, age)


bin_age1 <- read_excel(here("data/BINGO1y_ActualAge.xlsx")) %>%
  select(id = FamilyID, age = "1y-actual age_corrected in days") %>%
  filter(!is.na(age))

# we need to add the correct sample ids to the file
map_files <- c(
  "2019_11_06_Sequencing sample collection_2019_0084.xlsx",
  "2019_11_06_Sequencing sample collection_2019_0085.xlsx",
  "2020_01_14_Sequencing sample collection_2020_0005.xlsx"
)
bingo_map <- map_dfr(map_files, function(mapfile) {
  xl <- read_excel(glue(here("data/bingo_mapfiles/{mapfile}"))) %>%
    select(internal_sample_id, ProjectName, Seq_ID) %>%
    filter(ProjectName == "BINGO")
})

bingo_map_1 <- filter(bingo_map, str_detect(Seq_ID, "1_")) %>%
  mutate(id = str_pad(str_remove(Seq_ID, "1_"), 3, side = "left", "0")) %>%
  select(sample_id = internal_sample_id, id)
bin_age1 <- full_join(bin_age1, bingo_map_1, by = "id") %>%
  select(-id)

bin_age2 <- read_csv(here("data/bingo_age_1_3_years.csv")) %>%
  select(
    ID,
    birthdate = "Geboortedatum baby",
    coldate1 = "Poop collection tube arrives",
    coldate2 = "Poop collection date"
  ) %>%
  filter(!is.na(ID), ID != "Test") %>%
  mutate(
    coldate2 = str_remove(coldate2, "^\\w+,\\s"),
    coldate1 = mdy(coldate1),
    coldate2 = mdy(coldate2),
    birthdate = dmy(birthdate),
    age = ifelse(
      is.na(coldate2),
      coldate1 - birthdate,
      coldate2 - birthdate)
  ) %>%
  select(id = ID, age)

bingo_map_2 <- filter(bingo_map, str_detect(Seq_ID, "3_")) %>%
  mutate(id = str_pad(str_remove(Seq_ID, "3_"), 3, side = "left", "0")) %>%
  select(sample_id = internal_sample_id, id)
bin_age2 <- full_join(bin_age2, bingo_map_2, by = "id") %>%
  select(-id)

bin_age <- bind_rows(bin_age0, bin_age1, bin_age2) %>%
  na.omit()
bin_age



# BIBO
# first add proper sample ids
raw_path <- here::here("data/bibo_mapfiles/")
map_files <- list.files(raw_path, pattern = ".xlsx$")
lib_nums <- str_extract(map_files, "\\d\\d\\d\\d_00\\d\\d")
library_number_merged <- seq_along(map_files)
map_file_merged <- map2_dfr(
  map_files,
  library_number_merged,
  function(filename, lnm) {
  # extract library number
  lib_num <- str_extract(filename, "\\d\\d\\d\\d_00\\d\\d")

  # read provided map file
  xl <- read_excel(glue("{raw_path}/{filename}"))

  # create new map file for NGtax2.0
  map_file <- xl %>%
    filter(
      !grepl("L\\d\\d_NC_\\d", `Seq_ID`),
      Seq_ID != "empty",
      ProjectName == "BIBO"
    ) %>%
    select(sample_id = internal_sample_id, Seq_ID)

  return(map_file)
})

#---code book: a, 1m; b, 3m; c, 4m; d, 6y; e, 10y.
bibo_map <- map_file_merged %>%
  mutate(
    id = str_extract(Seq_ID, "\\d+"),
    time = str_extract(Seq_ID, "\\w{1}"),
    time = ifelse(time == "a", 28,
            ifelse(time == "b", 75, ifelse(time == "c", 105, NA))),
    id = as.numeric(id),
    time = as.numeric(time)
  ) %>%
  filter(time %in% c(28, 75, 105))
bibo_age_path <- "data/bibo_age_per_sample_infancy.xlsx"
bibo_age <- read_excel(bibo_age_path) %>%
  select(
    id = ID,
    day28 = "28 days",
    day75 = "CC -2 days",
    day105 = "CC+28 days"
  ) %>%
  mutate(id = as.integer(id)) %>%
  filter(!is.na(id)) %>%
  pivot_longer(contains("day"), names_to = "time", values_to = "age") %>%
  mutate(
    id = as.numeric(id),
    time = as.numeric(str_extract(time, "\\d+"))
  )

bibo_age <- dplyr::full_join(bibo_age, bibo_map, by = c("id", "time")) %>%
  select(sample_id, age) %>%
  filter(!is.na(sample_id))


bibo_age


# now we got all ages. We need to create a df that has sample ids as rownames
# and genera as colnames + age.

df <- full_join(tses[[1]], tses[[2]], by = "taxon") %>%
  full_join(tses[[3]], by = "taxon") %>%
  full_join(tses[[4]], by = "taxon") %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
  column_to_rownames("taxon")
colnames(df)
# first we do prevalence filtering.
prev <- 0.1
# for that level we will loose % of all genera -->
mean(rowSums(df[, -1]) <= 0.1)
df <- df[c(rowSums(df[, -1]) > 0.1), ]
df <- as.data.frame(t(df)) %>%
  rownames_to_column("sample_id")







sample_ages <- bind_rows(bin_age, bibo_age, skage) %>% na.omit()
filter(sample_ages) %>%
  ggplot(aes(age)) +
    geom_histogram(bins = 25)


d <- full_join(df, sample_ages, by = "sample_id") %>%
  na.omit()


# now we can train the model
colnames(d) <- make.names(colnames(d))
# in the paper they use ntree = 10k and mtry = p/3
p <- ((ncol(d) - 2) / 3) %>% round()
# i exclude the skippy samples for training
traindata <- filter(d, !sample_id %in% skage$sample_id)
# also leave out samples that are far above the age we work with here.
traindata <- filter(traindata, age <= 1000)
testdata <- filter(d, sample_id %in% skage$sample_id)
model <- ranger(
  formula = age ~ .,
  data = select(traindata, -sample_id),
  num.trees = 1e4,
  mtry = p,
  importance = "permutation"
)
model

# calculate correlation of predictions and actual values
pred <- predict(model, data = select(testdata, -sample_id))
testdata$pred <- pred$predictions
r <- cor.test(testdata$pred, testdata$age)$estimate
r
source(here::here("R/ml_helper.R"))

if (!file.exists(here::here("data/rvalues.Rds"))) {
  rvalues <- rf_null(
    y = "age",
    features = select(traindata, -sample_id, -age) %>% colnames(),
    train = traindata,
    test = testdata,
    n_perm = 500,
    ntree = 500
  )
  save(rvalues, file = here::here("data/rvalues.Rds"))
 } else {
  load(here::here("data/rvalues.Rds"))
}

# pvalue
mean(rvalues > r)


out <- select(testdata, internal_sample_id = sample_id, pred) %>%
  full_join(id_swap, by = "internal_sample_id") %>%
  filter(!is.na(pred)) %>%
  select(sample_id, pred)
out



# now we have the data to calculate the microbiota age and MAZ

# calculate microbiota for age z score
# first we need the median microbiota age of healthy children that are in the
# same age (same month). I.E. month 1 (2 weeks), month 2 (5 weeks) and month 12
# (52 weeks).

traindata$mbage <- predict(model, data = select(traindata, -sample_id))$predictions
median_mbage <- select(traindata, sample_id, age, mbage) %>%
  mutate(month = as.integer(age / (365 / 12))) %>%
  group_by(month) %>%
  summarise(md = median(mbage), sd = sd(mbage), n = n())
median_mbage
age_range <- ggplot(traindata, aes(age)) +
  geom_histogram(bins = 50, color = "white", fill = "#606060") +
  theme_bw(base_size = 25) +
  scale_x_continuous(breaks = seq(0, 600, 100)) +
  scale_y_continuous(breaks = seq(0, 110, 10)) +
  xlab("Age") + ylab("n")
age_range


# the formula is: microbiota age - median microbiota age of healthy children
# of same echronological age / sd of the healthy childrens mage
maz <- select(testdata, internal_sample_id = sample_id, age, mbage = pred) %>%
  left_join(id_swap, by = "internal_sample_id") %>%
  select(-internal_sample_id) %>%
  mutate(
    month = as.integer(age / (365 / 12))
  )


maz$maz <- NA

for (i in seq_along(maz$mbage)) {
  maz$maz[i] <- as.numeric((maz$mbage[i] -
            median_mbage[median_mbage$month == maz$month[i], "md"]) /
              median_mbage[median_mbage$month == maz$month[i], "sd"])
}

# for analyses it seems better to use microbiota age rather than MAZ because to calculate 
# the sd, more coverage in certain months would have been needed. However, including age in
# the regression model is comparable in meaning and works fine.







###############################################################################
#########################           1. ITT       ##############################
###############################################################################



####################### 1.1 Complete Case Analysis ############################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

# for analyses we apply prevalence fitlering and analyze at genus level
fvars <- c("siblings", "condition")
# add metadata to tse
d <- colData(tse) %>%
  as.data.frame() %>%
  select(-siblings) %>%
  left_join(select(d, age, sample_id, siblings), by = "sample_id") %>%
  mutate(across(all_of(fvars), function(x) as.factor(x))) 
d$age <- d$age + as.numeric(d$week) * 7
d$age <- scale(d$age)[, 1]

if(!file.exists(here::here("data/mbagedataa.Rds"))) {
  d_cc <- full_join(select(d, -age), maz, by = "sample_id") %>%
  mutate(
    condition_label = ifelse(condition == 0, "CAU", ifelse(
      condition == 1, "SSC", NA)),
    week_label = glue::glue("Week {week}")
  ) %>% 
  group_by(week) %>%
  mutate(
    maz_s = scale(maz)[, 1],
    maz_c = maz - median(maz, na.rm = TRUE)
    ) %>%
  ungroup()
  save(d_cc, file = here::here("data/mbagedataa.Rds"))
 } else {
  load(here::here("data/mbagedataa.Rds"))
}

# create plot
mba_plot <- ggplot(d_cc, aes(condition_label, mbage, fill = condition_label)) +
  geom_boxplot(outlier.alpha = 0) +
  #geom_beeswarm(alpha = 0.4) + 
  geom_jitter(width = 0.1, size = 2) + 
  #stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red") +
  facet_wrap(~week_label, strip.position = "bottom") +
  scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
  theme_bw(base_size = 25) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 20)
    ) +
  xlab("") + ylab("Microbiota Age")
mba_plot
save(mba_plot, age_range, file = here::here("data/mba_out.Rds"))

# how is bfexcl distruted
d_cc$bfexcl
dim(d_cc)
select(d_cc, skippy_id, bfexcl) %>% arrange(skippy_id) %>%
  distinct() %>%
  ggplot(aes(bfexcl)) +
  geom_histogram()



# check optimal model structure and then fit the models
coefs <- c(
  "csection", 
  "birthweight_s", 
  "siblings",
  "sex",
  "apgar_5_s",
  "ges_age_s",
  "edlevel_s"
)


loo_comp <- map_dfr(coefs, function(coef) {
  map_dfr(1:5, function(i) {
    map_dfr(list(c(2, 5), 52), function(n_week) {
      if (52 %in% n_week) {
        # base model 
        f1 <- bf(mbage ~ age * condition)
        m1 <- brm(
          data = filter(implist[[i]], week == 52),
          formula = f1,
          file = here::here(glue("data/m1_age_imp{i}_{n_week[1]}_mbage.Rds"))
        )
        loo_m1 <- add_criterion(
          m1,
          "loo",
          file = here::here(glue("data/loo_m1_age_imp{i}_{n_week[1]}_mbage")),
          moment_match = FALSE
        )

        f2 <- bf(glue("mbage ~ age * condition + {coef}"))
        m2 <- brm(
          data = filter(implist[[i]], week == 52),
          formula = f2,
          file = here::here(glue("data/m2_age_imp{i}_{coef}_{n_week[1]}_mbage.Rds"))
        )

        loo_m2 <- add_criterion(
          m2,
          "loo",
          file = here::here(glue("data/loo_m2_age_imp{i}_{coef}_{n_week[1]}_mbage")),
          moment_match = FALSE
        )
        lcomp <- loo_compare(loo_m2, loo_m1)
        score <- ifelse(rownames(lcomp)[1] == "loo_m1", 0, 1)
      } else {
        # base model 
        f1 <- bf(mbage ~ age * condition + (1|skippy_id))
        m1 <- brm(
          data = filter(implist[[i]], week != 52),
          formula = f1,
          file = here::here(glue("data/m1_age_imp{i}_{n_week[1]}_mbage.Rds"))
        )
        loo_m1 <- add_criterion(
          m1,
          "loo",
          file = here::here(glue("data/loo_m1_age_imp{i}_{n_week[1]}_mbage")),
          moment_match = FALSE
        )

        f2 <- bf(glue("mbage ~ age * condition + {coef} + (1|skippy_id)"))
        m2 <- brm(
          data = filter(implist[[i]], week != 52),
          formula = f2,
          file = here::here(glue("data/m2_age_imp{i}_{coef}_{n_week[1]}_mbage.Rds"))
        )

        loo_m2 <- add_criterion(
          m2,
          "loo",
          file = here::here(glue("data/loo_m2_age_imp{i}_{coef}_{n_week[1]}_mbage")),
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

group_by(loo_comp, coef, model) %>%
  summarise(ss = sum(score))










formula <- bf(
  mbage ~ condition * age + siblings + sex + csection + (1 | skippy_id)
)
model_infancy_cc <- brm(
  family = student(),
  formula = formula,
  data = filter(d_cc, week != 52),
  file = here::here("data/mbage_infancy_cc_cov")
)

formula <- bf(mbage ~ age + condition + siblings + sex + csection)
model_year1_cc <- brm(
  family = student(),
  formula = formula,
  data = filter(d_cc, week == 52),
  file = here::here("data/mbage_year1_cc_cov")
)

formula <- bf(mbage ~ condition * age + siblings + sex + csection)
model_all_cc <- brm(
  family = student(),
  formula = formula,
  data = d_cc,
  file = here::here("data/mbage_all_cc_cov")
)

summary(model_infancy_cc)
summary(model_year1_cc)
summary(model_all_cc)

# sensitivity analyses (no covs)
formula <- bf(mbage ~ age + condition)
model_year1_cc_bez <- brm(
  family = student(),
  formula = formula,
  data = filter(d_cc, week == 52),
  file = here::here("data/mbage_year1_cc_bezcov")
)
summary(model_year1_cc_bez)






######################## 1.2 Multiple imputation  #############################


implist <- map(implist, function(imp) {
  imp_new <- imp %>%
    left_join(
      select(d_cc, sample_id, mbage),
      by = "sample_id") 
  mice::complete(mice::mice(imp_new))
})


# fit model
formula <- bf(
  mbage ~ age * condition + siblings + sex + csection + (1 | skippy_id)
)
model_infancy <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week != 52)),
  file = here::here("data/mbage_infancy_cov")
)

formula <- bf(
  mbage ~ age * condition + siblings + sex + csection + bfexcl + (1 | skippy_id)
)
model_infancy_bf <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week != 52)),
  file = here::here("data/mbage_infancy_cov_bf")
)

formula <- bf(mbage ~ age + condition + siblings + sex + csection)
model_year1 <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_cov")
)

post <- posterior_samples(model_year1)
HDInterval::hdi(post$b_condition1, prob = 0.95)
mean(post$b_condition1 < 0)
summary(model_year1)

# sensitivity analyses without covs
formula <- bf(mbage ~ age + condition)
model_year1_bez <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_bez")
)

post <- posterior_samples(model_year1_bez)
HDInterval::hdi(post$b_condition1, prob = 0.95)
mean(post$b_condition1 < 0)


formula <- bf(mbage ~ age * condition + siblings + sex + csection + (1 | skippy_id))
model_all <- brm_multiple(
  family = student(),
  formula = formula,
  data = implist,
  file = here::here("data/mbage_all_cov")
)

summary(model_infancy)
summary(model_year1)
summary(model_all)

formula <- bf(mbage ~ age + condition + siblings + sex + csection + bfexcl)
model_year1_bf <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_cov_bf")
)
summary(model_year1_bf)
post <- posterior_samples(model_year1_bf)
HDInterval::hdi(post$b_condition1, prob = 0.95)
mean(post$b_condition1 < 0)
HDInterval::hdi(post$b_bfexcl, prob = 0.95)
mean(post$b_bfexcl < 0)
implist[[1]]$bfexcl


formula <- bf(mbage ~ age + condition + siblings + sex + csection + weaning)
model_year1_weaning <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_cov_w")
)
summary(model_year1_weaning)


implist_w <- map(implist, function(d) {
  d %>% mutate(months_bf_weaning = weaning - bfexcl)
})
formula <- bf(mbage ~ age + condition + siblings + sex + csection + bfexcl + months_bf_weaning)
model_year1_weaning <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_w, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_cov_months_bf_weaning")
)
summary(model_year1_weaning)


save(model_infancy, model_infancy_bf, model_year1, model_all, model_year1_bf, file = here::here("data/mbage_models_itt.Rds"))






###############################################################################
#########################           2. PP        ##############################
###############################################################################

fvars <- c("siblings", "condition", "csection", "sex")
# obtain ids that were selected for PP analyses 
pp_indicator <- foreign::read.spss(here::here("data/raw_data/kelly141022/Data_ITT_PP_ExploratoryDRselections.sav"), to.data.frame = TRUE)
pp_indicator <- select(pp_indicator, skippy_id = ID, pp = PP) %>%
  mutate(skippy_id = as.character(skippy_id))
# add pp info to existing data 
if (!"pp" %in% colnames(d_cc)) {
  d_cc <- left_join(d_cc, pp_indicator, by = "skippy_id")
}
d_pp <- filter(d_cc, pp == 1) %>%
    mutate(
    across(all_of(fvars), function(x) as.factor(x)),
    ges_age_s = scale(ges_age)[, 1],
    birthweight_s = scale(birthweight)[, 1],
    edlevel_s = scale(edlevel)[, 1]
    )
pp_indicator <- mutate(pp_indicator, skippy_id = as.integer(skippy_id))
implist_pp <- map(implist, function(dimp) {
  dimp_new <- left_join(dimp, pp_indicator, by = "skippy_id") %>%
                filter(pp == 1)
  dimp_new
})




####################### 2.1 Complete Case Analysis ############################


formula <- bf(
  mbage ~ age * condition + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex + (1 | skippy_id)
)
model_infancy_cc_pp <- brm(
  family = student(),
  formula = formula,
  data = filter(d_pp, week != 52),
  file = here::here("data/mbage_infancy_cc_pp")
)

formula <- bf(mbage ~ age + condition + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex)
model_year1_cc_pp <- brm(
  family = student(),
  formula = formula,
  data = filter(d_pp, week == 52),
  file = here::here("data/mbage_year1_cc_pp")
)

summary(model_infancy_cc_pp)
summary(model_year1_cc_pp)







######################## 2.2 Multiple imputation  #############################




# fit model
formula <- bf(
  mbage ~ age * condition + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex + (1 | skippy_id)
)
model_infancy <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_pp, ~filter(.x, week != 52)),
  file = here::here("data/mbage_infancy_pp_cov")
)

formula <- bf(
  mbage ~ age * condition + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex + bfexcl + (1 | skippy_id)
)
model_infancy_bf <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_pp, ~filter(.x, week != 52)),
  file = here::here("data/mbage_infancy_pp_cov_bf")
)

formula <- bf(mbage ~ age + condition + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex)
model_year1 <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_pp, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_pp_cov")
)

summary(model_infancy)
summary(model_year1)

# same results here. In conclusion, ssc might have effect on microbiota age
# such that infants microbiota is longer infant like, potentially because of
# longer breastfeeding but this does not explain effect entirely.


formula <- bf(mbage ~ age + condition + siblings + bfexcl + birthweight_s + ges_age_s + edlevel_s + csection + sex)
model_year1_bf <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist_pp, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_bf_pp_cov")
)
summary(model_year1_bf)

save(model_infancy, model_infancy_bf, model_year1, model_year1_bf, file = here::here("data/mbage_models_pp.Rds"))



# we can argue that there is evidence for an effect of SSC on MBA
# therefore we will look at the DR analysis as well here:

dr <- foreign::read.spss(
  here::here("data/kelly_documents/data_itt_pp_dr.sav"),
  to.data.frame = TRUE
  ) %>%
  select(skippy_id = ID, ITT, SSC = TotalSSCwk1wk5) %>%
  mutate(SSC_s = scale(SSC)[, 1])
head(dr)
ggplot(dr, aes(SSC)) +
  geom_histogram()


implist <- map(implist, function(imp) {
  imp_new <- imp %>%
    dplyr::left_join(
      select(dr, skippy_id, SSC_s),
      by = "skippy_id")
  mice::complete(mice::mice(imp_new))
})


dr$skippy_id <- as.character(dr$skippy_id)
df <- left_join(d_cc, select(dr, skippy_id, SSC_s, SSC), by = "skippy_id") 
implist[[1]] %>% colnames()

# fit model
formula <- bf(
  mbage ~ SSC_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex + age + (1 | skippy_id)
)
model_infancy <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week != 52)),
  file = here::here("data/mbage_infancy_dr_cov")
)

formula <- bf(mbage ~ age + SSC_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex)
model_year1 <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_dr_cov")
)



summary(model_infancy)
summary(model_year1)


formula <- bf(mbage ~ age + SSC_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex + bfexcl)
model_year1_bf <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/maz_year1_bf_dr_cov")
)
summary(model_year1_bf)


# only within SSC

formula <- bf(mbage ~  age + SSC_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex)
model_year1_within <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52, condition == 1)),
  file = here::here("data/mbage_year1_dr_within_cov")
)



summary(model_year1_within)




formula <- bf(mbage ~  age + SSC_s + siblings + birthweight_s + ges_age_s + edlevel_s + csection + sex)
model_year1_out <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52, condition == 0)),
  file = here::here("data/mbage_year1_dr_out_cov")
)

summary(model_year1_out)

save(model_year1, model_year1_within, model_year1_bf, file = here::here("data/mbage_models_dr.Rds"))


filter(d_cc, is.na(condition))
ssc_plot1 <- full_join(d_cc, dr, by = "skippy_id") %>% 
  select(skippy_id, condition, SSC) %>%
   mutate(
    condition_label = ifelse(condition == 0, "CAU", ifelse(
      condition == 1, "SSC", NA)),
    SSC_h = SSC/60) %>% 
  distinct(.keep_all = TRUE) %>%
  filter(skippy_id %in% d_cc$skippy_id) %>%
  ggplot(aes(condition_label, SSC_h, fill = condition_label)) +
  geom_boxplot(outlier.alpha = 0) +
  #geom_beeswarm() +
  geom_jitter(width = 0.1, size = 2) + 
  #scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
  scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
  theme_bw(base_size = 25) +
  theme(
  legend.position = "none",
  strip.placement = "outside",
  strip.background = element_blank()) +
xlab("") + ylab("SSC Total Hours")


ssc_plot2 <- full_join(d_cc, dr, by = "skippy_id") %>% 
  select(skippy_id, condition, SSC, mbage, week) %>%
   mutate(
    condition_label = ifelse(condition == 0, "CAU", ifelse(
      condition == 1, "SSC", NA)),
    SSC_h = SSC/60) %>% 
  filter(skippy_id %in% d_cc$skippy_id, week == 52) %>%
  ggplot(aes(SSC_h, mbage)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_vline(xintercept = 1455.561/60, linetype = "dashed") +
  scale_fill_manual(values = c("#ffffff", "#c0c1c2")) +
  #scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
  #scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
  theme_bw(base_size = 25) +
  theme(
  legend.position = "none",
  strip.placement = "outside",
  strip.background = element_blank()) +
xlab("SSC Total Hours") + ylab("Microbiota Age")
ssc_plot1
ssc_plot2

save(ssc_plot1, ssc_plot2, file = here::here("data/sscplot.Rds"))





implist <- map(implist, ~mutate(.x, upper = ifelse(SSC_s <= -0.5, "no", "yes")))


formula <- bf(mbage ~ age + upper + siblings + sex + csection)
model_year1_cat <- brm_multiple(
  family = student(),
  formula = formula,
  data = map(implist, ~filter(.x, week == 52)),
  file = here::here("data/mbage_year1_dr_cat_cov")
)

model_year1_cat
implist[[1]] %>% colnames()
test <- implist[[1]]
test$rs <- resid(lm(mbage ~ age, data = implist[[1]]))

ggplot(test, aes(SSC_s, rs)) +
  geom_point() +
  geom_smooth()










################################################################################
############################## make heatmap from DAA and RF ####################
################################################################################

load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))
source(here::here("R/ml_helper.R"))

files <- c(
  here::here("data/maaslin2_tables_itt.Rds"),
  # here::here("data/ancom_tables_itt.Rds"),
  here::here("data/tables_linda_itt.Rds")
)
for (file in files) load(file)

linda_identified
maaslin2_identified


plot_importance(model, top_n = 20)

rf_feats <- as.character(arrange(extract_importance(model, n = 20), desc(importance))$features)
plot_importance(model, top_n = 25) + theme_bw(base_size = 25)
# plot relabundance of the features but only for those who are nonzero in SKIPPY
rf_feats <- map_chr(rf_feats, function(feat) {
  out <- NA
  if(sum(testdata[[feat]]) > 0) {
    out <- feat
  }
  out}) %>%
  na.omit() %>%
  as.character()
rf_feats


msl <- make.names(maaslin2_identified$feature) %>% str_replace("genus.Eubacterium._hallii_group", "genus..Eubacterium._hallii_group")
lnd <- make.names(linda_identified$taxon) %>% str_replace("genus.Eubacterium._hallii_group", "genus..Eubacterium._hallii_group")
msl
feats <- union(union(rf_feats, lnd), msl)
feats

# now extract clr values for the heatmap
tse <- agglomerateByRank(tse, rank = "genus")
tse <- transformSamples(tse, method = "relabundance")
tse <- transformSamples(tse, abund_values = "relabundance", method = "clr", pseudocount = 1)
hmd <- as.data.frame(assay(tse, "clr")) %>% rownames_to_column("taxon")
hmd$taxon <- make.names(hmd$taxon)
hmd <- filter(hmd, taxon %in% feats) %>%
  mutate(
    RF = ifelse(taxon %in% rf_feats, 1, 0),
    Maaslin2 = ifelse(taxon %in% msl, 1, 0),
    LinDA = ifelse(taxon %in% lnd, 1, 0),
) 
side_vars <- select(hmd, taxon, RF, Maaslin2, LinDA) %>% 
  mutate(score = RF + Maaslin2 + LinDA) %>%
  arrange(desc(score))
  # the genus must be non-zero and identifiable in SKIPPY. 

side_vars
hmd <- select(hmd, -RF, -Maaslin2, -LinDA) %>% 
  column_to_rownames("taxon") %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(all_of(feats), function(x) scale(x)[, 1])) %>%
  rownames_to_column("sample_id")
temp <- colData(tse) %>% 
  as.data.frame() %>%
  mutate(Infancy = ifelse(week == 52, "Late Infancy", "Early Infancy")) %>%
  select(sample_id, condition, Infancy)
hmd <- left_join(hmd, temp, by = "sample_id")


side_vars
hmd <- arrange(
  hmd,
  Infancy,
  condition,
  genus.Faecalibacterium,
  genus.Megasphaera,
  genus.Bacteroides,
  genus.Flavonifractor,
  genus.Rothia,
  genus..Eubacterium._hallii_group,
  genus.Bifidobacterium,
  genus.Parabacteroides,
  genus.Enterococcus,
  genus.Erysipelatoclostridium,
  genus.Blautia,
  family.Lachnospiraceae,
  genus.Butyricicoccus,
  genus.Actinomyces,
  genus.Staphylococcus,
  genus.Anaerostipes,
  genus.Dialister,
  genus.Lacticaseibacillus,
  genus.Lachnospira,
  genus.Lachnospiraceae_UCG.004,
  genus.Monoglobus
)

sample_id <- hmd$sample_id
clu <- hmd$condition
age <- hmd$Infancy

method_label <- mutate(
  side_vars, 
  p1 = ifelse(RF, "R, ", ""), 
  p2 = ifelse(Maaslin2, "M, ", ""), 
  p3 = ifelse(LinDA, "L, ", ""),
  label = glue("{p1}{p2}{p3}"),
  label = str_replace(label, ",\\s+$", ""),
  rownames = glue('"{taxon}"^"{label}",'))

n <- length(method_label$label)
split <- method_label$label[n:1]
split
method_label
df_fin <- hmd[, side_vars$taxon]
df_fin <- t(df_fin)
# sort only based on small matrix 
colnames(df_fin) <- clu

mat <- df_fin
# Legend color
nvec <- seq(-1, 1, length.out = 8)
col_fun <- colorRamp2(nvec, 
                     c(
                       "#2166ac", 
                       "#4393c3", 
                       "#92c5de", 
                       "#d1e5f0", 
                       "#fddbc7", 
                       "#f4a582", 
                       "#d6604d", 
                       "#b2182b"
                   ))


# Annotation
ha <- HeatmapAnnotation(
  SSC = clu,
  col = list(
    SSC = c("0"="#B0B0B0", "1"="#000000")
    ),
  annotation_name_side = "left"
  )


# Heatmap
png(here::here("fig/rf_feats.png"),
  width = 50, 
  height = 40, 
  units = "cm", 
  res= 900)

  Heatmap(mat,
      width = unit(32, "cm"),
      height = unit(32, "cm"),
      top_annotation = ha,
      name = "RA", 
      col = col_fun,
      column_title = "",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_title = "",
      row_title_side = "right",
      row_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_names_gp = gpar(fontsize = 15, fontface = "italic"),
      row_names_side = "left",
      show_column_names = F,
      border = T,
      # rect_gp = gpar(col = "white", lwd = 0),
      cluster_rows = F,
      cluster_columns = F,
      show_heatmap_legend = TRUE,
      heatmap_legend_param = list(ncol =1),
      # row_split = split,
      # row_title_rot = 0,
      # row_order = side_vars$taxon,
      row_labels = expression(
        "Faecalibacterium"^"R, M, L",
        "Megasphaera"^"M, L",
        "Bacteroides"^"M, L",
        "Flavonifractor"^"M, L",
        "Rothia"^"M, L",
        "Eubacterium hallii group"^"R, M",
        "Bifidobacterium"^"R",
        "Parabacteroides"^"L",
        "Enterococcus"^"L",
        "Erysipelatoclostridium"^"L",
        "Blautia"^"R",
        "family:Lachnospiraceae"^"R",
        "Butyricicoccus"^"R",
        "Actinomyces"^"R",
        "Staphylococcus"^"R",
        "Anaerostipes"^"R",
        "Dialister"^"R",
        "Lacticaseibacillus"^"M",
        "Lachnospira"^"R",
        "Lachnospiraceae UCG 004"^"R",
        "Monoglobus"^"R"
      ),
      column_split = age,
    )
  dev.off()

pdf(here::here("fig/rf_feats.pdf"), width = 25, height = 20)
  Heatmap(mat,
      width = unit(32, "cm"),
      height = unit(32, "cm"),
      top_annotation = ha,
      name = "RA", 
      col = col_fun,
      column_title = "",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_title = "",
      row_title_side = "right",
      row_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_names_gp = gpar(fontsize = 15, fontface = "italic"),
      row_names_side = "left",
      show_column_names = F,
      border = T,
      # rect_gp = gpar(col = "white", lwd = 0),
      cluster_rows = F,
      cluster_columns = F,
      show_heatmap_legend = TRUE,
      heatmap_legend_param = list(ncol =1),
      # row_split = split,
      # row_title_rot = 0,
      # row_order = side_vars$taxon,
      row_labels = expression(
        "Faecalibacterium"^"R, M, L",
        "Megasphaera"^"M, L",
        "Bacteroides"^"M, L",
        "Flavonifractor"^"M, L",
        "Rothia"^"M, L",
        "Eubacterium hallii group"^"R, M",
        "Bifidobacterium"^"R",
        "Parabacteroides"^"L",
        "Enterococcus"^"L",
        "Erysipelatoclostridium"^"L",
        "Blautia"^"R",
        "family:Lachnospiraceae"^"R",
        "Butyricicoccus"^"R",
        "Actinomyces"^"R",
        "Staphylococcus"^"R",
        "Anaerostipes"^"R",
        "Dialister"^"R",
        "Lacticaseibacillus"^"M",
        "Lachnospira"^"R",
        "Lachnospiraceae UCG 004"^"R",
        "Monoglobus"^"R"
      ),
      column_split = age,
    )
  dev.off()

side_vars


# further data exploration for discussion
load(here::here("data/data.Rds"))
load(file = here::here("data/data_imp.Rds"))

tse <- agglomerateByRank(tse, rank = "genus")
tse2 <- tidySummarizedExperiment::filter(tse, week == 2)
tse5 <- tidySummarizedExperiment::filter(tse, week == 5)
tse52 <- tidySummarizedExperiment::filter(tse, week == 52)
prevs <- map(list(tse2, tse5, tse52), function(tse) {
  prev <- getPrevalence(
    tse, 
    detection = 1, 
    sort = TRUE, 
    assay.type = "counts",
    as_relative = FALSE
    )
    prev[names(prev) == "genus:Faecalibacterium"]
})
prevs

cd <- colData(tse) %>% as.data.frame() %>%
  select(skippy_id, sample_id, week, condition, bfexcl, siblings)
a <- as.data.frame(t(assay(tse, "counts"))) %>% 
  rownames_to_column("sample_id") %>%
  select(sample_id, "genus:Faecalibacterium")

df <- full_join(cd, a, by = "sample_id") %>%
  rename(faecalibacterium = `genus:Faecalibacterium`) %>%
  mutate(
    present = ifelse(faecalibacterium > 0, 1, 0),
    l_faecalibacterium = log((faecalibacterium + 1))
    )
group_by(df, week, condition) %>%
  summarise(m = mean(present))

summary(lm(l_faecalibacterium ~ siblings + bfexcl, data = df))
