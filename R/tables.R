# this code is from Jennifer but I (Henrik) checked it

library(table1)
library(dplyr)
library(tidyverse)
library(lubridate)

load(here::here("data/data_imp.Rds"))
load(here::here("data/data.Rds"))

bf <- readxl::read_excel(f, sheet = "Feeding", na = "\\") %>% 
  select(skippy_id = ID, bfexcl = ExclusiveBF) %>%
  mutate(skippy_id = as.character(skippy_id))

itt_var <-  read_csv2(here::here("data/skippy_stool_data_cleaned.csv")) %>%
  mutate_all(function(x) ifelse(x == 99999, NA, ifelse(x == 88888, NA, x))) %>%
  filter(id %in% sav$id) %>%
  mutate(across(contains("dat_"), function(x) {
    x_new <- str_replace(x, "okt", "october")
    x_new <- str_replace(x_new, "mrt", "march")
    x_new <- str_replace(x_new, "mei", "may")
    x_new <- str_replace(x_new, "15-jan-17/22-jan-17", "18/jan/17")
    x_new <- dmy(x_new)
    return(x_new)
  }),
  datbirth_infant = dmy(datbirth_infant),
  datbirth_mom = dmy(datbirth_mom),
  age_week2 = as.numeric(dat_week2 - datbirth_infant),
  age_week5 = as.numeric(dat_week5 - datbirth_infant),
  age_week52 = as.numeric(dat_1year - datbirth_infant),
  age_mom = as.numeric((dat_week2 - datbirth_mom)/365),
  stdat_week2 = as.numeric(dat_week2 - dmy('01-01-2016')),
  stdat_week5 = as.numeric(dat_week5 - dmy('01-01-2016')),
  stdat_week52 = as.numeric(dat_1year - dmy('01-01-2016')),
  id = as.character(id)) %>%
  select(
    skippy_id = id, matches("age_week\\d+"), condition, csection, birthweight, 
    siblings, sex, antibiotic_week2,
    antibiotic_week5, antibiotic_1year, apgar_5, ges_age, age_mom, smoking,
    drinking, weaning, parity, stdat_week2, stdat_week5, stdat_week52) %>%
  left_join(bf, by = "skippy_id")

# there is a typo in the maternal date of birth of ID 259, so her age will be removed 
itt_var$age_mom[52] <- NA

itt_var$condition <- factor(itt_var$condition,
                            levels = c(0, 1),
                            labels = c("CAU", "SSC"))
itt_var$sex <- factor(itt_var$sex,
                      levels = c(0, 1),
                      labels = c("Male", "Female"))
itt_var$csection <- factor(itt_var$csection,
                           levels = c("0", "1"),
                           labels = c("No", "Yes"))

itt_var$siblings <- as.integer(itt_var$siblings)
itt_var$parity <- as.integer(itt_var$parity)




itt_var$antibiotic_week2 <- factor(itt_var$antibiotic_week2,
                                     levels = c("0", "1"),
                                     labels = c("No", "Yes"))
itt_var$antibiotic_week5 <- factor(itt_var$antibiotic_week5,
                                     levels = c("0", "1"),
                                     labels = c("No", "Yes"))
itt_var$antibiotic_1year <- factor(itt_var$antibiotic_1year,
                                     levels = c("0", "1"),
                                     labels = c("No", "Yes"))

#change labels. Of course there is a better way, but due to time...
label(itt_var$sex) <- "Sex"
label(itt_var$ges_age) <- "Gestational age"
label(itt_var$csection) <- "C-section"
label(itt_var$siblings)<- "Siblings"
label(itt_var$parity) <- "Parity"
label(itt_var$birthweight) <- "Birth weight"
label(itt_var$condition) <- "Condition"
label(itt_var$bfexcl) <- "Exclusive breastfeeding duration"



label(itt_var$age_week2) <- "Age week 2"
label(itt_var$age_week5) <- "Age week 5"
label(itt_var$age_week52) <- "Age week 52"

label(itt_var$antibiotic_week2) <- "Antibiotics week 2"
label(itt_var$antibiotic_week5) <- "Antibiotics week 5"
label(itt_var$antibiotic_1year) <- "Antibiotics week 52"

units(itt_var$ges_age)       <- "weeks"
units(itt_var$birthweight) <- "grams"
units(itt_var$bfexcl) <- "months"
units(itt_var$age_week2) <- "days"
units(itt_var$age_week5) <- "days"
units(itt_var$age_week52) <- "days"


# produce table
tbl1 <- table1(~ sex + csection + ges_age + birthweight + siblings + bfexcl +
         age_week2 + age_week5 + age_week52 + antibiotic_week2 +
         antibiotic_week5 + antibiotic_1year | condition, data = itt_var,
         caption = "Demographics Table",
         )

save(tbl1, file = here::here("data/tbl1.Rds"))

tbl1
