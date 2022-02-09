# set up -----------------------------------------------------------------------
# load libraries
library(ggplot2)
library(lubridate)
library(tidyr)
library(stringr)
library(dplyr)

# set default working directory
setwd("~/GitHub/ggro_sharpie_prey/prey_preference")

# settings
options(scipen = 9999999)
# ------------------------------------------------------------------------------





# get and format datasets ------------------------------------------------------
# get diet sample data
d1 <- read.csv("sharpie1516long.csv") %>% 
  select(sample_id, date, age, sex, weight=ssha_weight, crop,
         prey_name, prey_mig=migratory, prey_mass) %>% # grab most impartant vars
  mutate(date=dmy(date), # fixed dates
         week_of_yr=week(date),
         year=year(date)) %>%
  filter(week_of_yr>=36, week_of_yr<=48) %>% # reconcile weeks
  mutate(age=tolower(age), # bin age groups
         age=ifelse(age=="sy", "ahy", age),
         age=ifelse(age=="asy", "ahy", age),
         sex=tolower(sex)) %>%
  select(sample_id, date, week_of_yr, year, age, sex, crop, prey_name) %>%
  mutate(use=1) %>%
  arrange(week_of_yr, sample_id) %>%
  ungroup()
summary(d1)
d1_taxa <- sort(unique(d1$prey_name)) # 69 taxa

# get ebird availability data
d2 <- read.csv("relative_predator_occurrence_prey_abundance_north_coast_CA_fall.csv") %>% 
  select(prey_name=prey_spp, 
         date=week,
         availability=sharp_wt_mean_prey_ai) %>%
  mutate(date=mdy(date),
         week_of_yr=week(date)) %>%
  filter(week_of_yr>=36, week_of_yr<=48) %>% # reconcile weeks
  filter(prey_name %in% d1_taxa) %>%
  select(prey_name, week_of_yr, prey_avail=availability) %>%
  arrange(prey_name, week_of_yr) %>%
  ungroup()
summary(d2)
d2_taxa <- sort(unique(d2$prey_name)) # 69 taxa

# check taxa overlap
all(d1_taxa %in% d2_taxa) # should be T
all(d2_taxa %in% d1_taxa) # should be T

# add prey availability and traits
d3 <- read.csv("traits.csv")
d3 <- left_join(d2, d3)
summary(d3)
d3_taxa <- sort(unique(d3$prey_name)) # 69 taxa

# check overlap all should be T
all(d1_taxa %in% d3_taxa)
all(d2_taxa %in% d3_taxa)
all(d1_taxa %in% d2_taxa)
all(d3_taxa %in% d2_taxa)

# all possible species by sharpie sample combinations. this is needed
# for discrete choice model analysis data format
d4 <- expand.grid(sample_id=unique(d1$sample_id),
                  prey_name=unique(d3$prey_name)) %>%
  mutate(sample_id=as.character(sample_id),
         prey_name=as.character(prey_name)) %>%
  left_join(d1) %>%
  group_by(sample_id) %>%
  select(-date) %>%
  mutate(week_of_yr=unique(week_of_yr)[which(unique(week_of_yr)!="NA")],
         year=unique(week_of_yr)[which(unique(week_of_yr)!="NA")],
         age=unique(age)[which(unique(age)!="NA")], 
         sex=unique(sex)[which(unique(sex)!="NA")],
         crop=unique(crop)[which(unique(crop)!="NA")],
         use=ifelse(is.na(use), 0, use)) %>%
  arrange(sample_id, prey_name) %>%
  ungroup()
summary(d4)
d4_taxa <- sort(unique(d4$prey_name)) # 69 taxa

# check taxon overlap
all(d1_taxa %in% d4_taxa)
all(d3_taxa %in% d4_taxa)
all(d2_taxa %in% d4_taxa)
all(d4_taxa %in% d3_taxa)
  
# join it all together into one long dataset
d5 <- d4 %>% left_join(d3) %>% 
  arrange(week_of_yr, sample_id, prey_name)
summary(d5)
d5_taxa <- sort(unique(d5$prey_name)) # 69 taxa

# check out abundance distributions for kicks
names(d5)
med_avails <- d5 %>% group_by(prey_name) %>% 
  summarise(median_avail=median(prey_avail)) %>%
  arrange(desc(median_avail)) %>% pull(median_avail)
hist(med_avails)
plot(med_avails)
d6 <- d5

# now break apart multiple choices per individual sharpie. again, long format
samp_ids <- sort(unique(d6$sample_id))
stacked_data <- c()
for(i in 1:length(samp_ids)){
  dat_i <- d6 %>% filter(sample_id==samp_ids[i])
  prey_j <- dat_i$prey_name[dat_i$use==1]
  new_i <- c()
  for(j in 1:length(prey_j)){
    dat_j <- dat_i
    dat_j$use <- 0
    dat_j$use[dat_j$prey_name==prey_j[j]] <- 1
    dat_j$choice <- j
    new_i <- rbind(new_i, dat_j)
  }
  stacked_data <- rbind(stacked_data, new_i)
  print(samp_ids[i])
}
d7 <- stacked_data %>% 
  mutate(sample_choice=paste(sample_id, choice, sep="-"),
         sample_choice_int=as.integer(as.factor(sample_choice))) %>%
  select(sample_id, choice_number=choice, sample_choice, 
         sample_choice_int, everything())
d7 <- d7 %>% mutate_if(is.character, as.factor)
summary(d7)

# more data formatting
head(as.data.frame(d7))
unique(d7$mig_strat)
d8 <- d7 %>% 
  dplyr::select(sharp_id=sample_id, sharp_sex=sex, trial_number=sample_choice_int, 
         week_of_yr, 
         prey_name, prey_use=use, prey_avail, 
         prey_mass=log10_mass, mig_strat,
         prey_feeder=feeder_bird, prey_urban=urban, 
         prey_habitat=habitat, prey_diet=diet, sharp_age=age,
         sharp_crop=crop) %>%
  mutate(sharp_id=as.integer(factor(sharp_id)),
         sharp_sex=as.integer(ifelse(sharp_sex=="m", 0, 1)), # make 0, 1
         prey_mig=factor(ifelse(mig_strat=="Resident", 0, 1)), # make dichotomous
         prey_wood=factor(ifelse(prey_habitat=="Woodland", 1, 0)), # make dichotomous
         prey_feeder=ifelse(prey_feeder=="Yes", 1, 0), # make dichotomous
         prey_urban=factor(ifelse(prey_urban=="Yes", 1, 0)), # make dichotomous
         prey_idx=as.integer(factor(prey_name)),
         trial_idx=as.integer(factor(trial_number))) %>% # indices for analysis
  select(-mig_strat)
head(as.data.frame(d8))
str(d8)

# split sexes, ad indiv random effect
id_male <- d8 %>% filter(sharp_sex==0) %>% 
  mutate(sharp_id=as.numeric(factor(sharp_id)),
         trial_number=as.numeric(factor(trial_number)),
         trial_idx=as.numeric(factor(trial_idx))) %>%
  mutate(sharp_re=ifelse(prey_name=="Hermit Thrush", sharp_id, NA))
id_fem <- d8 %>% filter(sharp_sex==1) %>% 
  mutate(sharp_id=as.numeric(factor(sharp_id)),
         trial_number=as.numeric(factor(trial_number)),
         trial_idx=as.numeric(factor(trial_idx))) %>%
  mutate(sharp_re=ifelse(prey_name=="Hermit Thrush", sharp_id, NA))
# ------------------------------------------------------------------------------





# maxlik conditional logit discrete choice models with mlogit ------------------
library(mlogit)
library(lmtest)
library(sandwich)

# male model mlogit
head(id_male)
md_male <- id_male %>% 
  dplyr::select(case=trial_idx, alt=prey_name, choice=prey_use, 
                prey_avail, prey_mass, id=sharp_id, 
                prey_mig, prey_habitat, prey_diet, prey_wood, prey_urban) 
md_male <- dfidx(md_male, choice="choice", idx = list(c("case", "id")),
                 idnames = c("chid", "alt")) # format for mlogit
mm1 <- mlogit(choice ~ 0 + prey_avail + prey_mass + 
                prey_mig + prey_habitat + prey_diet + prey_urban | 0 | 0, 
              data = md_male) # model run
summary(mm1) # model summary
coeftest(mm1,vcov=sandwich) # robust tests that account for repeat measures
data.frame(species=unique(id_male$prey_name), # raw proportions per species
           proportion=as.numeric(mm1$freq/sum(mm1$freq))) %>% 
  arrange(desc(proportion))
data.frame(species=unique(id_male$prey_name), # residuals per species
           avg_resid=as.numeric(apply(mm1$residuals, 2, mean))) %>% 
  arrange(desc(avg_resid))

# female model mlogit
head(id_fem)
md_fem <- id_fem %>% 
  dplyr::select(case=trial_idx, alt=prey_name, choice=prey_use, 
                prey_avail, prey_mass, id=sharp_id, 
                prey_mig, prey_habitat, prey_diet, prey_wood, prey_urban)
md_fem <- dfidx(md_fem, choice="choice", idx = list(c("case", "id")),
                idnames = c("chid", "alt")) # format for mlogit
mm2 <- mlogit(choice ~ 0 + prey_avail + prey_mass + I(prey_mass^2) +
                prey_mig + prey_habitat + prey_diet + prey_urban, 
              data = md_fem) # model run
summary(mm2) # model summary
coeftest(mm2,vcov=sandwich) # robust tests that account for repeat measures
data.frame(species=unique(id_fem$prey_name), # raw proportions per species
           proportion=as.numeric(mm2$freq/sum(mm2$freq))) %>% 
  arrange(desc(proportion))
data.frame(species=unique(id_fem$prey_name), # residuals per species
           avg_resid=as.numeric(apply(mm2$residuals, 2, mean))) %>% 
  arrange(desc(avg_resid))


# random model that takes forever
# mm3 <- mlogit(choice ~ prey_avail + prey_mass + I(prey_mass^2) +
#                 prey_mig + prey_urban + prey_habitat + prey_diet | - 1, md_fem,
#               panel = TRUE, rpar = c(prey_avail = "n", prey_mass = "n"), R = 100,
#               correlation = FALSE, halton = NA, method = "bhhh")
# ------------------------------------------------------------------------------





# bayesian conditional logit discrete choice models with inla ------------------
library(INLA)
library(brinla)

# male model with predictors for preference
head(as.data.frame(id_male))
form2 = prey_use ~ -1 + prey_avail + prey_mass + # no global intercept
  prey_mig + prey_habitat + prey_diet + prey_urban + # predictors as fixed effects
  f(prey_idx, fixed = T, constr = T) + # prey intercepts makes it conditional
  f(trial_idx, fixed = T, initial = -10) + # trial intercepts makes it conditional
  f(sharp_re, model="iid", constr=T) # random effect per bird to account for repeats
im2 = inla(form2, data = id_male, family = 'Poisson', verbose=T)
round(im2$summary.fixed[,c(1,3,5)], 3) # look at fixed effects section for importance of predictors
brinla::bri.hyperpar.summary(im2)
res2 <- im2$summary.random$prey_idx[c(5,4,6)] %>% 
  mutate(prey_species=c(unique(as.character(id_male$prey_name)))) %>%
  arrange(desc(`0.5quant`)) %>% 
  select(prey_species, estimate_median=`0.5quant`, 
         lcl=`0.025quant`, ucl=`0.975quant`)
View(res2) # these are preference indices, after accounting for predictors

# female with predictors 
form4 = prey_use ~ -1 + prey_avail + prey_mass + I(prey_mass^2) +
  prey_mig + prey_habitat + prey_diet + prey_urban +
  f(prey_idx, fixed = T, constr = T) +
  f(trial_idx, fixed = T, initial = -10) +
  f(sharp_re, model="iid", constr=T)
im4 = inla(form4, data = id_fem, family = 'Poisson', verbose=T)
round(im4$summary.fixed[,c(1,3,5)], 3)
brinla::bri.hyperpar.summary(im4)
res4 <- im4$summary.random$prey_idx[c(5,4,6)] %>% 
  mutate(prey_species=c(unique(as.character(id_fem$prey_name)))) %>%
  arrange(desc(`0.5quant`)) %>% 
  select(prey_species, estimate_median=`0.5quant`, 
         lcl=`0.025quant`, ucl=`0.975quant`)
View(res4)
# ------------------------------------------------------------------------------



