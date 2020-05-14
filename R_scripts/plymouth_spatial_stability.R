
# Title: Spatial stability in Plymouth rock pools

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)

# read in the data
ply_dat_raw <- read_delim( here("data/plymouth_rockpool_data.csv"), delim = ";" )

# duplicate this data
ply_dat <- ply_dat_raw

# check the variable structure
str(ply_dat)

# view the data
View(ply_dat)

# aggregate the data because some rows have multiple measurements

# set up a variable names vector
var_names <- names(ply_dat)

# set up a vector of id_vars
id_vars <- c("month", "shore", "pool", "composition", "removed")

# set up a vector of response aggregate variables
agg_vars <- var_names[!(names(ply_dat) %in% id_vars)]

# aggregate pools with multiple measurements by summing them
ply_dat <- ply_dat %>%
  group_by_at(id_vars) %>%
  summarise_at(vars(agg_vars), ~sum(., na.rm = TRUE)) %>%
  ungroup()

# matches with Lars and Bagchi's script results


# create a few extra variables for this analysis

# create a species richness variable
# select only relevant columns
tot_spp <- 4

ply_dat <- 
  ply_dat %>%
  mutate(species_richness = (tot_spp - removed) ) %>%
  select(shore, pool, grid_squares, month, composition, 
         removed, species_richness, totcover, cover_nocrusts)

# reorder these data so that it makes more sense
ply_dat <- 
  ply_dat %>% 
  arrange(shore, month, composition, removed, species_richness, pool)

# have a look at the data again
ply_dat %>% View()

# how many species were removed?
ply_dat$removed %>% unique()

# create a new mixture/monoculture column
# reorder the columns again
ply_dat <- 
  ply_dat %>%
  mutate(mono_mix = if_else(removed == "3", "monoculture",
                            if_else(removed == "4", "control", 
                                    if_else(removed == "0", "max_mixture", 
                                            "mixture")))) %>%
  select(shore, pool, grid_squares, month, composition, mono_mix, removed, species_richness, totcover, cover_nocrusts)

ply_dat %>% View()

# remove the composition = "None" treatment as this is a general control
ply_dat <- 
  ply_dat %>% filter(composition != "None")

ply_dat


# calculate the CV among replicates for different richnesses

# use a function for this
cv_func <- function(x) { (sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))*100 } 

test <- 
  ply_dat %>%
  group_by(shore, month, species_richness, composition) %>%
  summarise_at(vars(c("grid_squares", "totcover")), list(cv = ~cv_func(.)) ) %>%
  ungroup()

ggplot(data = test,
       mapping = aes(x = species_richness, y = totcover_cv, colour = month)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~shore, scales = "free") +
  theme_classic()






# Code for analyses of rock pool data from Plymouth field experiment 2007-2008
# Two sites: Challaborough and Kingsand, UK
# Three sampling times: March, July, September (all 2008)
#
# Code initiated 2012
#
# Code history
# 14-06-03: code modified by Bagchi in order to include Month as a continuous
# numerical variable. New code included to make predictions.
#
# 19-12-06: reviving the effort to write a manuscript.

rm(list=ls(all=TRUE))

# set working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Plymouth_rock_pools")

# install libraries
library(lme4)
library(car)
library(ggplot2)
library(plyr)
library(pbkrtest)
library(Hmisc)
library(vegan)
library(lattice)
library(tidyverse)
library(broom)

# read in the raw data
rdat<- read.csv(here("data/plymouth_rockpool_data.csv"), sep=";")

#dim(rdat)
summary(rdat)

#str(rdat)
#View(rdat)

plot(jitter(rdat$sargassum), jitter(rdat$ceramium))
plot(jitter(rdat$f_serratus), jitter(rdat$ceramium))
plot(jitter(rdat$sargassum), jitter(rdat$u_lactuca + rdat$u_intestinalis + rdat$u_linza))
plot(jitter(rdat$f_serratus), jitter(rdat$u_lactuca + rdat$u_intestinalis + rdat$u_linza))

ggplot(subset(rdat, f_serratus > 1), aes(x = f_serratus, y = ceramium)) +
  geom_point()
dim(rdat)
dim(subset(rdat, f_serratus > 1))
summary(subset(rdat, f_serratus > 1))


# reorder composition and month data, so the order of the levels make more sense
levels(rdat$composition)
rdat$composition <- factor(rdat$composition, levels=c('All', 'FLB', 'SFB', 'SFL', 'SLB', 'B', 'F', 'L', 'S', 'None'))
rdat$month <- factor(rdat$month, levels=c('March', 'July', 'September'))

### The code below is redundant code, in which the levels with the composition variable 
### were renamed. 
# Renaming the levels of the 'composition' variable. in the imported read.csv file, the 
# levels refer to the species *not* removed from the pool (e.g. FLB means S was removed). 
# The levels are here remaned to refer to the species removed (e.g. FLB means FLB was removed). 
# rdat$composition<- revalue(rdat$composition, c("All"="none", "FLB"="S", "SLB"="F", "SFB"="L", "SFL"="B", "B"="SFL", "L"="SFB", "F"="SLB", "S"="FLB", "None"="all"))
# checking that the 'revalue' function does what I want
# mydata_raw$composition[1:20]
# rdat$composition[1:20]

# Columns for the cover of some taxa are not numeric. Fixing that:
rdat$grid_squares<- as.numeric(rdat$grid_squares)
rdat$red_crusts<- as.numeric(rdat$red_crusts)
rdat$l_saccharina<- as.numeric(rdat$l_saccharina)
rdat$dictyota_dicotoma<- as.numeric(rdat$dictyota_dicotoma)
rdat$diatoms<- as.numeric(rdat$diatoms)

# Code to aggregate data, because some pools have data in two rows. 
sumfunction <- function(x){
  x <- x[!is.na(x)]
  if(length(x) ==0) return(0)  
  if(any(x =='p')) {
    s <- sum(as.numeric(as.character(x[x!='p'])))
    s <- ifelse(s==0, 'p', s)}
  else
    s <- sum(as.numeric(as.character(x)))
  return(s)} 

rdat2<- with(rdat, aggregate(rdat[, 6:48], list(month=month, shore=shore, 
                                                pool=pool, composition=composition, removed=removed), FUN = sumfunction))

rdat2 %>% head()

# Checking the new data frame
#dim(rdat2) # 249 rows, 48 columns (rdat: 319 rows, 48 columns)
#str(rdat2)

# Check unique values for the different variables
lapply(rdat2, function(x) { unique(x) })

# Have a look at the data
View(rdat2)

# What variables do we have?
colnames(rdat2)

# Code new variables for this analysis
tot_spp <- 4

ply_dat <- rdat2 %>%
  mutate(species_richness = (tot_spp - removed) ) %>%
  select(shore, pool, grid_squares, month, composition, 
         removed, species_richness, totcover, cover_nocrusts) %>%
  as_tibble()

# Convert factors to characters
ply_dat <- ply_dat %>% mutate_if(is.factor, as.character)

# Reorder these data so that it makes more sense
ply_dat <- ply_dat %>% arrange(shore, month, composition, removed, species_richness, pool)

# Have a look at the data again
ply_dat %>% View()

# How many species were removed?
ply_dat$removed %>% unique()

# We want a new column with either a mixture or monoculture variable
ply_dat <- ply_dat %>%
  mutate(mono_mix = if_else(removed == "3", "monoculture",
                            if_else(removed == "4", "control", 
                                    if_else(removed == "0", "max_mixture", "mixture")))) %>%
  select(shore, pool, grid_squares, month, composition, mono_mix, removed, species_richness, totcover, cover_nocrusts)

ply_dat %>% View()


### Test out some code

# Remove the composition = None data
test_dat <- ply_dat %>% filter(composition != "None")

# What compositions do we have?
test_dat$composition %>% unique()

# Make a unique environment column
test_dat <- test_dat %>%
  mutate(environment = as.factor( paste(shore, month, sep = "_") ) )

levels(test_dat$environment) <- c(1:length( unique(test_dat$environment) ) )

test_dat <- test_dat %>% 
  mutate(environment = as.integer(environment))

View(test_dat)

# Define the different environments
envs <- test_dat$environment %>%
  unique()
envs

# Define the number of different environments
n_env <- test_dat$environment %>%
  unique() %>%
  length()
n_env

# Define the different combinations of the environments

env_combs <- vector("list", n_env)
names(env_combs) <- seq_along(1:n_env)

for (i in seq_along(1:n_env)) {
  env_combs[[i]] <- combn(x = envs, m = i) %>%
    as_tibble() %>%
    lapply(function(x) {x})
}

env_combs <- env_combs %>% 
  unlist(recursive = FALSE) %>%
  unname()

# Now we have a list where each element is a different vector of environmental combinations
env_combs

# Replicate the dataframe the same number of times as the environment
replicate(length(env_combs), test_dat, simplify = FALSE) %>%
  head()

# Write a function to take the mean of each composition and then plot tot_cover versus richness
# for each combination of environments

bef_slope <- function(data, envs) {
  
  x <- data %>% 
    filter(environment %in% envs) %>%
    group_by(mono_mix, species_richness, composition) %>%
    summarise(ecosystem_function = mean(totcover)) %>%
    ungroup()
  
  bef_slope <- lm(ecosystem_function ~ species_richness, 
                  data = x) %>%
    tidy() %>%
    filter(term == "species_richness") %>%
    pull(estimate)
  
  ave_oy <- x %>%
    filter(mono_mix %in% c("max_mixture", "monoculture")) %>%
    group_by(mono_mix) %>%
    summarise(ecosystem_function_aoy = mean(ecosystem_function)) %>%
    ungroup() %>%
    spread(mono_mix, ecosystem_function_aoy) %>%
    mutate(average_overyielding = (max_mixture/monoculture) ) %>%
    select(-max_mixture, -monoculture) %>%
    pull(average_overyielding)
  
  t_oy <- x %>%
    filter(mono_mix %in% c("max_mixture", "monoculture")) %>%
    group_by(mono_mix) %>%
    summarise(ecosystem_function_aoy = max(ecosystem_function)) %>%
    ungroup() %>%
    spread(mono_mix, ecosystem_function_aoy) %>%
    mutate(transgressive_overyielding = (max_mixture/monoculture) ) %>%
    select(-max_mixture, -monoculture) %>%
    pull(transgressive_overyielding)
  
  tibble(environment = paste(envs, collapse = ""),
         BEF_slope = bef_slope,
         average_overyielding = ave_oy,
         transgressive_overyielding = t_oy)
  
}


# Use mapply to run bef_slope function on the test data for each environmental combination
out <- mapply(bef_slope, 
       replicate(length(env_combs), test_dat, simplify = FALSE),
       env_combs, SIMPLIFY = FALSE) %>%
  bind_rows()

out <- out %>% mutate(n_environments = nchar(environment)) %>%
  select(environment, n_environments, BEF_slope, average_overyielding, transgressive_overyielding)

gather(out, BEF_slope, average_overyielding, transgressive_overyielding,
       key = "eco_func", value = "value") %>%
  ggplot(mapping = aes(x = n_environments, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ eco_func, scales = "free") +
  theme_classic()
  









### Quick and dirty BEF-slope analysis
test <- ply_dat %>% filter(month == "July") %>%
  filter(composition != "None")

# Summarise data by shore and also by month
test <- test %>% group_by(shore, composition) %>%
  summarise(species_richness = mean(species_richness),
            totcover_m = mean(totcover)) %>%
  ungroup() %>%
  filter(composition %in% c("All", "B", "F", "L", "S")) %>%
  mutate(environment = if_else(shore == "Challaborough", 1, 2))

# Define the different environments
envs <- test$environment %>%
  unique()

# Define the number of different environments
n_env <- test$environment %>%
  unique() %>%
  length()

# Define the different combinations of the environments

env_combs <- vector("list", n_env)
names(env_combs) <- seq_along(1:n_env)

for (i in seq_along(1:n_env)) {
  env_combs[[i]] <- combn(x = envs, m = i) %>%
    as_tibble() %>%
    lapply(function(x) {x})
}

env_combs <- env_combs %>% 
  unlist(recursive = FALSE) %>%
  unname()

# Now we have a list where each element is a different vector of environmental combinations
env_combs

# Replicate the dataframe the same number of times as the environment
replicate(length(env_combs), test, simplify = FALSE)

# Write a function to take the mean of each composition and then plot tot_cover versus richness
# for each combination of environments

bef_slope <- function(data, envs) {
  x <- data %>% filter(environment %in% envs) %>%
    group_by(species_richness, composition) %>%
    summarise(ecosystem_function = mean(totcover_m)) %>%
    ungroup()
  
  bef_slope <- lm(ecosystem_function ~ species_richness, 
                  data = x) %>%
    tidy() %>%
    filter(term == "species_richness") %>%
    pull(estimate)
  
}

# Use mapply to run bef_slope function on the test data for each environmental combination
mapply(bef_slope, 
       replicate(length(env_combs), test, simplify = FALSE),
       env_combs)
