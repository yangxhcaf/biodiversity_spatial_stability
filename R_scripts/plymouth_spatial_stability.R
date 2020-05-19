
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

# aggregate the data because some rows have multiple measurements

# set up a variable names vector
var_names <- names(ply_dat)

# set up a vector of id_vars
id_vars <- c("month", "shore", "pool", "composition", "removed")

# set up a vector of response aggregate variables
agg_vars <- var_names[!(names(ply_dat) %in% id_vars)]

# aggregate pools with multiple measurements by summing them
ply_dat <- 
  ply_dat %>%
  group_by_at(id_vars) %>%
  summarise_at(vars(agg_vars), ~sum(., na.rm = TRUE)) %>%
  ungroup()

# matches with Lars and Bagchi's script results


# examine which 'species' classes have notable cover

# analyse only the focal species
foc_spp <- 
  c("sargassum", "bifurcaria", "f_serratus", "l_digitata")

# analyse all species (or categories) with more than 5 cover in any data point
foc_spp <-
  ply_dat %>%
  mutate_at(vars(agg_vars[agg_vars != "grid_squares"]), ~(./totcover)*100) %>%
  select(agg_vars[!(agg_vars %in% c("grid_squares", "totcover", "cover_nocrusts", "bare rock"))] ) %>%
  lapply(., function(x) { 
    if_else(x > 5, 1, 0) %>% sum()
    } ) %>%
  bind_rows() %>%
  gather() %>%
  filter(value > 0) %>%
  pull(key)


# create a few extra variables for this analysis

# create a species richness variable
# select only relevant columns
tot_spp <- 4

ply_dat <- 
  ply_dat %>%
  mutate(species_richness = (tot_spp - removed) ) %>%
  select(id_vars, grid_squares, species_richness, totcover, cover_nocrusts, foc_spp)

# create a new mixture/monoculture column
# reorder the columns again
ply_dat <- 
  ply_dat %>%
  mutate(mono_mix = if_else(removed == "3", "monoculture",
                            if_else(removed == "4", "control", 
                                    if_else(removed == "0", "max_mixture", 
                                            "mixture"))))


# remove the composition = "None" treatment as this is a general control
ply_dat <- 
  ply_dat %>% filter(composition != "None")


### stability across different levels of organisation

# create a total cover variable from the constituent species
ply_dat$total_cover <- 
  ply_dat %>%
  select(foc_spp) %>%
  rowSums()

# sort data to match with Lamy et al.'s function
ply_dat <-
  ply_dat %>%
  arrange(shore, composition, pool, month)


### calculate stability at different scales

# my own functions to calculate these different quantities

# create a test dataset to see if I can reproduce these metrics (Lamy et al. 2019)
stab_test <- 
  ply_dat %>%
  filter(composition == "All")

# metacommunity variability (CVmc)

stab_test %>%
  group_by(month) %>%
  summarise(total_cover = sum(total_cover, na.rm = TRUE)) %>%
  ungroup() %>%
  summarise(CVmc = (sd(total_cover, na.rm = TRUE)/mean(total_cover, na.rm = TRUE)))


# metapopulation variability (CVmp)
stab_test %>%
  gather(foc_spp, key = "species", value = "cover") %>%
  group_by(species, month) %>%
  summarise(cover = sum(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(sppcv = (sd(cover, na.rm = TRUE)/mean(cover, na.rm = TRUE)),
            mean_cover = mean(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  summarise(CVmp = weighted.mean(sppcv, w = mean_cover, na.rm = TRUE))

# local community variability (CVlc)
stab_test %>%
  group_by(pool) %>%
  summarise(pool_cv = (sd(total_cover, na.rm = TRUE)/mean(total_cover, na.rm = TRUE)),
            mean_cover = mean(total_cover, na.rm = TRUE) ) %>%
  ungroup() %>%
  summarise(CVlc = weighted.mean(pool_cv, w = mean_cover, na.rm = TRUE))

# local population variability (CVlp)
stab_test %>%
  gather(foc_spp, key = "species", value = "cover") %>%
  group_by(pool, species) %>%
  summarise(spp_cv = sd(cover, na.rm = TRUE)/mean(cover, na.rm = TRUE),
            mean_cover = mean(cover, na.rm = TRUE) ) %>%
  ungroup() %>%
  group_by(pool) %>%
  summarise(spp_pool_cv = weighted.mean(spp_cv, w = mean_cover, na.rm = TRUE),
            patch_cover = mean(mean_cover, na.rm = TRUE)) %>%
  ungroup() %>%
  summarise(CVlp = weighted.mean(spp_pool_cv, w = patch_cover))


### calculate synchrony values

# average local-scale species synchrony: phi_SC_L
full_join(stab_test %>%
            group_by(pool) %>%
            summarise(sd_patch = sd(total_cover, na.rm = TRUE)),
          stab_test %>%
            gather(foc_spp, key = "species", value = "cover") %>%
            group_by(pool, species) %>%
            summarise(sd_spp_patch = sd(cover, na.rm = TRUE)) %>%
            ungroup() %>%
            group_by(pool) %>%
            summarise(sd_spp_pool = sum(sd_spp_patch, na.rm = TRUE)),
          by = "pool") %>%
  mutate(phi_SC_k = sd_patch/sd_spp_pool) %>%
  summarise(phi_SC_L = weighted.mean(phi_SC_k, w = sd_spp_pool))
  
# spatial synchrony of species i: phi_i_LR
# phi_i_LR <- sd_iT/colSums(sd_ik)

stab_test %>%
  gather(foc_spp, key = "species", value = "cover") %>%
  group_by(species, month) %>%
  summarise(cover = sum(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(sppcv = (sd(cover, na.rm = TRUE)/mean(cover, na.rm = TRUE)),
            mean_cover = mean(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  pull(sppcv)

stab_test %>%
  gather(foc_spp, key = "species", value = "cover") %>%
  group_by(pool, species) %>%
  summarise(sd_spp_patch = sd(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(sd_spp_patch = sum(sd_spp_patch, na.rm = TRUE)) %>%
  pull(sd_spp_patch)

# this requires a weighted average...


### Lamy et al. (2019) functions

# these functions need measurements for all months for all pools
comp_pools <- 
  ply_dat %>%
  mutate(shore_comp_id = paste(shore, composition, sep = "_")) %>%
  group_by(shore_comp_id, pool) %>%
  summarise(comp = unique(month) %>% length()) %>%
  filter(comp < 3) %>%
  pull(pool)

comp_pools

# create a data set with only complete pools  
stab_dat <- 
  ply_dat %>%
  filter( !(pool %in% comp_pools) ) 

# at what scale should calculations be done?
# within shore or between shore?

scale <- c("between")

stab_dat <- 
  if (scale == "within") {
    stab_dat %>%
      mutate(shore_comp_id = paste(shore, composition, sep = "_")) %>%
      split(., .$shore_comp_id)
    } else { 
      split(stab_dat, stab_dat$composition) 
      }

stab_dat[[1]]

stab_out <- 
  lapply(stab_dat, function(x) {
  
  y <- part_stab_comp(Y = select(x, foc_spp),
                 s = length(unique(x$pool)),
                 t = length(unique(x$month)))
  y$part
  
} )

# join these stability metrics to the richness data

stab_out <- 
  stab_out %>%
  bind_rows(.id = "comp_id") %>% 
  as_tibble()

# create moderator variables
mod_vars <- 
  if (scale == "within") {
    ply_dat %>%
      select(shore, composition, species_richness) %>%
      mutate(comp_id = paste(shore, composition, sep = "_")) %>%
      distinct()
  } else {
    ply_dat %>%
      select(composition, species_richness) %>%
      distinct() %>%
      rename(comp_id = "composition")
  }
  
stab_full <- 
  full_join(mod_vars, stab_out, c("comp_id"))


stab_full %>%
  gather(-comp_id, -species_richness,
         key = "stability", value = "metric") %>%
  ggplot(data = .,
         mapping = aes(x = species_richness, y = metric)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~stability, scales = "free") +
  theme_classic()














### spatial stability

# calculate the CV among replicates for different richnesses

# write a function to calculate the coefficient of variation
cv_func <- function(x) { (sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))*100 } 


### low heterogeneity among replicates pools

# calculate cv among replicate pools within each month and shore
with_month <- 
  ply_dat %>%
  group_by(shore, month, species_richness, composition) %>%
  summarise_at(vars(c("grid_squares", "totcover")), list(cv = ~cv_func(.), n = ~n()) ) %>%
  ungroup()

ggplot(data = with_month,
       mapping = aes(x = species_richness, y = totcover_cv, colour = month)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~shore, scales = "free") +
  theme_classic()

ggplot(data = with_month,
       mapping = aes(x = month, y = grid_squares_cv)) +
  geom_point() +
  facet_wrap(~shore, scales = "free") +
  theme_classic()

# calculate mean cv among replicate pools across months within each shore
mean_month <- 
  ply_dat %>%
  group_by(shore, month, species_richness, composition) %>%
  summarise(totcover_cv = cv_func(totcover)) %>%
  ungroup() %>%
  group_by(shore, species_richness, composition) %>%
  summarise(totcover_cv = mean(totcover_cv, na.rm = TRUE) ) %>%
  ungroup()

ggplot(data = mean_month,
       mapping = aes(x = species_richness, y = totcover_cv)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~shore, scales = "free") +
  theme_classic()


### high heterogeneity among replicate pools from different shores

mean_shore <- 
  ply_dat %>%
  group_by(shore, month, species_richness, composition) %>%
  summarise_at(vars(c("totcover")), ~mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(shore, species_richness, composition) %>%
  summarise(totcover = mean(totcover, na.rm = TRUE) ) %>%
  ungroup() %>%
  group_by(species_richness, composition) %>%
  summarise(totcover_cv = cv_func(totcover)) %>%
  ungroup()

ggplot(data = mean_shore,
       mapping = aes(x = species_richness, y = totcover_cv)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()



### Test out some code for scaling up as previously:

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
