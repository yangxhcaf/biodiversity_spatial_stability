
# Title: simulation of spatial stability in heterogenous environments

# load relevant libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)

# define number of species and environments
specnum <- 3
envnum <- 3

# set up species matrix
spec <- paste("strain_", c(1:specnum), sep = "")
spec_comb <- sapply(c(1:specnum), function(x) combn(spec, x)) 
names(spec_comb) <- paste("richness", c(1:specnum))

spec_comb

# set up environmental matrix
env <- paste("Environment_", c(1:envnum), sep = "")
env_comb <- sapply(c(1:envnum), function(x) combn(env, x)) 
names(env_comb) <- paste("Heterogeneity", c(1:envnum))

env_comb

set.seed(1346)
Traits <- matrix(
  nrow = length(env),
  ncol = length(spec),
  dimnames = list(env, spec))
Traits[1:nrow(Traits), 1:ncol(Traits)] <-
  runif( n = prod( dim( Traits)), 0, 1)
Traits


plot_values <- sapply(env, #for each environment
                      function(x) {
                        lapply(spec_comb, #for each richness level
                               function(y) {
                                 apply(y, 2, function(z){ #for each species combination
                                     mean(sum(Traits[x,z]^2) / abs(sum(Traits[x,z])))
                                 }
                                 )
                               }
                        )
                      }
)

plot_values_df <- 
  tibble(
    richness = rep( rep(1:specnum, unlist(lapply(spec_comb, ncol))), envnum),
    spec_comb = rep(lapply(spec_comb, 
                           function(x) apply(x, 2, function(y) 
                             paste(y, collapse = " "))) %>% unlist(),
                    envnum),
    environment = rep(env, each = sum(choose(specnum,1:specnum))),
    functioning = unlist(plot_values)
  )

plot_values_df

### landscapes

env <- paste("env_", c(1:envnum), sep = "")
env_comb <- rep(list(env), envnum) %>%
  expand.grid() %>%
  as_tibble() %>%
  mutate_all(~as.character(.))
names(env_comb) <- paste("patch", c(1:envnum), sep = "_")

env_comb <- env_comb %>%
  mutate(landscape = seq_along(1:nrow(env_comb)))

env_comb <- env_comb %>%
  gather(-landscape, key = "patch", value = "environment") %>%
  arrange(landscape) %>%
  group_by(landscape) %>%
  mutate(heterogeneity = unique(environment) %>% length() ) %>%
  ungroup()








