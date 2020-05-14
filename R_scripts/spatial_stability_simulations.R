
# Title: simulation of spatial stability in heterogenous environments

# Code adapted from: https://github.com/FabianRoger/Bacteria_BEF_Scale/blob/master/Simulations.Rmd

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
specnum <- 5
envnum <- 5

# set up species matrix
spec <- paste("species_", c(1:specnum), sep = "")
spec_comb <- sapply(c(1:specnum), function(x) combn(spec, x)) 
names(spec_comb) <- paste("richness", c(1:specnum))

spec_comb

# set up environmental matrix
env <- paste("env_", c(1:envnum), sep = "")
env_comb <- sapply(c(1:envnum), function(x) combn(env, x)) 
names(env_comb) <- paste("heterogeneity", c(1:envnum))

env_comb

spp_traits <- matrix(
  nrow = length(env),
  ncol = length(spec),
  dimnames = list(env, spec))
spp_traits[1:nrow(spp_traits), 1:ncol(spp_traits)] <-
  runif( n = prod( dim( spp_traits)), 0, 1)
spp_traits


# 1 = Dominance (# pure selection effect and perfect environmental sorting)
# 2 = Weighted mean (# partial selection effect)
# 3 = Random dominance (# dominance that is random with respect to function)

Scenario <- 2

plot_values <- sapply(env, #for each environment
                      function(x) {
                        lapply(spec_comb, #for each richness level
                               function(y) {
                                 apply(y, 2, function(z){ #for each species combination
                                   if(Scenario == 1){
                                     max(spp_traits[x,z])
                                   } else if(Scenario == 2){
                                     mean(sum(spp_traits[x,z]^2) / abs(sum(spp_traits[x,z])))
                                   } else {
                                     sample(spp_traits[x,z], size = 1)
                                   } 
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

land_comb <- paste("env_", c(1:envnum), sep = "")

land_comb <- 
  rep(list(land_comb), envnum) %>%
  expand.grid() %>%
  as_tibble() %>%
  mutate_all(~as.character(.))
names(land_comb) <- paste("patch", c(1:envnum), sep = "_")

land_comb <- 
  land_comb %>%
  mutate(landscape = seq_along(1:nrow(land_comb))) %>%
  gather(-landscape, key = "patch", value = "environment") %>%
  arrange(landscape) %>%
  group_by(landscape) %>%
  mutate(heterogeneity = unique(environment) %>% length() ) %>%
  ungroup()

land_comb <- 
  split(select(land_comb, -heterogeneity), land_comb$heterogeneity) %>%
  lapply(function(x) {x %>% 
    mutate(n = rep(c(1:envnum), (n()/envnum) )) %>%
    select(-patch) %>%
    spread(key = "landscape", value = "environment") %>%
    select(-n) %>%
    as.matrix() }
    )  
names(land_comb) <- paste0("het_", c(1:envnum))

# we have a list with different levels of heterogeneity all combinations of environments in columns
plot_values_wide <- plot_values_df %>%
  spread(environment, functioning)


func_values <- 
  lapply(land_comb, function(x) { #for all scales (1:5)
    apply(x, 2, function(y) { #for all environmental combinations
      
      df <- data.frame( #average functioning for this environment combination
        env_comb = paste0(y, collapse = " "),
        envnum = unique(y) %>% length(),
        richness = plot_values_wide$richness,
        functioning_mean = rowMeans(plot_values_wide[, y]),
        functioning_sd = apply(plot_values_wide[, y], 1, sd, na.rm = TRUE)
      )
    })
  })

# bind these data into a single dataframe
# remove highest richness level as this is unreplicated (combinations wise)
func_values <- 
  bind_rows( unlist(func_values, recursive = FALSE), .id = "landscape_id" ) %>%
  as_tibble() %>%
  filter(richness < specnum)

# add small random errors to both the functioning_mean and functioning_sd
# calculate coefficient of variation
func_values <- 
  func_values %>%
  mutate(functioning_mean = functioning_mean + rnorm(n(), 0, 0.01),
         functioning_sd = functioning_sd + rnorm(n(), 0, 0.01)) %>%
  mutate(functioning_cv = (functioning_sd/functioning_mean)*100 )


# plot some of these data

ggplot(data = func_values[sample(seq(1:nrow(func_values)), size = 10000),  ],
       mapping = aes(x = richness, y = functioning_cv, colour = as.character(envnum))) +
  geom_jitter(alpha = 0.1, width = 0.1) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_d() +
  theme_classic()







