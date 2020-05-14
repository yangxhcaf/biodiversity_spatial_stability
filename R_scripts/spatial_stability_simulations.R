
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


plot_values <- sapply(env, #for each environment
                      function(x) {
                        lapply(spec_comb, #for each richness level
                               function(y) {
                                 apply(y, 2, function(z){ #for each species combination
                                     mean(sum(spp_traits[x,z]^2) / abs(sum(spp_traits[x,z])))
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

land <- paste("env_", c(1:envnum), sep = "")
land_comb <- rep(list(land), envnum) %>%
  expand.grid() %>%
  as_tibble() %>%
  mutate_all(~as.character(.))
names(land_comb) <- paste("patch", c(1:envnum), sep = "_")

land_dat <- land_comb %>%
  mutate(landscape = seq_along(1:nrow(land_comb)))

land_dat <- 
  land_dat %>%
  gather(-landscape, key = "patch", value = "environment") %>%
  arrange(landscape) %>%
  group_by(landscape) %>%
  mutate(heterogeneity = unique(environment) %>% length() ) %>%
  ungroup()

land_dat

plot_values_df

as.vector(DF[1,])

# make a list to output functioning to
land_dat <- split(select(land_dat, -landscape), land_dat$landscape)
land_dat

# make landscape combination dataframe
as.vector(land_comb[1,])

unname( unlist(land_comb[1,]) )

rowMeans(plot_values_wide[,y])


lapply(land_comb, function(x) { pull(x, environment) } )

env_id <- pull(land_comb[[5]], environment)
env_id

land_comb[[5]]

funcs <- plot_values_df %>% filter( environment %in% env_id )

funcs

plot_values_wide <- plot_values_df %>%
  spread(environment, functioning)


env_values <- 
  lapply(land_comb, function(x) { #for all scales (1:5)
    apply(x, 2, function(y) { #for all environmental combinations
      
      df <- data.frame( #average functioning for this environment combination
        env_comb = paste0(y, collapse = " "),
        envnum = length(y),
        richness = plot_values_wide$richness,
        functioning = rowMeans(plot_values_wide[,y])
      )
    })
  })

apply(x, 2, function(y) { #for all environmental combinations
  
  df <- data.frame( #average functioning for this environment combination
    env_comb = paste0(y, collapse = " "),
    envnum = length(y),
    richness = plot_values_wide$richness,
    functioning = rowMeans(plot_values_wide[,y])
  )
})




sample_n(size = envnum, replace = FALSE, weight = NULL, .env = NULL)





