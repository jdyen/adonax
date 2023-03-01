# Analysis of the effects of Arundo donax on fish assemblages in
#   north-eastern Spain
#
# Authors: Alberto Maceda Veiga, Ralph Mac Nally, Jian Yen
# Analysis: Jian Yen
#
# last updated: 24 February 2023

# load packages
library(qs)
library(readxl)
library(dplyr)
library(tidyr)
library(rstan)
library(ggplot2)
library(patchwork)
library(tidybayes)
library(ROCR)

# set a flag to re-run models
refit_models <- FALSE

# load some helper functions
source("R/utils.R")

# load species info
spp_info <- read_xlsx(
  "data/BD_Fish_Arundo.xlsx", sheet = "Metadata and other info",
  range = "A1:U30"
)

# load fish guild info
guilds <- read_xlsx(
  "data/FISH GUILDS.xlsx",
  skip = 1,
  col_names = c("species", "habitat", "diet", "breeding", "migration", "habitat_use")
)

# rename one species that has an updated name
guilds <- guilds %>%
  mutate(species = gsub("LUGR", "BAGR", species))

# remove second (incorrect) ebro basin classifications
#   and rename other one
spp_info <- spp_info %>% 
  select(-Ebro...14) %>% 
  rename(Ebro = Ebro...8)

# load fish cpue data
cpue <- read_xlsx(
  "data/BD_Fish_Arundo.xlsx", 
  sheet = "BD_env_CPUE_richness"
)

# rename some predictors to avoid naming issues
cpue <- cpue %>%
  rename(
    ammonia_no2 = `[NH3+NO2]`,
    nitrate_phosphate = `[NO3+PO4]`,
    Channel_conservation =  `Channel conservation`,
    Habitat_diversity = `Habitat diversity`
  )

# and remove NA values
cpue <- cpue %>% filter(!apply(cpue, 1, anyNA))

# and load fish individual SMI data
smi <- read_xlsx(
  "data/BD_Fish_Arundo.xlsx", 
  sheet = "BD_SMI per ind_sp_site"
)

# load list of paired sites for model 1
paired_sites <- read_xlsx("data/BD_Arundo.xlsx")

# prepare data for model 0.1 and model 0.2
model1_predictors <- cpue %>%
  select(
    Stream_order, Elevation, Water_temperature,
    Water_depth, Water_velocity, Conductivity,
    pH, ammonia_no2, nitrate_phosphate,
    Channel_conservation, Habitat_diversity
  )
model1_data <- list(
  N = nrow(cpue),
  K = ncol(model1_predictors),
  y = cpue %>% pull(Arundo),
  X = scale(model1_predictors),
  basin = cpue %>% pull(Basin) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_random = 2.
)
model1_data$nbasin <- max(model1_data$basin)

# filter data to paired sites
paired_fish <- paired_sites %>%
  left_join(cpue %>% select(-Arundo), by = c("CODE1", "CODE2"))

# extract CPUE info for model 1
paired_cpue <- paired_fish %>%
  select(
    ACAR, ALAL, AMME, ANAN, BAQU, BAGR, BAHA, BAME,
    CAAU, COCA, COPA, CYCA, ESLU, GAHO, GOLO, LEGI,
    MISA, ONMY, PAMY, PHBI, PSPA, RURU, SAFL, SATR,
    SALU, SCER, SIGL, SQLA, GAGY
  )

# drop species that are never observed
paired_cpue <- paired_cpue %>% select(which(colSums(paired_cpue) > 0))

# attach guild info to species-level CPUE measurements
cpue_guild <- paired_cpue %>%
  pivot_longer(cols = everything(), names_to = "species") %>%
  mutate(id = rep(seq_len(nrow(paired_cpue)), each = ncol(paired_cpue))) %>%
  left_join(guilds, by = "species")

# extract origin info for guild groupings
origin_info <- spp_info %>% 
  select(-Order, -Family, -Genus, -`TROPHIC LEVEL`) %>%
  pivot_longer(
    cols = c(
      Besos, Daro, Ebro, Fluvia, Francoli, Foix, Gaia,
      Garona, Llobregat, Muga, Riudecanyes, Ridaura,
      Senia, Ter, Tordera
    ),
    names_to = "Basin",
    values_to = "Origin"
  )

# and predictor variables
paired_predictors <- paired_fish %>%
  select(
    Stream_order, Elevation, Water_temperature,
    Water_depth, Water_velocity, Conductivity, pH, 
    ammonia_no2, nitrate_phosphate, Channel_conservation,
    Habitat_diversity
  )

# reformat origin info for species model
origin_matrix <- origin_info %>%
  filter(
    SPECIES %in% colnames(paired_cpue),
    Basin %in% unique(paired_fish$Basin)
  ) %>%
  mutate(origin_level = rebase_factor(Origin)) %>%
  select(-Origin) %>%
  pivot_wider(id_cols = SPECIES, names_from = Basin, values_from = origin_level) %>%
  select(-SPECIES)

# dump this all together in a list for all species and for each guild
model2a_data <- list(
  N = nrow(paired_cpue),
  Q = ncol(paired_cpue),
  K = ncol(paired_predictors),
  X = t(scale(paired_predictors)),
  arundo = paired_fish %>% pull(Arundo),
  origin = origin_matrix,
  basin = paired_fish %>% pull(Basin) %>% rebase_factor(),
  block_id = paired_fish %>% pull(BLOCK) %>% rebase_factor(),
  site = paired_fish %>% pull(CODE1) %>% rebase_factor(),
  scale_factor = 10.,
  sigma_fixed = 5.,
  sigma_covar = 2.,
  sigma_resid = 2.,
  sigma_random = 2.
)

# add random effect counts
model2a_data$norigin <- max(model2a_data$origin)
model2a_data$nbasin <- max(model2a_data$basin)
model2a_data$nblock <- max(model2a_data$block_id)
model2a_data$nsite <- max(model2a_data$site)

# and add response variable and zero identifiers
model2a_data$yflat <- round(c(t(paired_cpue)) * model2a_data$scale_factor)
model2a_data$nflat <- length(model2a_data$yflat)

# grab summary stats for model 2
paired_richness <- paired_fish %>%
  select(contains("Sp_Rich_"))

# and calculate abundance-weighted trophic level for model 2
wtl <- paired_cpue %>% 
  sweep(1, rowSums(paired_cpue), "/") %>% 
  sweep(2, spp_info$`TROPHIC LEVEL`[match(colnames(paired_cpue), spp_info$SPECIES)], "*") %>%
  apply(1, sum)
wtl[is.na(wtl)] <- median(wtl, na.rm = TRUE)

# combine into a list
model2b_data <- list(
  N = nrow(paired_richness),
  Q = ncol(paired_richness),
  K = ncol(paired_predictors),
  X = t(scale(paired_predictors)),
  arundo = paired_fish %>% pull(Arundo),
  basin = paired_fish %>% pull(Basin) %>% rebase_factor(),
  block_id = paired_fish %>% pull(BLOCK) %>% rebase_factor(),
  site = paired_fish %>% pull(CODE1) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_covar = 2.,
  sigma_resid = 2.,
  sigma_random = 2.,
  scale_factor = 1.
)

# add random effect counts
model2b_data$nbasin <- max(model2b_data$basin)
model2b_data$nblock <- max(model2b_data$block_id)
model2b_data$nsite <- max(model2b_data$site)

# and add response variable and zero identifiers
model2b_data$yflat <- c(t(paired_richness))
model2b_data$nflat <- length(model2b_data$yflat)

# and moving onto data for model 3 (individual sizes)
smi_spp_info <- spp_info %>% 
  rename(trophic_level = `TROPHIC LEVEL`) %>%
  pivot_longer(
    cols = c(
      Besos, Daro, Ebro, Fluvia, Francoli, Foix, Gaia,
      Garona, Llobregat, Muga, Riudecanyes, Ridaura,
      Senia, Ter, Tordera
    ),
    names_to = "Basin",
    values_to = "Origin"
  )
smi_full <- smi %>%
  filter(as.character(CODE1) %in% as.character(paired_sites$CODE1)) %>%
  left_join(paired_sites %>% select(CODE1, BLOCK, Arundo), by = "CODE1") %>%
  left_join(
    cpue %>% 
      select(
        CODE1, Basin, Stream_order, Elevation, Water_temperature,
        Water_depth, Water_velocity, Conductivity, pH, 
        ammonia_no2, nitrate_phosphate, Channel_conservation,
        Habitat_diversity
      ),
    by = "CODE1"
  ) %>%
  left_join(smi_spp_info, by = c("SPECIES", "Basin"))
smi_full <- smi_full %>%
  filter(!apply(smi_full, 1, anyNA))

# extract predictor variables
smi_predictors <- smi_full %>%
  select(
    Stream_order, Elevation, Water_temperature,
    Water_depth, Water_velocity, Conductivity, pH, 
    ammonia_no2, nitrate_phosphate, Channel_conservation,
    Habitat_diversity
  )

# put this all in a list
model2c_data <- list(
  N = nrow(smi_full),
  K = ncol(smi_predictors),
  y = smi_full %>% pull(SMI),
  X = scale(smi_predictors),
  arundo = smi_full %>% pull(Arundo),
  origin = smi_full %>% select(SPECIES, Basin) %>% left_join(smi_spp_info, by = c("SPECIES", "Basin")) %>% select(Origin) %>% unlist() %>% rebase_factor(),
  basin = smi_full %>% pull(Basin) %>% rebase_factor(),
  site = smi_full %>% pull(CODE1) %>% rebase_factor(),
  block_id = smi_full %>% pull(BLOCK) %>% rebase_factor(),
  species = smi_full %>% pull(SPECIES) %>% rebase_factor(),
  order = smi_full %>% pull(Order) %>% rebase_factor(),
  family = smi_full %>% pull(Family) %>% rebase_factor(),
  genus = smi_full %>% pull(Genus) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_random = 2.
)

# add counts for all random effects
model2c_data$norigin <- max(model2c_data$origin)
model2c_data$nbasin <- max(model2c_data$basin)
model2c_data$nsite <- max(model2c_data$site)
model2c_data$nblock <- max(model2c_data$block_id)
model2c_data$norder <- max(model2c_data$order)
model2c_data$nfamily <- max(model2c_data$family)
model2c_data$ngenus <- max(model2c_data$genus)
model2c_data$Q <- max(model2c_data$species)

# settings for all MCMC models
seed <- 352124142
iter <- 10000
warmup <- floor(iter / 2)
thin <- 4
chains <- 4
cores <- 4

# fit all models if required
if (refit_models) {
  
  # Fit model 0: Bayesian logistic regression, where does A. donax occur?
  # compile model
  model1_stan <- stan_model("src/bern.stan")
  
  # sample from model
  draws_model1 <- sampling(
    model1_stan,
    data = model1_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model1, file = paste0("outputs/fitted/draws-model1-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model1)

  # Fit model 1a: model of fish CPUE, how does A. donax affect fish species?
  # compile model
  model2a_stan <- stan_model("src/nb-simple.stan")
  
  # sample from model
  draws_model2a <- sampling(
    model2a_stan,
    data = model2a_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 20),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model2a, file = paste0("outputs/fitted/draws-model2a-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model2a)
  
  # sample from model for species richness by classification (native or otherwise)
  # compile model
  model2b_stan <- stan_model("src/nb-simple-no-origin.stan")

  # sample from model
  draws_model2b <- sampling(
    model2b_stan,
    data = model2b_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model2b, file = paste0("outputs/fitted/draws-model2b-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model2b)
  
  # compile model
  model2c_stan <- stan_model("src/lnorm-simple.stan")
  
  # sample from model
  draws_model2c <- sampling(
    model2c_stan,
    data = model2c_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.82, max_treedepth = 11),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model2c, file = paste0("outputs/fitted/draws-model2c-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model2c)
  
}

# calculate sample sizes with and without A. donax for models 2a and 2c
#    (species-level analyses)
model2a_sample_size <- paired_cpue %>% 
  group_by(paired_fish %>% pull(Arundo)) %>%
  summarise_all(~sum(.x > 0))
model2a_sample_size <- tibble(
  GROUP = colnames(model2a_sample_size)[-1],
  with_adonax = unlist(model2a_sample_size[2, -1]),
  without_adonax = unlist(model2a_sample_size[1, -1])
)
model2b_sample_size <- paired_richness %>% 
  group_by(paired_fish %>% pull(Arundo)) %>%
  summarise_all(~sum(.x > 0))
model2b_sample_size <- tibble(
  GROUP = colnames(model2b_sample_size)[-1],
  with_adonax = unlist(model2b_sample_size[2, -1]),
  without_adonax = unlist(model2b_sample_size[1, -1])
)
model2c_sample_size <- smi_full %>% 
  group_by(SPECIES, Arundo) %>%
  summarise(count = n()) %>%
  mutate(Arundo = ifelse(Arundo == 0, "without_adonax", "with_adonax")) %>%
  pivot_wider(
    id_cols = SPECIES,
    names_from = Arundo, 
    values_from = count
  ) %>%
  mutate(without_adonax = ifelse(is.na(without_adonax), 0, without_adonax)) %>%
  rename(GROUP = SPECIES)
sample_sizes <- rbind(
  c(GROUP = "Model 2a", with_adonax = "", without_adonax = ""),
  model2a_sample_size,
  c(GROUP = "Model 2b", with_adonax = "", without_adonax = ""),
  model2b_sample_size,
  c(GROUP = "Model 2c", with_adonax = "", without_adonax = ""),
  model2c_sample_size
)
write.csv(sample_sizes, file = "outputs/tables/Table3.csv")

# load all fitted models
file_names <- c(
  "draws-model1",
  "draws-model2a",
  "draws-model2b",
  "draws-model2c"
)

# loop over each and load most recent version
all_files <- dir("outputs/fitted/")
out <- vector("list", length = length(file_names))
for (i in seq_along(file_names)) {
  file_sub <- all_files[grepl(paste0(file_names[i], "-"), all_files)]
  file_sub <- sort(file_sub, decreasing = TRUE)[1]
  print(file_sub)
  assign(gsub("-", "_", file_names[i]), qread(paste0("outputs/fitted/", file_sub)))
}

# plot diagnostics
model_names <- c(
  "draws_model1",
  "draws_model2a",
  "draws_model2b", 
  "draws_model2c"
)
model_summaries <- vector("list", length = length(model_names))
for (i in seq_along(model_names))
  model_summaries[[i]] <- summary(get(model_names[i]))

get_diagnostics <- function(x) {
  x$summary[, c("Rhat", "n_eff")]
}
diagnostics <- lapply(model_summaries, get_diagnostics)
diagnostics <- data.frame(
  model = rep(model_names, times = sapply(diagnostics, nrow)),
  par = do.call(c, lapply(diagnostics, function(x) rownames(x))),
  Rhat = do.call(c, lapply(diagnostics, function(x) x[, "Rhat"])),
  n_eff = do.call(c, lapply(diagnostics, function(x) x[, "n_eff"]))
)
diagnostic_plot <- diagnostics %>%
  filter(!grepl("_term\\[", par)) %>%
  pivot_longer(
    cols = c(Rhat, n_eff),
    names_to = "statistic"
  ) %>%
  mutate(
    model = gsub("draws_", "", model),
    model = gsub("model", "Model ", model),
    model = gsub("Model 2", "Model 2b", model),
    model = gsub("Model 1a", "Model 2a", model),
    model = gsub("Model 3", "Model 2c", model),
    model = gsub("Model 0", "Model 1", model),
    statistic = gsub("n_eff", "Effective sample size", statistic)
  ) %>%
  ggplot() +
  geom_histogram(aes(x = value), bins = 50) +
  xlab("Value") +
  ylab("Count") +
  facet_grid(model ~ statistic, scales = "free")
ggsave(
  diagnostic_plot,
  filename = "outputs/figures/FigS1.png",
  device = png,
  width = 5,
  height = 8,
  units = "in", 
  dpi = 600
)

# create lookup table for predictor variables
predictor_names <- c(
  "Stream_order" = "Stream order",
  "Elevation" = "Elevation",
  "Water_temperature" = "Water temperature",
  "Water_depth" = "Water depth",
  "Water_velocity" = "Water velocity",
  "Conductivity" = "Conductivity",
  "pH" = "pH",
  "ammonia_no2" = "[NH3+NO2]",
  "nitrate_phosphate" = "[NO3+PO4]",
  "Channel_conservation" = "Channel conservation",
  "Habitat_diversity" = "Habitat diversity"
)

# extract parameters from model 0.a: Bayesian logistic regression, where does A. donax occur?
model1_effects <- draws_model1 %>% 
  spread_draws(
    beta[Predictor],
  ) %>% 
  median_qi(
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Predictor = colnames(model1_data$X)[Predictor],
    Predictor = predictor_names[Predictor]
  )
model1_beta <- model1_effects %>%
  ggplot(aes(y = Predictor, x = beta, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Parameter estimate")
ggsave(
  model1_beta,
  filename = "outputs/figures/Fig2.png",
  device = png,
  width = 5,
  height = 5,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 1a: model of fish CPUE, how does
#    A. donax affect fish species?
model2a_effects <- draws_model2a %>% 
  spread_draws(
    theta[species],
    beta[species, predictor],
    delta[species, predictor]
  ) %>% 
  median_qi(
    theta, beta, delta = theta + delta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = colnames(paired_cpue)[species],
    predictor = rownames(model2a_data$X)[predictor],
    Predictor = predictor_names[predictor]
  )
model2a_theta <- model2a_effects %>%
  select(contains("theta"), Species, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model2a_beta <- model2a_effects %>%
  select(contains("beta"), contains("delta"), Species, Predictor, .width, .point, .interval) %>%
  pivot_longer(cols = c(contains("beta"), contains("delta"))) %>%
  mutate(
    level = "mid",
    level = ifelse(grepl("lower", name), "lower", level),
    level = ifelse(grepl("upper", name), "upper", level),
    name = gsub("\\.lower|\\.upper", "", name)
  ) %>%
  pivot_wider(
    id_cols = c(species, Species, Predictor, .width, .point, .interval, name),
    names_from = level,
    values_from = value
  ) %>%
  ggplot(aes(y = Species, x = mid, xmin = lower, xmax = upper, col = name), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap( ~ Predictor, ncol = 4, scales = "free_x") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")

# save to file
ggsave(
  model2a_beta,
  filename = "outputs/figures/Fig4.png",
  device = png,
  width = 9,
  height = 12,
  units = "in",
  dpi = 600
)

# extract effects from model 2: how does A. donax affect species richness?
origin_list <- c("Exotic", "Native", "Translocated")
model2b_effects <- draws_model2b %>% 
  spread_draws(
    theta[origin],
    beta[origin, predictor],
    delta[origin, predictor]
  ) %>% 
  median_qi(
    theta, beta, delta = beta + delta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Origin = origin_list[origin],
    predictor = rownames(model2b_data$X)[predictor],
    Predictor = predictor_names[predictor]
  )
model2b_theta <- model2b_effects %>%
  select(contains("theta"), Origin, .width, .point, .interval) %>%
  ggplot(aes(y = Origin, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model2b_beta <- model2b_effects %>%
  select(contains("beta"), contains("delta"), Origin, Predictor, .width, .point, .interval) %>%
  pivot_longer(cols = c(contains("beta"), contains("delta"))) %>%
  mutate(
    level = "mid",
    level = ifelse(grepl("lower", name), "lower", level),
    level = ifelse(grepl("upper", name), "upper", level),
    name = gsub("\\.lower|\\.upper", "", name)
  ) %>%
  pivot_wider(
    id_cols = c(origin, Origin, Predictor, .width, .point, .interval, name),
    names_from = level,
    values_from = value
  ) %>%
  ggplot(aes(y = Predictor, x = mid, xmin = lower, xmax = upper, col = name), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap( ~ Origin, ncol = 3, scales = "free_x") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")

# save plots to file
ggsave(
  model2b_beta,
  filename = "outputs/figures/Fig5.png",
  device = png,
  width = 8,
  height = 6,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 3: SMI model
model2c_effects <- draws_model2c %>% 
  spread_draws(
    theta_main,
    theta[species],
    sigma_beta[predictor],
    zbeta_main[predictor],
    zbeta[species, predictor],
    sigma_delta[predictor],
    zdelta_main[predictor],
    zdelta[species, predictor]
  ) %>% 
  median_qi(
    theta = theta_main + theta,
    beta = model2c_data$sigma_fixed * zbeta_main +
      model2c_data$sigma_random * sigma_beta * zbeta,
    delta = beta + model2c_data$sigma_fixed * zdelta_main +
      model2c_data$sigma_random * sigma_delta * zdelta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = levels(factor(smi_full$SPECIES))[species],
    predictor = colnames(model2c_data$X)[predictor],
    Predictor = predictor_names[predictor]
  )
model2c_theta <- model2c_effects %>%
  select(contains("theta"), theta.lower, theta.upper, Species, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Parameter estimate")
model2c_beta <- model2c_effects %>%
  select(contains("beta"), contains("delta"), Species, Predictor, .width, .point, .interval) %>%
  pivot_longer(cols = c(contains("beta"), contains("delta"))) %>%
  mutate(
    level = "mid",
    level = ifelse(grepl("lower", name), "lower", level),
    level = ifelse(grepl("upper", name), "upper", level),
    name = gsub("\\.lower|\\.upper", "", name)
  ) %>%
  pivot_wider(
    id_cols = c(species, Species, Predictor, .width, .point, .interval, name),
    names_from = level,
    values_from = value
  ) %>%
  ggplot(aes(y = Species, x = mid, xmin = lower, xmax = upper, col = name), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  facet_wrap( ~ Predictor, ncol = 3) +
  xlab("Parameter estimate")

# combine theta plots into a single figure
theta_plot <- (model2a_theta + ggtitle("Model 2a")) | 
  ((model2b_theta + ggtitle("Model 2b")) / (model2c_theta + ggtitle("Model 2c")))

# save plots to file
ggsave(
  theta_plot,
  filename = "outputs/figures/Fig3.png",
  device = png,
  width = 5,
  height = 6,
  units = "in", 
  dpi = 600
)
ggsave(
  model2c_beta,
  filename = "outputs/figures/Fig6.png",
  device = png,
  width = 8,
  height = 8,
  units = "in", 
  dpi = 600
)

# extract fitted values
model1_fitted <- draws_model1 %>% 
  gather_draws(mu[obs]) %>% 
  median_qi() %>%
  mutate(
    fitted = plogis(.value),
    lower = plogis(.lower),
    upper = plogis(.upper),
    observed = model1_data$y
  )
model2a_fitted <- draws_model2a %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model2b_fitted <- draws_model2b %>%
  gather_draws(mu[origin, obs]) %>%
  median_qi()
model2c_fitted <- draws_model2c %>%
  gather_draws(mu[obs]) %>%
  median_qi()

# calculate fit stats for model 0a and 0b
model1_pred <- prediction(model1_fitted$fitted, model1_fitted$observed)
model1_auc <- performance(model1_pred, "auc")@y.values[[1]]
model1_curve <- performance(model1_pred, "tpr", "fpr")

# collate fitted values
fitted_values <- data.frame(
  model = rep(
    model_names,
    times = c(
      length(model1_curve@y.values[[1]]),
      nrow(model2a_fitted),
      nrow(model2b_fitted),
      nrow(model2c_fitted)
    )
  ),
  fitted = c(
    model1_curve@y.values[[1]],
    exp(model2a_fitted$.value - log(model2a_data$scale_factor)),
    exp(model2b_fitted$.value),
    exp(model2c_fitted$.value)
  ),
  lower = c(
    rep(NA, length(model1_curve@y.values[[1]])),
    exp(model2a_fitted$.lower - log(model2a_data$scale_factor)),
    exp(model2b_fitted$.lower),
    exp(model2c_fitted$.lower)
  ),
  upper = c(
    rep(NA, length(model1_curve@y.values[[1]])),
    exp(model2a_fitted$.upper - log(model2a_data$scale_factor)),
    exp(model2b_fitted$.upper),
    exp(model2c_fitted$.upper)
  ),
  observed = c(
    model1_curve@x.values[[1]],
    unlist(paired_cpue),
    unlist(paired_richness),
    model2c_data$y
  )
)

# plot fitted values or TPR/FPR plots
fitted_plot <- fitted_values %>%
  mutate(
    model = gsub("draws_", "", model),
    model = gsub("model", "Model ", model),
    model = gsub("Model 2", "Model 2b", model),
    model = gsub("Model 3", "Model 2c", model),
    model = gsub("Model 1a", "Model 2a", model),
    model = gsub("Model 0", "Model 1", model),
    upper = ifelse(upper > 10 * fitted, 10 * fitted, upper)
  ) %>%
  ggplot() +
  geom_point(aes(y = fitted, x = observed)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = observed)) +
  facet_wrap( ~ model, scales = "free") +
  xlab("False positive rate (Model 1) or observed value (Models 2a, 2b, or 2c)") +
  ylab("True positive rate (Model 1) or fitted value (Models 2a, 2b, or 2c)") +
  theme(axis.title = element_text(size = 8))

# save to file
ggsave(
  fitted_plot,
  filename = "outputs/figures/FigS2.png",
  device = png,
  width = 6,
  height = 4,
  units = "in", 
  dpi = 600
)

# grab fit stats
fit_stats <- c(
  model1_auc,
  cor(exp(model2a_fitted$.value - log(model2a_data$scale_factor)), c(unlist(paired_cpue))),
  cor(exp(model2b_fitted$.value), unlist(paired_richness)),
  cor(exp(model2c_fitted$.value), model2c_data$y)
)
fit_stats <- data.frame(
  model = model_names,
  statistic = c("AUC", rep("r", 3)),
  values = fit_stats
)
fit_stats <- fit_stats %>%
  mutate(model = gsub("draws_m", "M", model))
write.csv(fit_stats, file = "outputs/tables/Table2.csv")

# calculate posterior predictive distributions and plot against observed data
model1_pp <- pp_check(draws_model1, model1_data$y, breaks = seq(-0.5, 1.5, by = 1), xlim = c(-0.5, 1.5)) + 
  ggtitle("Model 1")
model2a_pp <- pp_check(draws_model2a, model2a_data$yflat, breaks = seq(-0.5, 101, by = 1), xlim = c(-0.5, 100)) + 
  ggtitle("Model 2a")
model2b_pp <- pp_check(draws_model2b, model2b_data$yflat, breaks = seq(-0.5, 225, by = 1), xlim = c(-0.5, 220)) + 
  ggtitle("Model 2b")
model2c_pp <- pp_check(draws_model2c, model2c_data$y, breaks = seq(-0.5, 285, by = 1), xlim = c(-0.5, 280)) + 
  ggtitle("Model 2c")

# save pp check plot to file
ggsave(
  (model1_pp / model2b_pp) | (model2a_pp / model2c_pp),
  filename = "outputs/figures/FigS3.png",
  device = png,
  width = 7,
  height = 5,
  units = "in", 
  dpi = 600
)
