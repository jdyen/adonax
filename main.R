# Analysis of the effects of Arundo donax on fish assemblages in
#   north-eastern Spain
#
# Authors: Alberto Maceda Veiga, Ralph Mac Nally, Jian Yen
# 
# last updated: 26 August 2022 

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
refit_models <- TRUE
refit_model3 <- FALSE

# load some helper functions
source("R/data-helpers.R")
source("R/summary-helpers.R")

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
model0_predictors <- cpue %>%
  select(
    Stream_order, Elevation, Water_temperature,
    Water_depth, Water_velocity, Conductivity,
    pH, ammonia_no2, nitrate_phosphate,
    Channel_conservation, Habitat_diversity
  )
model0_data <- list(
  N = nrow(cpue),
  K = ncol(model0_predictors),
  y = cpue %>% pull(Arundo),
  X = scale(model0_predictors),
  basin = cpue %>% pull(Basin) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_random = 2.
)
model0_data$nbasin <- max(model0_data$basin)

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

# diet guild
cpue_diet <- cpue_guild %>%
  group_by(id, diet) %>%
  summarise(cpue = sum(value)) %>%
  pivot_wider(id_cols = id, values_from = cpue, names_from = diet) %>%
  ungroup() %>%
  select(-id)

# migration guild
cpue_migration <- cpue_guild %>%
  group_by(id, migration) %>%
  summarise(cpue = sum(value)) %>%
  pivot_wider(id_cols = id, values_from = cpue, names_from = migration) %>%
  ungroup() %>%
  select(-id)

# habitat use guild
cpue_habitat_use <- cpue_guild %>%
  group_by(id, habitat_use) %>%
  summarise(cpue = sum(value)) %>%
  pivot_wider(id_cols = id, values_from = cpue, names_from = habitat_use) %>%
  ungroup() %>%
  select(-id)

# and predictor variables
paired_predictors <- paired_fish %>%
  select(
    Stream_order, Elevation, Water_temperature,
    Water_depth, Water_velocity, Conductivity, pH, 
    ammonia_no2, nitrate_phosphate, Channel_conservation,
    Habitat_diversity
  )

# dump this all together in a list for all species and for each guild
model1a_data <- list(
  N = nrow(paired_cpue),
  Q = ncol(paired_cpue),
  K = ncol(paired_predictors),
  X = t(scale(paired_predictors)),
  arundo = paired_fish %>% pull(Arundo),
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
model1a_data$nbasin <- max(model1a_data$basin)
model1a_data$nblock <- max(model1a_data$block_id)
model1a_data$nsite <- max(model1a_data$site)

# and add response variable and zero identifiers
model1a_data$yflat <- round(c(t(paired_cpue)) * model1a_data$scale_factor)
model1a_data$nflat <- length(model1a_data$yflat)

# create similar data sets for the diet guild
model1b_data <- model1a_data
model1b_data$Q <- ncol(cpue_diet)
model1b_data$yflat <- round(c(t(cpue_diet)) * model1b_data$scale_factor)
model1b_data$nflat <- length(model1b_data$yflat)

# create similar data sets for the migration guild
model1c_data <- model1a_data
model1c_data$Q <- ncol(cpue_migration)
model1c_data$yflat <- round(c(t(cpue_migration)) * model1c_data$scale_factor)
model1c_data$nflat <- length(model1c_data$yflat)

# create similar data sets for the habitat use guild
model1d_data <- model1a_data
model1d_data$Q <- ncol(cpue_habitat_use)
model1d_data$yflat <- round(c(t(cpue_habitat_use)) * model1d_data$scale_factor)
model1d_data$nflat <- length(model1d_data$yflat)

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
model2_data <- list(
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
model2_data$nbasin <- max(model2_data$basin)
model2_data$nblock <- max(model2_data$block_id)
model2_data$nsite <- max(model2_data$site)

# and add response variable and zero identifiers
model2_data$yflat <- c(t(paired_richness))
model2_data$nflat <- length(model2_data$yflat)

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
model3_data <- list(
  N = nrow(smi_full),
  K = ncol(smi_predictors),
  y = smi_full %>% pull(SMI),
  X = scale(smi_predictors),
  arundo = smi_full %>% pull(Arundo),
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
model3_data$nbasin <- max(model3_data$basin)
model3_data$nsite <- max(model3_data$site)
model3_data$nblock <- max(model3_data$block_id)
model3_data$norder <- max(model3_data$order)
model3_data$nfamily <- max(model3_data$family)
model3_data$ngenus <- max(model3_data$genus)
model3_data$Q <- max(model3_data$species)

# settings for all MCMC models
seed <- 352124142
iter <- 1000
warmup <- floor(iter / 2)
thin <- 1 #4
chains <- 4
cores <- 4

# fit all models if required
if (refit_models) {
  
  # Fit model 0: Bayesian logistic regression, where does A. donax occur?
  # compile model
  model0_stan <- stan_model("src/bern.stan")
  
  # sample from model
  draws_model0 <- sampling(
    model0_stan,
    data = model0_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model0, file = paste0("outputs/fitted/draws-model0-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model0)

  # Fit model 1a: model of fish CPUE, how does A. donax affect fish species?
  # compile model
  model1_stan <- stan_model("src/nb-simple.stan")
  
  # sample from model
  draws_model1a <- sampling(
    model1_stan,
    data = model1a_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 20),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model1a, file = paste0("outputs/fitted/draws-model1a-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model1a)

  # sample from model for diet guilds
  draws_model1b <- sampling(
    model1_stan,
    data = model1b_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 20),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model1b, file = paste0("outputs/fitted/draws-model1b-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model1b)
  
  # sample from model for migration guilds
  draws_model1c <- sampling(
    model1_stan,
    data = model1c_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 20),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model1c, file = paste0("outputs/fitted/draws-model1c-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model1c)
  
  # sample from model for habitat use guilds
  draws_model1d <- sampling(
    model1_stan,
    data = model1d_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 20),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model1d, file = paste0("outputs/fitted/draws-model1d-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model1d)
  
  # sample from model for species richness by classification (native or otherwise)
  draws_model2 <- sampling(
    model1_stan,
    data = model2_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model2, file = paste0("outputs/fitted/draws-model2-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model2)
  
  # only refit model3 if required (very slow and model is stable)
  if (refit_model3) {
    
    # compile model
    model3_stan <- stan_model("src/lnorm-simple.stan")
    
    # sample from model
    draws_model3 <- sampling(
      model3_stan,
      data = model3_data,
      iter = iter,
      warmup = warmup,
      chains = chains,
      thin = thin,
      cores = cores,
      control = list(adapt_delta = 0.95, max_treedepth = 18),
      seed = seed
    )
    
    # save fitted model
    qsave(draws_model3, file = paste0("outputs/fitted/draws-model3-", Sys.Date(), ".qs"))
    
    # free up some space
    rm(draws_model3)
    
  }
  
}

# load all fitted models
file_names <- c(
  "draws-model0",
  "draws-model1a",
  "draws-model1b",
  "draws-model1c",
  "draws-model1d",
  "draws-model2",
  "draws-model3"
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
  "draws_model0",
  "draws_model1a",
  "draws_model1b",
  "draws_model1c",
  "draws_model1d",
  "draws_model2", 
  "draws_model3"
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
  pivot_longer(
    cols = c(Rhat, n_eff),
    names_to = "statistic"
  ) %>%
  mutate(
    model = gsub("draws_", "", model),
    model = gsub("model", "Model ", model),
    statistic = gsub("n_eff", "Effective sample size", statistic)
  ) %>%
  ggplot() +
  geom_histogram(aes(x = value), bins = 50) +
  facet_grid(model ~ statistic, scales = "free")
ggsave(
  diagnostic_plot,
  filename = "outputs/figures/diagnostics.png",
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
model0_effects <- draws_model0 %>% 
  spread_draws(
    beta[Predictor],
  ) %>% 
  median_qi(
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Predictor = colnames(model0_data$X)[Predictor],
    Predictor = predictor_names[Predictor]
  )
model0_beta <- model0_effects %>%
  ggplot(aes(y = Predictor, x = beta, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Parameter estimate")
ggsave(
  model0_beta,
  filename = "outputs/figures/model0_beta.png",
  device = png,
  width = 5,
  height = 5,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 1a: model of fish CPUE, how does
#    A. donax affect fish species?
model1a_effects <- draws_model1a %>% 
  spread_draws(
    theta[species]#,
    # beta[species, predictor]
  ) %>% 
  median_qi(
    theta, #beta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = colnames(paired_cpue)[species]#,
    # predictor = rownames(model1a_data$X)[predictor],
    # Predictor = predictor_names[predictor]
  )
model1a_theta <- model1a_effects %>%
  select(contains("theta"), Species, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
# model1a_beta <- model1a_effects %>%
#   select(contains("beta"), Species, Predictor, .width, .point, .interval) %>%
#   ggplot(aes(y = Species, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
#   geom_pointinterval(position = position_dodge(0.4)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   facet_wrap( ~ Predictor, ncol = 4, scales = "free_x") +
#   scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
#   xlab("Parameter estimate") +
#   theme(legend.position = "bottom")

# extract and plot parameters from model 1b: model of CPUE of fish diet guilds, how does
#    A. donax affect fish species?
model1b_effects <- draws_model1b %>% 
  spread_draws(
    theta[guild] #,
    # beta[guild, predictor]
  ) %>% 
  median_qi(
    theta, #beta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Guild = colnames(cpue_diet)[guild]#,
    # predictor = rownames(model1b_data$X)[predictor],
    # Predictor = predictor_names[predictor]
  )
model1b_theta <- model1b_effects %>%
  select(contains("theta"), Guild, .width, .point, .interval) %>%
  ggplot(aes(y = Guild, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
# model1b_beta <- model1b_effects %>%
#   select(contains("beta"), Guild, Predictor, .width, .point, .interval) %>%
#   ggplot(aes(y = Guild, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
#   geom_pointinterval(position = position_dodge(0.4)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   facet_wrap( ~ Predictor, ncol = 1, scales = "free_x") +
#   scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
#   xlab("Parameter estimate") +
#   theme(legend.position = "bottom")

# extract and plot parameters from model 1c: model of CPUE of fish migration guilds, how does
#    A. donax affect fish species?
model1c_effects <- draws_model1c %>% 
  spread_draws(
    theta[guild] #,
    # beta[guild, predictor]
  ) %>% 
  median_qi(
    theta, #beta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Guild = colnames(cpue_migration)[guild]#,
    # predictor = rownames(model1b_data$X)[predictor],
    # Predictor = predictor_names[predictor]
  )
model1c_theta <- model1c_effects %>%
  select(contains("theta"), Guild, .width, .point, .interval) %>%
  ggplot(aes(y = Guild, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
# model1c_beta <- model1c_effects %>%
#   select(contains("beta"), Guild, Predictor, .width, .point, .interval) %>%
#   ggplot(aes(y = Guild, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
#   geom_pointinterval(position = position_dodge(0.4)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   facet_wrap( ~ Predictor, ncol = 1, scales = "free_x") +
#   scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
#   xlab("Parameter estimate") +
#   theme(legend.position = "bottom")

# extract and plot parameters from model 1d: model of CPUE of fish habitat-use guilds, how does
#    A. donax affect fish species?
model1d_effects <- draws_model1d %>% 
  spread_draws(
    theta[guild]#,
    # beta[guild, predictor]
  ) %>% 
  median_qi(
    theta, #beta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Guild = colnames(cpue_habitat_use)[guild]#,
    # predictor = rownames(model1d_data$X)[predictor],
    # Predictor = predictor_names[predictor]
  )
model1d_theta <- model1d_effects %>%
  select(contains("theta"), Guild, .width, .point, .interval) %>%
  ggplot(aes(y = Guild, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
# model1d_beta <- model1d_effects %>%
#   select(contains("beta"), Guild, Predictor, .width, .point, .interval) %>%
#   ggplot(aes(y = Guild, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
#   geom_pointinterval(position = position_dodge(0.4)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   facet_wrap( ~ Predictor, ncol = 1, scales = "free_x") +
#   scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
#   xlab("Parameter estimate") +
#   theme(legend.position = "bottom")

# combine theta plots and save to file
model1_theta <- model1a_theta | (model1b_theta / model1c_theta / model1d_theta)
ggsave(
  model1_theta,
  filename = "outputs/figures/model1_theta.png",
  device = png,
  width = 5,
  height = 7,
  units = "in", 
  dpi = 600
)

# combine beta plots and save to file
# model1_beta <- model1b_beta | model1c_beta | model1d_beta
# ggsave(
#   model1_beta,
#   filename = "outputs/figures/model1_beta.png",
#   device = png,
#   width = 10,
#   height = 14,
#   units = "in", 
#   dpi = 600
# )
# ggsave(
#   model1a_beta,
#   filename = "outputs/figures/model1a_beta.png",
#   device = png,
#   width = 7,
#   height = 10,
#   units = "in", 
#   dpi = 600
# )

# extract effects from model 2: how does A. donax affect species richness?
origin_list <- c("Exotic", "Native", "Translocated")
model2_effects <- draws_model2 %>% 
  spread_draws(
    theta[origin],
    beta[origin, predictor]
  ) %>% 
  median_qi(
    theta, beta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Origin = origin_list[origin],
    predictor = rownames(model2_data$X)[predictor],
    Predictor = predictor_names[predictor]
  )
model2_theta <- model2_effects %>%
  select(contains("theta"), Origin, .width, .point, .interval) %>%
  ggplot(aes(y = Origin, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model2_beta <- model2_effects %>%
  select(contains("beta"), Origin, Predictor, .width, .point, .interval) %>%
  ggplot(aes(y = Predictor, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap( ~ Origin, ncol = 3, scales = "free_x") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")

# save plots to file
ggsave(
  model2_theta,
  filename = "outputs/figures/model2_theta.png",
  device = png,
  width = 4,
  height = 4,
  units = "in", 
  dpi = 600
)
ggsave(
  model2_beta,
  filename = "outputs/figures/model2_beta.png",
  device = png,
  width = 8,
  height = 6,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 3: SMI model
model3_effects <- draws_model3 %>% 
  spread_draws(
    theta_main,
    theta[species],
    sigma_beta[predictor],
    zbeta_main[predictor],
    zbeta[species, predictor]
  ) %>% 
  median_qi(
    theta = theta_main + theta,
    beta = model3_data$sigma_fixed * zbeta_main +
      model3_data$sigma_random * sigma_beta * zbeta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = levels(factor(smi_full$SPECIES))[species],
    predictor = colnames(model3_data$X)[predictor],
    Predictor = predictor_names[predictor]
  )
model3_theta <- model3_effects %>%
  select(contains("theta"), theta.lower, theta.upper, Species, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Parameter estimate")
model3_beta <- model3_effects %>%
  select(contains("beta"), Species, Predictor, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap( ~ Predictor, ncol = 3) +
  xlab("Parameter estimate")

# save plots to file
ggsave(
  model3_theta,
  filename = "outputs/figures/model3_theta.png",
  device = png,
  width = 5,
  height = 6,
  units = "in", 
  dpi = 600
)
ggsave(
  model3_beta,
  filename = "outputs/figures/model3_beta.png",
  device = png,
  width = 8,
  height = 8,
  units = "in", 
  dpi = 600
)

# extract fitted values
model0_fitted <- draws_model0 %>% 
  gather_draws(mu[obs]) %>% 
  median_qi() %>%
  mutate(
    fitted = plogis(.value),
    lower = plogis(.lower),
    upper = plogis(.upper),
    observed = model0_data$y
  )
model1a_fitted <- draws_model1a %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model1b_fitted <- draws_model1b %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model1c_fitted <- draws_model1c %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model1d_fitted <- draws_model1d %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model2_fitted <- draws_model2 %>%
  gather_draws(mu[origin, obs]) %>%
  median_qi()
model3_fitted <- draws_model3 %>%
  gather_draws(mu[obs]) %>%
  median_qi()

# calculate fit stats for model 0a and 0b
model0_pred <- prediction(model0_fitted$fitted, model0_fitted$observed)
model0_auc <- performance(model0_pred, "auc")@y.values[[1]]
model0_curve <- performance(model0_pred, "tpr", "fpr")

# collate fitted values
fitted_values <- data.frame(
  model = rep(
    model_names,
    times = c(
      length(model0_curve@y.values[[1]]),
      nrow(model1a_fitted),
      nrow(model1b_fitted),
      nrow(model1c_fitted),
      nrow(model1d_fitted),
      nrow(model2_fitted),
      nrow(model3_fitted)
    )
  ),
  fitted = c(
    model0_curve@y.values[[1]],
    exp(model1a_fitted$.value - log(model1a_data$scale_factor)),
    exp(model1b_fitted$.value - log(model1b_data$scale_factor)),
    exp(model1c_fitted$.value - log(model1c_data$scale_factor)),
    exp(model1d_fitted$.value - log(model1d_data$scale_factor)),
    exp(model2_fitted$.value),
    exp(model3_fitted$.value)
  ),
  lower = c(
    rep(NA, length(model0_curve@y.values[[1]])),
    exp(model1a_fitted$.lower - log(model1a_data$scale_factor)),
    exp(model1b_fitted$.lower - log(model1b_data$scale_factor)),
    exp(model1c_fitted$.lower - log(model1c_data$scale_factor)),
    exp(model1d_fitted$.lower - log(model1d_data$scale_factor)),
    exp(model2_fitted$.lower),
    exp(model3_fitted$.lower)
  ),
  upper = c(
    rep(NA, length(model0_curve@y.values[[1]])),
    exp(model1a_fitted$.upper - log(model1a_data$scale_factor)),
    exp(model1b_fitted$.upper - log(model1b_data$scale_factor)),
    exp(model1c_fitted$.upper - log(model1c_data$scale_factor)),
    exp(model1d_fitted$.upper - log(model1d_data$scale_factor)),
    exp(model2_fitted$.upper),
    exp(model3_fitted$.upper)
  ),
  observed = c(
    model0_curve@x.values[[1]],
    unlist(paired_cpue),
    unlist(cpue_diet),
    unlist(cpue_migration),
    unlist(cpue_habitat_use),
    unlist(paired_richness),
    model3_data$y
  )
)

# plot fitted values or TPR/FPR plots
fitted_plot <- fitted_values %>%
  mutate(
    model = gsub("draws_", "", model),
    model = gsub("model", "Model ", model),
    upper = ifelse(upper > 10 * fitted, 10 * fitted, upper)
  ) %>%
  ggplot() +
  geom_point(aes(y = fitted, x = observed)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = observed)) +
  facet_wrap( ~ model, scales = "free") +
  xlab("Observed value (Models 1-3) or false positive rate (Model 0)") +
  ylab("Fitted value (Models 1-3) or true positive rate (Model 0)") +
  theme(axis.title = element_text(size = 8))

# save to file
ggsave(
  fitted_plot,
  filename = "outputs/figures/fitted_vs_observed.png",
  device = png,
  width = 6,
  height = 4,
  units = "in", 
  dpi = 600
)

# grab fit stats
fit_stats <- c(
  model0_auc,
  cor(exp(model1a_fitted$.value - log(model1a_data$scale_factor)), c(unlist(paired_cpue))),
  cor(exp(model1b_fitted$.value - log(model1b_data$scale_factor)), c(unlist(cpue_diet))),
  cor(exp(model1c_fitted$.value - log(model1c_data$scale_factor)), c(unlist(cpue_migration))),
  cor(exp(model1d_fitted$.value - log(model1d_data$scale_factor)), c(unlist(cpue_habitat_use))),
  cor(exp(model2_fitted$.value), unlist(paired_richness)),
  cor(exp(model3_fitted$.value), model3_data$y)
)
fit_stats <- data.frame(
  model = model_names,
  statistic = c("AUC", rep("r", 6)),
  values = fit_stats
)
fit_stats <- fit_stats %>%
  mutate(model = gsub("draws_m", "M", model))
write.csv(fit_stats, file = "outputs/tables/fit-statistics.csv")

# calculate posterior predictive distributions and plot against observed data
model0_pp <- pp_check(draws_model0, model0_data$y, breaks = seq(-0.5, 1.5, by = 1), xlim = c(-0.5, 1.5))
model1a_pp <- pp_check(draws_model1a, model1a_data$yflat, breaks = seq(-0.5, 101, by = 1), xlim = c(-0.5, 100))
model1b_pp <- pp_check(draws_model1b, model1b_data$yflat, breaks = seq(-0.5, 101, by = 1), xlim = c(-0.5, 100))
model1c_pp <- pp_check(draws_model1c, model1c_data$yflat, breaks = seq(-0.5, 101, by = 1), xlim = c(-0.5, 100))
model1d_pp <- pp_check(draws_model1d, model1d_data$yflat, breaks = seq(-0.5, 101, by = 1), xlim = c(-0.5, 100))
model2_pp <- pp_check(draws_model2, model2_data$yflat, breaks = seq(-0.5, 225, by = 1), xlim = c(-0.5, 220))
model3_pp <- pp_check(draws_model3, model3_data$y, breaks = seq(-0.5, 285, by = 1), xlim = c(-0.5, 280))

# save pp check plot to file
ggsave(
  (model1a_pp / model1b_pp) | (model1c_pp / model1d_pp),
  filename = "outputs/figures/model1_pp.png",
  device = png,
  width = 7,
  height = 5,
  units = "in", 
  dpi = 600
)
