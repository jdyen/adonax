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
library(tidybayes)
library(ROCR)

# set a flag to re-run models
refit_models <- FALSE

# load some helper functions
source("R/data-helpers.R")

# load species info
spp_info <- read_xlsx(
  "data/BD_Fish_Arundo.xlsx", sheet = "Metadata and other info",
  range = "A1:U30"
)

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
model0a_data <- list(
  N = nrow(cpue),
  K = ncol(model0_predictors),
  y = cpue %>% pull(Arundo),
  X = scale(model0_predictors),
  basin = cpue %>% pull(Basin) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_random = 2.
)
model0a_data$nbasin <- max(model0a_data$basin)

# create same again but with fish as predictors
model0_fish <- cpue %>%
  select(
    ACAR, ALAL, AMME, ANAN, BAQU, BAGR, BAHA,
    BAME, CAAU, COCA, COPA, CYCA, ESLU, GAHO, GOLO,
    LEGI, MISA, ONMY, PAMY, PHBI, PSPA, RURU, SAFL,
    SATR, SALU, SCER, SIGL, SQLA, GAGY
  ) 
model0b_data <- list(
  N = nrow(cpue),
  K = ncol(model0_fish),
  y = cpue %>% pull(Arundo),
  X = scale(model0_fish),
  basin = cpue %>% pull(Basin) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_random = 2.
)
model0b_data$nbasin <- max(model0b_data$basin)

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

# and predictor variables
paired_predictors <- paired_fish %>%
  select(
    Stream_order, Elevation, Water_temperature,
    Water_depth, Water_velocity, Conductivity, pH, 
    ammonia_no2, nitrate_phosphate, Channel_conservation,
    Habitat_diversity
  )

# and grab species info
paired_origin <- spp_info %>% 
  pivot_longer(
    cols = c(
      Besos, Daro, Ebro, Fluvia, Francoli, Foix, Gaia,
      Garona, Llobregat, Muga, Riudecanyes, Ridaura,
      Senia, Ter, Tordera
    ),
    names_to = "Basin",
    values_to = "Origin"
  ) %>% 
  mutate(Origin = rebase_factor(Origin)) %>%
  select(SPECIES, Basin, Origin) %>%
  pivot_wider(
    id_cols = Basin,
    values_from = Origin,
    names_from = SPECIES
  ) %>%
  right_join(paired_fish %>% select(Basin), by = "Basin") %>%
  select(-Basin) %>%
  select(all_of(colnames(paired_cpue)))

# dump this all together in a list
model1_data <- list(
  N = nrow(paired_cpue),
  Q = ncol(paired_cpue),
  K = ncol(paired_predictors),
  nflat = length(unlist(paired_cpue)),
  X = t(scale(paired_predictors)),
  arundo = paired_fish %>% pull(Arundo),
  ftle = c(scale(spp_info$`TROPHIC LEVEL`[match(colnames(paired_cpue), spp_info$SPECIES)])),
  norigin = max(paired_origin),
  origin = paired_origin,
  basin = paired_fish %>% pull(Basin) %>% rebase_factor(),
  block_id = paired_fish %>% pull(BLOCK) %>% rebase_factor(),
  site = paired_fish %>% pull(CODE1) %>% rebase_factor(),
  scale_factor = 10.,
  sigma_fixed = 5.,
  sigma_covar = 2.,
  sigma_random = 2.
)

# add random effect counts
model1_data$nbasin <- max(model1_data$basin)
model1_data$nblock <- max(model1_data$block_id)
model1_data$nsite <- max(model1_data$site)

# and add response variable and zero identifiers
model1_data$yflat <- round(c(t(paired_cpue)) * model1_data$scale_factor)
model1_data$nflat <- length(model1_data$yflat)
model1_data$zero_idx <- which(model1_data$yflat == 0)
model1_data$nzero <- length(model1_data$zero_idx)
model1_data$nonzero_idx <- which(model1_data$yflat > 0)
model1_data$notzero <- length(model1_data$nonzero_idx)

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
  nflat = length(unlist(paired_richness)),
  X = t(scale(paired_predictors)),
  arundo = paired_fish %>% pull(Arundo),
  wtl = c(scale(wtl)),
  basin = paired_fish %>% pull(Basin) %>% rebase_factor(),
  block_id = paired_fish %>% pull(BLOCK) %>% rebase_factor(),
  site = paired_fish %>% pull(CODE1) %>% rebase_factor(),
  sigma_fixed = 5.,
  sigma_covar = 2.,
  sigma_random = 2.
)

# add random effect counts
model2_data$nbasin <- max(model2_data$basin)
model2_data$nblock <- max(model2_data$block_id)
model2_data$nsite <- max(model2_data$site)

# and add response variable and zero identifiers
model2_data$yflat <- c(t(paired_richness))
model2_data$nflat <- length(model2_data$yflat)
model2_data$zero_idx <- which(model2_data$yflat == 0)
model2_data$nzero <- length(model2_data$zero_idx)
model2_data$nonzero_idx <- which(model2_data$yflat > 0)
model2_data$notzero <- length(model2_data$nonzero_idx)

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
  ftle = c(scale(smi_full %>% pull(trophic_level))),
  origin = smi_full %>% pull(Origin) %>% rebase_factor(),
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
model3_data$norigin <- max(model3_data$origin)
model3_data$nbasin <- max(model3_data$basin)
model3_data$nsite <- max(model3_data$site)
model3_data$nblock <- max(model3_data$block_id)
model3_data$norder <- max(model3_data$order)
model3_data$nfamily <- max(model3_data$family)
model3_data$ngenus <- max(model3_data$genus)
model3_data$Q <- max(model3_data$species)

# settings for all MCMC models
seed <- 3524142
iter <- 5000
warmup <- floor(iter / 2)
thin <- 2
chains <- 4
cores <- 4

# fit all models if required
if (refit_models) {
  
  # Fit model 0.a: Bayesian logistic regression, where does A. donax occur?
  # compile model
  model0_stan <- stan_model("src/bern.stan")
  
  # sample from model
  draws_model0a <- sampling(
    model0_stan,
    data = model0a_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model0a, file = paste0("outputs/fitted/draws-model0a-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model0a)
  
  # Fit model 0.b: Bayesian logistic regression, which fish assemblages are associated with A. donax?
  # sample from model
  draws_model0b <- sampling(
    model0_stan,
    data = model0b_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model0b, file = paste0("outputs/fitted/draws-model0b-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model0b)
  
  # Fit model 1: multivariate model of fish CPUE, how does A. donax affect fish species?
  
  # compile model
  model1_stan <- stan_model("src/mvn.stan")
  
  # define initial values
  empirical_corr <- cor(paired_cpue)
  init <- lapply(
    seq_len(chains),
    function(x) list(
      L = t(chol(empirical_corr))
    )
  )
  
  # sample from model
  draws_model1 <- sampling(
    model1_stan,
    data = model1_data,
    init = init,
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
  
  # compile model
  model1a_stan <- stan_model("src/zip-simple.stan")
  
  # create data set from model1_data
  model1a_data <- model1_data
  model1a_data$sigma_resid <- model1a_data$sigma_covar

  # sample from model
  draws_model1a <- sampling(
    model1a_stan,
    data = model1a_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model1a, file = paste0("outputs/fitted/draws-model1a-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model1a)
  
  # Fit model 2: multivariate model of species richness by origin, how does A. donax affect fish species?
  # compile model
  model2_stan <- stan_model("src/zip.stan")
  
  # define initial values
  empirical_corr <- cor(paired_richness)
  init <- lapply(
    seq_len(chains),
    function(x) list(
      L = t(chol(empirical_corr))
    )
  )
  
  # sample from model
  draws_model2 <- sampling(
    model2_stan,
    data = model2_data,
    init = init,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model2, file = paste0("outputs/fitted/draws-model2-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model2)
  
  # create data set from model1_data
  model2a_data <- model2_data
  model2a_data$sigma_resid <- model2a_data$sigma_covar
  model2a_data$scale_factor <- 1
  
  # sample from model
  draws_model2a <- sampling(
    model1a_stan,  # uses model 1a (identical)
    data = model2a_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model2a, file = paste0("outputs/fitted/draws-model2a-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model2a)
  
  # Fit model 3: SMI model
  # compile model
  model3_stan <- stan_model("src/lnorm.stan")
  
  # sample from model
  draws_model3 <- sampling(
    model3_stan,
    data = model3_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model3, file = paste0("outputs/fitted/draws-model3-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model3)
  
  # compile model
  model3a_stan <- stan_model("src/lnorm-simple.stan")
  
  # sample from model
  draws_model3a <- sampling(
    model3a_stan,
    data = model3_data,  # uses same data as model3
    iter = iter,
    warmup = warmup,
    chains = chains,
    thin = thin,
    cores = cores,
    control = list(adapt_delta = 0.9, max_treedepth = 18),
    seed = seed
  )
  
  # save fitted model
  qsave(draws_model3a, file = paste0("outputs/fitted/draws-model3a-", Sys.Date(), ".qs"))
  
  # free up some space
  rm(draws_model3a)
  
} else {
  
  # list all file names
  file_names <- c(
    "draws-model0a", "draws-model0b",
    "draws-model1", "draws-model1a",
    "draws-model2", "draws-model2a",
    "draws-model3", "draws-model3a"
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
  
}

# plot diagnostics
model_names <- c(
  "draws_model0a", "draws_model0b",
  "draws_model1", "draws_model1a",
  "draws_model2", "draws_model2a", 
  "draws_model3", "draws_model3a"
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
model0a_effects <- draws_model0a %>% 
  spread_draws(
    beta[Predictor],
  ) %>% 
  median_qi(
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Predictor = colnames(model0a_data$X)[Predictor],
    Predictor = predictor_names[Predictor]
  )
model0a_beta <- model0a_effects %>%
  ggplot(aes(y = Predictor, x = beta, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  xlab("Parameter estimate")
ggsave(
  model0a_beta,
  filename = "outputs/figures/model0a_beta.png",
  device = png,
  width = 5,
  height = 5,
  units = "in", 
  dpi = 600
)

# extract parameters from model 0.b: Bayesian logistic regression,
#     which fish assemblages are associated with A. donax?
model0b_effects <- draws_model0b %>% 
  spread_draws(beta[Species]) %>% 
  median_qi(.width = c(0.95, 0.66)) %>%
  mutate(Species = colnames(model0b_data$X)[Species])
model0b_beta <- model0b_effects %>%
  ggplot(aes(y = Species, x = beta, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  xlab("Parameter estimate")
ggsave(
  model0b_beta,
  filename = "outputs/figures/model0b_beta.png",
  device = png,
  width = 5,
  height = 8,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 1: multivariate model of fish CPUE, how does
#    A. donax affect fish species?
origin_list <- c("Exotic", "Native", "Translocated")
model1_effects <- draws_model1 %>% 
  spread_draws(
    beta[species, predictor],
    beta_arundo[species, predictor],
    theta[species],
    gamma,
    gamma_arundo,
    delta[origin, species]
  ) %>% 
  median_qi(
    beta,
    beta_arundo = beta + beta_arundo,
    gamma,
    gamma_arundo = gamma + gamma_arundo,
    theta,
    delta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = colnames(paired_cpue)[species],
    predictor = rownames(model1_data$X)[predictor],
    Predictor = predictor_names[predictor],
    Origin = origin_list[origin]
  )
model1_beta <- model1_effects %>%
  select(contains("beta"), Species, Predictor, .width, .point, .interval) %>%
  pivot_longer(cols = c(beta, beta_arundo), names_to = "Parameter") %>%
  mutate(
    .lower = ifelse(Parameter == "beta", beta.lower, beta_arundo.lower),
    .upper = ifelse(Parameter == "beta", beta.upper, beta_arundo.upper)
  ) %>%
  select(!contains("beta")) %>%
  ggplot(aes(y = Predictor, x = value, xmin = .lower, xmax = .upper, col = Parameter), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  facet_wrap( ~ Species, ncol = 5, scales = "free_x") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")
model1_delta <- model1_effects %>%
  select(contains("delta"), Species, Origin, .width, .point, .interval) %>%
  ungroup() %>%
  select(-predictor) %>%
  distinct() %>%
  ggplot(aes(y = Origin, x = delta, xmin = delta.lower, xmax = delta.upper)) +
  geom_pointinterval() +
  facet_wrap( ~ Species) +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model1_gamma <- model1_effects %>%
  select(contains("gamma"), Predictor, .width, .point, .interval) %>%
  pivot_longer(cols = c(gamma, gamma_arundo), names_to = "Parameter") %>%
  mutate(
    .lower = ifelse(Parameter == "gamma", gamma.lower, gamma_arundo.lower),
    .upper = ifelse(Parameter == "gamma", gamma.upper, gamma_arundo.upper)
  ) %>%
  ungroup() %>%
  select(!contains("gamma")) %>% 
  select(-species, -predictor) %>%
  distinct() %>%
  mutate(Parameter = ifelse(Parameter == "gamma_arundo", "A. donax present", "A. donax absent")) %>%
  ggplot(aes(y = Parameter, x = value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model1_theta <- model1_effects %>%
  select(contains("theta"), Species, .width, .point, .interval) %>%
  ungroup() %>%
  select(-predictor) %>%
  distinct() %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")

# save plots to file
ggsave(
  model1_beta,
  filename = "outputs/figures/model1_beta.png",
  device = png,
  width = 12,
  height = 12,
  units = "in", 
  dpi = 600
)
ggsave(
  model1_delta,
  filename = "outputs/figures/model1_delta.png",
  device = png,
  width = 7,
  height = 7,
  units = "in", 
  dpi = 600
)
ggsave(
  model1_gamma,
  filename = "outputs/figures/model1_gamma.png",
  device = png,
  width = 4,
  height = 4,
  units = "in", 
  dpi = 600
)
ggsave(
  model1_theta,
  filename = "outputs/figures/model1_theta.png",
  device = png,
  width = 5,
  height = 7,
  units = "in", 
  dpi = 600
)

model1a_effects <- draws_model1a %>% 
  spread_draws(
    theta[species],
    beta[species, predictor]
  ) %>% 
  median_qi(
    theta, beta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = colnames(paired_cpue)[species],
    predictor = rownames(model1_data$X)[predictor],
    Predictor = predictor_names[predictor]
  )
model1a_theta <- model1a_effects %>%
  select(contains("theta"), Species, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model1a_beta <- model1a_effects %>%
  select(contains("beta"), Species, Predictor, .width, .point, .interval) %>%
  ggplot(aes(y = Predictor, x = beta, xmin = beta.lower, xmax = beta.upper), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  facet_wrap( ~ Species, ncol = 4, scales = "free_x") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")

# save plots to file
ggsave(
  model1a_theta,
  filename = "outputs/figures/model1a_theta.png",
  device = png,
  width = 5,
  height = 7,
  units = "in", 
  dpi = 600
)
ggsave(
  model1a_beta,
  filename = "outputs/figures/model1a_beta.png",
  device = png,
  width = 10,
  height = 12,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 2: multivariate model of species
#   richness by origin, how does A. donax affect fish species?
model2_effects <- draws_model2 %>% 
  spread_draws(
    beta[origin, predictor],
    beta_arundo[origin, predictor],
    theta[origin],
    gamma[origin],
    gamma_arundo[origin]
  ) %>% 
  median_qi(
    beta,
    beta_arundo = beta + beta_arundo,
    gamma,
    gamma_arundo = gamma + gamma_arundo,
    theta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    predictor = rownames(model2_data$X)[predictor],
    Predictor = predictor_names[predictor],
    Origin = origin_list[origin]
  )
model2_beta <- model2_effects %>%
  select(contains("beta"), Predictor, Origin, .width, .point, .interval) %>%
  pivot_longer(cols = c(beta, beta_arundo), names_to = "Parameter") %>%
  mutate(
    .lower = ifelse(Parameter == "beta", beta.lower, beta_arundo.lower),
    .upper = ifelse(Parameter == "beta", beta.upper, beta_arundo.upper)
  ) %>%
  select(!contains("beta")) %>%
  ggplot(aes(y = Predictor, x = value, xmin = .lower, xmax = .upper, col = Parameter), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  facet_wrap( ~ Origin, scales = "free_x") +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")
model2_gamma <- model2_effects %>%
  select(contains("gamma"), Origin, .width, .point, .interval) %>%
  pivot_longer(cols = c(gamma, gamma_arundo), names_to = "Parameter") %>%
  mutate(
    .lower = ifelse(Parameter == "gamma", gamma.lower, gamma_arundo.lower),
    .upper = ifelse(Parameter == "gamma", gamma.upper, gamma_arundo.upper)
  ) %>%
  ungroup() %>%
  select(!contains("gamma")) %>% 
  distinct() %>%
  mutate(Parameter = ifelse(Parameter == "gamma_arundo", "A. donax present", "A. donax absent")) %>%
  ggplot(aes(y = Origin, x = value, xmin = .lower, xmax = .upper, col = Parameter), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  scale_color_brewer(type = "qual", palette = "Set2") +
  xlab("Parameter estimate")
model2_theta <- model2_effects %>%
  select(contains("theta"), Origin, .width, .point, .interval) %>%
  ungroup() %>%
  distinct() %>%
  ggplot(aes(y = Origin, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")

# save plots to file
ggsave(
  model2_beta,
  filename = "outputs/figures/model2_beta.png",
  device = png,
  width = 6,
  height = 6,
  units = "in", 
  dpi = 600
)
ggsave(
  model2_gamma,
  filename = "outputs/figures/model2_gamma.png",
  device = png,
  width = 5,
  height = 5,
  units = "in", 
  dpi = 600
)
ggsave(
  model2_theta,
  filename = "outputs/figures/model2_theta.png",
  device = png,
  width = 4,
  height = 4,
  units = "in", 
  dpi = 600
)

model2a_effects <- draws_model2a %>% 
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
model2a_theta <- model2a_effects %>%
  select(contains("theta"), Origin, .width, .point, .interval) %>%
  ggplot(aes(y = Origin, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model2a_beta <- model2a_effects %>%
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
  model2a_theta,
  filename = "outputs/figures/model2a_theta.png",
  device = png,
  width = 4,
  height = 4,
  units = "in", 
  dpi = 600
)
ggsave(
  model2a_beta,
  filename = "outputs/figures/model2a_beta.png",
  device = png,
  width = 8,
  height = 6,
  units = "in", 
  dpi = 600
)

# extract and plot parameters from model 3: SMI model
smi_origin <- levels(factor(smi_full$Origin))
model3_effects <- draws_model3 %>% 
  spread_draws(
    beta_main[predictor],
    beta[species, predictor],
    beta_main_arundo[predictor],
    beta_arundo[species, predictor],
    theta_main,
    theta[species],
    gamma,
    gamma_arundo,
    delta_main[origin],
    delta[species, origin]
  ) %>% 
  median_qi(
    beta_full = beta_main + beta,
    beta_arundo = beta_main + beta + beta_main_arundo + beta_arundo,
    gamma,
    gamma_arundo = gamma + gamma_arundo,
    theta = theta_main + theta,
    delta = delta_main + delta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(
    Species = levels(factor(smi_full$SPECIES))[species],
    predictor = colnames(model3_data$X)[predictor],
    Predictor = predictor_names[predictor],
    Origin = smi_origin[origin]
  ) %>% 
  rename(
    beta = beta_full, 
    beta.lower = beta_full.lower,
    beta.upper = beta_full.upper
  )
model3_beta <- model3_effects %>%
  select(contains("beta"), Species, Predictor, .width, .point, .interval) %>%
  pivot_longer(cols = c(beta, beta_arundo), names_to = "Parameter") %>%
  mutate(
    .lower = ifelse(Parameter == "beta", beta.lower, beta_arundo.lower),
    .upper = ifelse(Parameter == "beta", beta.upper, beta_arundo.upper)
  ) %>%
  select(!contains("beta")) %>%
  ggplot(aes(y = Predictor, x = value, xmin = .lower, xmax = .upper, col = Parameter), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  facet_wrap( ~ Species, ncol = 3) +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("A. donax absent", "A. donax present")) +
  xlab("Parameter estimate") +
  theme(legend.position = "bottom")
model3_gamma <- model3_effects %>%
  select(contains("gamma"), .width, .point, .interval) %>%
  pivot_longer(cols = c(gamma, gamma_arundo), names_to = "Parameter") %>%
  mutate(
    .lower = ifelse(Parameter == "gamma", gamma.lower, gamma_arundo.lower),
    .upper = ifelse(Parameter == "gamma", gamma.upper, gamma_arundo.upper)
  ) %>%
  ungroup() %>%
  select(!contains("gamma")) %>%
  select(-predictor, -species) %>%
  distinct() %>%
  mutate(Parameter = ifelse(Parameter == "gamma_arundo", "A. donax present", "A. donax absent")) %>%
  ggplot(aes(y = Parameter, x = value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  xlab("Parameter estimate")
model3_theta <- model3_effects %>%
  select(contains("theta"), Species, .width, .point, .interval) %>%
  ungroup() %>%
  select(-predictor) %>%
  distinct() %>%
  ggplot(aes(y = Species, x = theta, xmin = theta.lower, xmax = theta.upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")
model3_delta <- model3_effects %>%
  select(contains("delta"), Species, Origin, .width, .point, .interval) %>%
  ungroup() %>%
  select(-predictor) %>%
  distinct() %>%
  ggplot(aes(y = Species, x = delta, xmin = delta.lower, xmax = delta.upper, col = Origin), position = position_dodge(0.4)) +
  geom_pointinterval(position = position_dodge(0.4)) +
  scale_color_brewer(type = "qual", palette = "Set2", labels = c("Native", "Translocated")) +
  xlab("Parameter estimate")

# save plots to file
ggsave(
  model3_beta,
  filename = "outputs/figures/model3_beta.png",
  device = png,
  width = 8,
  height = 8,
  units = "in", 
  dpi = 600
)
ggsave(
  model3_delta,
  filename = "outputs/figures/model3_delta.png",
  device = png,
  width = 5,
  height = 6,
  units = "in", 
  dpi = 600
)
ggsave(
  model3_gamma,
  filename = "outputs/figures/model3_gamma.png",
  device = png,
  width = 4,
  height = 4,
  units = "in", 
  dpi = 600
)
ggsave(
  model3_theta,
  filename = "outputs/figures/model3_theta.png",
  device = png,
  width = 5,
  height = 6,
  units = "in", 
  dpi = 600
)

model3a_effects <- draws_model3a %>% 
  spread_draws(theta[species]) %>% 
  median_qi(
    theta,
    .width = c(0.95, 0.66)
  ) %>%
  mutate(Species = levels(factor(smi_full$SPECIES))[species])
model3a_theta <- model3a_effects %>%
  select(contains("theta"), .lower, .upper, Species, .width, .point, .interval) %>%
  ggplot(aes(y = Species, x = theta, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme(legend.position = "none") +
  xlab("Parameter estimate")

# save plots to file
ggsave(
  model3a_theta,
  filename = "outputs/figures/model3a_theta.png",
  device = png,
  width = 5,
  height = 6,
  units = "in", 
  dpi = 600
)

# extract fitted values
model0a_fitted <- draws_model0a %>% 
  gather_draws(mu[obs]) %>% 
  median_qi() %>%
  mutate(
    fitted = plogis(.value),
    lower = plogis(.lower),
    upper = plogis(.upper),
    observed = model0a_data$y
  )
model0b_fitted <- draws_model0b %>% 
  gather_draws(mu[obs]) %>% 
  median_qi() %>%
  mutate(
    fitted = plogis(.value),
    lower = plogis(.lower),
    upper = plogis(.upper),
    observed = model0a_data$y
  )
model1_fitted <- draws_model1 %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model1a_fitted <- draws_model1a %>%
  gather_draws(mu[species, obs]) %>%
  median_qi()
model2_fitted <- draws_model2 %>%
  gather_draws(mu[origin, obs]) %>%
  median_qi()
model2a_fitted <- draws_model2a %>%
  gather_draws(mu[origin, obs]) %>%
  median_qi()
model3_fitted <- draws_model3 %>%
  gather_draws(mu[obs]) %>%
  median_qi()
model3a_fitted <- draws_model3a %>%
  gather_draws(mu[obs]) %>%
  median_qi()

# calculate fit stats for model 0a and 0b
model0a_pred <- prediction(model0a_fitted$fitted, model0a_fitted$observed)
model0a_auc <- performance(model0a_pred, "auc")@y.values[[1]]
model0a_curve <- performance(model0a_pred, "tpr", "fpr")
model0b_pred <- prediction(model0b_fitted$fitted, model0b_fitted$observed)
model0b_auc <- performance(model0b_pred, "auc")@y.values[[1]]
model0b_curve <- performance(model0b_pred, "tpr", "fpr")

# collate fitted values
fitted_values <- data.frame(
  model = rep(
    model_names,
    times = c(
      length(model0a_curve@y.values[[1]]),
      length(model0b_curve@y.values[[1]]),
      nrow(model1_fitted),
      nrow(model1a_fitted),
      nrow(model2_fitted),
      nrow(model2a_fitted),
      nrow(model3_fitted),
      nrow(model3a_fitted)
    )
  ),
  fitted = c(
    model0a_curve@y.values[[1]],
    model0b_curve@y.values[[1]],
    exp(model1_fitted$.value - log(model1_data$scale_factor)),
    exp(model1a_fitted$.value),
    exp(model2_fitted$.value),
    exp(model2a_fitted$.value),
    exp(model3_fitted$.value),
    exp(model3a_fitted$.value)
  ),
  lower = c(
    rep(NA, length(model0a_curve@y.values[[1]])),
    rep(NA, length(model0b_curve@y.values[[1]])),
    exp(model1_fitted$.lower - log(model1_data$scale_factor)),
    exp(model1a_fitted$.lower),
    exp(model2_fitted$.lower),
    exp(model2a_fitted$.lower),
    exp(model3_fitted$.lower),
    exp(model3a_fitted$.lower)
  ),
  upper = c(
    rep(NA, length(model0a_curve@y.values[[1]])),
    rep(NA, length(model0b_curve@y.values[[1]])),
    exp(model1_fitted$.upper - log(model1_data$scale_factor)),
    exp(model1a_fitted$.upper),
    exp(model2_fitted$.upper),
    exp(model2a_fitted$.upper),
    exp(model3_fitted$.upper),
    exp(model3a_fitted$.upper)
  ),
  observed = c(
    model0a_curve@x.values[[1]],
    model0b_curve@x.values[[1]],
    unlist(paired_cpue),
    unlist(paired_cpue),
    unlist(paired_richness),
    unlist(paired_richness),
    model3_data$y,
    model3_data$y
  )
)

# plot fitted values or TPR/FPR plots
fitted_plot <- fitted_values %>%
  mutate(
    model = gsub("draws_", "", model),
    model = gsub("model", "Model ", model),
    fitted = ifelse(fitted > 1000 & model == "Model 1", 1000, fitted),
    upper = ifelse(upper > 2 * fitted, 2 * fitted, upper)
  ) %>%
  ggplot() +
  geom_point(aes(y = fitted, x = observed)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = observed)) +
  facet_wrap( ~ model, scales = "free") +
  xlab("Observed value (Models 1-3) or false positive rate (Model 0a and 0b)") +
  ylab("Fitted value (Models 1-3) or true positive rate (Model 0a and 0b)")

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
model1_fitted_val <- exp(model1_fitted$.value - log(model1_data$scale_factor))
model1_fitted_val <- ifelse(model1_fitted_val > 1000, 1000, model1_fitted_val)
fit_stats <- c(
  model0a_auc,
  model0b_auc,
  cor(model1_fitted_val, c(unlist(paired_cpue))),
  cor(exp(model1a_fitted$.value), c(unlist(paired_cpue))),
  cor(exp(model2_fitted$.value), unlist(paired_richness)),
  cor(exp(model2a_fitted$.value), unlist(paired_richness)),
  cor(exp(model3_fitted$.value), model3_data$y),
  cor(exp(model3a_fitted$.value), model3_data$y)
)
fit_stats <- data.frame(
  model = model_names,
  statistic = c(rep("AUC", 2), rep("r", 6)),
  values = fit_stats
)
fit_stats <- fit_stats %>%
  mutate(model = gsub("draws_m", "M", model))
write.csv(fit_stats, file = "outputs/tables/fit-statistics.csv")

# extract and plot correlations for models 1 and 2
corr_model1 <- draws_model1 %>% 
  spread_draws(Omega[species_a, species_b]) %>% 
  median_qi(.width = c(0.95, 0.66))
corr_model2 <- draws_model2 %>% 
  spread_draws(Omega[origin_a, origin_b]) %>% 
  median_qi(.width = c(0.95, 0.66))
corr_plot1 <- corr_model1 %>% 
  mutate(Omega = ifelse(species_a == species_b, NA, Omega)) %>%
  ggplot(aes(x = species_a, y = species_b)) +
  geom_tile(aes(fill = Omega)) +
  scale_fill_viridis_c(limits = c(-0.2, 0.2))
corr_plot2 <- corr_model2 %>% 
  mutate(Omega = ifelse(origin_a == origin_b, NA, Omega)) %>%
  ggplot(aes(x = origin_a, y = origin_b)) +
  geom_tile(aes(fill = Omega)) +
  scale_fill_viridis_c(limits = c(-0.2, 0.2))

# save to file
ggsave(
  corr_plot1,
  filename = "outputs/figures/model1_corr.png",
  device = png,
  width = 6,
  height = 6,
  units = "in", 
  dpi = 600
)
ggsave(
  corr_plot2,
  filename = "outputs/figures/model2_corr.png",
  device = png,
  width = 4,
  height = 4,
  units = "in", 
  dpi = 600
)
