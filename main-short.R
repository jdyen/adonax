

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
model2_fitted <- draws_model2 %>%
  gather_draws(mu[origin, obs]) %>%
  median_qi()
model3_fitted <- draws_model3 %>%
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
      nrow(model2_fitted),
      nrow(model3_fitted)
    )
  ),
  fitted = c(
    model0a_curve@y.values[[1]],
    model0b_curve@y.values[[1]],
    exp(model1_fitted$.value - log(model1_data$scale_factor)),
    exp(model2_fitted$.value),
    exp(model3_fitted$.value)
  ),
  lower = c(
    rep(NA, length(model0a_curve@y.values[[1]])),
    rep(NA, length(model0b_curve@y.values[[1]])),
    exp(model1_fitted$.lower - log(model1_data$scale_factor)),
    exp(model2_fitted$.lower),
    exp(model3_fitted$.lower)
  ),
  upper = c(
    rep(NA, length(model0a_curve@y.values[[1]])),
    rep(NA, length(model0b_curve@y.values[[1]])),
    exp(model1_fitted$.upper - log(model1_data$scale_factor)),
    exp(model2_fitted$.upper),
    exp(model3_fitted$.upper)
  ),
  observed = c(
    model0a_curve@x.values[[1]],
    model0b_curve@x.values[[1]],
    unlist(paired_cpue),
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
  xlab("Observed value (Models 1-3) or false positive rate (Model 0a and 0b)") +
  ylab("Fitted value (Models 1-3) or true positive rate (Model 0a and 0b)") +
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
  model0a_auc,
  model0b_auc,
  cor(exp(model1_fitted$.value - log(model1_data$scale_factor)), c(unlist(paired_cpue))),
  cor(exp(model2_fitted$.value), unlist(paired_richness)),
  cor(exp(model3_fitted$.value), model3_data$y)
)
fit_stats <- data.frame(
  model = model_names,
  statistic = c(rep("AUC", 2), rep("r", 3)),
  values = fit_stats
)
fit_stats <- fit_stats %>%
  mutate(model = gsub("draws_m", "M", model))
write.csv(fit_stats, file = "outputs/tables/fit-statistics.csv")
