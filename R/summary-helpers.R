# function to calculate posterior predictive plot from stan model draws
#   and stan data list
# - breaks and xlim must be specified. xlim is open at hte lower bound and 
#   closed at the upper bound
pp_check <- function(draws, obs, breaks, xlim) {
  
  # pull out predicted values
  pred <- draws %>%
    gather_draws(ypred[idx])

  # calculate count frequencies for both
  density_observed <- hist(obs[obs > xlim[1] & obs <= xlim[2]], breaks = breaks, plot = FALSE)
  density_pred <- pred %>%
    filter(.value > xlim[1], .value <= xlim[2]) %>%
    group_by(.chain, .iteration, .draw) %>%
    summarise(
      idx  = hist(.value, breaks = breaks, plot = FALSE)$mids,
      val = hist(.value, breaks = breaks, plot = FALSE)$count
    ) %>%
    ungroup() %>%
    group_by(idx) %>%
    summarise(
      lower = quantile(val, p = 0.1),
      upper = quantile(val, p = 0.9),
      val = median(val)
    )
  
  # combine these into a single dataframe for plotting
  density_all <- data.frame(
    value = c(
      density_observed$counts, 
      density_pred$val
    ),
    upper = c(
      density_observed$counts, 
      density_pred$upper
    ),
    lower = c(
      density_observed$counts, 
      density_pred$lower
    ),
    x = c(density_observed$mids, density_pred$idx),
    category = c(
      rep("Observed", length(density_observed$counts)),
      rep("Simulated", nrow(density_pred))
    )
  )
  
  # plot and return
  density_all %>%
    filter(!is.infinite(x), x > xlim[1], x <= xlim[2]) %>%
    ggplot(aes(x = x)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, col = category)) +
    geom_point(aes(y = value, col = category)) +
    ylab("Count") +
    xlab("Recruit abundance") +
    scale_colour_brewer(type = "qual", palette = "Set2", name = "Type") +
    theme(legend.position = "none")
  
}
