label_alpha = label_bquote(cols = {alpha~"="~.(alpha)})

extract_intervals <- function(fit) {
  tibble(
    lower = fit$intervals[,1],
    upper = fit$intervals[,2],
    pred  = fit$predictions[,1]
  )
}

aci_theme <- theme(axis.title.y = element_text(size = 8))

basic_theme <- theme(
  axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(), 
  legend.position = "none"
)

simulation_one_plot <- function(results) {
  labels = ggplot2::label_bquote(
    rows = {alpha~"="~.(alpha)}, 
    cols = {psi~"="~theta~"="~.(param)}
  )
  
  plot_coverr <- results %>%
    ggplot(aes(x = method, y = coverage - alpha, color = method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
    facet_grid(alpha ~ param, labeller = labels) +
    labs(x = "", y = expression(CovErr(T)), subtitle = "Coverage Error", title = "Simulation study: ARMA errors") +
    basic_theme +
    aci_theme
  
  plot_piwidth <- results %>%
    ggplot(aes(x = method, y = mean_width, color = method)) +
    geom_boxplot() +
    facet_grid(alpha ~ param, scale = "free_y", labeller = labels) +
    labs(x = "", y = expression(MeanWidth(T)), subtitle = "Mean Interval Width") +
    basic_theme +
    aci_theme
  
  plot_pathlength <- results %>%
    ggplot(aes(x = method, y = path_length, color = method)) +
    geom_boxplot() +
    facet_grid(alpha ~ param, scale = "free_y", labeller = labels) +
    labs(x = "", y = expression(PathLength(T)), subtitle = "Path Length") +
    basic_theme +
    aci_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  plot_coverr / plot_piwidth / plot_pathlength
}

simulation_two_plot <- function(results) {
  results <- results %>% mutate(
    distribution_shift = ifelse(distribution_shift == 0.5, "Shift", "No shift")
  )
  
  labels = ggplot2::label_bquote(
    rows = .(distribution_shift),
    cols = {alpha~"="~.(alpha)}
  )
  
  plot_coverr <- results %>%
    ggplot(aes(x = method,
               y = coverage - alpha, color = method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
    facet_grid(distribution_shift ~ alpha, labeller = labels) +
    labs(x = "", y = expression(CovErr(T)), subtitle = "Coverage Error", 
         title = "Simulation Study: Distribution Shift") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    aci_theme
  
  plot_piwidth <- results %>%
    ggplot(aes(x = method,
               y = mean_width, color = method)) +
    geom_boxplot() +
    facet_grid(distribution_shift ~ alpha, labeller = labels, scales = "free_y") +
    labs(x = "", y = expression(MeanWidth(T)), subtitle = "Mean Interval Width") +
    aci_theme +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
  
  plot_pathlength <- results %>%
    ggplot(aes(x = method,
               y = path_length, color = method)) +
    geom_boxplot() +
    facet_grid(distribution_shift ~ alpha, labeller = labels, scales = "free_y") +
    labs(x = "ACI Method", y = expression(PathLength(T)), subtitle = "Path Length") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    aci_theme
  
  plot_coverr / plot_piwidth / plot_pathlength
}

case_study_plot <- function(results) {
  plot_coverr <- results %>%
    ggplot(aes(x = method, y = coverage - alpha, color = method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
    facet_wrap(~alpha, labeller = label_alpha) +
    labs(x = "", y = expression(CovErr(T))) +
    labs(subtitle = "Coverage Error", title = "Case Study Results") +
    aci_theme +
    basic_theme
  
  plot_piwidth <- results %>%
    ggplot(aes(x = method, y = mean_width, color = method)) +
    geom_boxplot() +
    facet_wrap(~alpha, labeller = label_alpha, scales = "free_y") +
    labs(x = "", y = expression(MeanWidth(T))) +
    labs(subtitle = "Interval Width") +
    aci_theme +
    basic_theme
  
  plot_pathlength <- results %>%
    ggplot(aes(x = method, y = path_length, color = method)) +
    geom_boxplot() +
    facet_wrap(~alpha, labeller = label_alpha) +
    labs(x = "ACI Method", y = expression(PathLength(T))) +
    labs(subtitle = "Path Length") +
    aci_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  
  plot_coverr / plot_piwidth / plot_pathlength
}


run_simulation_study1 <- function(data, setup, estimate_model, fit, workers = 8) {
  cache <- "~/aci-paper/cache/simulation_study1.rds"
  if(file.exists(cache)) return(read_rds(cache))
  
  plan(multisession, workers = workers)
  with_progress({
    p <- progressor(steps = nrow(data))
    data <- data %>%
      mutate(data = pmap(list(index, psi, xi, N), simulate),
             preds = future_map(data, estimate_model, p = p, 
                                .options = furrr_options(seed = TRUE)))
    
    setup <- expand_grid(
      data,
      setup
    )
    
    p <- progressor(steps = nrow(setup))
    simulation_study <- setup %>%
      mutate(
        fit = future_pmap(list(data, preds, method, alpha), fit, p = p, 
                          .options = furrr_options(seed = TRUE)),
        metrics = future_map(fit, metrics, .options = furrr_options(seed = TRUE))
      ) %>%
      select(-data)
  })
  
  results <- list(
    example_fits = simulation_study %>% 
      filter(param == 0.8, index == 1, alpha == 0.9),
    results =  simulation_study1 %>%
      mutate(metrics = map(metrics, as_tibble)) %>%
      select(-fit) %>%
      unnest(c(metrics)) 
  )
  
  write_rds(results, cache)
  
  results
}

run_simulation_study2 <- function(setup, fit, workers = 8) {
  cache <- "~/aci-paper/cache/simulation_study2.rds"
  if(file.exists(cache)) return(read_rds(cache))
  
  plan(multisession, workers = workers)
  
  with_progress({
    p <- progressor(steps = nrow(setup))
    
    simulation_study <- setup %>%
      mutate(
        fit = future_pmap(list(data, method, alpha), fit, p = p, 
                          .options = furrr_options(seed = TRUE)),
        metrics = map(fit, metrics)
      )
  })
  
  results <- list(
    example_fits = simulation_study %>% filter(distribution_shift == 0.5, index == 1, alpha == 0.9),
    results = simulation_study %>%
      select(-data) %>%
      mutate(metrics = map(metrics, as_tibble)) %>%
      select(-fit) %>%
      unnest(c(metrics)) 
  )
  
  write_rds(results, cache)
  
  results
}
