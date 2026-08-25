# Response-function Monte Carlo
# Mean LP, Median LP, Modal LP, and VAR
# Nested samples across sample sizes

source("1_estimators.R")
source("2_dgps.R")


# Monte Carlo settings ----

sample_sizes <- c(
  250,
  500,
  1000
)

H <- 20
horizons <- seq_len(H)

n_replications <- 100

dgp_names <- c(
  "gaussian",
  "skewed",
  "disaster",
  "rich_state",
  "downside_risk"
)

estimator_names <- c(
  "Mean LP",
  "Median LP",
  "Modal LP",
  "VAR"
)

cheap_estimators <- c(
  "Mean LP",
  "Median LP",
  "VAR"
)

max_sample_size <- max(sample_sizes)


# preallocate Modal LP results ----

n_modal_results <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(horizons)

modal_results <- data.frame(
  dgp = rep(NA_character_, n_modal_results),
  estimator = rep("Modal LP", n_modal_results),
  sample_size = rep(NA_integer_, n_modal_results),
  replication = rep(NA_integer_, n_modal_results),
  horizon = rep(NA_integer_, n_modal_results),
  estimate = rep(NA_real_, n_modal_results),
  truth = rep(NA_real_, n_modal_results),
  error = rep(NA_real_, n_modal_results),
  converged = rep(NA, n_modal_results),
  iterations = rep(NA_integer_, n_modal_results),
  selected_start = rep(NA_character_, n_modal_results),
  n_converged_starts = rep(NA_integer_, n_modal_results),
  failure = rep(NA_character_, n_modal_results)
)

modal_row <- 1


# preallocate Mean LP, Median LP, and VAR results ----

n_other_results <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(horizons) *
  length(cheap_estimators)

other_results <- data.frame(
  dgp = rep(NA_character_, n_other_results),
  estimator = rep(NA_character_, n_other_results),
  sample_size = rep(NA_integer_, n_other_results),
  replication = rep(NA_integer_, n_other_results),
  horizon = rep(NA_integer_, n_other_results),
  estimate = rep(NA_real_, n_other_results),
  truth = rep(NA_real_, n_other_results),
  error = rep(NA_real_, n_other_results),
  failure = rep(NA_character_, n_other_results)
)

other_row <- 1


# response Monte Carlo ----

set.seed(12345)

for (replication in seq_len(n_replications)) {
  
  for (dgp_name in dgp_names) {
    
    # Simulate once at the largest sample size
    simulated_dgp <- switch(
      dgp_name,
      gaussian = simulate_gaussian_dgp(
        T = max_sample_size,
        horizon = H
      ),
      skewed = simulate_skewed_dgp(
        T = max_sample_size,
        horizon = H
      ),
      disaster = simulate_disaster_dgp(
        T = max_sample_size,
        horizon = H
      ),
      rich_state = simulate_rich_state_dgp(
        T = max_sample_size,
        horizon = H
      ),
      downside_risk = simulate_downside_risk_dgp(
        T = max_sample_size,
        horizon = H
      )
    )
    
    truth_index <- match(
      horizons,
      simulated_dgp$true_response$horizon
    )
    
    # Population response assigned to each estimator
    if (dgp_name == "downside_risk") {
      
      population_truth <- list(
        `Mean LP` =
          simulated_dgp$true_response$
          mean_response[truth_index],
        `Median LP` =
          simulated_dgp$true_response$
          median_response[truth_index],
        `Modal LP` =
          simulated_dgp$true_response$
          mode_response[truth_index],
        VAR =
          simulated_dgp$true_response$
          mean_response[truth_index]
      )
      
    } else {
      
      common_truth <- simulated_dgp$true_response$
        response[truth_index]
      
      population_truth <- list(
        `Mean LP` = common_truth,
        `Median LP` = common_truth,
        `Modal LP` = common_truth,
        VAR = common_truth
      )
    }
    
    if (anyNA(unlist(population_truth))) {
      stop(
        paste(
          "Population responses are not aligned for",
          dgp_name
        )
      )
    }
    
    for (sample_size in sample_sizes) {
      
      # Nested post-burn sample
      sample_data <- simulated_dgp$data[
        seq_len(sample_size),
        ,
        drop = FALSE
      ]
      
      y <- sample_data$y
      x <- sample_data$x
      
      
      # LP conditioning vector ----
      
      if (dgp_name == "rich_state") {
        
        z <- sample_data[, c(
          "y",
          "y_lag",
          "s"
        )]
        
        reference_z <- c(0, 0, 0)
        
      } else if (dgp_name == "downside_risk") {
        
        z <- sample_data$r
        reference_z <- 0
        
      } else {
        
        z <- sample_data$y
        reference_z <- 0
      }
      
      
      # VAR system and lag order ----
      
      if (dgp_name == "rich_state") {
        
        var_data <- sample_data[, c(
          "x",
          "s",
          "y"
        )]
        
        var_lags <- 2
        
      } else if (dgp_name == "downside_risk") {
        
        var_data <- sample_data[, c(
          "x",
          "r",
          "y"
        )]
        
        var_lags <- 1
        
      } else {
        
        var_data <- sample_data[, c(
          "x",
          "y"
        )]
        
        var_lags <- 1
      }
      
      
      # Mean LP: all horizons in one call ----
      
      mean_attempt <- tryCatch(
        {
          mean_fit <- fit_mean_lp(
            y = y,
            x = x,
            z = z,
            horizons = horizons
          )
          
          mean_response_fit <- lp_response(
            fitted_lp = mean_fit,
            x = 0,
            z = reference_z,
            delta = 1
          )
          
          data.frame(
            horizon = mean_response_fit$horizon,
            estimate = mean_response_fit$response
          )
        },
        error = function(e) e
      )
      
      
      # Median LP: all horizons in one call ----
      
      median_attempt <- tryCatch(
        {
          median_fit <- fit_median_lp(
            y = y,
            x = x,
            z = z,
            horizons = horizons
          )
          
          median_response_fit <- lp_response(
            fitted_lp = median_fit,
            x = 0,
            z = reference_z,
            delta = 1
          )
          
          data.frame(
            horizon = median_response_fit$horizon,
            estimate = median_response_fit$response
          )
        },
        error = function(e) e
      )
      
      
      # VAR: one fit and all response horizons ----
      
      var_attempt <- tryCatch(
        {
          var_fit <- fit_var(
            data = var_data,
            lags = var_lags
          )
          
          var_response_fit <- var_response(
            fitted_var = var_fit,
            shock = "x",
            response = "y",
            horizon = H,
            delta = 1
          )
          
          data.frame(
            horizon = var_response_fit$horizon,
            estimate = var_response_fit$response
          )
        },
        error = function(e) e
      )
      
      
      # Modal LP: isolate each horizon ----
      
      for (h_index in seq_along(horizons)) {
        
        h <- horizons[h_index]
        truth_h <- population_truth[["Modal LP"]][h_index]
        
        modal_attempt <- tryCatch(
          {
            modal_fit <- fit_modal_lp(
              y = y,
              x = x,
              z = z,
              horizons = h
            )
            
            modal_response_fit <- lp_response(
              fitted_lp = modal_fit,
              x = 0,
              z = reference_z,
              delta = 1
            )
            
            list(
              fit = modal_fit,
              estimate =
                modal_response_fit$response[1]
            )
          },
          error = function(e) e
        )
        
        if (inherits(modal_attempt, "error")) {
          
          modal_results[modal_row, ] <- list(
            dgp_name,
            "Modal LP",
            sample_size,
            replication,
            h,
            NA_real_,
            truth_h,
            NA_real_,
            FALSE,
            NA_integer_,
            NA_character_,
            NA_integer_,
            conditionMessage(modal_attempt)
          )
          
        } else {
          
          estimate_h <- modal_attempt$estimate
          error_h <- estimate_h - truth_h
          
          modal_results[modal_row, ] <- list(
            dgp_name,
            "Modal LP",
            sample_size,
            replication,
            h,
            estimate_h,
            truth_h,
            error_h,
            modal_attempt$fit$converged[1],
            modal_attempt$fit$iterations[1],
            modal_attempt$fit$selected_start[1],
            modal_attempt$fit$
              n_converged_starts[1],
            NA_character_
          )
        }
        
        modal_row <- modal_row + 1
      }
      
      
      # Store Mean LP, Median LP, and VAR results ----
      
      cheap_attempts <- list(
        `Mean LP` = mean_attempt,
        `Median LP` = median_attempt,
        VAR = var_attempt
      )
      
      for (estimator in cheap_estimators) {
        
        estimator_attempt <-
          cheap_attempts[[estimator]]
        
        estimator_truth <-
          population_truth[[estimator]]
        
        rows <- other_row:(
          other_row + length(horizons) - 1
        )
        
        if (inherits(estimator_attempt, "error")) {
          
          other_results[rows, ] <- data.frame(
            dgp = rep(
              dgp_name,
              length(horizons)
            ),
            estimator = rep(
              estimator,
              length(horizons)
            ),
            sample_size = rep(
              sample_size,
              length(horizons)
            ),
            replication = rep(
              replication,
              length(horizons)
            ),
            horizon = horizons,
            estimate = rep(
              NA_real_,
              length(horizons)
            ),
            truth = estimator_truth,
            error = rep(
              NA_real_,
              length(horizons)
            ),
            failure = rep(
              conditionMessage(estimator_attempt),
              length(horizons)
            )
          )
          
        } else {
          
          estimator_estimate <-
            estimator_attempt$estimate[
              match(
                horizons,
                estimator_attempt$horizon
              )
            ]
          
          estimator_error <-
            estimator_estimate -
            estimator_truth
          
          other_results[rows, ] <- data.frame(
            dgp = rep(
              dgp_name,
              length(horizons)
            ),
            estimator = rep(
              estimator,
              length(horizons)
            ),
            sample_size = rep(
              sample_size,
              length(horizons)
            ),
            replication = rep(
              replication,
              length(horizons)
            ),
            horizon = horizons,
            estimate = estimator_estimate,
            truth = estimator_truth,
            error = estimator_error,
            failure = rep(
              NA_character_,
              length(horizons)
            )
          )
        }
        
        other_row <-
          other_row + length(horizons)
      }
    }
  }
  
  if (replication %% 25 == 0) {
    message(
      "Completed replication ",
      replication,
      " of ",
      n_replications
    )
  }
}


# Monte Carlo summaries ----

raw_result_tables <- list(
  modal = modal_results,
  other = other_results
)

summary_tables <- vector(
  "list",
  length(raw_result_tables)
)

names(summary_tables) <- names(raw_result_tables)

for (table_name in names(raw_result_tables)) {
  
  results <- raw_result_tables[[table_name]]
  
  result_groups <- split(
    results,
    interaction(
      results$dgp,
      results$estimator,
      results$sample_size,
      results$horizon,
      drop = TRUE
    )
  )
  
  summary_tables[[table_name]] <- do.call(
    rbind,
    lapply(
      result_groups,
      function(group) {
        
        successful <- is.na(group$failure) &
          is.finite(group$estimate)
        
        n_success <- sum(successful)
        
        data.frame(
          dgp = group$dgp[1],
          estimator = group$estimator[1],
          sample_size = group$sample_size[1],
          horizon = group$horizon[1],
          bias = if (n_success > 0) {
            mean(group$error[successful])
          } else {
            NA_real_
          },
          variance = if (n_success > 1) {
            var(group$estimate[successful])
          } else {
            NA_real_
          },
          rmse = if (n_success > 0) {
            sqrt(mean(
              group$error[successful]^2
            ))
          } else {
            NA_real_
          },
          n_replications = nrow(group),
          n_success = n_success,
          success_rate =
            n_success / nrow(group)
        )
      }
    )
  )
  
  rownames(
    summary_tables[[table_name]]
  ) <- NULL
}

modal_summary <- summary_tables$modal
other_summary <- summary_tables$other

response_comparison_summary <- rbind(
  modal_summary,
  other_summary
)

rownames(response_comparison_summary) <- NULL

response_comparison_summary <-
  response_comparison_summary[
    order(
      match(
        response_comparison_summary$dgp,
        dgp_names
      ),
      match(
        response_comparison_summary$estimator,
        estimator_names
      ),
      response_comparison_summary$sample_size,
      response_comparison_summary$horizon
    ),
  ]


# inspect results ----

head(modal_results)
head(other_results)
head(response_comparison_summary)

response_comparison_summary[
  response_comparison_summary$success_rate < 1,
  c(
    "dgp",
    "estimator",
    "sample_size",
    "horizon",
    "n_success",
    "success_rate"
  )
]

# Save output

saveRDS(response_comparison_summary, "response_comparison_summary.rds")

saveRDS(modal_results, "modal_response_results.rds")

saveRDS(other_results, "other_response_results.rds")

dir.create(
  "results",
  showWarnings = FALSE
)


# extract population truths ----

response_truth <- bind_rows(
  modal_results |>
    dplyr::select(
      dgp,
      estimator,
      horizon,
      truth
    ),
  other_results |>
    dplyr::select(
      dgp,
      estimator,
      horizon,
      truth
    )
) |>
  distinct()


# data for plotting ----

response_plot_data <-
  response_comparison_summary |>
  left_join(
    response_truth,
    by = c(
      "dgp",
      "estimator",
      "horizon"
    )
  ) |>
  mutate(
    estimator = factor(
      estimator,
      levels = c(
        "Mean LP",
        "Median LP",
        "Modal LP",
        "VAR"
      )
    ),
    average_estimate = truth + bias,
    mcse = sqrt(
      variance / n_success
    ),
    mc_lower =
      average_estimate - 1.96 * mcse,
    mc_upper =
      average_estimate + 1.96 * mcse
  )


# plot 1. RMSE overview ----

plot_rmse <- ggplot(
  response_plot_data,
  aes(
    horizon,
    rmse,
    colour = estimator,
    linetype = estimator
  )
) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    dgp ~ sample_size,
    scales = "free_y"
  ) +
  scale_colour_brewer(
    palette = "Dark2"
  ) +
  labs(
    x = "Horizon",
    y = "RMSE",
    colour = "Estimator",
    linetype = "Estimator"
  ) +
  theme_minimal()


# plot 2. bias overview ----

plot_bias <- ggplot(
  response_plot_data,
  aes(
    horizon,
    bias,
    colour = estimator,
    linetype = estimator
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    dgp ~ sample_size,
    scales = "free_y"
  ) +
  scale_colour_brewer(
    palette = "Dark2"
  ) +
  labs(
    x = "Horizon",
    y = "Bias",
    colour = "Estimator",
    linetype = "Estimator"
  ) +
  theme_minimal()


# plot 3. DGP 5 response separation ----

downside_plot_data <-
  response_plot_data |>
  filter(
    dgp == "downside_risk",
    estimator != "VAR"
  ) |>
  dplyr::select(
    estimator,
    sample_size,
    horizon,
    truth,
    average_estimate
  ) |>
  pivot_longer(
    c(
      truth,
      average_estimate
    ),
    names_to = "series",
    values_to = "response"
  )

plot_downside <- ggplot(
  downside_plot_data,
  aes(
    horizon,
    response,
    colour = estimator,
    linetype = series
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~ sample_size,
    nrow = 1
  ) +
  scale_colour_brewer(
    palette = "Dark2"
  ) +
  scale_linetype_manual(
    values = c(
      truth = "dashed",
      average_estimate = "solid"
    ),
    labels = c(
      truth = "Population truth",
      average_estimate =
        "Monte Carlo average"
    )
  ) +
  labs(
    x = "Horizon",
    y = "Response",
    colour = "Estimator",
    linetype = NULL
  ) +
  theme_minimal()


# plot 4. Gaussian common-target responses ----

gaussian_plot_data <-
  response_plot_data |>
  filter(
    dgp == "gaussian"
  )

plot_gaussian <- ggplot(
  gaussian_plot_data,
  aes(
    horizon,
    average_estimate,
    colour = estimator
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.9) +
  geom_line(
    data = gaussian_plot_data |>
      filter(
        estimator == "Mean LP"
      ),
    aes(
      x = horizon,
      y = truth
    ),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ sample_size,
    nrow = 1
  ) +
  scale_colour_brewer(
    palette = "Dark2"
  ) +
  labs(
    x = "Horizon",
    y = "Response",
    colour = "Estimator"
  ) +
  theme_minimal()


# plot 5. DGP 5 modal LP Monte Carlo uncertainty ----

plot_downside_modal_mcse <- ggplot(
  response_plot_data |>
    filter(
      dgp == "downside_risk",
      estimator == "Modal LP"
    ),
  aes(
    horizon,
    average_estimate
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_ribbon(
    aes(
      ymin = mc_lower,
      ymax = mc_upper
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      y = truth
    ),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ sample_size,
    nrow = 1
  ) +
  labs(
    x = "Horizon",
    y = "Modal response"
  ) +
  theme_minimal()


# store plots together in list ----

response_plots <- list(
  rmse = plot_rmse,
  bias = plot_bias,
  downside_response = plot_downside,
  gaussian_response = plot_gaussian,
  downside_modal_mcse =
    plot_downside_modal_mcse
)


# save everything .RDS----

run_label <- paste0(
  "R",
  n_replications
)

saveRDS(
  modal_results,
  file.path(
    "results",
    paste0(
      "modal_results_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  other_results,
  file.path(
    "results",
    paste0(
      "other_results_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  response_comparison_summary,
  file.path(
    "results",
    paste0(
      "response_summary_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  response_truth,
  file.path(
    "results",
    paste0(
      "response_truth_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  response_plot_data,
  file.path(
    "results",
    paste0(
      "response_plot_data_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  response_plots,
  file.path(
    "results",
    paste0(
      "response_plots_",
      run_label,
      ".rds"
    )
  )
)