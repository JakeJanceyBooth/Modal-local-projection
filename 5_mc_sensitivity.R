# all the modal LP sensitivity experiments
# different bandwidth constants
# optionally, later: different bandwidth rules entirely
# sensitivity of MEM to different starts
# MEM convergence/failure rates

# Modal LP sensitivity Monte Carlo
# Bandwidth constants and deterministic starting values
# Nested samples across sample sizes

source("1_estimators.R")
source("2_dgps.R")


# Sensitivity settings ----

dgp_names <- c(
  "gaussian",
  "downside_risk"
)

sample_sizes <- c(
  250,
  500,
  1000
)

horizons <- c(
  1,
  5,
  10,
  20
)

bw_constants <- c(
  0.8,
  1.2,
  1.6,
  2.0,
  2.4
)

baseline_bw_constant <- 1.6

start_schemes <- list(
  two_starts = 0.5,
  six_starts = c(
    0.10,
    0.25,
    0.50,
    0.75,
    0.90
  )
)

n_replications <- 5

max_sample_size <- max(sample_sizes)
H <- max(horizons)


# One Modal LP fit with isolated failure handling ----

fit_modal_once <- function(y, x, z, horizon,
                           bw_constant,
                           start_quantiles,
                           reference_z,
                           return_start_diagnostics = FALSE) {
  
  tryCatch(
    {
      modal_fit <- fit_modal_lp(
        y = y,
        x = x,
        z = z,
        horizons = horizon,
        bw_constant = bw_constant,
        start_quantiles = start_quantiles,
        return_start_diagnostics = return_start_diagnostics
      )
      
      modal_response <- lp_response(
        fitted_lp = modal_fit,
        x = 0,
        z = reference_z,
        delta = 1
      )
      
      list(
        fit = modal_fit,
        estimate = modal_response$response[1]
      )
    },
    error = function(e) e
  )
}


# preallocate bandwidth-sensitivity results ----

n_bandwidth_results <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(horizons) *
  length(bw_constants)

bandwidth_sensitivity_raw <- data.frame(
  dgp = rep(NA_character_, n_bandwidth_results),
  sample_size = rep(NA_integer_, n_bandwidth_results),
  replication = rep(NA_integer_, n_bandwidth_results),
  horizon = rep(NA_integer_, n_bandwidth_results),
  bw_constant = rep(NA_real_, n_bandwidth_results),
  bandwidth = rep(NA_real_, n_bandwidth_results),
  estimate = rep(NA_real_, n_bandwidth_results),
  truth = rep(NA_real_, n_bandwidth_results),
  error = rep(NA_real_, n_bandwidth_results),
  objective = rep(NA_real_, n_bandwidth_results),
  converged = rep(NA, n_bandwidth_results),
  iterations = rep(NA_integer_, n_bandwidth_results),
  selected_start = rep(NA_character_, n_bandwidth_results),
  n_converged_starts = rep(NA_integer_, n_bandwidth_results),
  failure = rep(NA_character_, n_bandwidth_results)
)

bandwidth_row <- 1


# preallocate starting-value-sensitivity results ----

n_start_results <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(horizons) *
  length(start_schemes)

start_sensitivity_raw <- data.frame(
  dgp = rep(NA_character_, n_start_results),
  sample_size = rep(NA_integer_, n_start_results),
  replication = rep(NA_integer_, n_start_results),
  horizon = rep(NA_integer_, n_start_results),
  start_scheme = rep(NA_character_, n_start_results),
  bw_constant = rep(NA_real_, n_start_results),
  bandwidth = rep(NA_real_, n_start_results),
  estimate = rep(NA_real_, n_start_results),
  truth = rep(NA_real_, n_start_results),
  error = rep(NA_real_, n_start_results),
  objective = rep(NA_real_, n_start_results),
  converged = rep(NA, n_start_results),
  iterations = rep(NA_integer_, n_start_results),
  selected_start = rep(NA_character_, n_start_results),
  n_converged_starts = rep(NA_integer_, n_start_results),
  failure = rep(NA_character_, n_start_results)
)

start_row <- 1


# paired sensitivity experiments ----

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
      downside_risk = simulate_downside_risk_dgp(
        T = max_sample_size,
        horizon = H
      )
    )
    
    truth_index <- match(
      horizons,
      simulated_dgp$true_response$horizon
    )
    
    # Population modal-response truth
    if (dgp_name == "gaussian") {
      
      modal_truth <- simulated_dgp$true_response$
        response[truth_index]
      
    } else {
      
      modal_truth <- simulated_dgp$true_response$
        mode_response[truth_index]
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
      
      if (dgp_name == "gaussian") {
        
        z <- sample_data$y
        reference_z <- 0
        
      } else {
        
        z <- sample_data$r
        reference_z <- 0
      }
      
      for (h_index in seq_along(horizons)) {
        
        horizon <- horizons[h_index]
        truth_h <- modal_truth[h_index]
        
        # Experiment A: bandwidth sensitivity ----
        
        baseline_two_start_attempt <- NULL
        
        for (bw_constant in bw_constants) {
          
          bandwidth_attempt <- fit_modal_once(
            y = y,
            x = x,
            z = z,
            horizon = horizon,
            bw_constant = bw_constant,
            start_quantiles =
              start_schemes$two_starts,
            reference_z = reference_z
          )
          
          # Reuse this exact fit in Experiment B
          if (bw_constant == baseline_bw_constant) {
            baseline_two_start_attempt <-
              bandwidth_attempt
          }
          
          if (inherits(bandwidth_attempt, "error")) {
            
            bandwidth_sensitivity_raw[
              bandwidth_row,
            ] <- list(
              dgp_name,
              sample_size,
              replication,
              horizon,
              bw_constant,
              NA_real_,
              NA_real_,
              truth_h,
              NA_real_,
              NA_real_,
              FALSE,
              NA_integer_,
              NA_character_,
              NA_integer_,
              conditionMessage(bandwidth_attempt)
            )
            
          } else {
            
            fit <- bandwidth_attempt$fit
            estimate_h <- bandwidth_attempt$estimate
            
            bandwidth_sensitivity_raw[
              bandwidth_row,
            ] <- list(
              dgp_name,
              sample_size,
              replication,
              horizon,
              bw_constant,
              fit$bandwidth[1],
              estimate_h,
              truth_h,
              estimate_h - truth_h,
              fit$objective[1],
              fit$converged[1],
              fit$iterations[1],
              fit$selected_start[1],
              fit$n_converged_starts[1],
              NA_character_
            )
          }
          
          bandwidth_row <- bandwidth_row + 1
        }
        
        
        # Experiment B: starting-value sensitivity ----
        
        six_start_attempt <- fit_modal_once(
          y = y,
          x = x,
          z = z,
          horizon = horizon,
          bw_constant = baseline_bw_constant,
          start_quantiles =
            start_schemes$six_starts,
          reference_z = reference_z,
          return_start_diagnostics = TRUE
        )
        
        start_attempts <- list(
          two_starts =
            baseline_two_start_attempt,
          six_starts =
            six_start_attempt
        )
        
        for (start_scheme in names(start_attempts)) {
          
          start_attempt <-
            start_attempts[[start_scheme]]
          
          if (inherits(start_attempt, "error")) {
            
            start_sensitivity_raw[
              start_row,
            ] <- list(
              dgp_name,
              sample_size,
              replication,
              horizon,
              start_scheme,
              baseline_bw_constant,
              NA_real_,
              NA_real_,
              truth_h,
              NA_real_,
              NA_real_,
              FALSE,
              NA_integer_,
              NA_character_,
              NA_integer_,
              conditionMessage(start_attempt)
            )
            
          } else {
            
            fit <- start_attempt$fit
            estimate_h <- start_attempt$estimate
            
            start_sensitivity_raw[
              start_row,
            ] <- list(
              dgp_name,
              sample_size,
              replication,
              horizon,
              start_scheme,
              baseline_bw_constant,
              fit$bandwidth[1],
              estimate_h,
              truth_h,
              estimate_h - truth_h,
              fit$objective[1],
              fit$converged[1],
              fit$iterations[1],
              fit$selected_start[1],
              fit$n_converged_starts[1],
              fit$start_diagnostics[[1]],
              NA_character_
            )
          }
          
          start_row <- start_row + 1
        }
      }
    }
  }
  
  if (
    replication %% 25 == 0 ||
    replication == n_replications
  ) {
    message(
      "Completed replication ",
      replication,
      " of ",
      n_replications
    )
  }
}


# bandwidth-sensitivity summary ----

bandwidth_sensitivity_summary <-
  bandwidth_sensitivity_raw |>
  dplyr::mutate(
    success =
      is.na(failure) &
      converged %in% TRUE &
      is.finite(estimate)
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    horizon,
    bw_constant
  ) |>
  dplyr::summarise(
    bias = if (any(success)) {
      mean(error[success])
    } else {
      NA_real_
    },
    variance = if (sum(success) > 1) {
      var(estimate[success])
    } else {
      NA_real_
    },
    rmse = if (any(success)) {
      sqrt(mean(error[success]^2))
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_success = sum(success),
    success_rate = mean(success),
    mean_bandwidth = if (any(success)) {
      mean(bandwidth[success])
    } else {
      NA_real_
    },
    mean_iterations = if (any(success)) {
      mean(iterations[success])
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    horizon,
    bw_constant
  )


# starting-value summary by scheme ----

start_scheme_summary <-
  start_sensitivity_raw |>
  dplyr::mutate(
    success =
      is.na(failure) &
      converged %in% TRUE &
      is.finite(estimate)
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    horizon,
    start_scheme
  ) |>
  dplyr::summarise(
    bias = if (any(success)) {
      mean(error[success])
    } else {
      NA_real_
    },
    variance = if (sum(success) > 1) {
      var(estimate[success])
    } else {
      NA_real_
    },
    rmse = if (any(success)) {
      sqrt(mean(error[success]^2))
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_success = sum(success),
    success_rate = mean(success),
    mean_bandwidth = if (any(success)) {
      mean(bandwidth[success])
    } else {
      NA_real_
    },
    mean_iterations = if (any(success)) {
      mean(iterations[success])
    } else {
      NA_real_
    },
    mean_objective = if (any(success)) {
      mean(objective[success])
    } else {
      NA_real_
    },
    .groups = "drop"
  )


# paired two-start versus six-start diagnostics ----

pair_keys <- c(
  "dgp",
  "sample_size",
  "replication",
  "horizon"
)

two_start_results <- start_sensitivity_raw[
  start_sensitivity_raw$start_scheme ==
    "two_starts",
  c(
    pair_keys,
    "estimate",
    "objective",
    "converged",
    "failure"
  ),
  drop = FALSE
]

names(two_start_results)[
  names(two_start_results) == "estimate"
] <- "estimate_two"

names(two_start_results)[
  names(two_start_results) == "objective"
] <- "objective_two"

names(two_start_results)[
  names(two_start_results) == "converged"
] <- "converged_two"

names(two_start_results)[
  names(two_start_results) == "failure"
] <- "failure_two"


six_start_results <- start_sensitivity_raw[
  start_sensitivity_raw$start_scheme ==
    "six_starts",
  c(
    pair_keys,
    "estimate",
    "objective",
    "converged",
    "selected_start",
    "failure"
  ),
  drop = FALSE
]

names(six_start_results)[
  names(six_start_results) == "estimate"
] <- "estimate_six"

names(six_start_results)[
  names(six_start_results) == "objective"
] <- "objective_six"

names(six_start_results)[
  names(six_start_results) == "converged"
] <- "converged_six"

names(six_start_results)[
  names(six_start_results) == "selected_start"
] <- "selected_start_six"

names(six_start_results)[
  names(six_start_results) == "failure"
] <- "failure_six"


start_pair_comparison <- merge(
  two_start_results,
  six_start_results,
  by = pair_keys,
  all = TRUE
)

start_pair_comparison$paired_success <-
  is.na(start_pair_comparison$failure_two) &
  is.na(start_pair_comparison$failure_six) &
  start_pair_comparison$converged_two %in% TRUE &
  start_pair_comparison$converged_six %in% TRUE &
  is.finite(start_pair_comparison$estimate_two) &
  is.finite(start_pair_comparison$estimate_six)

start_pair_comparison$estimate_difference <-
  start_pair_comparison$estimate_six -
  start_pair_comparison$estimate_two

start_pair_comparison$objective_gain <-
  start_pair_comparison$objective_six -
  start_pair_comparison$objective_two

start_pair_comparison$new_quantile_selected <-
  start_pair_comparison$selected_start_six %in%
  c(
    "q10",
    "q25",
    "q75",
    "q90"
  )


start_pair_summary <-
  start_pair_comparison |>
  dplyr::group_by(
    dgp,
    sample_size,
    horizon
  ) |>
  dplyr::summarise(
    n_paired_success =
      sum(paired_success),
    mean_absolute_estimate_difference =
      if (any(paired_success)) {
        mean(abs(
          estimate_difference[paired_success]
        ))
      } else {
        NA_real_
      },
    material_estimate_difference_rate =
      if (any(paired_success)) {
        mean(
          abs(
            estimate_difference[paired_success]
          ) > 1e-6
        )
      } else {
        NA_real_
      },
    mean_objective_gain =
      if (any(paired_success)) {
        mean(
          objective_gain[paired_success]
        )
      } else {
        NA_real_
      },
    higher_objective_rate =
      if (any(paired_success)) {
        mean(
          objective_gain[paired_success] >
            1e-10
        )
      } else {
        NA_real_
      },
    new_quantile_selection_rate =
      if (any(paired_success)) {
        mean(
          new_quantile_selected[
            paired_success
          ]
        )
      } else {
        NA_real_
      },
    new_solution_rate =
      if (any(paired_success)) {
        mean(
          new_quantile_selected[
            paired_success
          ] &
            objective_gain[paired_success] >
            1e-10
        )
      } else {
        NA_real_
      },
    .groups = "drop"
  )


# Pairwise columns on the six-start row compare it with two starts
start_sensitivity_summary <- merge(
  start_scheme_summary,
  start_pair_summary,
  by = c(
    "dgp",
    "sample_size",
    "horizon"
  ),
  all.x = TRUE,
  sort = FALSE
)

paired_columns <- setdiff(
  names(start_pair_summary),
  c(
    "dgp",
    "sample_size",
    "horizon"
  )
)

start_sensitivity_summary[
  start_sensitivity_summary$start_scheme ==
    "two_starts",
  paired_columns
] <- NA

start_sensitivity_summary <-
  start_sensitivity_summary[
    order(
      match(
        start_sensitivity_summary$dgp,
        dgp_names
      ),
      start_sensitivity_summary$sample_size,
      start_sensitivity_summary$horizon,
      match(
        start_sensitivity_summary$start_scheme,
        names(start_schemes)
      )
    ),
  ]

rownames(start_sensitivity_summary) <- NULL


# bandwidth diagnostic plots ----

plot_rmse <- ggplot(
  bandwidth_sensitivity_summary,
  aes(
    x = bw_constant,
    y = rmse,
    colour = factor(sample_size),
    group = factor(sample_size)
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(
    dgp ~ horizon,
    scales = "free_y",
    labeller = label_both
  ) +
  scale_x_continuous(
    breaks = bw_constants
  ) +
  labs(
    x = "Bandwidth constant",
    y = "RMSE",
    colour = "Sample size"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


plot_variance <- ggplot(
  bandwidth_sensitivity_summary,
  aes(
    x = bw_constant,
    y = variance,
    colour = factor(sample_size),
    group = factor(sample_size)
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(
    dgp ~ horizon,
    scales = "free_y",
    labeller = label_both
  ) +
  scale_x_continuous(
    breaks = bw_constants
  ) +
  labs(
    x = "Bandwidth constant",
    y = "Empirical variance",
    colour = "Sample size"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


plot_bias <- ggplot(
  bandwidth_sensitivity_summary,
  aes(
    x = bw_constant,
    y = bias,
    colour = factor(sample_size),
    group = factor(sample_size)
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(
    dgp ~ horizon,
    scales = "free_y",
    labeller = label_both
  ) +
  scale_x_continuous(
    breaks = bw_constants
  ) +
  labs(
    x = "Bandwidth constant",
    y = "Bias",
    colour = "Sample size"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


plot_success_rate <- ggplot(
  bandwidth_sensitivity_summary,
  aes(
    x = bw_constant,
    y = success_rate,
    colour = factor(sample_size),
    group = factor(sample_size)
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(
    dgp ~ horizon,
    labeller = label_both
  ) +
  scale_x_continuous(
    breaks = bw_constants
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    x = "Bandwidth constant",
    y = "Success rate",
    colour = "Sample size"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


sensitivity_plots <- list(
  rmse = plot_rmse,
  variance = plot_variance,
  bias = plot_bias,
  success_rate = plot_success_rate
)


# save results ----

dir.create(
  "results",
  showWarnings = FALSE
)

run_label <- paste0(
  "R",
  n_replications
)

saveRDS(
  bandwidth_sensitivity_raw,
  file.path(
    "results",
    paste0(
      "bandwidth_sensitivity_raw_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  bandwidth_sensitivity_summary,
  file.path(
    "results",
    paste0(
      "bandwidth_sensitivity_summary_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  start_sensitivity_raw,
  file.path(
    "results",
    paste0(
      "start_sensitivity_raw_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  start_sensitivity_summary,
  file.path(
    "results",
    paste0(
      "start_sensitivity_summary_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  sensitivity_plots,
  file.path(
    "results",
    paste0(
      "sensitivity_plots_",
      run_label,
      ".rds"
    )
  )
)