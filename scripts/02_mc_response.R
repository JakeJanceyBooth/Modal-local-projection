# Response-function Monte Carlo
# Mean LP, Median LP, Modal LP, and VAR
# Nested samples across sample sizes

library(dplyr)

source(file.path("R", "estimators.R"))
source(file.path("R", "dgps.R"))

# Monte Carlo settings ----
sample_sizes <- c(250, 500, 1000)
mc_seed <- 12345
bw_constant <- 2.4
start_quantiles <- 0.5
H <- 20
horizons <- seq_len(H)
n_replications <- 1000
dgp_names <- c("gaussian", "skewed", "disaster", "rich_state", "downside_risk")
estimator_names <- c("Mean LP", "Median LP", "Modal LP", "VAR")
cheap_estimators <- c("Mean LP", "Median LP", "VAR")
max_sample_size <- max(sample_sizes)
results_directory <- file.path("results", "responses")

dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

run_label <- paste0("R", n_replications)
dgp_parameters <- vector("list", length(dgp_names))
names(dgp_parameters) <- dgp_names

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
set.seed(mc_seed)

for (replication in seq_len(n_replications)) {
  for (dgp_name in dgp_names) {
    # Simulate once at the largest sample size
    simulated_dgp <- switch(
      dgp_name,
      gaussian = simulate_gaussian_dgp(T = max_sample_size, horizon = H),
      skewed = simulate_skewed_dgp(T = max_sample_size, horizon = H),
      disaster = simulate_disaster_dgp(T = max_sample_size, horizon = H),
      rich_state = simulate_rich_state_dgp(T = max_sample_size, horizon = H),
      downside_risk = simulate_downside_risk_dgp(
        T = max_sample_size,
        horizon = H
      )
    )

    if (replication == 1) {
      dgp_parameters[[dgp_name]] <- simulated_dgp$parameters
    }

    truth_index <- match(horizons, simulated_dgp$true_response$horizon)

    # Population response assigned to each estimator
    # For DGPs 1-4 these are exact population response targets.
    # For DGP 5 they are pointwise nonlinear population benchmarks.
    # The pseudo-true targets of the linear estimators under DGP 5
    # are not calculated here.
    if (dgp_name == "downside_risk") {
      population_truth <- list(
        `Mean LP` = simulated_dgp$true_response$mean_response[truth_index],
        `Median LP` = simulated_dgp$true_response$median_response[truth_index],
        `Modal LP` = simulated_dgp$true_response$mode_response[truth_index],
        VAR = simulated_dgp$true_response$mean_response[truth_index]
      )
    } else {
      common_truth <- simulated_dgp$true_response$response[truth_index]

      population_truth <- list(
        `Mean LP` = common_truth,
        `Median LP` = common_truth,
        `Modal LP` = common_truth,
        VAR = common_truth
      )
    }

    if (anyNA(unlist(population_truth))) {
      stop(paste("Population responses are not aligned for", dgp_name))
    }

    for (sample_size in sample_sizes) {
      # Nested post-burn sample
      sample_data <- simulated_dgp$data[seq_len(sample_size), , drop = FALSE]
      y <- sample_data$y
      x <- sample_data$x

      # LP conditioning vector ----
      if (dgp_name == "rich_state") {
        z <- sample_data[, c("y", "y_lag", "s")]
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
        var_data <- sample_data[, c("x", "s", "y")]
        var_lags <- 2
      } else if (dgp_name == "downside_risk") {
        var_data <- sample_data[, c("x", "r", "y")]
        var_lags <- 1
      } else {
        var_data <- sample_data[, c("x", "y")]
        var_lags <- 1
      }

      # Mean LP: all horizons in one call ----
      mean_attempt <- tryCatch(
        {
          mean_fit <- fit_mean_lp(y = y, x = x, z = z, horizons = horizons)
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
          median_fit <- fit_median_lp(y = y, x = x, z = z, horizons = horizons)
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
          var_fit <- fit_var(data = var_data, lags = var_lags)
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
              horizons = h,
              bw_constant = bw_constant,
              start_quantiles = start_quantiles
            )

            modal_response_fit <- lp_response(
              fitted_lp = modal_fit,
              x = 0,
              z = reference_z,
              delta = 1
            )

            list(fit = modal_fit, estimate = modal_response_fit$response[1])
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
            modal_attempt$fit$n_converged_starts[1],
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
        estimator_attempt <- cheap_attempts[[estimator]]
        estimator_truth <- population_truth[[estimator]]
        rows <- other_row:(other_row + length(horizons) - 1)

        if (inherits(estimator_attempt, "error")) {
          other_results[rows, ] <- data.frame(
            dgp = rep(dgp_name, length(horizons)),
            estimator = rep(estimator, length(horizons)),
            sample_size = rep(sample_size, length(horizons)),
            replication = rep(replication, length(horizons)),
            horizon = horizons,
            estimate = rep(NA_real_, length(horizons)),
            truth = estimator_truth,
            error = rep(NA_real_, length(horizons)),
            failure = rep(conditionMessage(estimator_attempt), length(horizons))
          )
        } else {
          estimator_estimate <- estimator_attempt$estimate[
            match(horizons, estimator_attempt$horizon)
          ]

          estimator_error <- estimator_estimate - estimator_truth
          other_results[rows, ] <- data.frame(
            dgp = rep(dgp_name, length(horizons)),
            estimator = rep(estimator, length(horizons)),
            sample_size = rep(sample_size, length(horizons)),
            replication = rep(replication, length(horizons)),
            horizon = horizons,
            estimate = estimator_estimate,
            truth = estimator_truth,
            error = estimator_error,
            failure = rep(NA_character_, length(horizons))
          )
        }

        other_row <- other_row + length(horizons)
      }
    }
  }

  if (replication %% 25 == 0) {
    message("Completed replication ", replication, " of ", n_replications)
  }
}

# Save expensive raw output immediately ----
saveRDS(
  modal_results,
  file.path(results_directory, paste0("modal_results_", run_label, ".rds"))
)

saveRDS(
  other_results,
  file.path(results_directory, paste0("other_results_", run_label, ".rds"))
)

response_run_settings <- list(
  dgp_names = dgp_names,
  dgp_parameters = dgp_parameters,
  sample_sizes = sample_sizes,
  horizons = horizons,
  estimator_names = estimator_names,
  n_replications = n_replications,
  mc_seed = mc_seed,
  bw_constant = bw_constant,
  start_quantiles = start_quantiles,
  session_info = utils::sessionInfo(),
)

saveRDS(
  response_run_settings,
  file.path(
    results_directory,
    paste0("response_run_settings_", run_label, ".rds")
  )
)

# Monte Carlo summaries ----
response_comparison_summary <- dplyr::bind_rows(
  modal_results,
  other_results
) |>
  dplyr::mutate(
    success = is.na(failure) &
      is.finite(estimate)
  ) |>
  dplyr::group_by(
    dgp,
    estimator,
    sample_size,
    horizon
  ) |>
  dplyr::summarise(
    mean_estimate = if (any(success)) {
      mean(estimate[success])
    } else {
      NA_real_
    },
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
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    match(estimator, estimator_names),
    sample_size,
    horizon
  )

# extract population truths ----
response_truth <- bind_rows(
  modal_results |>
    dplyr::select(dgp, estimator, horizon, truth),
  other_results |>
    dplyr::select(dgp, estimator, horizon, truth)
) |>
  distinct()

# save everything .RDS----
saveRDS(
  modal_results,
  file.path(results_directory, paste0("modal_results_", run_label, ".rds"))
)

saveRDS(
  other_results,
  file.path(results_directory, paste0("other_results_", run_label, ".rds"))
)

saveRDS(
  response_comparison_summary,
  file.path(results_directory, paste0("response_summary_", run_label, ".rds"))
)

saveRDS(
  response_truth,
  file.path(results_directory, paste0("response_truth_", run_label, ".rds"))
)

saveRDS(
  response_plot_data,
  file.path(results_directory, paste0("response_plot_data_", run_label, ".rds"))
)
