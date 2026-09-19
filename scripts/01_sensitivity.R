# Modal LP sensitivity Monte Carlo
# Bandwidth constants and deterministic starting values
# Nested samples across sample sizes

source(file.path("R", "estimators.R"))
source(file.path("R", "dgps.R"))

# Sensitivity settings ----
mc_seed <- 12345
set.seed(mc_seed)

dgp_names <- c("gaussian", "downside_risk")
sample_sizes <- c(250, 500, 1000)
horizons <- c(1, 5, 10, 20)
bw_constants <- c(1.6, 2.0, 2.4, 3.0, 3.6, 4.2, 6.4)
start_experiment_bw_constant <- 1.6
start_schemes <- list(
  two_starts = 0.5,
  six_starts = c(0.10, 0.25, 0.50, 0.75, 0.90)
)

n_replications <- 10
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
        estimate = modal_response$response[1],
        start_diagnostics = if (return_start_diagnostics) {
          modal_fit$start_diagnostics[[1]]
        } else {
          NULL
        }
      )
    },
    error = function(e) e
  )
}

# preallocate results ----
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
  start_diagnostics = I(vector("list", n_start_results)),
  failure = rep(NA_character_, n_start_results)
)

start_scalar_columns <- setdiff(
  names(start_sensitivity_raw),
  "start_diagnostics"
)

start_row <- 1

# bandwidth + number-of-MEM-starts sensitivity experiments ----
for (replication in seq_len(n_replications)) {
  for (dgp_name in dgp_names) {
    # Simulate once at the largest sample size
    simulated_dgp <- switch(
      dgp_name,
      gaussian = simulate_gaussian_dgp(T = max_sample_size, horizon = H),
      downside_risk = simulate_downside_risk_dgp(
        T = max_sample_size,
        horizon = H
      )
    )

    truth_index <- match(horizons, simulated_dgp$true_response$horizon)

    # Population modal-response truth
    if (dgp_name == "gaussian") {
      modal_truth <- simulated_dgp$true_response$response[truth_index]
    } else {
      modal_truth <- simulated_dgp$true_response$mode_response[truth_index]
    }

    for (sample_size in sample_sizes) {
      # Nested post-burn sample
      sample_data <- simulated_dgp$data[seq_len(sample_size), , drop = FALSE]
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
            start_quantiles = start_schemes$two_starts,
            reference_z = reference_z
          )

          # save to reuse this exact fit in Experiment B
          if (bw_constant == start_experiment_bw_constant) {
            baseline_two_start_attempt <- bandwidth_attempt
          }

          if (inherits(bandwidth_attempt, "error")) {
            bandwidth_sensitivity_raw[bandwidth_row, ] <- list(
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
            bandwidth_sensitivity_raw[bandwidth_row, ] <- list(
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
          bw_constant = start_experiment_bw_constant,
          start_quantiles = start_schemes$six_starts,
          reference_z = reference_z,
          return_start_diagnostics = TRUE
        )

        start_attempts <- list(
          two_starts = baseline_two_start_attempt,
          six_starts = six_start_attempt
        )

        for (start_scheme in names(start_attempts)) {
          start_attempt <- start_attempts[[start_scheme]]

          if (inherits(start_attempt, "error")) {
            start_sensitivity_raw[start_row, start_scalar_columns] <- list(
              dgp_name,
              sample_size,
              replication,
              horizon,
              start_scheme,
              start_experiment_bw_constant,
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

            start_sensitivity_raw$start_diagnostics[start_row] <- list(NULL)
          } else {
            fit <- start_attempt$fit
            estimate_h <- start_attempt$estimate
            start_sensitivity_raw[start_row, start_scalar_columns] <- list(
              dgp_name,
              sample_size,
              replication,
              horizon,
              start_scheme,
              start_experiment_bw_constant,
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

            start_sensitivity_raw$start_diagnostics[start_row] <-
              list(start_attempt$start_diagnostics)
          }

          start_row <- start_row + 1
        }
      }
    }
  }

  if (replication %% 5 == 0 || replication == n_replications) {
    message("Completed replication ", replication, " of ", n_replications)
  }
}

# bandwidth-sensitivity summary ----
bandwidth_sensitivity_summary <- bandwidth_sensitivity_raw |>
  dplyr::mutate(
    success = is.na(failure) &
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
  dplyr::arrange(match(dgp, dgp_names), sample_size, horizon, bw_constant)

# starting-value summary by scheme ----
start_scheme_summary <- start_sensitivity_raw |>
  dplyr::mutate(
    success = is.na(failure) &
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
pair_keys <- c("dgp", "sample_size", "replication", "horizon")
objective_tolerance <- 1e-8
coefficient_tolerance <- 1e-6

# Extract two-start results ----
two_start_results <- start_sensitivity_raw |>
  dplyr::filter(
    start_scheme == "two_starts"
  ) |>
  dplyr::select(
    dplyr::all_of(pair_keys),
    estimate,
    objective,
    converged,
    failure
  ) |>
  dplyr::rename(
    estimate_two = estimate,
    objective_two = objective,
    converged_two = converged,
    failure_two = failure
  )

# Extract six-start results ----
six_start_results <- start_sensitivity_raw |>
  dplyr::filter(
    start_scheme == "six_starts"
  ) |>
  dplyr::select(
    dplyr::all_of(pair_keys),
    estimate,
    objective,
    converged,
    selected_start,
    start_diagnostics,
    failure
  ) |>
  dplyr::rename(
    estimate_six = estimate,
    objective_six = objective,
    converged_six = converged,
    selected_start_six = selected_start,
    failure_six = failure
  )

# Two-start versus six-start comparison ----
start_pair_comparison <- dplyr::full_join(
  two_start_results,
  six_start_results |>
    dplyr::select(-start_diagnostics),
  by = pair_keys
) |>
  dplyr::mutate(
    paired_success = is.na(failure_two) &
      is.na(failure_six) &
      converged_two %in% TRUE &
      converged_six %in% TRUE &
      is.finite(estimate_two) &
      is.finite(estimate_six),

    estimate_difference = estimate_six - estimate_two,
    objective_gain = objective_six - objective_two,
    new_quantile_selected = selected_start_six %in%
      c("q10", "q25", "q75", "q90")
  )

start_pair_summary <- start_pair_comparison |>
  dplyr::group_by(
    dgp,
    sample_size,
    horizon
  ) |>
  dplyr::summarise(
    n_paired_success = sum(paired_success),
    mean_absolute_estimate_difference = if (any(paired_success)) {
      mean(abs(estimate_difference[paired_success]))
    } else {
      NA_real_
    },
    material_estimate_difference_rate = if (any(paired_success)) {
      mean(abs(estimate_difference[paired_success]) > coefficient_tolerance)
    } else {
      NA_real_
    },
    mean_objective_gain = if (any(paired_success)) {
      mean(objective_gain[paired_success])
    } else {
      NA_real_
    },
    higher_objective_rate = if (any(paired_success)) {
      mean(objective_gain[paired_success] > objective_tolerance)
    } else {
      NA_real_
    },
    negative_objective_gain_rate = if (any(paired_success)) {
      mean(objective_gain[paired_success] < -objective_tolerance)
    } else {
      NA_real_
    },
    new_quantile_selection_rate = if (any(paired_success)) {
      mean(new_quantile_selected[paired_success])
    } else {
      NA_real_
    },
    new_solution_rate = if (any(paired_success)) {
      mean(
        new_quantile_selected[paired_success] &
          objective_gain[paired_success] > objective_tolerance
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Individual-start diagnostics from six-start fits ----
diagnostic_rows <- lapply(
  seq_len(nrow(six_start_results)),
  function(i) {
    diagnostics_i <- six_start_results$start_diagnostics[[i]]

    if (is.null(diagnostics_i) ||
        !is.data.frame(diagnostics_i) ||
        nrow(diagnostics_i) == 0) {
      return(NULL)
    }

    successful_starts <- diagnostics_i$converged %in% TRUE &
      is.finite(diagnostics_i$objective)

    if (!any(successful_starts)) {
      return(NULL)
    }

    best_objective <- max(diagnostics_i$objective[successful_starts])

    diagnostics_i$objective_gap_from_best <- best_objective -
      diagnostics_i$objective
    diagnostics_i$tied_best <- successful_starts &
      diagnostics_i$objective_gap_from_best <= objective_tolerance

    diagnostics_i$unique_best <- diagnostics_i$tied_best &
      sum(diagnostics_i$tied_best) == 1

    diagnostics_i$selected <- diagnostics_i$start ==
      six_start_results$selected_start_six[i]

    selected_index <- match(
      six_start_results$selected_start_six[i],
      diagnostics_i$start
    )

    selected_coefficients <- diagnostics_i$coefficients[[selected_index]]

    diagnostics_i$coefficient_distance_from_selected <- vapply(
      diagnostics_i$coefficients,
      function(theta) {
        sqrt(sum((theta - selected_coefficients)^2))
      },
      numeric(1)
    )

    diagnostics_i$dgp <- six_start_results$dgp[i]
    diagnostics_i$sample_size <- six_start_results$sample_size[i]
    diagnostics_i$replication <- six_start_results$replication[i]
    diagnostics_i$horizon <- six_start_results$horizon[i]

    diagnostics_i
  }
)

start_diagnostics_long <- dplyr::bind_rows(diagnostic_rows)

# Summary by individual starting value ----
start_diagnostics_summary <- start_diagnostics_long |>
  dplyr::group_by(
    dgp,
    sample_size,
    horizon,
    start
  ) |>
  dplyr::summarise(
    n_runs = dplyr::n(),
    convergence_rate = mean(converged %in% TRUE),
    selection_rate = mean(selected %in% TRUE),
    tied_best_rate = mean(tied_best %in% TRUE),
    unique_best_rate = mean(unique_best %in% TRUE),
    mean_objective_gap = if (any(converged %in% TRUE)) {
      mean(objective_gap_from_best[converged %in% TRUE], na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_coefficient_distance = if (any(converged %in% TRUE)) {
      mean(
        coefficient_distance_from_selected[converged %in% TRUE],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    mean_iterations = if (any(converged %in% TRUE)) {
      mean(iterations[converged %in% TRUE], na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Add paired diagnostics----
start_sensitivity_summary <- start_scheme_summary |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    horizon,
    match(start_scheme, names(start_schemes))
  )

# save results ----
results_directory <- file.path("results", "sensitivity")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)
run_label <- paste0("R", n_replications, "_bw_extended")
sensitivity_run_settings <- list(
  dgp_names = dgp_names,
  sample_sizes = sample_sizes,
  horizons = horizons,
  bw_constants = bw_constants,
  start_experiment_bw_constant = start_experiment_bw_constant,
  start_schemes = start_schemes,
  n_replications = n_replications,
  mc_seed = mc_seed,
  session_info = utils::sessionInfo()
)

saveRDS(
  sensitivity_run_settings,
  file.path(
    results_directory,
    paste0("sensitivity_run_settings_", run_label, ".rds")
  )
)

saveRDS(
  bandwidth_sensitivity_raw,
  file.path(
    results_directory,
    paste0("bandwidth_sensitivity_raw_", run_label, ".rds")
  )
)

saveRDS(
  bandwidth_sensitivity_summary,
  file.path(
    results_directory,
    paste0("bandwidth_sensitivity_summary_", run_label, ".rds")
  )
)

saveRDS(
  start_sensitivity_raw,
  file.path(
    results_directory,
    paste0("start_sensitivity_raw_", run_label, ".rds")
  )
)

saveRDS(
  start_sensitivity_summary,
  file.path(
    results_directory,
    paste0("start_sensitivity_summary_", run_label, ".rds")
  )
)

saveRDS(
  start_pair_summary,
  file.path(results_directory, paste0("start_pair_summary_", run_label, ".rds"))
)

saveRDS(
  start_diagnostics_long,
  file.path(
    results_directory,
    paste0("start_diagnostics_long_", run_label, ".rds")
  )
)

saveRDS(
  start_diagnostics_summary,
  file.path(
    results_directory,
    paste0("start_diagnostics_summary_", run_label, ".rds")
  )
)
