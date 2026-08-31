# Realized-outcome forecasting Monte Carlo
# Fixed estimation with multiple out-of-sample forecast origins
# Mean LP, Median LP, Modal LP, and VAR

source(file.path("R", "estimators.R"))
source(file.path("R", "dgps.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required. 
    Install it with install.packages('data.table')."
  )
}

# Forecasting settings ----
dgp_names <- c("gaussian", "skewed", "disaster", "rich_state")
sample_sizes <- c(250, 500, 1000)
horizons <- c(1, 2, 5, 10)
estimator_names <- c("Mean LP", "Median LP", "Modal LP", "VAR")
n_forecast_origins <- 20
n_replications <- 5000
bw_constant <- 2.4
start_quantiles <- 0.5
disaster_probability <- 0.05
disaster_size <- 7.5
mc_seed <- 12345
max_sample_size <- max(sample_sizes)
H <- max(horizons)
validate_var_recursion <- TRUE

# Retain exact DGP parameters
dgp_parameters <- vector("list", length(dgp_names))
names(dgp_parameters) <- dgp_names

# The first forecast origin is the common estimation cutoff.
# The last required outcome is: cutoff + (n_forecast_origins - 1) + H
estimation_cutoff <- max_sample_size
simulation_size <- max_sample_size + n_forecast_origins - 1 + H
forecast_origins <- estimation_cutoff + 0:(n_forecast_origins - 1)
half_widths <- c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00)
forecast_design_label <- "fixed_width_coverage"
results_directory <- file.path("results", "forecasting")

dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

run_label <- paste0(
  forecast_design_label, "_R", n_replications,
  "_O", n_forecast_origins, "_h",
  paste(horizons, collapse = "-")
)

save_run_object <- function(object, file_stem) {
  saveRDS(
    object,
    file.path(results_directory, paste0(file_stem, "_", run_label, ".rds"))
  )

  invisible(NULL)
}

run_started_at <- Sys.time()

# helper functions ----
# Fixed-coefficient VAR forecast from an arbitrary origin
var_predict_from_origin <- function(fitted_var, history, variable, horizon) {
  variables <- fitted_var$variables
  p <- fitted_var$lags
  K <- length(variables)
  history <- as.matrix(history)
  history <- history[, variables, drop = FALSE]

  if (nrow(history) < p) {
    stop("history must contain at least p observed rows.")
  }

  if (!(variable %in% variables)) {
    stop("variable must match a VAR variable name.")
  }

  # Estimated VAR lag matrices
  A <- vars::Acoef(fitted_var$fit)

  # Bcoef contains the lag coefficients followed by
  # the deterministic coefficients
  B <- vars::Bcoef(fitted_var$fit)

  if ("const" %in% colnames(B)) {
    intercept <- as.numeric(B[, "const"])
  } else {
    intercept <- numeric(K)
  }

  # Retain the p observed lag vectors ending at the origin
  observed_lags <- history[
    (nrow(history) - p + 1):nrow(history), , drop = FALSE
  ]

  forecast_path <- matrix(
    NA_real_,
    nrow = p + horizon,
    ncol = K,
    dimnames = list(NULL, variables)
  )

  forecast_path[seq_len(p), ] <- observed_lags

  # Recursively forecast all variables while holding the
  # estimated coefficients fixed
  for (step in seq_len(horizon)) {
    next_forecast <- intercept

    for (lag in seq_len(p)) {
      lagged_state <- forecast_path[p + step - lag, , drop = TRUE]
      next_forecast <- next_forecast + as.numeric(A[[lag]] %*% lagged_state)
    }

    forecast_path[p + step, ] <- next_forecast
  }

  variable_index <- match(variable, variables)

  data.frame(
    horizon = seq_len(horizon),
    variable = variable,
    prediction = forecast_path[p + seq_len(horizon), variable_index]
  )
}

# Construct one raw forecast row
forecast_result_row <- function(dgp_name,
                                sample_size,
                                replication,
                                estimator,
                                forecast_origin,
                                target_time,
                                horizon,
                                forecast,
                                realized,
                                n_disasters_to_target,
                                failure) {
  if (is.na(failure) && !is.finite(forecast)) {
    failure <- "Non-finite forecast."
  }

  successful <- is.na(failure) && is.finite(forecast) && is.finite(realized)

  if (successful) {
    # Note: forecast error is defined as forecast minus realized.
    forecast_error <- forecast - realized
    absolute_error <- abs(forecast_error)
    squared_error <- forecast_error^2
  } else {
    forecast_error <- NA_real_
    absolute_error <- NA_real_
    squared_error <- NA_real_
  }

  list(
    dgp = dgp_name,
    sample_size = sample_size,
    replication = replication,
    estimator = estimator,
    forecast_origin = forecast_origin,
    target_time = target_time,
    horizon = horizon,
    forecast = forecast,
    realized = realized,
    forecast_error = forecast_error,
    absolute_error = absolute_error,
    squared_error = squared_error,
    n_disasters_to_target = n_disasters_to_target,
    failure = failure
  )
}

# Store one forecast row without copying the full data frame
store_forecast_row <- function(results, row_index, row_values) {
  data.table::set(
    results,
    i = row_index,
    j = names(row_values),
    value = unname(row_values)
  )

  invisible(NULL)
}

# Count disasters to target (for DGP 3)
count_disasters_to_target <- function(data, dgp_name,
                                      forecast_origin, target_time) {
  if (dgp_name != "disaster") {
    return(NA_integer_)
  }

  sum(data$disaster[(forecast_origin + 1):target_time])
}

# Preallocate results ----
#forecast results
n_forecast_results <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(estimator_names) *
  length(horizons) *
  n_forecast_origins

forecast_results_raw <- data.frame(
  dgp = rep(NA_character_, n_forecast_results),
  sample_size = rep(NA_integer_, n_forecast_results),
  replication = rep(NA_integer_, n_forecast_results),
  estimator = rep(NA_character_, n_forecast_results),
  forecast_origin = rep(NA_integer_, n_forecast_results),
  target_time = rep(NA_integer_, n_forecast_results),
  horizon = rep(NA_integer_, n_forecast_results),
  forecast = rep(NA_real_, n_forecast_results),
  realized = rep(NA_real_, n_forecast_results),
  forecast_error = rep(NA_real_, n_forecast_results),
  absolute_error = rep(NA_real_, n_forecast_results),
  squared_error = rep(NA_real_, n_forecast_results),
  n_disasters_to_target = rep(NA_integer_, n_forecast_results),
  failure = rep(NA_character_, n_forecast_results)
)

forecast_row <- 1

#Modal-LP fit diagnostics
n_modal_diagnostics <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(horizons)

modal_forecast_fit_diagnostics <- data.frame(
  dgp = rep(NA_character_, n_modal_diagnostics),
  sample_size = rep(NA_integer_, n_modal_diagnostics),
  replication = rep(NA_integer_, n_modal_diagnostics),
  horizon = rep(NA_integer_, n_modal_diagnostics),
  bandwidth = rep(NA_real_, n_modal_diagnostics),
  objective = rep(NA_real_, n_modal_diagnostics),
  converged = rep(NA, n_modal_diagnostics),
  iterations = rep(NA_integer_, n_modal_diagnostics),
  selected_start = rep(NA_character_, n_modal_diagnostics),
  n_converged_starts = rep(NA_integer_, n_modal_diagnostics),
  failure = rep(NA_character_, n_modal_diagnostics)
)

modal_diagnostic_row <- 1

# Forecasting Monte Carlo ----
set.seed(mc_seed)

for (replication in seq_len(n_replications)) {
  for (dgp_name in dgp_names) {
    # Simulate once per replication and DGP, then reuse the
    # same path across sample sizes and estimators
    simulated_dgp <- switch(
      dgp_name,
      gaussian = simulate_gaussian_dgp(T = simulation_size, horizon = H),
      skewed = simulate_skewed_dgp(T = simulation_size, horizon = H),
      disaster = simulate_disaster_dgp(
        T = simulation_size,
        horizon = H,
        disaster_probability = disaster_probability,
        disaster_size = disaster_size
      ),
      rich_state = simulate_rich_state_dgp(T = simulation_size, horizon = H)
    )

    full_data <- simulated_dgp$data

    if (replication == 1) {
      dgp_parameters[[dgp_name]] <- simulated_dgp$parameters
    }

    for (sample_size in sample_sizes) {
      # Nested trailing training sample ending at the
      # common estimation cutoff
      training_rows <- (estimation_cutoff - sample_size + 1):estimation_cutoff
      training_data <- full_data[training_rows, , drop = FALSE]
      y_training <- training_data$y
      x_training <- training_data$x

      # LP conditioning vector ----
      if (dgp_name == "rich_state") {
        z_training <- training_data[, c("y", "y_lag", "s"), drop = FALSE]
      } else {
        z_training <- training_data$y
      }

      # VAR system and lag order ----
      if (dgp_name == "rich_state") {
        var_columns <- c("x", "s", "y")
        var_lags <- 2
      } else {
        var_columns <- c("x", "y")
        var_lags <- 1
      }

      var_data_full <- full_data[, var_columns, drop = FALSE]
      var_training_data <- var_data_full[training_rows, , drop = FALSE]

      # Direct-LP estimation alignment:
      #
      # Only observations through the estimation cutoff are
      # supplied to the estimators. At horizon h,
      # .prepare_lp_data() uses origins 1, ..., T-h and
      # outcomes 1+h, ..., T. Therefore, in calendar time,
      # the latest usable estimation origin is cutoff-h.
      # No outcome after the cutoff can enter estimation.
      # Mean LP: all requested horizons in one call ----
      mean_attempt <- tryCatch(
        fit_mean_lp(
          y = y_training,
          x = x_training,
          z = z_training,
          horizons = horizons
        ),
        error = function(e) e
      )

      # Median LP: all requested horizons in one call ----
      median_attempt <- tryCatch(
        fit_median_lp(
          y = y_training,
          x = x_training,
          z = z_training,
          horizons = horizons
        ),
        error = function(e) e
      )

      # VAR: estimate once and hold coefficients fixed ----
      var_attempt <- tryCatch(
        {
          var_fit <- fit_var(data = var_training_data, lags = var_lags)

          # Validate the custom recursion against the package forecast
          # once, using the first successful VAR fit.
          if (!var_recursion_validated) {
            package_forecast <- var_predict(
              fitted_var = var_fit,
              variable = "y",
              horizon = H
            )

            fixed_forecast <- var_predict_from_origin(
              fitted_var = var_fit,
              history = var_training_data,
              variable = "y",
              horizon = H
            )

            forecasts_match <- isTRUE(
              all.equal(
                package_forecast$prediction,
                fixed_forecast$prediction,
                tolerance = 1e-8,
                check.attributes = FALSE
              )
            )

            if (!forecasts_match) {
              stop(
                paste(
                  "Fixed-coefficient VAR recursion does not",
                  "match var_predict() at the estimation cutoff."
                )
              )
            }

            var_recursion_validated <- TRUE
          }

          var_fit
        },
        error = function(e) e
      )

      # Store Mean-LP and Median-LP forecasts ----
      lp_attempts <- list(
        `Mean LP` = mean_attempt,
        `Median LP` = median_attempt
      )

      for (estimator in names(lp_attempts)) {
        fitted_attempt <- lp_attempts[[estimator]]

        for (forecast_origin in forecast_origins) {
          x_origin <- full_data$x[forecast_origin]
          y_origin <- full_data$y[forecast_origin]
          y_lag1_origin <- full_data$y_lag[forecast_origin]

          if (dgp_name == "rich_state") {
            s_origin <- full_data$s[forecast_origin]
            z_origin <- c(y_origin, y_lag1_origin, s_origin)
          } else {
            s_origin <- NA_real_
            z_origin <- y_origin
          }

          if (inherits(fitted_attempt, "error")) {
            forecast_values <- rep(NA_real_, length(horizons))
            failure_message <- conditionMessage(fitted_attempt)
          } else {
            prediction_attempt <- tryCatch(
              lp_predict(
                fitted_lp = fitted_attempt,
                x = x_origin,
                z = z_origin
              ),
              error = function(e) e
            )

            if (inherits(prediction_attempt, "error")) {
              forecast_values <- rep(NA_real_, length(horizons))
              failure_message <- conditionMessage(prediction_attempt)
            } else {
              forecast_values <- prediction_attempt$prediction[
                match(horizons, prediction_attempt$horizon)
              ]

              failure_message <- NA_character_
            }
          }

          for (h_index in seq_along(horizons)) {
            horizon <- horizons[h_index]
            target_time <- forecast_origin + horizon
            realized <- full_data$y[target_time]
            n_disasters_to_target <- count_disasters_to_target(
              data = full_data,
              dgp_name = dgp_name,
              forecast_origin = forecast_origin,
              target_time = target_time
            )

            forecast_row_values <- forecast_result_row(
              dgp_name = dgp_name,
              sample_size = sample_size,
              replication = replication,
              estimator = estimator,
              forecast_origin = forecast_origin,
              target_time = target_time,
              horizon = horizon,
              forecast = forecast_values[h_index],
              realized = realized,
              n_disasters_to_target = n_disasters_to_target,
              failure = failure_message
            )

            store_forecast_row(
              results = forecast_results_raw,
              row_index = forecast_row,
              row_values = forecast_row_values
            )

            forecast_row <- forecast_row + 1
          }
        }
      }

      # Modal LP: fit and forecast one horizon at a time ----
      for (h_index in seq_along(horizons)) {
        horizon <- horizons[h_index]
        modal_attempt <- tryCatch(
          fit_modal_lp(
            y = y_training,
            x = x_training,
            z = z_training,
            horizons = horizon,
            bw_constant = bw_constant,
            start_quantiles = start_quantiles
          ),
          error = function(e) e
        )

        if (inherits(modal_attempt, "error")) {
          modal_forecast_fit_diagnostics[modal_diagnostic_row, ] <- list(
            dgp_name,
            sample_size,
            replication,
            horizon,
            NA_real_,
            NA_real_,
            FALSE,
            NA_integer_,
            NA_character_,
            NA_integer_,
            conditionMessage(modal_attempt)
          )
        } else {
          modal_forecast_fit_diagnostics[modal_diagnostic_row, ] <- list(
            dgp_name,
            sample_size,
            replication,
            horizon,
            modal_attempt$bandwidth[1],
            modal_attempt$objective[1],
            modal_attempt$converged[1],
            modal_attempt$iterations[1],
            modal_attempt$selected_start[1],
            modal_attempt$n_converged_starts[1],
            NA_character_
          )
        }

        modal_diagnostic_row <- modal_diagnostic_row + 1

        for (forecast_origin in forecast_origins) {
          x_origin <- full_data$x[forecast_origin]
          y_origin <- full_data$y[forecast_origin]
          y_lag1_origin <- full_data$y_lag[forecast_origin]

          if (dgp_name == "rich_state") {
            s_origin <- full_data$s[forecast_origin]
            z_origin <- c(y_origin, y_lag1_origin, s_origin)
          } else {
            s_origin <- NA_real_
            z_origin <- y_origin
          }

          if (inherits(modal_attempt, "error")) {
            forecast_value <- NA_real_
            failure_message <- conditionMessage(modal_attempt)
          } else {
            prediction_attempt <- tryCatch(
              lp_predict(fitted_lp = modal_attempt, x = x_origin, z = z_origin),
              error = function(e) e
            )

            if (inherits(prediction_attempt, "error")) {
              forecast_value <- NA_real_
              failure_message <- conditionMessage(prediction_attempt)
            } else {
              forecast_value <- prediction_attempt$prediction[1]
              failure_message <- NA_character_
            }
          }

          target_time <- forecast_origin + horizon
          realized <- full_data$y[target_time]
          n_disasters_to_target <- count_disasters_to_target(
            data = full_data,
            dgp_name = dgp_name,
            forecast_origin = forecast_origin,
            target_time = target_time
          )

          forecast_row_values <- forecast_result_row(
            dgp_name = dgp_name,
            sample_size = sample_size,
            replication = replication,
            estimator = "Modal LP",
            forecast_origin = forecast_origin,
            target_time = target_time,
            horizon = horizon,
            forecast = forecast_value,
            realized = realized,
            n_disasters_to_target = n_disasters_to_target,
            failure = failure_message
          )

          store_forecast_row(
            results = forecast_results_raw,
            row_index = forecast_row,
            row_values = forecast_row_values
          )

          forecast_row <- forecast_row + 1
        }
      }

      # Fixed-coefficient VAR forecasts ----
      for (forecast_origin in forecast_origins) {
        x_origin <- full_data$x[forecast_origin]
        y_origin <- full_data$y[forecast_origin]
        y_lag1_origin <- full_data$y_lag[forecast_origin]

        if (dgp_name == "rich_state") {
          s_origin <- full_data$s[forecast_origin]
        } else {
          s_origin <- NA_real_
        }

        if (inherits(var_attempt, "error")) {
          forecast_values <- rep(NA_real_, length(horizons))
          failure_message <- conditionMessage(var_attempt)
        } else {
          # Use the actually observed p lag vectors ending at
          # this forecast origin. Future states are then
          # forecast recursively using fixed coefficients.
          history_rows <- (forecast_origin - var_lags + 1):forecast_origin
          observed_history <- var_data_full[history_rows, , drop = FALSE]
          prediction_attempt <- tryCatch(
            var_predict_from_origin(
              fitted_var = var_attempt,
              history = observed_history,
              variable = "y",
              horizon = H
            ),
            error = function(e) e
          )

          if (inherits(prediction_attempt, "error")) {
            forecast_values <- rep(NA_real_, length(horizons))
            failure_message <- conditionMessage(prediction_attempt)
          } else {
            forecast_values <- prediction_attempt$prediction[
              match(horizons, prediction_attempt$horizon)
            ]

            failure_message <- NA_character_
          }
        }

        for (h_index in seq_along(horizons)) {
          horizon <- horizons[h_index]
          target_time <- forecast_origin + horizon
          realized <- full_data$y[target_time]
          n_disasters_to_target <- count_disasters_to_target(
            data = full_data,
            dgp_name = dgp_name,
            forecast_origin = forecast_origin,
            target_time = target_time
          )

          forecast_row_values <- forecast_result_row(
            dgp_name = dgp_name,
            sample_size = sample_size,
            replication = replication,
            estimator = "VAR",
            forecast_origin = forecast_origin,
            target_time = target_time,
            horizon = horizon,
            forecast = forecast_values[h_index],
            realized = realized,
            n_disasters_to_target = n_disasters_to_target,
            failure = failure_message
          )

          store_forecast_row(
            results = forecast_results_raw,
            row_index = forecast_row,
            row_values = forecast_row_values
          )

          forecast_row <- forecast_row + 1
        }
      }
    }
  }

  if (replication %% 50 == 0 || replication == n_replications) {
    message("Completed replication ", replication, " of ", n_replications)
  }
}

run_finished_at <- Sys.time()

# Check that every preallocated forecast row was filled
if (forecast_row != n_forecast_results + 1) {
  stop(
    paste(
      "Forecast row count mismatch:",
      forecast_row - 1,
      "rows filled but",
      n_forecast_results,
      "were expected."
    )
  )
}

# Fixed-width forecast-coverage settings ----
# Every interval is centered on an estimator's point forecast:
# [forecast - half_width, forecast + half_width].
# This is a fixed-width coverage comparison, not a
# loss function that relates specifically to the conditional mode.
mcse_from_replications <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) <= 1) {
    return(NA_real_)
  }

  stats::sd(x) / sqrt(length(x))
}

forecast_run_settings <- list(
  forecast_design_label = forecast_design_label,
  dgp_names = dgp_names,
  dgp_parameters = dgp_parameters,
  sample_sizes = sample_sizes,
  horizons = horizons,
  estimator_names = estimator_names,
  half_widths = half_widths,
  total_widths = 2 * half_widths,
  coverage_rule = "absolute_error <= half_width",
  n_forecast_origins = n_forecast_origins,
  n_replications = n_replications,
  bw_constant = bw_constant,
  start_quantiles = start_quantiles,
  disaster_probability = disaster_probability,
  disaster_size = disaster_size,
  estimation_cutoff = estimation_cutoff,
  simulation_size = simulation_size,
  mc_seed = mc_seed,
  run_started_at = run_started_at,
  run_finished_at = run_finished_at,
  session_info = utils::sessionInfo()
)

# Checkpoint expensive output before analysis ----
#
# Save these objects before substantial summarization or plotting so
# that a later analysis error cannot erase the expensive simulation.
save_run_object(forecast_results_raw, "forecast_results_raw")

save_run_object(
  modal_forecast_fit_diagnostics,
  "modal_forecast_fit_diagnostics"
)

save_run_object(forecast_run_settings, "forecast_run_settings")

# Common success indicator used by the derived analyses ----
forecast_evaluation_data <- forecast_results_raw |>
  dplyr::mutate(
    success = is.na(failure) &
      is.finite(forecast) &
      is.finite(realized)
  )

# Secondary conventional forecast scores ----
#
# Squared-error loss naturally favors the conditional mean.
# Absolute-error loss naturally favors the conditional median.
# Neither criterion is inherently targeted to the conditional mode.
# Origins are averaged within replication before the final summary.
forecast_scores_by_replication <- forecast_evaluation_data |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_forecast_error = if (any(success)) {
      mean(forecast_error[success])
    } else {
      NA_real_
    },
    mean_absolute_error = if (any(success)) {
      mean(absolute_error[success])
    } else {
      NA_real_
    },
    mean_squared_error = if (any(success)) {
      mean(squared_error[success])
    } else {
      NA_real_
    },
    n_origins = dplyr::n(),
    n_success = sum(success),
    success_rate = mean(success),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    replication,
    match(estimator, estimator_names),
    horizon
  )

forecast_summary <- forecast_scores_by_replication |>
  dplyr::mutate(
    score_available = n_success > 0 &
      is.finite(mean_squared_error)
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_forecast_error = if (any(score_available)) {
      mean(mean_forecast_error[score_available])
    } else {
      NA_real_
    },
    mae = if (any(score_available)) {
      mean(mean_absolute_error[score_available])
    } else {
      NA_real_
    },
    rmse = if (any(score_available)) {
      sqrt(mean(mean_squared_error[score_available]))
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_successful_replications = sum(score_available),
    n_forecasts = sum(n_origins),
    n_successful_forecasts = sum(n_success),
    success_rate = n_successful_forecasts / n_forecasts,
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(estimator, estimator_names),
    horizon
  )

# Average issued forecasts ----
average_forecast_by_replication <- forecast_evaluation_data |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    average_forecast = if (any(success)) {
      mean(forecast[success])
    } else {
      NA_real_
    },
    n_origins = dplyr::n(),
    n_success = sum(success),
    success_rate = mean(success),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    replication,
    match(estimator, estimator_names),
    horizon
  )

average_forecast_summary <- average_forecast_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mcse = mcse_from_replications(average_forecast),
    mean_average_forecast = if (any(is.finite(average_forecast))) {
      mean(average_forecast[is.finite(average_forecast)])
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_replications_available = sum(is.finite(average_forecast)),
    n_forecasts = sum(n_origins),
    n_successful_forecasts = sum(n_success),
    success_rate = n_successful_forecasts / n_forecasts,
    .groups = "drop"
  ) |>
  dplyr::rename(
    average_forecast = mean_average_forecast
  ) |>
  dplyr::mutate(
    mc_lower = average_forecast - 1.96 * mcse,
    mc_upper = average_forecast + 1.96 * mcse
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(estimator, estimator_names),
    horizon
  )

# Absolute fixed-width coverage ----
#
# The raw forecast table is never expanded by width. Each width is
# evaluated and immediately aggregated over origins within replication.
coverage_replication_list <- vector("list", length(half_widths))

for (width_index in seq_along(half_widths)) {
  current_half_width <- half_widths[width_index]
  coverage_replication_list[[width_index]] <- forecast_evaluation_data |>
    dplyr::group_by(
      dgp,
      sample_size,
      replication,
      estimator,
      horizon
    ) |>
    dplyr::summarise(
      coverage = if (any(success)) {
        mean(absolute_error[success] <= current_half_width)
      } else {
        NA_real_
      },
      n_origins = dplyr::n(),
      n_success = sum(success),
      success_rate = mean(success),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      half_width = current_half_width,
      total_width = 2 * current_half_width
    )
}

forecast_coverage_by_replication <-
  dplyr::bind_rows(coverage_replication_list) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    replication,
    match(estimator, estimator_names),
    horizon,
    half_width
  )

forecast_coverage_summary <- forecast_coverage_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon,
    half_width,
    total_width
  ) |>
  dplyr::summarise(
    mcse = mcse_from_replications(coverage),
    mean_coverage = if (any(is.finite(coverage))) {
      mean(coverage[is.finite(coverage)])
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_replications_available = sum(is.finite(coverage)),
    n_forecasts = sum(n_origins),
    n_successful_forecasts = sum(n_success),
    success_rate = n_successful_forecasts / n_forecasts,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    # Monte Carlo uncertainty for the simulated average coverage;
    # these are not forecast intervals or estimator confidence bands.
    mc_lower = mean_coverage - 1.96 * mcse,
    mc_upper = mean_coverage + 1.96 * mcse
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(estimator, estimator_names),
    horizon,
    half_width
  )

# Paired Modal-LP coverage comparisons ----
#
# Pair forecasts at the event level. Each comparison uses only events
# for which Modal LP and its comparator both produced valid forecasts.
paired_event_data <- forecast_evaluation_data |>
  dplyr::filter(estimator %in% c("Mean LP", "Median LP", "Modal LP")) |>
  dplyr::mutate(
    estimator_key = dplyr::case_when(
      estimator == "Mean LP" ~ "mean",
      estimator == "Median LP" ~ "median",
      estimator == "Modal LP" ~ "modal"
    )
  ) |>
  dplyr::select(
    dgp,
    sample_size,
    replication,
    forecast_origin,
    target_time,
    horizon,
    estimator_key,
    absolute_error,
    success
  ) |>
  tidyr::pivot_wider(
    names_from = estimator_key,
    values_from = c(absolute_error, success),
    names_sep = "__"
  )

paired_event_comparisons <- dplyr::bind_rows(
  paired_event_data |>
    dplyr::transmute(
      dgp,
      sample_size,
      replication,
      forecast_origin,
      target_time,
      horizon,
      comparison = "Modal - Mean",
      modal_absolute_error = absolute_error__modal,
      comparator_absolute_error = absolute_error__mean,
      common_success = dplyr::coalesce(success__modal, FALSE) &
        dplyr::coalesce(success__mean, FALSE)
    ),
  paired_event_data |>
    dplyr::transmute(
      dgp,
      sample_size,
      replication,
      forecast_origin,
      target_time,
      horizon,
      comparison = "Modal - Median",
      modal_absolute_error = absolute_error__modal,
      comparator_absolute_error = absolute_error__median,
      common_success = dplyr::coalesce(success__modal, FALSE) &
        dplyr::coalesce(success__median, FALSE)
    )
)

paired_common_success_counts <- paired_event_comparisons |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    horizon,
    comparison
  ) |>
  dplyr::summarise(
    expected_n_common_success = sum(common_success),
    .groups = "drop"
  )

paired_replication_list <- vector("list", length(half_widths))

for (width_index in seq_along(half_widths)) {
  current_half_width <- half_widths[width_index]
  paired_replication_list[[width_index]] <- paired_event_comparisons |>
    dplyr::mutate(
      modal_hit = common_success &
        modal_absolute_error <= current_half_width,
      comparator_hit = common_success &
        comparator_absolute_error <= current_half_width
    ) |>
    dplyr::group_by(
      dgp,
      sample_size,
      replication,
      horizon,
      comparison
    ) |>
    dplyr::summarise(
      mean_coverage_difference = if (any(common_success)) {
        mean(
          as.numeric(modal_hit[common_success]) -
            as.numeric(comparator_hit[common_success])
        )
      } else {
        NA_real_
      },
      modal_coverage = if (any(common_success)) {
        mean(modal_hit[common_success])
      } else {
        NA_real_
      },
      comparator_coverage = if (any(common_success)) {
        mean(comparator_hit[common_success])
      } else {
        NA_real_
      },
      n_events = dplyr::n(),
      n_common_success = sum(common_success),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      half_width = current_half_width,
      total_width = 2 * current_half_width
    )
}

paired_coverage_by_replication <- dplyr::bind_rows(paired_replication_list) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    replication,
    match(comparison, c("Modal - Mean", "Modal - Median")),
    horizon,
    half_width
  )

paired_coverage_summary <- paired_coverage_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    comparison,
    horizon,
    half_width,
    total_width
  ) |>
  dplyr::summarise(
    mcse = mcse_from_replications(mean_coverage_difference),
    mean_difference = if (any(is.finite(mean_coverage_difference))) {
      mean(mean_coverage_difference[is.finite(mean_coverage_difference)])
    } else {
      NA_real_
    },
    mean_modal_coverage = if (any(is.finite(modal_coverage))) {
      mean(modal_coverage[is.finite(modal_coverage)])
    } else {
      NA_real_
    },
    mean_comparator_coverage = if (any(is.finite(comparator_coverage))) {
      mean(comparator_coverage[is.finite(comparator_coverage)])
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_replications_available = sum(is.finite(mean_coverage_difference)),
    n_paired_forecasts = sum(n_common_success),
    .groups = "drop"
  ) |>
  dplyr::rename(
    mean_coverage_difference = mean_difference
  ) |>
  dplyr::mutate(
    # Monte Carlo uncertainty for the paired simulated average.
    mc_lower = mean_coverage_difference - 1.96 * mcse,
    mc_upper = mean_coverage_difference + 1.96 * mcse
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(comparison, c("Modal - Mean", "Modal - Median")),
    horizon,
    half_width
  )

# DGP-3 disaster-path coverage diagnostic ----
disaster_forecast_data <- forecast_evaluation_data |>
  dplyr::mutate(
    disaster_path = dplyr::case_when(
      n_disasters_to_target == 0 ~ "No disaster",
      n_disasters_to_target > 0 ~ "At least one disaster",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(dgp == "disaster", !is.na(disaster_path))

disaster_coverage_replication_list <- vector("list", length(half_widths))

for (width_index in seq_along(half_widths)) {
  current_half_width <- half_widths[width_index]
  disaster_coverage_replication_list[[width_index]] <- disaster_forecast_data |>
    dplyr::group_by(
      dgp,
      sample_size,
      replication,
      estimator,
      horizon,
      disaster_path
    ) |>
    dplyr::summarise(
      coverage = if (any(success)) {
        mean(absolute_error[success] <= current_half_width)
      } else {
        NA_real_
      },
      n_event_forecasts = dplyr::n(),
      n_success = sum(success),
      success_rate = mean(success),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      half_width = current_half_width,
      total_width = 2 * current_half_width
    )
}

disaster_coverage_by_replication <- dplyr::bind_rows(
  disaster_coverage_replication_list
) |>
  dplyr::arrange(
    sample_size,
    replication,
    match(estimator, estimator_names),
    horizon,
    match(disaster_path, c("No disaster", "At least one disaster")),
    half_width
  )

disaster_coverage_summary <- disaster_coverage_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon,
    disaster_path,
    half_width,
    total_width
  ) |>
  dplyr::summarise(
    mcse = mcse_from_replications(coverage),
    mean_coverage = if (any(is.finite(coverage))) {
      mean(coverage[is.finite(coverage)])
    } else {
      NA_real_
    },
    n_replications_with_event = dplyr::n(),
    n_replications_available = sum(is.finite(coverage)),
    n_event_forecasts = sum(n_event_forecasts),
    n_successful_event_forecasts = sum(n_success),
    success_rate = n_successful_event_forecasts / n_event_forecasts,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mc_lower = mean_coverage - 1.96 * mcse,
    mc_upper = mean_coverage + 1.96 * mcse
  ) |>
  dplyr::arrange(
    sample_size,
    match(estimator, estimator_names),
    horizon,
    match(disaster_path, c("No disaster", "At least one disaster")),
    half_width
  )

# Final checks before saving derived output ----
# 1. Coverage must be weakly increasing in half-width within every
# replication-level evaluation cell.
coverage_monotonicity_failures <- forecast_coverage_by_replication |>
  dplyr::filter(is.finite(coverage)) |>
  dplyr::arrange(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon,
    half_width
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    monotone = all(diff(coverage) >= -1e-12),
    .groups = "drop"
  ) |>
  dplyr::filter(!monotone)

disaster_monotonicity_failures <- disaster_coverage_by_replication |>
  dplyr::filter(is.finite(coverage)) |>
  dplyr::arrange(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon,
    disaster_path,
    half_width
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon,
    disaster_path
  ) |>
  dplyr::summarise(
    monotone = all(diff(coverage) >= -1e-12),
    .groups = "drop"
  ) |>
  dplyr::filter(!monotone)

if (nrow(coverage_monotonicity_failures) > 0 ||
    nrow(disaster_monotonicity_failures) > 0) {
  stop(
    paste(
      "Coverage is not weakly increasing in",
      "half-width in at least one cell."
    )
  )
}

# 2. Stored paired counts must be consistent across widths
# with the pair-specific common-success counts.
paired_count_check <- paired_coverage_by_replication |>
  dplyr::select(
    dgp,
    sample_size,
    replication,
    horizon,
    comparison,
    half_width,
    n_common_success
  ) |>
  dplyr::left_join(
    paired_common_success_counts,
    by = c("dgp", "sample_size", "replication", "horizon", "comparison")
  ) |>
  dplyr::filter(
    is.na(expected_n_common_success) |
      n_common_success != expected_n_common_success
  )

if (nrow(paired_count_check) > 0) {
  stop(
    paste(
      "At least one paired coverage cell does not",
      "use the correct common-success observations."
    )
  )
}

# 3. The raw table must retain one row per original issued forecast
# and must not acquire width-indexed analysis columns.
if (nrow(forecast_results_raw) != n_forecast_results) {
  stop("The raw forecast table changed size during analysis.")
}

if (any(c("half_width", "total_width", "coverage", "hit") %in%
        names(forecast_results_raw))) {
  stop("Width-indexed columns were added to forecast_results_raw.")
}

# 4. Required outputs must be nonempty and width-indexed outputs must
# contain exactly the six requested half-widths.
required_output_objects <- list(
  average_forecast_by_replication = average_forecast_by_replication,
  average_forecast_summary = average_forecast_summary,
  forecast_coverage_by_replication = forecast_coverage_by_replication,
  forecast_coverage_summary = forecast_coverage_summary,
  paired_coverage_by_replication = paired_coverage_by_replication,
  paired_coverage_summary = paired_coverage_summary,
  disaster_coverage_by_replication = disaster_coverage_by_replication,
  disaster_coverage_summary = disaster_coverage_summary
)

empty_output_objects <- names(required_output_objects)[
  vapply(
    required_output_objects,
    function(x) {
      !is.data.frame(x) || nrow(x) == 0
    },
    logical(1)
  )
]

if (length(empty_output_objects) > 0) {
  stop(
    paste(
      "Required output objects are empty:",
      paste(empty_output_objects, collapse = ", ")
    )
  )
}

width_indexed_objects <- list(
  forecast_coverage_by_replication,
  forecast_coverage_summary,
  paired_coverage_by_replication,
  paired_coverage_summary,
  disaster_coverage_by_replication,
  disaster_coverage_summary
)

valid_width_sets <- vapply(
  width_indexed_objects,
  function(x) {
    setequal(sort(unique(x$half_width)), half_widths)
  },
  logical(1)
)

if (!all(valid_width_sets)) {
  stop(
    paste(
      "At least one width-indexed output does not",
      "contain exactly the six requested half-widths."
    )
  )
}

if (!setequal(
    unique(paired_coverage_summary$comparison),
    c("Modal - Mean", "Modal - Median")
)) {
  stop("The paired output does not contain the requested comparisons.")
}

# Save analysis objects ----
save_run_object(
  average_forecast_by_replication,
  "average_forecast_by_replication"
)

save_run_object(average_forecast_summary, "average_forecast_summary")

save_run_object(
  forecast_coverage_by_replication,
  "forecast_coverage_by_replication"
)

save_run_object(forecast_coverage_summary, "forecast_coverage_summary")

save_run_object(
  paired_coverage_by_replication,
  "paired_coverage_by_replication"
)

save_run_object(paired_coverage_summary, "paired_coverage_summary")

save_run_object(
  disaster_coverage_by_replication,
  "disaster_coverage_by_replication"
)

save_run_object(disaster_coverage_summary, "disaster_coverage_summary")

save_run_object(forecast_summary, "forecast_summary")

message("Forecast analysis complete. Run label: ", run_label)
