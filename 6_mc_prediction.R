# now compare simulation performance across all the estimators
# mode LP, mean LP, median LP, and VARs
# "common-target vs different-target" tests, maybe

# Realized-outcome forecasting Monte Carlo
# Fixed estimation with multiple forecast origins
# Mean LP, Median LP, Modal LP, and VAR

source("1_estimators.R")
source("2_dgps.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required. 
    Install it with install.packages('data.table')."
  )
}

# Forecasting settings ----

dgp_names <- c(
  "gaussian",
  "skewed",
  "disaster",
  "rich_state"
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

estimator_names <- c(
  "Mean LP",
  "Median LP",
  "Modal LP",
  "VAR"
)

n_forecast_origins <- 20
n_replications <- 100

bw_constant <- 2.4
start_quantiles <- 0.5

mc_seed <- 12345

max_sample_size <- max(sample_sizes)
H <- max(horizons)

# The first forecast origin is the common estimation cutoff.
# The last required outcome is:
# cutoff + (n_forecast_origins - 1) + H
estimation_cutoff <- max_sample_size

simulation_size <- max_sample_size +
  n_forecast_origins - 1 +
  H

forecast_origins <- estimation_cutoff +
  0:(n_forecast_origins - 1)


# Fixed-coefficient VAR forecast from an arbitrary origin ----

var_predict_from_origin <- function(fitted_var,
                                    history,
                                    variable,
                                    horizon) {
  
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
    (nrow(history) - p + 1):nrow(history),
    ,
    drop = FALSE
  ]
  
  forecast_path <- matrix(
    NA_real_,
    nrow = p + horizon,
    ncol = K,
    dimnames = list(
      NULL,
      variables
    )
  )
  
  forecast_path[seq_len(p), ] <- observed_lags
  
  # Recursively forecast all variables while holding the
  # estimated coefficients fixed
  for (step in seq_len(horizon)) {
    
    next_forecast <- intercept
    
    for (lag in seq_len(p)) {
      
      lagged_state <- forecast_path[
        p + step - lag,
        ,
        drop = TRUE
      ]
      
      next_forecast <- next_forecast +
        as.numeric(
          A[[lag]] %*% lagged_state
        )
    }
    
    forecast_path[p + step, ] <-
      next_forecast
  }
  
  variable_index <- match(
    variable,
    variables
  )
  
  data.frame(
    horizon = seq_len(horizon),
    variable = variable,
    prediction = forecast_path[
      p + seq_len(horizon),
      variable_index
    ]
  )
}


# Construct one raw forecast row ----

forecast_result_row <- function(dgp_name,
                                sample_size,
                                replication,
                                estimator,
                                forecast_origin,
                                target_time,
                                horizon,
                                forecast,
                                realized,
                                x_origin,
                                y_origin,
                                y_lag1_origin,
                                s_origin,
                                disaster_target,
                                n_disasters_to_target,
                                failure) {
  
  if (
    is.na(failure) &&
    !is.finite(forecast)
  ) {
    failure <- "Non-finite forecast."
  }
  
  successful <- is.na(failure) &&
    is.finite(forecast) &&
    is.finite(realized)
  
  if (successful) {
    
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
    x_origin = x_origin,
    y_origin = y_origin,
    y_lag1_origin = y_lag1_origin,
    s_origin = s_origin,
    disaster_target = disaster_target,
    n_disasters_to_target =
      n_disasters_to_target,
    failure = failure
  )
}

# Store one forecast row without copying the full data frame ----

store_forecast_row <- function(results,
                               row_index,
                               row_values) {
  
  for (column_name in names(row_values)) {
    
    data.table::set(
      results,
      i = row_index,
      j = column_name,
      value = row_values[[column_name]]
    )
  }
  
  invisible(NULL)
}


# Preallocate raw forecast results ----

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
  x_origin = rep(NA_real_, n_forecast_results),
  y_origin = rep(NA_real_, n_forecast_results),
  y_lag1_origin = rep(NA_real_, n_forecast_results),
  s_origin = rep(NA_real_, n_forecast_results),
  disaster_target = rep(NA_integer_, n_forecast_results),
  n_disasters_to_target = rep(NA_integer_, n_forecast_results),
  failure = rep(NA_character_, n_forecast_results)
)

forecast_row <- 1


# Preallocate Modal-LP fit diagnostics ----

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
  n_converged_starts = rep(
    NA_integer_,
    n_modal_diagnostics
  ),
  failure = rep(NA_character_, n_modal_diagnostics)
)

modal_diagnostic_row <- 1


# Retain the exact DGP parameters used in the run ----

dgp_parameters <- vector(
  "list",
  length(dgp_names)
)

names(dgp_parameters) <- dgp_names


# Forecasting Monte Carlo ----

set.seed(mc_seed)

for (replication in seq_len(n_replications)) {
  
  for (dgp_name in dgp_names) {
    
    # Simulate once for every sample size and estimator
    simulated_dgp <- switch(
      dgp_name,
      gaussian = simulate_gaussian_dgp(
        T = simulation_size,
        horizon = H
      ),
      skewed = simulate_skewed_dgp(
        T = simulation_size,
        horizon = H
      ),
      disaster = simulate_disaster_dgp(
        T = simulation_size,
        horizon = H
      ),
      rich_state = simulate_rich_state_dgp(
        T = simulation_size,
        horizon = H
      )
    )
    
    full_data <- simulated_dgp$data
    
    if (replication == 1) {
      dgp_parameters[[dgp_name]] <-
        simulated_dgp$parameters
    }
    
    for (sample_size in sample_sizes) {
      
      # Nested trailing training sample ending at the
      # common estimation cutoff
      training_rows <- (
        estimation_cutoff - sample_size + 1
      ):estimation_cutoff
      
      training_data <- full_data[
        training_rows,
        ,
        drop = FALSE
      ]
      
      y_training <- training_data$y
      x_training <- training_data$x
      
      
      # LP conditioning vector ----
      
      if (dgp_name == "rich_state") {
        
        z_training <- training_data[
          ,
          c(
            "y",
            "y_lag",
            "s"
          ),
          drop = FALSE
        ]
        
      } else {
        
        z_training <- training_data$y
      }
      
      
      # VAR system and lag order ----
      
      if (dgp_name == "rich_state") {
        
        var_columns <- c(
          "x",
          "s",
          "y"
        )
        
        var_lags <- 2
        
      } else {
        
        var_columns <- c(
          "x",
          "y"
        )
        
        var_lags <- 1
      }
      
      var_data_full <- full_data[
        ,
        var_columns,
        drop = FALSE
      ]
      
      var_training_data <- var_data_full[
        training_rows,
        ,
        drop = FALSE
      ]
      
      
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
          var_fit <- fit_var(
            data = var_training_data,
            lags = var_lags
          )
          
          # Validate the local recursion at the original
          # estimation cutoff against the package forecast
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
          y_lag1_origin <-
            full_data$y_lag[forecast_origin]
          
          if (dgp_name == "rich_state") {
            
            s_origin <- full_data$s[forecast_origin]
            
            z_origin <- c(
              y_origin,
              y_lag1_origin,
              s_origin
            )
            
          } else {
            
            s_origin <- NA_real_
            z_origin <- y_origin
          }
          
          if (inherits(fitted_attempt, "error")) {
            
            forecast_values <- rep(
              NA_real_,
              length(horizons)
            )
            
            failure_message <-
              conditionMessage(fitted_attempt)
            
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
              
              forecast_values <- rep(
                NA_real_,
                length(horizons)
              )
              
              failure_message <-
                conditionMessage(prediction_attempt)
              
            } else {
              
              forecast_values <-
                prediction_attempt$prediction[
                  match(
                    horizons,
                    prediction_attempt$horizon
                  )
                ]
              
              failure_message <- NA_character_
            }
          }
          
          for (h_index in seq_along(horizons)) {
            
            horizon <- horizons[h_index]
            target_time <- forecast_origin + horizon
            realized <- full_data$y[target_time]
            
            if (dgp_name == "disaster") {
              
              disaster_target <-
                full_data$disaster[target_time]
              
              n_disasters_to_target <- sum(
                full_data$disaster[
                  (forecast_origin + 1):target_time
                ]
              )
              
            } else {
              
              disaster_target <- NA_integer_
              n_disasters_to_target <- NA_integer_
            }
            
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
              x_origin = x_origin,
              y_origin = y_origin,
              y_lag1_origin =
                y_lag1_origin,
              s_origin = s_origin,
              disaster_target = disaster_target,
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
          
          modal_forecast_fit_diagnostics[
            modal_diagnostic_row,
          ] <- list(
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
          
          modal_forecast_fit_diagnostics[
            modal_diagnostic_row,
          ] <- list(
            dgp_name,
            sample_size,
            replication,
            horizon,
            modal_attempt$bandwidth[1],
            modal_attempt$objective[1],
            modal_attempt$converged[1],
            modal_attempt$iterations[1],
            modal_attempt$selected_start[1],
            modal_attempt$
              n_converged_starts[1],
            NA_character_
          )
        }
        
        modal_diagnostic_row <-
          modal_diagnostic_row + 1
        
        for (forecast_origin in forecast_origins) {
          
          x_origin <- full_data$x[forecast_origin]
          y_origin <- full_data$y[forecast_origin]
          y_lag1_origin <-
            full_data$y_lag[forecast_origin]
          
          if (dgp_name == "rich_state") {
            
            s_origin <- full_data$s[forecast_origin]
            
            z_origin <- c(
              y_origin,
              y_lag1_origin,
              s_origin
            )
            
          } else {
            
            s_origin <- NA_real_
            z_origin <- y_origin
          }
          
          if (inherits(modal_attempt, "error")) {
            
            forecast_value <- NA_real_
            
            failure_message <-
              conditionMessage(modal_attempt)
            
          } else {
            
            prediction_attempt <- tryCatch(
              lp_predict(
                fitted_lp = modal_attempt,
                x = x_origin,
                z = z_origin
              ),
              error = function(e) e
            )
            
            if (inherits(prediction_attempt, "error")) {
              
              forecast_value <- NA_real_
              
              failure_message <-
                conditionMessage(prediction_attempt)
              
            } else {
              
              forecast_value <-
                prediction_attempt$prediction[1]
              
              failure_message <- NA_character_
            }
          }
          
          target_time <- forecast_origin + horizon
          realized <- full_data$y[target_time]
          
          if (dgp_name == "disaster") {
            
            disaster_target <-
              full_data$disaster[target_time]
            
            n_disasters_to_target <- sum(
              full_data$disaster[
                (forecast_origin + 1):target_time
              ]
            )
            
          } else {
            
            disaster_target <- NA_integer_
            n_disasters_to_target <- NA_integer_
          }
          
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
            x_origin = x_origin,
            y_origin = y_origin,
            y_lag1_origin =
              y_lag1_origin,
            s_origin = s_origin,
            disaster_target = disaster_target,
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
        y_lag1_origin <-
          full_data$y_lag[forecast_origin]
        
        if (dgp_name == "rich_state") {
          s_origin <- full_data$s[forecast_origin]
        } else {
          s_origin <- NA_real_
        }
        
        if (inherits(var_attempt, "error")) {
          
          forecast_values <- rep(
            NA_real_,
            length(horizons)
          )
          
          failure_message <-
            conditionMessage(var_attempt)
          
        } else {
          
          # Use the actually observed p lag vectors ending at
          # this forecast origin. Future states are then
          # forecast recursively using fixed coefficients.
          history_rows <- (
            forecast_origin - var_lags + 1
          ):forecast_origin
          
          observed_history <- var_data_full[
            history_rows,
            ,
            drop = FALSE
          ]
          
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
            
            forecast_values <- rep(
              NA_real_,
              length(horizons)
            )
            
            failure_message <-
              conditionMessage(prediction_attempt)
            
          } else {
            
            forecast_values <-
              prediction_attempt$prediction[
                match(
                  horizons,
                  prediction_attempt$horizon
                )
              ]
            
            failure_message <- NA_character_
          }
        }
        
        for (h_index in seq_along(horizons)) {
          
          horizon <- horizons[h_index]
          target_time <- forecast_origin + horizon
          realized <- full_data$y[target_time]
          
          if (dgp_name == "disaster") {
            
            disaster_target <-
              full_data$disaster[target_time]
            
            n_disasters_to_target <- sum(
              full_data$disaster[
                (forecast_origin + 1):target_time
              ]
            )
            
          } else {
            
            disaster_target <- NA_integer_
            n_disasters_to_target <- NA_integer_
          }
          
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
            x_origin = x_origin,
            y_origin = y_origin,
            y_lag1_origin =
              y_lag1_origin,
            s_origin = s_origin,
            disaster_target = disaster_target,
            n_disasters_to_target =
              n_disasters_to_target,
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
  
  if (
    replication %% 5 == 0 ||
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

# Replication-level forecast scores ----
#
# Squared-error loss naturally favors the conditional mean.
# Absolute-error loss naturally favors the conditional median.
# Neither criterion is inherently targeted to the conditional mode.
#
# Forecast origins within a replication are dependent, so
# replication-level averages are formed before the final summary.

forecast_scores_by_replication <-
  forecast_results_raw |>
  dplyr::mutate(
    success =
      is.na(failure) &
      is.finite(forecast) &
      is.finite(realized)
  ) |>
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


# Overall forecast summary ----

forecast_summary <-
  forecast_scores_by_replication |>
  dplyr::mutate(
    score_available =
      n_success > 0 &
      is.finite(mean_squared_error)
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_forecast_error =
      if (any(score_available)) {
        mean(
          mean_forecast_error[
            score_available
          ]
        )
      } else {
        NA_real_
      },
    mae = if (any(score_available)) {
      mean(
        mean_absolute_error[
          score_available
        ]
      )
    } else {
      NA_real_
    },
    rmse = if (any(score_available)) {
      sqrt(
        mean(
          mean_squared_error[
            score_available
          ]
        )
      )
    } else {
      NA_real_
    },
    n_replications = dplyr::n(),
    n_successful_replications =
      sum(score_available),
    n_forecasts = sum(n_origins),
    n_successful_forecasts =
      sum(n_success),
    success_rate =
      n_successful_forecasts /
      n_forecasts,
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(estimator, estimator_names),
    horizon
  )


# Forecast plots ----

forecast_plot_data <-
  forecast_summary |>
  dplyr::mutate(
    dgp = factor(
      dgp,
      levels = dgp_names
    ),
    sample_size = factor(
      sample_size,
      levels = sample_sizes
    ),
    estimator = factor(
      estimator,
      levels = estimator_names
    )
  )


# RMSE 

forecast_rmse_plot <-
  ggplot2::ggplot(
    forecast_plot_data,
    ggplot2::aes(
      x = horizon,
      y = rmse,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size)
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "RMSE",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


# MAE 

forecast_mae_plot <-
  ggplot2::ggplot(
    forecast_plot_data,
    ggplot2::aes(
      x = horizon,
      y = mae,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size)
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "MAE",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


# Mean forecast error

forecast_error_plot <-
  ggplot2::ggplot(
    forecast_plot_data,
    ggplot2::aes(
      x = horizon,
      y = mean_forecast_error,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Mean forecast error",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

# Forecast success rates

forecast_success_plot <-
  ggplot2::ggplot(
    forecast_plot_data,
    ggplot2::aes(
      x = horizon,
      y = success_rate,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size)
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    labels = function(x) {
      paste0(
        round(100 * x),
        "%"
      )
    }
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Forecast success rate",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

# Average issued forecasts

average_forecast_data <-
  forecast_results_raw |>
  dplyr::filter(
    is.na(failure),
    is.finite(forecast)
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    average_forecast = mean(forecast),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    dgp = factor(
      dgp,
      levels = dgp_names
    ),
    sample_size = factor(
      sample_size,
      levels = sample_sizes
    ),
    estimator = factor(
      estimator,
      levels = estimator_names
    )
  )


average_forecast_plot <-
  ggplot2::ggplot(
    average_forecast_data,
    ggplot2::aes(
      x = horizon,
      y = average_forecast,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Average forecast",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


# Paired forecast-loss diagnostics ----
#
# Because all estimators forecast the same realization at the same
# origin, compare losses within the same forecast observation.
#
# Squared loss benchmark: Mean LP
# Absolute loss benchmark: Median LP


forecast_pair_data <-
  forecast_results_raw |>
  dplyr::mutate(
    success =
      is.na(failure) &
      is.finite(forecast) &
      is.finite(realized),
    estimator_key = dplyr::case_when(
      estimator == "Mean LP" ~ "mean_lp",
      estimator == "Median LP" ~ "median_lp",
      estimator == "Modal LP" ~ "modal_lp",
      estimator == "VAR" ~ "var"
    )
  ) |>
  dplyr::filter(success) |>
  dplyr::select(
    dgp,
    sample_size,
    replication,
    forecast_origin,
    target_time,
    horizon,
    estimator_key,
    squared_error,
    absolute_error
  ) |>
  tidyr::pivot_wider(
    names_from = estimator_key,
    values_from = c(
      squared_error,
      absolute_error
    ),
    names_sep = "__"
  )


# Paired squared-loss differences relative to Mean LP ----

paired_mse_by_replication <-
  dplyr::bind_rows(
    
    forecast_pair_data |>
      dplyr::transmute(
        dgp,
        sample_size,
        replication,
        forecast_origin,
        horizon,
        estimator = "Median LP",
        loss_difference =
          squared_error__median_lp -
          squared_error__mean_lp
      ),
    
    forecast_pair_data |>
      dplyr::transmute(
        dgp,
        sample_size,
        replication,
        forecast_origin,
        horizon,
        estimator = "Modal LP",
        loss_difference =
          squared_error__modal_lp -
          squared_error__mean_lp
      ),
    
    forecast_pair_data |>
      dplyr::transmute(
        dgp,
        sample_size,
        replication,
        forecast_origin,
        horizon,
        estimator = "VAR",
        loss_difference =
          squared_error__var -
          squared_error__mean_lp
      )
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_loss_difference =
      mean(loss_difference),
    .groups = "drop"
  )


paired_mse_summary <-
  paired_mse_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_loss_difference =
      mean(mean_loss_difference),
    mcse =
      stats::sd(mean_loss_difference) /
      sqrt(dplyr::n()),
    n_replications = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    lower =
      mean_loss_difference -
      1.96 * mcse,
    upper =
      mean_loss_difference +
      1.96 * mcse,
    dgp = factor(
      dgp,
      levels = dgp_names
    ),
    sample_size = factor(
      sample_size,
      levels = sample_sizes
    )
  )


paired_mse_plot <-
  ggplot2::ggplot(
    paired_mse_summary,
    ggplot2::aes(
      x = horizon,
      y = mean_loss_difference,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.4,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Squared-loss difference vs Mean LP",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


# Paired absolute-loss differences relative to Median LP ----

paired_mae_by_replication <-
  dplyr::bind_rows(
    
    forecast_pair_data |>
      dplyr::transmute(
        dgp,
        sample_size,
        replication,
        forecast_origin,
        horizon,
        estimator = "Mean LP",
        loss_difference =
          absolute_error__mean_lp -
          absolute_error__median_lp
      ),
    
    forecast_pair_data |>
      dplyr::transmute(
        dgp,
        sample_size,
        replication,
        forecast_origin,
        horizon,
        estimator = "Modal LP",
        loss_difference =
          absolute_error__modal_lp -
          absolute_error__median_lp
      ),
    
    forecast_pair_data |>
      dplyr::transmute(
        dgp,
        sample_size,
        replication,
        forecast_origin,
        horizon,
        estimator = "VAR",
        loss_difference =
          absolute_error__var -
          absolute_error__median_lp
      )
  ) |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_loss_difference =
      mean(loss_difference),
    .groups = "drop"
  )


paired_mae_summary <-
  paired_mae_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mean_loss_difference =
      mean(mean_loss_difference),
    mcse =
      stats::sd(mean_loss_difference) /
      sqrt(dplyr::n()),
    n_replications = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    lower =
      mean_loss_difference -
      1.96 * mcse,
    upper =
      mean_loss_difference +
      1.96 * mcse,
    dgp = factor(
      dgp,
      levels = dgp_names
    ),
    sample_size = factor(
      sample_size,
      levels = sample_sizes
    )
  )


paired_mae_plot <-
  ggplot2::ggplot(
    paired_mae_summary,
    ggplot2::aes(
      x = horizon,
      y = mean_loss_difference,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.4,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Absolute-loss difference vs Median LP",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


# Disaster-path diagnostics ----
#
# Split forecasts according to whether at least one disaster
# occurs after the forecast origin and through the target date.

disaster_path_data <-
  forecast_results_raw |>
  dplyr::mutate(
    success =
      is.na(failure) &
      is.finite(forecast) &
      is.finite(realized)
  ) |>
  dplyr::filter(
    dgp == "disaster",
    success,
    !is.na(n_disasters_to_target)
  ) |>
  dplyr::mutate(
    disaster_path = dplyr::if_else(
      n_disasters_to_target > 0,
      "At least one disaster",
      "No disaster"
    )
  )


# Counts use Mean LP only so the same forecast event is not
# counted four times simply because four estimators forecast it.

disaster_path_counts <-
  disaster_path_data |>
  dplyr::filter(
    estimator == "Mean LP"
  ) |>
  dplyr::count(
    sample_size,
    horizon,
    disaster_path,
    name = "n_forecasts"
  )


disaster_path_by_replication <-
  disaster_path_data |>
  dplyr::group_by(
    sample_size,
    replication,
    estimator,
    horizon,
    disaster_path
  ) |>
  dplyr::summarise(
    mse = mean(squared_error),
    mae = mean(absolute_error),
    n_origins = dplyr::n(),
    .groups = "drop"
  )


disaster_path_summary <-
  disaster_path_by_replication |>
  dplyr::group_by(
    sample_size,
    estimator,
    horizon,
    disaster_path
  ) |>
  dplyr::summarise(
    rmse = sqrt(
      mean(mse)
    ),
    mae = mean(mae),
    n_replications_with_event =
      dplyr::n(),
    n_event_forecasts =
      sum(n_origins),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    sample_size = factor(
      sample_size,
      levels = sample_sizes
    ),
    disaster_path = factor(
      disaster_path,
      levels = c(
        "No disaster",
        "At least one disaster"
      )
    ),
    estimator = factor(
      estimator,
      levels = estimator_names
    )
  )


disaster_path_rmse_plot <-
  ggplot2::ggplot(
    disaster_path_summary,
    ggplot2::aes(
      x = horizon,
      y = rmse,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(disaster_path),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Conditional RMSE",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


disaster_path_mae_plot <-
  ggplot2::ggplot(
    disaster_path_summary,
    ggplot2::aes(
      x = horizon,
      y = mae,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(disaster_path),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Conditional MAE",
    color = "Estimator",
    linetype = "Estimator"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )


forecast_plots <- list(
  rmse = forecast_rmse_plot,
  mae = forecast_mae_plot,
  mean_forecast_error = forecast_error_plot,
  success_rate = forecast_success_plot,
  average_forecast = average_forecast_plot,
  paired_mse_difference = paired_mse_plot,
  paired_mae_difference = paired_mae_plot,
  disaster_path_rmse = disaster_path_rmse_plot,
  disaster_path_mae = disaster_path_mae_plot
)
# Save results ----

dir.create(
  "results",
  showWarnings = FALSE
)

run_label <- paste0(
  "R",
  n_replications,
  "_O",
  n_forecast_origins,
  "_h",
  paste(
    horizons,
    collapse = "-"
  )
)

forecast_run_settings <- list(
  dgp_names = dgp_names,
  dgp_parameters = dgp_parameters,
  sample_sizes = sample_sizes,
  horizons = horizons,
  n_forecast_origins =
    n_forecast_origins,
  n_replications = n_replications,
  bw_constant = bw_constant,
  start_quantiles = start_quantiles,
  estimation_cutoff =
    estimation_cutoff,
  simulation_size = simulation_size,
  mc_seed = mc_seed
)

saveRDS(
  forecast_results_raw,
  file.path(
    "results",
    paste0(
      "forecast_results_raw_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  modal_forecast_fit_diagnostics,
  file.path(
    "results",
    paste0(
      "modal_forecast_fit_diagnostics_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  forecast_scores_by_replication,
  file.path(
    "results",
    paste0(
      "forecast_scores_by_replication_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  forecast_summary,
  file.path(
    "results",
    paste0(
      "forecast_summary_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  forecast_run_settings,
  file.path(
    "results",
    paste0(
      "forecast_run_settings_",
      run_label,
      ".rds"
    )
  )
)

saveRDS(
  forecast_plots,
  file.path(
    "results",
    paste0(
      "forecast_plots_",
      run_label,
      ".rds"
    )
  )
)