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
  2,
  5,
  10
)

estimator_names <- c(
  "Mean LP",
  "Median LP",
  "Modal LP",
  "VAR"
)

n_forecast_origins <- 20
n_replications <- 5000

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

# Fixed-width forecast-coverage settings ----

# Every interval is centered on an estimator's point forecast:
# [forecast - half_width, forecast + half_width].
# This is a fixed-width coverage comparison, not a newly proposed
# loss function that uniquely elicits the conditional mode.

half_widths <- c(
  0.25,
  0.50,
  0.75,
  1.00,
  1.50,
  2.00
)

results_directory <- file.path(
  "results",
  "forecasting"
)

dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

run_label <- paste0(
  forecast_design_label,
  "_R",
  n_replications,
  "_O",
  n_forecast_origins,
  "_h",
  paste(horizons, collapse = "-")
)

save_run_object <- function(object, file_stem) {
  saveRDS(
    object,
    file.path(
      results_directory,
      paste0(file_stem, "_", run_label, ".rds")
    )
  )
  
  invisible(NULL)
}

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
  coverage_rule =
    "absolute_error <= half_width",
  n_forecast_origins = n_forecast_origins,
  n_replications = n_replications,
  bw_constant = bw_constant,
  start_quantiles = start_quantiles,
  estimation_cutoff = estimation_cutoff,
  simulation_size = simulation_size,
  mc_seed = mc_seed
)


# Checkpoint expensive output before analysis ----
#
# Save these objects before substantial summarization or plotting so
# that a later analysis error cannot erase the expensive simulation.

save_run_object(
  forecast_results_raw,
  "forecast_results_raw"
)

save_run_object(
  modal_forecast_fit_diagnostics,
  "modal_forecast_fit_diagnostics"
)

save_run_object(
  forecast_run_settings,
  "forecast_run_settings"
)


# Common success indicator used by the derived analyses ----

forecast_evaluation_data <-
  forecast_results_raw |>
  dplyr::mutate(
    success =
      is.na(failure) &
      is.finite(forecast) &
      is.finite(realized)
  )


# Secondary conventional forecast scores ----
#
# Squared-error loss naturally favors the conditional mean.
# Absolute-error loss naturally favors the conditional median.
# Neither criterion is inherently targeted to the conditional mode.
# Origins are averaged within replication before the final summary.

forecast_scores_by_replication <-
  forecast_evaluation_data |>
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
    success_rate =
      n_successful_forecasts / n_forecasts,
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(estimator, estimator_names),
    horizon
  )


# Average issued forecasts ----

average_forecast_by_replication <-
  forecast_evaluation_data |>
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

average_forecast_summary <-
  average_forecast_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    estimator,
    horizon
  ) |>
  dplyr::summarise(
    mcse =
      mcse_from_replications(average_forecast),
    mean_average_forecast =
      if (any(is.finite(average_forecast))) {
        mean(
          average_forecast[
            is.finite(average_forecast)
          ]
        )
      } else {
        NA_real_
      },
    n_replications = dplyr::n(),
    n_replications_available =
      sum(is.finite(average_forecast)),
    n_forecasts = sum(n_origins),
    n_successful_forecasts = sum(n_success),
    success_rate =
      n_successful_forecasts / n_forecasts,
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

coverage_replication_list <- vector(
  "list",
  length(half_widths)
)

for (width_index in seq_along(half_widths)) {
  current_half_width <- half_widths[width_index]
  
  coverage_replication_list[[width_index]] <-
    forecast_evaluation_data |>
    dplyr::group_by(
      dgp,
      sample_size,
      replication,
      estimator,
      horizon
    ) |>
    dplyr::summarise(
      coverage = if (any(success)) {
        mean(
          absolute_error[success] <=
            current_half_width
        )
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

forecast_coverage_summary <-
  forecast_coverage_by_replication |>
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
    n_replications_available =
      sum(is.finite(coverage)),
    n_forecasts = sum(n_origins),
    n_successful_forecasts = sum(n_success),
    success_rate =
      n_successful_forecasts / n_forecasts,
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

paired_event_data <-
  forecast_evaluation_data |>
  dplyr::filter(
    estimator %in% c(
      "Mean LP",
      "Median LP",
      "Modal LP"
    )
  ) |>
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

paired_event_comparisons <-
  dplyr::bind_rows(
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
        comparator_absolute_error =
          absolute_error__mean,
        common_success =
          dplyr::coalesce(success__modal, FALSE) &
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
        comparator_absolute_error =
          absolute_error__median,
        common_success =
          dplyr::coalesce(success__modal, FALSE) &
          dplyr::coalesce(success__median, FALSE)
      )
  )

paired_common_success_counts <-
  paired_event_comparisons |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    horizon,
    comparison
  ) |>
  dplyr::summarise(
    expected_n_common_success =
      sum(common_success),
    .groups = "drop"
  )

paired_replication_list <- vector(
  "list",
  length(half_widths)
)

for (width_index in seq_along(half_widths)) {
  current_half_width <- half_widths[width_index]
  
  paired_replication_list[[width_index]] <-
    paired_event_comparisons |>
    dplyr::mutate(
      modal_hit =
        common_success &
        modal_absolute_error <= current_half_width,
      comparator_hit =
        common_success &
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
      mean_coverage_difference =
        if (any(common_success)) {
          mean(
            as.numeric(modal_hit[common_success]) -
              as.numeric(
                comparator_hit[common_success]
              )
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

paired_coverage_by_replication <-
  dplyr::bind_rows(paired_replication_list) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    replication,
    match(
      comparison,
      c("Modal - Mean", "Modal - Median")
    ),
    horizon,
    half_width
  )

paired_coverage_summary <-
  paired_coverage_by_replication |>
  dplyr::group_by(
    dgp,
    sample_size,
    comparison,
    horizon,
    half_width,
    total_width
  ) |>
  dplyr::summarise(
    mcse = mcse_from_replications(
      mean_coverage_difference
    ),
    mean_difference =
      if (any(is.finite(mean_coverage_difference))) {
        mean(
          mean_coverage_difference[
            is.finite(mean_coverage_difference)
          ]
        )
      } else {
        NA_real_
      },
    mean_modal_coverage =
      if (any(is.finite(modal_coverage))) {
        mean(modal_coverage[is.finite(modal_coverage)])
      } else {
        NA_real_
      },
    mean_comparator_coverage =
      if (any(is.finite(comparator_coverage))) {
        mean(
          comparator_coverage[
            is.finite(comparator_coverage)
          ]
        )
      } else {
        NA_real_
      },
    n_replications = dplyr::n(),
    n_replications_available =
      sum(is.finite(mean_coverage_difference)),
    n_paired_forecasts = sum(n_common_success),
    .groups = "drop"
  ) |>
  dplyr::rename(
    mean_coverage_difference = mean_difference
  ) |>
  dplyr::mutate(
    # Monte Carlo uncertainty for the paired simulated average.
    mc_lower =
      mean_coverage_difference - 1.96 * mcse,
    mc_upper =
      mean_coverage_difference + 1.96 * mcse
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    match(
      comparison,
      c("Modal - Mean", "Modal - Median")
    ),
    horizon,
    half_width
  )


# DGP-3 disaster-path coverage diagnostic ----

disaster_forecast_data <-
  forecast_evaluation_data |>
  dplyr::mutate(
    disaster_path = dplyr::case_when(
      n_disasters_to_target == 0 ~
        "No disaster",
      n_disasters_to_target > 0 ~
        "At least one disaster",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(
    dgp == "disaster",
    !is.na(disaster_path)
  )

disaster_coverage_replication_list <- vector(
  "list",
  length(half_widths)
)

for (width_index in seq_along(half_widths)) {
  current_half_width <- half_widths[width_index]
  
  disaster_coverage_replication_list[[width_index]] <-
    disaster_forecast_data |>
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
        mean(
          absolute_error[success] <=
            current_half_width
        )
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

disaster_coverage_by_replication <-
  dplyr::bind_rows(
    disaster_coverage_replication_list
  ) |>
  dplyr::arrange(
    sample_size,
    replication,
    match(estimator, estimator_names),
    horizon,
    match(
      disaster_path,
      c("No disaster", "At least one disaster")
    ),
    half_width
  )

disaster_coverage_summary <-
  disaster_coverage_by_replication |>
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
    n_replications_available =
      sum(is.finite(coverage)),
    n_event_forecasts = sum(n_event_forecasts),
    n_successful_event_forecasts = sum(n_success),
    success_rate =
      n_successful_event_forecasts /
      n_event_forecasts,
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
    match(
      disaster_path,
      c("No disaster", "At least one disaster")
    ),
    half_width
  )


# Plot data ----

dgp_labels <- c(
  gaussian = "Gaussian",
  skewed = "Skewed lognormal",
  disaster = "Disaster mixture",
  rich_state = "Rich state"
)

average_forecast_plot_data <-
  average_forecast_summary |>
  dplyr::mutate(
    dgp = factor(
      dgp,
      levels = dgp_names,
      labels = unname(dgp_labels[dgp_names])
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

forecast_coverage_plot_data <-
  forecast_coverage_summary |>
  dplyr::mutate(
    dgp = factor(
      dgp,
      levels = dgp_names,
      labels = unname(dgp_labels[dgp_names])
    ),
    horizon = factor(horizon, levels = horizons),
    estimator = factor(
      estimator,
      levels = estimator_names
    )
  )

paired_coverage_plot_data <-
  paired_coverage_summary |>
  dplyr::mutate(
    dgp = factor(
      dgp,
      levels = dgp_names,
      labels = unname(dgp_labels[dgp_names])
    ),
    horizon = factor(horizon, levels = horizons),
    comparison = factor(
      comparison,
      levels = c("Modal - Mean", "Modal - Median")
    )
  )

disaster_coverage_plot_data <-
  disaster_coverage_summary |>
  dplyr::mutate(
    horizon = factor(horizon, levels = horizons),
    estimator = factor(
      estimator,
      levels = estimator_names
    ),
    disaster_path = factor(
      disaster_path,
      levels = c(
        "No disaster",
        "At least one disaster"
      )
    )
  )


# Figure 1: average issued forecasts ----

average_forecast_plot <-
  ggplot2::ggplot(
    average_forecast_plot_data,
    ggplot2::aes(
      x = horizon,
      y = average_forecast,
      color = estimator,
      linetype = estimator,
      group = estimator
    )
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mc_lower, ymax = mc_upper),
    width = 0.35,
    linewidth = 0.35,
    alpha = 0.65,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    size = 1.8,
    na.rm = TRUE
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(dgp),
    cols = ggplot2::vars(sample_size),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(breaks = horizons) +
  ggplot2::labs(
    x = "Horizon",
    y = "Average issued forecast",
    color = "Estimator",
    linetype = "Estimator",
    caption = paste(
      "Vertical bars show +/- 1.96 Monte Carlo SE",
      "based on replication-level averages."
    )
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")


# Figure 2: absolute fixed-width coverage ----
#
# Create one faceted plot per sample size. This retains all sample
# sizes without placing twelve estimator/sample-size curves together.

forecast_coverage_plots <- lapply(
  sample_sizes,
  function(current_sample_size) {
    current_plot_data <-
      forecast_coverage_plot_data |>
      dplyr::filter(
        sample_size == current_sample_size
      )
    
    ggplot2::ggplot(
      current_plot_data,
      ggplot2::aes(
        x = half_width,
        y = mean_coverage,
        color = estimator,
        linetype = estimator,
        group = estimator
      )
    ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mc_lower, ymax = mc_upper),
        width = 0.03,
        linewidth = 0.3,
        alpha = 0.60,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        linewidth = 0.8,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        size = 1.6,
        na.rm = TRUE
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(dgp),
        cols = ggplot2::vars(horizon)
      ) +
      ggplot2::scale_x_continuous(breaks = half_widths) +
      ggplot2::scale_y_continuous(
        labels = function(x) {
          paste0(round(100 * x), "%")
        }
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste0(
          "Fixed-width forecast coverage: T = ",
          current_sample_size
        ),
        subtitle =
          "Half-width a; total interval width is 2a",
        x = "Interval half-width",
        y = "Coverage",
        color = "Estimator",
        linetype = "Estimator",
        caption = paste(
          "Vertical bars show +/- 1.96 Monte Carlo SE",
          "based on replication-level coverage rates."
        )
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
      )
  }
)

names(forecast_coverage_plots) <-
  paste0("T", sample_sizes)


# Figure 3: paired Modal coverage differences ----

paired_coverage_plots <- lapply(
  sample_sizes,
  function(current_sample_size) {
    current_plot_data <-
      paired_coverage_plot_data |>
      dplyr::filter(
        sample_size == current_sample_size
      )
    
    ggplot2::ggplot(
      current_plot_data,
      ggplot2::aes(
        x = half_width,
        y = mean_coverage_difference,
        color = comparison,
        linetype = comparison,
        group = comparison
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.4
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mc_lower, ymax = mc_upper),
        width = 0.03,
        linewidth = 0.3,
        alpha = 0.60,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        linewidth = 0.8,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        size = 1.6,
        na.rm = TRUE
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(dgp),
        cols = ggplot2::vars(horizon),
        scales = "free_y"
      ) +
      ggplot2::scale_x_continuous(breaks = half_widths) +
      ggplot2::scale_y_continuous(
        labels = function(x) {
          paste0(round(100 * x, 1), " pp")
        }
      ) +
      ggplot2::labs(
        title = paste0(
          "Paired Modal-LP coverage differences: T = ",
          current_sample_size
        ),
        subtitle = paste(
          "Positive values favor Modal LP;",
          "the ranking is an empirical question"
        ),
        x = "Interval half-width",
        y = "Paired coverage difference",
        color = "Comparison",
        linetype = "Comparison",
        caption = paste(
          "Comparisons use pair-specific common-success events.",
          "Vertical bars show +/- 1.96 Monte Carlo SE."
        )
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
      )
  }
)

names(paired_coverage_plots) <-
  paste0("T", sample_sizes)


# Small DGP-3 diagnostic: disaster versus no-disaster paths ----

disaster_coverage_plots <- lapply(
  sample_sizes,
  function(current_sample_size) {
    current_plot_data <-
      disaster_coverage_plot_data |>
      dplyr::filter(
        sample_size == current_sample_size
      )
    
    ggplot2::ggplot(
      current_plot_data,
      ggplot2::aes(
        x = half_width,
        y = mean_coverage,
        color = estimator,
        linetype = estimator,
        group = estimator
      )
    ) +
      ggplot2::geom_line(
        linewidth = 0.8,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        size = 1.6,
        na.rm = TRUE
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(disaster_path),
        cols = ggplot2::vars(horizon)
      ) +
      ggplot2::scale_x_continuous(breaks = half_widths) +
      ggplot2::scale_y_continuous(
        labels = function(x) {
          paste0(round(100 * x), "%")
        }
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste0(
          "Disaster-path coverage diagnostic: T = ",
          current_sample_size
        ),
        subtitle =
          "See disaster_coverage_summary for event counts",
        x = "Interval half-width",
        y = "Conditional coverage",
        color = "Estimator",
        linetype = "Estimator"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
      )
  }
)

names(disaster_coverage_plots) <-
  paste0("T", sample_sizes)

forecast_plots <- list(
  average_issued_forecasts = average_forecast_plot,
  fixed_width_coverage = forecast_coverage_plots,
  paired_modal_coverage_differences =
    paired_coverage_plots,
  disaster_coverage_diagnostic =
    disaster_coverage_plots
)


# Final smoke checks before saving derived output ----

# 1. Coverage must be weakly increasing in half-width within every
# replication-level evaluation cell.

coverage_monotonicity_failures <-
  forecast_coverage_by_replication |>
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

disaster_monotonicity_failures <-
  disaster_coverage_by_replication |>
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

if (
  nrow(coverage_monotonicity_failures) > 0 ||
  nrow(disaster_monotonicity_failures) > 0
) {
  stop(
    paste(
      "Coverage is not weakly increasing in",
      "half-width in at least one cell."
    )
  )
}


# 2. Stored paired counts must equal independently constructed
# pair-specific common-success counts.

paired_count_check <-
  paired_coverage_by_replication |>
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
    by = c(
      "dgp",
      "sample_size",
      "replication",
      "horizon",
      "comparison"
    )
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
  stop(
    "The raw forecast table changed size during analysis."
  )
}

if (
  any(
    c(
      "half_width",
      "total_width",
      "coverage",
      "hit"
    ) %in% names(forecast_results_raw)
  )
) {
  stop(
    "Width-indexed columns were added to forecast_results_raw."
  )
}


# 4. All estimators must face the same realized outcome for each
# forecast event.

realized_outcome_check <-
  forecast_results_raw |>
  dplyr::group_by(
    dgp,
    sample_size,
    replication,
    forecast_origin,
    target_time,
    horizon
  ) |>
  dplyr::summarise(
    n_realized_values = dplyr::n_distinct(realized),
    n_estimators = dplyr::n_distinct(estimator),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_realized_values != 1 |
      n_estimators != length(estimator_names)
  )

if (nrow(realized_outcome_check) > 0) {
  stop(
    paste(
      "Realized outcomes or estimator counts disagree",
      "within at least one forecast event."
    )
  )
}


# 5. DGP 3 must use the redesigned defaults.

disaster_parameters <- dgp_parameters[["disaster"]]

valid_disaster_probability <-
  !is.null(disaster_parameters$disaster_probability) &&
  isTRUE(
    all.equal(
      as.numeric(
        disaster_parameters$disaster_probability
      ),
      0.05,
      tolerance = 1e-12
    )
  )

valid_disaster_size <-
  !is.null(disaster_parameters$disaster_size) &&
  isTRUE(
    all.equal(
      as.numeric(disaster_parameters$disaster_size),
      7.5,
      tolerance = 1e-12
    )
  )

if (!valid_disaster_probability || !valid_disaster_size) {
  stop(
    paste(
      "Stored DGP-3 parameters do not equal",
      "disaster_probability = 0.05 and",
      "disaster_size = 7.5."
    )
  )
}


# 6. Required outputs must be nonempty and width-indexed outputs must
# contain exactly the six requested half-widths.

required_output_objects <- list(
  average_forecast_by_replication =
    average_forecast_by_replication,
  average_forecast_summary = average_forecast_summary,
  forecast_coverage_by_replication =
    forecast_coverage_by_replication,
  forecast_coverage_summary = forecast_coverage_summary,
  paired_coverage_by_replication =
    paired_coverage_by_replication,
  paired_coverage_summary = paired_coverage_summary,
  disaster_coverage_by_replication =
    disaster_coverage_by_replication,
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
    setequal(
      sort(unique(x$half_width)),
      half_widths
    )
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

if (
  !setequal(
    unique(paired_coverage_summary$comparison),
    c("Modal - Mean", "Modal - Median")
  )
) {
  stop(
    "The paired output does not contain the requested comparisons."
  )
}


# Save derived analysis objects ----

save_run_object(
  average_forecast_by_replication,
  "average_forecast_by_replication"
)

save_run_object(
  average_forecast_summary,
  "average_forecast_summary"
)

save_run_object(
  forecast_coverage_by_replication,
  "forecast_coverage_by_replication"
)

save_run_object(
  forecast_coverage_summary,
  "forecast_coverage_summary"
)

save_run_object(
  paired_coverage_by_replication,
  "paired_coverage_by_replication"
)

save_run_object(
  paired_coverage_summary,
  "paired_coverage_summary"
)

save_run_object(
  disaster_coverage_by_replication,
  "disaster_coverage_by_replication"
)

save_run_object(
  disaster_coverage_summary,
  "disaster_coverage_summary"
)

save_run_object(
  forecast_scores_by_replication,
  "forecast_scores_by_replication"
)

save_run_object(
  forecast_summary,
  "forecast_summary"
)

save_run_object(
  forecast_plots,
  "forecast_plots"
)

message(
  "Forecast analysis complete. Run label: ",
  run_label
)