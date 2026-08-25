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


# Preallocate Modal LP results ----

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


# Preallocate Mean LP, Median LP, and VAR results ----

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


# Response Monte Carlo ----

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


# Inspect results ----

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