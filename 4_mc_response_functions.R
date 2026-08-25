# full LMP monte carlo simulations
# multiple T's, horizons, and all the dgps
# bias/variance/RMSE of LMP/MIRF estimator for all dgps
# MRF plots

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

max_sample_size <- max(sample_sizes)


# Preallocate raw results ----

n_results <- length(dgp_names) *
  length(sample_sizes) *
  n_replications *
  length(horizons)

mc_results <- data.frame(
  dgp = rep(NA_character_, n_results),
  sample_size = rep(NA_integer_, n_results),
  replication = rep(NA_integer_, n_results),
  horizon = rep(NA_integer_, n_results),
  estimate = rep(NA_real_, n_results),
  truth = rep(NA_real_, n_results),
  error = rep(NA_real_, n_results),
  converged = rep(NA, n_results),
  iterations = rep(NA_integer_, n_results),
  selected_start = rep(NA_character_, n_results),
  n_converged_starts = rep(NA_integer_, n_results),
  failure = rep(NA_character_, n_results)
)

result_row <- 1


# Baseline Monte Carlo ----

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
    
    # Select the appropriate modal population response
    truth_column <- if (dgp_name == "downside_risk") {
      "mode_response"
    } else {
      "response"
    }
    
    dgp_truth <- simulated_dgp$true_response[
      [
        truth_column
      ]
    ][
      match(
        horizons,
        simulated_dgp$true_response$horizon
      )
    ]
    
    if (anyNA(dgp_truth)) {
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
      
      # Correct conditioning vector
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
        
        # The true modal function is mildly nonlinear in x.
        # Here the linear LMP is a deliberate substantive approximation,
        # not the exactly specified efficiency benchmark supplied by DGP 1.
        
      } else {
        
        z <- sample_data$y
        reference_z <- 0
      }
      
      for (h_index in seq_along(horizons)) {
        
        h <- horizons[h_index]
        truth_h <- dgp_truth[h_index]
        
        attempt <- tryCatch(
          {
            fitted_lp <- fit_modal_lp(
              y = y,
              x = x,
              z = z,
              horizons = h
            )
            
            fitted_response <- lp_response(
              fitted_lp = fitted_lp,
              x = 0,
              z = reference_z,
              delta = 1
            )
            
            list(
              fit = fitted_lp,
              estimate = fitted_response$response[1]
            )
          },
          error = function(e) e
        )
        
        if (inherits(attempt, "error")) {
          
          mc_results[result_row, ] <- list(
            dgp_name,
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
            conditionMessage(attempt)
          )
          
        } else {
          
          estimate_h <- attempt$estimate
          error_h <- estimate_h - truth_h
          
          mc_results[result_row, ] <- list(
            dgp_name,
            sample_size,
            replication,
            h,
            estimate_h,
            truth_h,
            error_h,
            attempt$fit$converged[1],
            attempt$fit$iterations[1],
            attempt$fit$selected_start[1],
            attempt$fit$n_converged_starts[1],
            NA_character_
          )
        }
        
        result_row <- result_row + 1
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


# Monte Carlo summary ----

result_groups <- split(
  mc_results,
  interaction(
    mc_results$dgp,
    mc_results$sample_size,
    mc_results$horizon,
    drop = TRUE
  )
)

mc_summary <- do.call(
  rbind,
  lapply(
    result_groups,
    function(group) {
      
      successful <- group$converged &
        is.finite(group$estimate)
      
      n_success <- sum(successful)
      
      data.frame(
        dgp = group$dgp[1],
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
        success_rate = n_success / nrow(group)
      )
    }
  )
)

rownames(mc_summary) <- NULL

mc_summary <- mc_summary[
  order(
    match(mc_summary$dgp, dgp_names),
    mc_summary$sample_size,
    mc_summary$horizon
  ),
]

# Inspect results ----

head(mc_results)
head(mc_summary)

mc_summary[
  mc_summary$success_rate < 1,
  c(
    "dgp",
    "sample_size",
    "horizon",
    "n_success",
    "success_rate"
  )
]