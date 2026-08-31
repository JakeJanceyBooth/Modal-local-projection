# all the tests used to determine the proper bandwidth for approximating
# the modes of DGP 5. We find b = 1.5 bw.nrd0(x) is optimal, and we should
# use S = 2,000,000

# One-step population-mode bandwidth calibration function ----

source("2_dgps.R")

calibrate_one_step_population_mode <- function(
    S_values,
    y_initial,
    x_baseline = 0,
    delta = 1,
    intercept = 0,
    rho = 0.6,
    beta = 1,
    sigma_y = 1,
    sdlog = 0.5,
    nrd_multipliers = c(0.5, 0.75, 1, 1.25, 1.5, 2),
    mode_rate_constants = c(0.5, 0.75, 1, 1.25, 1.5, 2),
    n = 8192) {
  
  S_values <- sort(unique(S_values))
  max_S <- max(S_values)
  
  lognormal_mode <- exp(-sdlog^2)
  
  lognormal_sd <- sqrt(
    (exp(sdlog^2) - 1) * exp(sdlog^2)
  )
  
  u <- sigma_y * (
    rlnorm(
      max_S,
      meanlog = 0,
      sdlog = sdlog
    ) -
      lognormal_mode
  ) / lognormal_sd
  
  baseline_population <- intercept +
    rho * y_initial +
    beta * x_baseline +
    exp(x_baseline) * u
  
  shocked_population <- intercept +
    rho * y_initial +
    beta * (x_baseline + delta) +
    exp(x_baseline + delta) * u
  
  exact_baseline_mode <- intercept +
    rho * y_initial +
    beta * x_baseline
  
  exact_shocked_mode <- intercept +
    rho * y_initial +
    beta * (x_baseline + delta)
  
  exact_response <- beta * delta
  
  results <- vector(
    "list",
    length(S_values)
  )
  
  for (j in seq_along(S_values)) {
    
    S <- S_values[j]
    keep <- seq_len(S)
    
    baseline <- baseline_population[keep]
    shocked <- shocked_population[keep]
    
    baseline_limits <- quantile(
      baseline,
      c(0.001, 0.999),
      names = FALSE
    )
    
    shocked_limits <- quantile(
      shocked,
      c(0.001, 0.999),
      names = FALSE
    )
    
    nrd_baseline <- bw.nrd0(baseline)
    nrd_shocked <- bw.nrd0(shocked)
    
    scale_baseline <- IQR(baseline) / 1.349
    scale_shocked <- IQR(shocked) / 1.349
    
    bandwidths <- rbind(
      data.frame(
        bandwidth_family = "nrd0_multiple",
        tuning_constant = nrd_multipliers,
        baseline_bandwidth =
          nrd_multipliers * nrd_baseline,
        shocked_bandwidth =
          nrd_multipliers * nrd_shocked
      ),
      data.frame(
        bandwidth_family = "mode_rate",
        tuning_constant = mode_rate_constants,
        baseline_bandwidth =
          mode_rate_constants *
          scale_baseline *
          S^(-1 / 7),
        shocked_bandwidth =
          mode_rate_constants *
          scale_shocked *
          S^(-1 / 7)
      )
    )
    
    estimated_baseline_mode <- vapply(
      seq_len(nrow(bandwidths)),
      function(i) {
        population_mode(
          baseline,
          bw = bandwidths$baseline_bandwidth[i],
          n = n
        )
      },
      numeric(1)
    )
    
    estimated_shocked_mode <- vapply(
      seq_len(nrow(bandwidths)),
      function(i) {
        population_mode(
          shocked,
          bw = bandwidths$shocked_bandwidth[i],
          n = n
        )
      },
      numeric(1)
    )
    
    estimated_response <-
      estimated_shocked_mode -
      estimated_baseline_mode
    
    baseline_grid_position <- (
      estimated_baseline_mode -
        baseline_limits[1]
    ) / diff(baseline_limits)
    
    shocked_grid_position <- (
      estimated_shocked_mode -
        shocked_limits[1]
    ) / diff(shocked_limits)
    
    results[[j]] <- data.frame(
      S = S,
      bandwidth_family =
        bandwidths$bandwidth_family,
      tuning_constant =
        bandwidths$tuning_constant,
      baseline_bandwidth =
        bandwidths$baseline_bandwidth,
      estimated_baseline_mode =
        estimated_baseline_mode,
      exact_baseline_mode =
        exact_baseline_mode,
      baseline_absolute_error =
        abs(
          estimated_baseline_mode -
            exact_baseline_mode
        ),
      baseline_boundary_fraction =
        pmin(
          baseline_grid_position,
          1 - baseline_grid_position
        ),
      shocked_bandwidth =
        bandwidths$shocked_bandwidth,
      estimated_shocked_mode =
        estimated_shocked_mode,
      exact_shocked_mode =
        exact_shocked_mode,
      shocked_absolute_error =
        abs(
          estimated_shocked_mode -
            exact_shocked_mode
        ),
      shocked_boundary_fraction =
        pmin(
          shocked_grid_position,
          1 - shocked_grid_position
        ),
      estimated_response =
        estimated_response,
      exact_response =
        exact_response,
      response_absolute_error =
        abs(
          estimated_response -
            exact_response
        )
    )
  }
  
  output <- do.call(
    rbind,
    results
  )
  
  rownames(output) <- NULL
  
  output
}

# tests to find best bandwidth for population mode estimation for DGP5 ----

calibration_seeds <- 1:5

calibration_runs <- vector(
  "list",
  length(calibration_seeds)
)

for (i in seq_along(calibration_seeds)) {
  
  set.seed(calibration_seeds[i])
  
  calibration_runs[[i]] <-
    calibrate_one_step_population_mode(
      S_values = c(
        50000,
        100000,
        250000,
        500000
      ),
      y_initial = 0,
      x_baseline = 0,
      delta = 1
    )
  
  calibration_runs[[i]]$seed <-
    calibration_seeds[i]
}

calibration_results <- do.call(
  rbind,
  calibration_runs
)

calibration_errors <- aggregate(
  cbind(
    baseline_absolute_error,
    shocked_absolute_error,
    response_absolute_error
  ) ~
    bandwidth_family +
    tuning_constant,
  data = calibration_results,
  FUN = function(x) {
    c(
      mean = mean(x),
      maximum = max(x)
    )
  }
)

calibration_summary <- data.frame(
  bandwidth_family =
    calibration_errors$bandwidth_family,
  tuning_constant =
    calibration_errors$tuning_constant,
  baseline_mean_error =
    calibration_errors$
    baseline_absolute_error[, "mean"],
  baseline_maximum_error =
    calibration_errors$
    baseline_absolute_error[, "maximum"],
  shocked_mean_error =
    calibration_errors$
    shocked_absolute_error[, "mean"],
  shocked_maximum_error =
    calibration_errors$
    shocked_absolute_error[, "maximum"],
  response_mean_error =
    calibration_errors$
    response_absolute_error[, "mean"],
  response_maximum_error =
    calibration_errors$
    response_absolute_error[, "maximum"]
)

calibration_summary <- calibration_summary[
  order(
    calibration_summary$response_mean_error,
    calibration_summary$response_maximum_error
  ),
]

rownames(calibration_summary) <- NULL

calibration_minimum_boundary_fraction <- min(
  calibration_results$baseline_boundary_fraction,
  calibration_results$shocked_boundary_fraction
)

calibration_summary
calibration_minimum_boundary_fraction



#stability tests with b = 1.25nrd ----


stability_S_values <- c(
  500000,
  1000000
)

stability_adjustments <- c(
  0.75,
  1,
  1.25
)

stability_horizon <- 20
selected_bw_multiplier <- 1.25

set.seed(123)

stability_paths <- simulate_shock_response(
  S = max(stability_S_values),
  horizon = stability_horizon,
  y_initial = 0,
  x_baseline = 0,
  delta = 1,
  mode_function = NULL
)

stability_results <- vector(
  "list",
  length(stability_S_values) *
    length(stability_adjustments) *
    stability_horizon
)

row <- 1

for (S in stability_S_values) {
  
  keep <- seq_len(S)
  
  for (adjustment in stability_adjustments) {
    
    for (h in 1:stability_horizon) {
      
      baseline <- stability_paths$
        baseline_paths[keep, h + 1]
      
      shocked <- stability_paths$
        shocked_paths[keep, h + 1]
      
      baseline_bw <- adjustment *
        selected_bw_multiplier *
        bw.nrd0(baseline)
      
      shocked_bw <- adjustment *
        selected_bw_multiplier *
        bw.nrd0(shocked)
      
      baseline_mode <- population_mode(
        baseline,
        bw = baseline_bw
      )
      
      shocked_mode <- population_mode(
        shocked,
        bw = shocked_bw
      )
      
      baseline_limits <- quantile(
        baseline,
        c(0.001, 0.999),
        names = FALSE
      )
      
      shocked_limits <- quantile(
        shocked,
        c(0.001, 0.999),
        names = FALSE
      )
      
      baseline_position <- (
        baseline_mode -
          baseline_limits[1]
      ) / diff(baseline_limits)
      
      shocked_position <- (
        shocked_mode -
          shocked_limits[1]
      ) / diff(shocked_limits)
      
      stability_results[[row]] <- data.frame(
        S = S,
        horizon = h,
        bandwidth_adjustment = adjustment,
        baseline_bandwidth = baseline_bw,
        shocked_bandwidth = shocked_bw,
        baseline_mode = baseline_mode,
        shocked_mode = shocked_mode,
        modal_response =
          shocked_mode - baseline_mode,
        boundary_fraction = min(
          baseline_position,
          1 - baseline_position,
          shocked_position,
          1 - shocked_position
        )
      )
      
      row <- row + 1
    }
  }
}

stability_results <- do.call(
  rbind,
  stability_results
)

# summary of tests for b = 1.25 nrd ---- 

central_curve <- subset(
  stability_results,
  bandwidth_adjustment == 1,
  select = c(
    S,
    horizon,
    modal_response
  )
)

lower_curve <- subset(
  stability_results,
  bandwidth_adjustment == 0.75,
  select = c(
    S,
    horizon,
    modal_response
  )
)

upper_curve <- subset(
  stability_results,
  bandwidth_adjustment == 1.25,
  select = c(
    S,
    horizon,
    modal_response
  )
)

names(central_curve)[3] <- "central_response"
names(lower_curve)[3] <- "lower_response"
names(upper_curve)[3] <- "upper_response"

bandwidth_comparison <- merge(
  central_curve,
  lower_curve,
  by = c("S", "horizon")
)

bandwidth_comparison <- merge(
  bandwidth_comparison,
  upper_curve,
  by = c("S", "horizon")
)

bandwidth_stability <- data.frame(
  comparison = c(
    "0.75 versus 1",
    "1.25 versus 1"
  ),
  maximum_absolute_difference = c(
    max(abs(
      bandwidth_comparison$lower_response -
        bandwidth_comparison$central_response
    )),
    max(abs(
      bandwidth_comparison$upper_response -
        bandwidth_comparison$central_response
    ))
  )
)

small_S_curve <- subset(
  central_curve,
  S == min(stability_S_values),
  select = c(horizon, central_response)
)

large_S_curve <- subset(
  central_curve,
  S == max(stability_S_values),
  select = c(horizon, central_response)
)

names(small_S_curve)[2] <- "small_S_response"
names(large_S_curve)[2] <- "large_S_response"

sample_size_comparison <- merge(
  small_S_curve,
  large_S_curve,
  by = "horizon"
)

sample_size_stability <- max(abs(
  sample_size_comparison$large_S_response -
    sample_size_comparison$small_S_response
))

stability_minimum_boundary_fraction <- min(
  stability_results$boundary_fraction
)

bandwidth_stability
sample_size_stability
stability_minimum_boundary_fraction

# Where is bandwidth sensitivity largest? ----

bandwidth_comparison$lower_difference <- abs(
  bandwidth_comparison$lower_response -
    bandwidth_comparison$central_response
)

bandwidth_comparison$upper_difference <- abs(
  bandwidth_comparison$upper_response -
    bandwidth_comparison$central_response
)

bandwidth_comparison[
  order(
    -pmax(
      bandwidth_comparison$lower_difference,
      bandwidth_comparison$upper_difference
    )
  ),
][1:10, ]


# Where is sample-size sensitivity largest? ----

sample_size_comparison$absolute_difference <- abs(
  sample_size_comparison$large_S_response -
    sample_size_comparison$small_S_response
)

sample_size_comparison[
  order(-sample_size_comparison$absolute_difference),
][1:10, ]

# finally testing across even larger sample sizes and the two
# best bandwidth rules so far ----

final_S_values <- c(
  500000,
  1000000,
  2000000
)

final_bw_multipliers <- c(
  1.25,
  1.5
)

set.seed(123)

final_paths <- simulate_shock_response(
  S = max(final_S_values),
  horizon = 5,
  y_initial = 0,
  x_baseline = 0,
  delta = 1,
  mode_function = NULL
)

final_results <- list()
row <- 1

for (S in final_S_values) {
  
  keep <- seq_len(S)
  
  for (bw_multiplier in final_bw_multipliers) {
    
    for (h in 1:5) {
      
      baseline <- final_paths$
        baseline_paths[keep, h + 1]
      
      shocked <- final_paths$
        shocked_paths[keep, h + 1]
      
      baseline_mode <- population_mode(
        baseline,
        bw = bw_multiplier * bw.nrd0(baseline)
      )
      
      shocked_mode <- population_mode(
        shocked,
        bw = bw_multiplier * bw.nrd0(shocked)
      )
      
      final_results[[row]] <- data.frame(
        S = S,
        horizon = h,
        bw_multiplier = bw_multiplier,
        modal_response =
          shocked_mode - baseline_mode
      )
      
      row <- row + 1
    }
  }
}

final_results <- do.call(
  rbind,
  final_results
)

final_results
