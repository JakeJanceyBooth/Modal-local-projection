# Quick demonstration of the estimators and Gaussian benchmark DGP.

source(file.path("R", "estimators.R"))git status
source(file.path("R", "dgps.R"))

demo_seed <- 12345
sample_size <- 250L
horizons <- c(1L, 5L, 10L)

set.seed(demo_seed)

# Simulate a Gaussian DGP. Mean, median, and modal responses share the
# same population target in this design.
simulation <- simulate_gaussian_dgp(
  T = sample_size,
  horizon = max(horizons)
)

data <- simulation$data

# LPs condition on the current outcome y_t.
controls <- data$y

mean_fit <- fit_mean_lp(
  y = data$y,
  x = data$x,
  z = controls,
  horizons = horizons
)

median_fit <- fit_median_lp(
  y = data$y,
  x = data$x,
  z = controls,
  horizons = horizons
)

modal_fit <- fit_modal_lp(
  y = data$y,
  x = data$x,
  z = controls,
  horizons = horizons,
  bw_constant = 2.4,
  start_quantiles = 0.5
)

# VAR(1) in the observed shock and outcome.
var_fit <- fit_var(
  data = data[, c("x", "y")],
  lags = 1L
)

var_irf <- var_response(
  fitted_var = var_fit,
  shock = "x",
  response = "y",
  horizon = max(horizons)
)

benchmark <- simulation$true_response$response[
  match(horizons, simulation$true_response$horizon)
]

# "estimate" denotes the estimated horizon-h response 
# to a \delta = 1 shock in x_t

demo_results <- rbind(
  data.frame(
    horizon = horizons,
    estimator = "Mean LP",
    estimate = mean_fit$beta,
    benchmark = benchmark,
    bandwidth = NA_real_,
    iterations = NA_integer_,
    selected_start = NA_character_,
    n_converged_starts = NA_integer_
  ),
  data.frame(
    horizon = horizons,
    estimator = "Median LP",
    estimate = median_fit$beta,
    benchmark = benchmark,
    bandwidth = NA_real_,
    iterations = NA_integer_,
    selected_start = NA_character_,
    n_converged_starts = NA_integer_
  ),
  data.frame(
    horizon = horizons,
    estimator = "Modal LP",
    estimate = modal_fit$beta,
    benchmark = benchmark,
    bandwidth = modal_fit$bandwidth,
    iterations = modal_fit$iterations,
    selected_start = modal_fit$selected_start,
    n_converged_starts = modal_fit$n_converged_starts
  ),
  data.frame(
    horizon = horizons,
    estimator = "VAR",
    estimate = var_irf$response[
      match(horizons, var_irf$horizon)
    ],
    benchmark = benchmark,
    bandwidth = NA_real_,
    iterations = NA_integer_,
    selected_start = NA_character_,
    n_converged_starts = NA_integer_
  )
)

demo_results <- demo_results[
  order(
    demo_results$horizon,
    match(
      demo_results$estimator,
      c("Mean LP", "Median LP", "Modal LP", "VAR")
    )
  ),
]

rownames(demo_results) <- NULL

print(demo_results, digits = 3)