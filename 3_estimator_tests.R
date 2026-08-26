# one simulated data set
# check that MEM is non-decreasing
# check multi-start MEM/stopping criterion works
# estimator gets kind of close to the actual LMP coefficient and MIRF

source("1_estimators.R")
source("2_dgps.R")

# Gaussian common mean-median-mode target ----

set.seed(12345)

sample_size <- 3000
test_horizon <- 8
horizons <- seq_len(test_horizon)

gaussian_dgp <- simulate_gaussian_dgp(
  T = sample_size,
  horizon = test_horizon
)

y <- gaussian_dgp$data$y
x <- gaussian_dgp$data$x

z <- gaussian_dgp$data$y


# Estimate the three LPs

mean_fit <- fit_mean_lp(
  y = y,
  x = x,
  z = z,
  horizons = horizons
)

median_fit <- fit_median_lp(
  y = y,
  x = x,
  z = z,
  horizons = horizons
)

modal_fit <- fit_modal_lp(
  y = y,
  x = x,
  z = z,
  horizons = horizons
)


# Estimate the correctly specified VAR(1)

var_fit <- fit_var(
  data = gaussian_dgp$data[, c("x", "y")],
  lags = 1
)

# Counterfactual responses

reference_z <- mean(z)

mean_response <- lp_response(
  fitted_lp = mean_fit,
  x = 0,
  z = reference_z,
  delta = 1
)

median_response <- lp_response(
  fitted_lp = median_fit,
  x = 0,
  z = reference_z,
  delta = 1
)

modal_response <- lp_response(
  fitted_lp = modal_fit,
  x = 0,
  z = reference_z,
  delta = 1
)

var_response_fit <- var_response(
  fitted_var = var_fit,
  shock = "x",
  response = "y",
  horizon = test_horizon,
  delta = 1
)

# Compare estimated and population responses

response_comparison <- data.frame(
  horizon = horizons,
  true = gaussian_dgp$true_response$response,
  mean_lp = mean_response$response,
  median_lp = median_response$response,
  modal_lp = modal_response$response,
  var = var_response_fit$response
)

response_comparison

# Response-function accuracy (MAE and RMSE) ----

response_accuracy <- data.frame(
  method = c(
    "Mean LP",
    "Median LP",
    "Modal LP",
    "VAR"
  ),
  MAE = c(
    mean(abs(
      response_comparison$mean_lp -
        response_comparison$true
    )),
    mean(abs(
      response_comparison$median_lp -
        response_comparison$true
    )),
    mean(abs(
      response_comparison$modal_lp -
        response_comparison$true
    )),
    mean(abs(
      response_comparison$var -
        response_comparison$true
    ))
  ),
  RMSE = c(
    sqrt(mean(
      (response_comparison$mean_lp -
         response_comparison$true)^2
    )),
    sqrt(mean(
      (response_comparison$median_lp -
         response_comparison$true)^2
    )),
    sqrt(mean(
      (response_comparison$modal_lp -
         response_comparison$true)^2
    )),
    sqrt(mean(
      (response_comparison$var -
         response_comparison$true)^2
    ))
  )
)

response_accuracy

# In-sample LP prediction accuracy (MAE and RMSE)----

lp_fits <- list(
  "Mean LP" = mean_fit,
  "Median LP" = median_fit,
  "Modal LP" = modal_fit
)

forecast_accuracy <- list()
row <- 1

for (method in names(lp_fits)) {
  
  fitted_lp <- lp_fits[[method]]
  
  for (h in horizons) {
    
    dat <- .prepare_lp_data(
      y = y,
      x = x,
      z = z,
      h = h
    )
    
    coefficients <- fitted_lp$coefficients[
      match(h, fitted_lp$horizons),
    ]
    
    prediction <- as.numeric(
      dat$D_h %*%
        coefficients[colnames(dat$D_h)]
    )
    
    forecast_error <- dat$y_h - prediction
    
    forecast_accuracy[[row]] <- data.frame(
      method = method,
      horizon = h,
      n_forecasts = length(forecast_error),
      MAE = mean(abs(forecast_error)),
      RMSE = sqrt(mean(forecast_error^2))
    )
    
    row <- row + 1
  }
}

forecast_accuracy <- do.call(
  rbind,
  forecast_accuracy
)

rownames(forecast_accuracy) <- NULL

forecast_accuracy

# ----
# We find modal LP performs the worst for response functions.
# All estimators have basically similar in-sample forecast performance.

# Quick plot of estimated response functions to see how they compare ----

response_plot_data <- data.frame(
  horizon = rep(
    response_comparison$horizon,
    times = 5
  ),
  method = rep(
    c(
      "True",
      "Mean LP",
      "Median LP",
      "Modal LP",
      "VAR"
    ),
    each = nrow(response_comparison)
  ),
  response = c(
    response_comparison$true,
    response_comparison$mean_lp,
    response_comparison$median_lp,
    response_comparison$modal_lp,
    response_comparison$var
  )
)

response_plot <- ggplot(
  response_plot_data,
  aes(
    x = horizon,
    y = response,
    color = method,
    linetype = method
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Horizon",
    y = "Response",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

response_plot
ggsave(
  "gaussian_response.pdf",
  plot = response_plot,
  width = 6.5,
  height = 4.5
)