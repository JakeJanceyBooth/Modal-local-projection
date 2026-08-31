# All DGPs return T post-burn observations and the population
# responses for h = 1, ..., horizon.
# Gaussian AR(1) DGP.
# y_t = intercept + rho * y_{t-1} + beta * x_{t-1} + epsilon_t.
simulate_gaussian_dgp <- function(T, horizon,
                                  intercept = 0,
                                  rho = 0.8,
                                  beta = 1,
                                  sigma_x = 1,
                                  sigma_y = 1,
                                  burn = 200) {
  if (horizon < 1) {
    stop("horizon must be at least 1.")
  }

  n_total <- T + burn

  # x_t is the observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)
  epsilon_y <- rnorm(n_total, mean = 0, sd = sigma_y)
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)

  # y_t = rho y_{t-1} + beta x_{t-1} + epsilon_t
  for (t in 2:n_total) {
    y[t] <- intercept + rho * y[t - 1] + beta * x[t - 1] + epsilon_y[t]
  }

  # Remove the burn-in observations
  keep <- (burn + 1):n_total
  data <- data.frame(
    x = x[keep],
    y = y[keep],
    y_lag = y[keep - 1],
    x_lag = x[keep - 1]
  )

  # R_0 = 0 and R_h = beta * rho^(h - 1) for h >= 1
  h <- seq_len(horizon)
  response <- beta * rho^(h - 1)
  true_response <- data.frame(horizon = h, response = response)

  list(
    data = data,
    true_response = true_response,
    parameters = list(
      intercept = intercept,
      rho = rho,
      beta = beta,
      sigma_x = sigma_x,
      sigma_y = sigma_y
    )
  )
}

# Additive-skewness DGP.
# same as DGP 1 except epsilon_t follows skewed distribution
simulate_skewed_dgp <- function(T, horizon,
                                intercept = 0,
                                rho = 0.8,
                                beta = 1,
                                sigma_x = 1,
                                sigma_y = 1,
                                sdlog = 0.5,
                                burn = 200) {
  if (horizon < 1) {
    stop("horizon must be at least 1.")
  }

  n_total <- T + burn

  # x_t is the observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)

  # Shift the lognormal innovation so its mode is zero,
  # then scale it so its variance is sigma_y^2
  lognormal_draw <- rlnorm(n_total, meanlog = 0, sdlog = sdlog)
  lognormal_mode <- exp(-sdlog^2)
  lognormal_sd <- sqrt((exp(sdlog^2) - 1) * exp(sdlog^2))
  epsilon_y <- sigma_y * (lognormal_draw - lognormal_mode) / lognormal_sd
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)

  for (t in 2:n_total) {
    y[t] <- intercept + rho * y[t - 1] + beta * x[t - 1] + epsilon_y[t]
  }

  # Remove the burn-in observations
  keep <- (burn + 1):n_total
  data <- data.frame(
    x = x[keep],
    y = y[keep],
    y_lag = y[keep - 1],
    x_lag = x[keep - 1]
  )

  # The response slope is unchanged by additive skewness
  h <- seq_len(horizon)
  response <- beta * rho^(h - 1)
  true_response <- data.frame(horizon = h, response = response)

  list(
    data = data,
    true_response = true_response,
    parameters = list(
      intercept = intercept,
      rho = rho,
      beta = beta,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      sdlog = sdlog
    )
  )
}

# Rare-disaster mixture DGP.
# same as DGP 1 except epsilon_t follows a "rare disaster" distribution
simulate_disaster_dgp <- function(T, horizon,
                                  intercept = 0,
                                  rho = 0.8,
                                  beta = 1,
                                  sigma_x = 1,
                                  sigma_y = 1,
                                  disaster_probability = 0.05,
                                  disaster_size = 7.5,
                                  burn = 200) {
  if (horizon < 1) {
    stop("horizon must be at least 1.")
  }

  n_total <- T + burn

  # x_t is the observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)

  # D_t equals one during a rare-disaster realization
  disaster <- rbinom(n_total, size = 1, prob = disaster_probability)
  epsilon_y <- sigma_y * (rnorm(n_total) - disaster_size * disaster)
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)

  for (t in 2:n_total) {
    y[t] <- intercept + rho * y[t - 1] + beta * x[t - 1] + epsilon_y[t]
  }

  # Remove the burn-in observations
  keep <- (burn + 1):n_total
  data <- data.frame(
    x = x[keep],
    y = y[keep],
    y_lag = y[keep - 1],
    x_lag = x[keep - 1],
    disaster = disaster[keep]
  )

  # Additive disasters do not change the response slope
  h <- seq_len(horizon)
  response <- beta * rho^(h - 1)
  true_response <- data.frame(horizon = h, response = response)

  list(
    data = data,
    true_response = true_response,
    parameters = list(
      intercept = intercept,
      rho = rho,
      beta = beta,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      disaster_probability = disaster_probability,
      disaster_size = disaster_size
    )
  )
}

# Rich-state AR(2) DGP.
# y_t = intercept + rho1*y_{t-1} + rho2*y_{t-2}
#       + beta*x_{t-1} + gamma_s*s_{t-1} + epsilon_t.
# The state s_t follows an AR(1).
simulate_rich_state_dgp <- function(T, horizon,
                                    intercept = 0,
                                    rho1 = 0.5,
                                    rho2 = 0.2,
                                    beta = 1,
                                    phi_s = 0.8,
                                    gamma_s = 0.4,
                                    sigma_x = 1,
                                    sigma_y = 1,
                                    sigma_s = 1,
                                    sdlog = 0.5,
                                    burn = 200) {
  if (horizon < 1) {
    stop("horizon must be at least 1.")
  }

  n_total <- T + burn

  # Observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)

  # Persistent state variable
  epsilon_s <- rnorm(n_total, mean = 0, sd = sigma_s)
  s <- numeric(n_total)

  for (t in 2:n_total) {
    s[t] <- phi_s * s[t - 1] + epsilon_s[t]
  }

  # Shifted and rescaled lognormal innovation
  lognormal_draw <- rlnorm(n_total, meanlog = 0, sdlog = sdlog)
  lognormal_mode <- exp(-sdlog^2)
  lognormal_sd <- sqrt((exp(sdlog^2) - 1) * exp(sdlog^2))
  epsilon_y <- sigma_y * (lognormal_draw - lognormal_mode) / lognormal_sd

  # Outcome dynamics
  y <- numeric(n_total)
  y[1:2] <- intercept / (1 - rho1 - rho2)

  for (t in 3:n_total) {
    y[t] <- intercept +
      rho1 * y[t - 1] +
      rho2 * y[t - 2] +
      beta * x[t - 1] +
      gamma_s * s[t - 1] +
      epsilon_y[t]
  }

  # Remove burn-in and align predetermined variables
  keep <- (burn + 1):n_total
  data <- data.frame(
    x = x[keep],
    y = y[keep],
    s = s[keep],
    y_lag = y[keep - 1],
    y_lag2 = y[keep - 2],
    x_lag = x[keep - 1],
    s_lag = s[keep - 1]
  )

  # AR(2) response recursion
  h <- seq_len(horizon)
  response <- numeric(horizon)

  # R_1 = beta
  response[1] <- beta

  # R_2 = rho1 * R_1 + rho2 * R_0, where R_0 = 0
  if (horizon >= 2) {
    response[2] <- rho1 * beta
  }

  # R_h = rho1 R_{h-1} + rho2 R_{h-2}
  if (horizon >= 3) {
    for (j in 3:horizon) {
      response[j] <- rho1 * response[j - 1] + rho2 * response[j - 2]
    }
  }

  true_response <- data.frame(horizon = h, response = response)

  list(
    data = data,
    true_response = true_response,
    parameters = list(
      intercept = intercept,
      rho1 = rho1,
      rho2 = rho2,
      beta = beta,
      phi_s = phi_s,
      gamma_s = gamma_s,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      sigma_s = sigma_s,
      sdlog = sdlog
    )
  )
}

# Shock-dependent downside-risk DGP.
# x_t affects the persistent risk state r_{t+1}; r_t determines the
# probability of a downside realization in y_t. The shock therefore
# changes the shape of the horizon-h conditional distribution.
simulate_downside_risk_dgp <- function(T, horizon,
                                       phi = 0.8,
                                       kappa = 0.6,
                                       sigma_r = 0.3,
                                       sigma_x = 1,
                                       p_max = 0.30,
                                       alpha = -1,
                                       mu = 0,
                                       sigma_y = 1,
                                       disaster_size = 2.5,
                                       reference_r = 0,
                                       x_baseline = 0,
                                       delta = 1,
                                       burn = 200) {
  if (horizon < 1) {
    stop("horizon must be at least 1.")
  }

  n_total <- T + burn

  # Observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)

  # Persistent downside-risk state
  eta_r <- rnorm(n_total, mean = 0, sd = sigma_r)
  r <- numeric(n_total)

  for (t in 2:n_total) {
    r[t] <- phi * r[t - 1] + kappa * x[t - 1] + eta_r[t]
  }

  # The shock changes the probability, but not the locations,
  # of the normal and disaster components
  disaster_probability <- p_max * pnorm(alpha + r)
  disaster <- rbinom(n_total, size = 1, prob = disaster_probability)
  y <- mu - disaster_size * disaster + rnorm(n_total, mean = 0, sd = sigma_y)

  # Remove the burn-in observations
  keep <- (burn + 1):n_total
  data <- data.frame(
    x = x[keep],
    y = y[keep],
    r = r[keep],
    disaster_probability = disaster_probability[keep],
    disaster = disaster[keep]
  )

  # Conditional distribution of r_{t+h}
  h <- seq_len(horizon)
  A_h <- (1 - phi^(2 * h)) / (1 - phi^2)
  A_h_minus_1 <- (1 - phi^(2 * (h - 1))) / (1 - phi^2)
  risk_variance <- sigma_r^2 * A_h + kappa^2 * sigma_x^2 * A_h_minus_1
  risk_scale <- sqrt(1 + risk_variance)
  risk_mean_baseline <- phi^h * reference_r + kappa * phi^(h - 1) * x_baseline
  risk_mean_shocked <- phi^h * reference_r +
    kappa * phi^(h - 1) * (x_baseline + delta)

  # Probit-normal identity gives the mixture weights analytically
  probability_baseline <- p_max * pnorm(
    (alpha + risk_mean_baseline) /
      risk_scale
  )

  probability_shocked <- p_max * pnorm((alpha + risk_mean_shocked) / risk_scale)
  mean_baseline <- mu - disaster_size * probability_baseline
  mean_shocked <- mu - disaster_size * probability_shocked

  # Deterministic Gaussian-mixture median
  mixture_median <- function(probability) {
    uniroot(
      function(q) {
        (1 - probability) *
          pnorm((q - mu) / sigma_y) +
          probability *
          pnorm((q - (mu - disaster_size)) / sigma_y) -
          0.5
      },
      lower = mu - disaster_size - 8 * sigma_y,
      upper = mu + 8 * sigma_y,
      tol = 1e-12
    )$root
  }

  median_baseline <- vapply(probability_baseline, mixture_median, numeric(1))
  median_shocked <- vapply(probability_shocked, mixture_median, numeric(1))

  # Deterministic Gaussian-mixture mode. With the baseline
  # calibration, the density is unimodal and this is its unique mode,
  # obtained by finding the unique root of the derivative of the PDF
  mixture_mode <- function(probability) {
    uniroot(
      function(q) {
        (1 - probability) *
          (q - mu) *
          dnorm((q - mu) / sigma_y) +
          probability *
          (q - (mu - disaster_size)) *
          dnorm((q - (mu - disaster_size)) / sigma_y)
      },
      lower = mu - disaster_size / 2,
      upper = mu,
      tol = 1e-12
    )$root
  }

  mode_baseline <- vapply(probability_baseline, mixture_mode, numeric(1))
  mode_shocked <- vapply(probability_shocked, mixture_mode, numeric(1))
  true_response <- data.frame(
    horizon = h,
    probability_baseline = probability_baseline,
    probability_shocked = probability_shocked,
    mean_baseline = mean_baseline,
    mean_shocked = mean_shocked,
    mean_response = mean_shocked - mean_baseline,
    median_baseline = median_baseline,
    median_shocked = median_shocked,
    median_response = median_shocked - median_baseline,
    mode_baseline = mode_baseline,
    mode_shocked = mode_shocked,
    mode_response = mode_shocked - mode_baseline
  )

  # Unlike DGPs 1--4, the true modal function is mildly nonlinear
  # in x. The linear LMP is therefore a substantive approximation,
  # not an exactly specified estimator-efficiency benchmark.
  list(
    data = data,
    true_response = true_response,
    parameters = list(
      phi = phi,
      kappa = kappa,
      sigma_r = sigma_r,
      sigma_x = sigma_x,
      p_max = p_max,
      alpha = alpha,
      mu = mu,
      sigma_y = sigma_y,
      disaster_size = disaster_size,
      reference_r = reference_r,
      x_baseline = x_baseline,
      delta = delta
    )
  )
}
