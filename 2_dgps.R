# simple one-horizon test DGP
# symmetric DGP (Gaussian)
# skewed/assymetric DGP (lognormal)
# rare-disaster-type DGP (Gaussian mixture)
# misspecification ()
# functions that give the true mean/mode/median/IRF/MIRF if known

# Gaussian common-target benchmark
simulate_gaussian_dgp <- function(T, horizon,
                                  intercept = 0,
                                  rho = 0.6,
                                  beta = 1,
                                  sigma_x = 1,
                                  sigma_y = 1,
                                  burn = 200) {
  
  n_total <- T + burn
  
  # x_t is the observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)
  epsilon_y <- rnorm(n_total, mean = 0, sd = sigma_y)
  
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)
  
  # y_t = rho y_{t-1} + beta x_{t-1} + epsilon_t
  for (t in 2:n_total) {
    y[t] <- intercept +
      rho * y[t - 1] +
      beta * x[t - 1] +
      epsilon_y[t]
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
  
  true_response <- data.frame(
    horizon = h,
    response = response
  )
  
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

# Additive-skewness benchmark
simulate_skewed_dgp <- function(T, horizon,
                                intercept = 0,
                                rho = 0.6,
                                beta = 1,
                                sigma_x = 1,
                                sigma_y = 1,
                                sdlog = 0.5,
                                burn = 200) {
  
  n_total <- T + burn
  
  # x_t is the observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)
  
  # Shift the lognormal innovation so its mode is zero,
  # then scale it so its variance is sigma_y^2
  lognormal_draw <- rlnorm(
    n_total,
    meanlog = 0,
    sdlog = sdlog
  )
  
  lognormal_mode <- exp(-sdlog^2)
  
  lognormal_sd <- sqrt(
    (exp(sdlog^2) - 1) * exp(sdlog^2)
  )
  
  epsilon_y <- sigma_y *
    (lognormal_draw - lognormal_mode) /
    lognormal_sd
  
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)
  
  for (t in 2:n_total) {
    y[t] <- intercept +
      rho * y[t - 1] +
      beta * x[t - 1] +
      epsilon_y[t]
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
  
  true_response <- data.frame(
    horizon = h,
    response = response
  )
  
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

# Rare-disaster mixture benchmark
simulate_disaster_dgp <- function(T, horizon,
                                  intercept = 0,
                                  rho = 0.6,
                                  beta = 1,
                                  sigma_x = 1,
                                  sigma_y = 1,
                                  disaster_probability = 0.01,
                                  disaster_size = 10,
                                  burn = 200) {
  
  n_total <- T + burn
  
  # x_t is the observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)
  
  # D_t equals one during a rare-disaster realization
  disaster <- rbinom(
    n_total,
    size = 1,
    prob = disaster_probability
  )
  
  epsilon_y <- sigma_y * (
    rnorm(n_total) -
      disaster_size * disaster
  )
  
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)
  
  for (t in 2:n_total) {
    y[t] <- intercept +
      rho * y[t - 1] +
      beta * x[t - 1] +
      epsilon_y[t]
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
  
  true_response <- data.frame(
    horizon = h,
    response = response
  )
  
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

# Rich dynamic-state benchmark
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
  lognormal_draw <- rlnorm(
    n_total,
    meanlog = 0,
    sdlog = sdlog
  )
  
  lognormal_mode <- exp(-sdlog^2)
  
  lognormal_sd <- sqrt(
    (exp(sdlog^2) - 1) * exp(sdlog^2)
  )
  
  epsilon_y <- sigma_y *
    (lognormal_draw - lognormal_mode) /
    lognormal_sd
  
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
      response[j] <-
        rho1 * response[j - 1] +
        rho2 * response[j - 2]
    }
  }
  
  true_response <- data.frame(
    horizon = h,
    response = response
  )
  
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

# Shock-dependent distributional shock_dependent_scale

simulate_shock_dependent_scale_dgp <- function(T,
                               intercept = 0,
                               rho = 0.6,
                               beta = 1,
                               sigma_x = 1,
                               sigma_y = 1,
                               sdlog = 0.5,
                               burn = 200) {
  
  n_total <- T + burn
  
  # Observed structural shock
  x <- rnorm(n_total, mean = 0, sd = sigma_x)
  
  # Shifted and rescaled lognormal innovation
  lognormal_draw <- rlnorm(
    n_total,
    meanlog = 0,
    sdlog = sdlog
  )
  
  lognormal_mode <- exp(-sdlog^2)
  
  lognormal_sd <- sqrt(
    (exp(sdlog^2) - 1) * exp(sdlog^2)
  )
  
  epsilon_y <- sigma_y *
    (lognormal_draw - lognormal_mode) /
    lognormal_sd
  
  y <- numeric(n_total)
  y[1] <- intercept / (1 - rho)
  
  for (t in 2:n_total) {
    y[t] <- intercept +
      rho * y[t - 1] +
      beta * x[t - 1] +
      exp(x[t - 1]) * epsilon_y[t]
  }
  
  # Remove burn-in and align lagged variables
  keep <- (burn + 1):n_total
  
  data <- data.frame(
    x = x[keep],
    y = y[keep],
    y_lag = y[keep - 1],
    x_lag = x[keep - 1]
  )
  
  list(
    data = data,
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

# Population responses under paired counterfactual paths

simulate_shock_dependent_scale_response <- function(S = 2000000, 
                                    horizon = 20, 
                                    y_initial = 0,
                                    x_baseline = 0,
                                    delta = 1,
                                    intercept = 0,
                                    rho = 0.6,
                                    beta = 1,
                                    sigma_x = 1,
                                    sigma_y = 1,
                                    sdlog = 0.5,
                                    mode_function = NULL) {
  
  lognormal_mode <- exp(-sdlog^2)
  
  lognormal_sd <- sqrt(
    (exp(sdlog^2) - 1) * exp(sdlog^2)
  )
  
  baseline_paths <- matrix(
    NA_real_,
    nrow = S,
    ncol = horizon + 1
  )
  
  shocked_paths <- matrix(
    NA_real_,
    nrow = S,
    ncol = horizon + 1
  )
  
  colnames(baseline_paths) <- paste0("h", 0:horizon)
  colnames(shocked_paths) <- paste0("h", 0:horizon)
  
  y_baseline <- rep(y_initial, S)
  y_shocked <- rep(y_initial, S)
  
  baseline_paths[, 1] <- y_baseline
  shocked_paths[, 1] <- y_shocked
  
  if (horizon >= 1) {
    
    for (h in 1:horizon) {
      
      # Only the time-t shock differs across experiments
      if (h == 1) {
        x_baseline_h <- x_baseline
        x_shocked_h <- x_baseline + delta
      } else {
        x_future <- rnorm(
          S,
          mean = 0,
          sd = sigma_x
        )
        
        x_baseline_h <- x_future
        x_shocked_h <- x_future
      }
      
      # Same innovation draw in both counterfactual paths
      u <- sigma_y * (
        rlnorm(
          S,
          meanlog = 0,
          sdlog = sdlog
        ) -
          lognormal_mode
      ) / lognormal_sd
      
      y_baseline <- intercept +
        rho * y_baseline +
        beta * x_baseline_h +
        exp(x_baseline_h) * u
      
      y_shocked <- intercept +
        rho * y_shocked +
        beta * x_shocked_h +
        exp(x_shocked_h) * u
      
      baseline_paths[, h + 1] <- y_baseline
      shocked_paths[, h + 1] <- y_shocked
    }
  }
  
  mean_baseline <- colMeans(baseline_paths)
  mean_shocked <- colMeans(shocked_paths)
  
  median_baseline <- apply(
    baseline_paths,
    2,
    median
  )
  
  median_shocked <- apply(
    shocked_paths,
    2,
    median
  )
  
  # The initial state has a known degenerate mode.
  # Remaining modes are filled once a mode function is supplied.
  mode_baseline <- rep(NA_real_, horizon + 1)
  mode_shocked <- rep(NA_real_, horizon + 1)
  
  mode_baseline[1] <- y_initial
  mode_shocked[1] <- y_initial
  
  if (!is.null(mode_function) && horizon >= 1) {
    
    columns <- 2:(horizon + 1)
    
    mode_baseline[columns] <- apply(
      baseline_paths[, columns, drop = FALSE],
      2,
      mode_function
    )
    
    mode_shocked[columns] <- apply(
      shocked_paths[, columns, drop = FALSE],
      2,
      mode_function
    )
  }
  
  # Exact mean-response benchmark
  u_mean <- sigma_y *
    (exp(sdlog^2 / 2) - exp(-sdlog^2)) /
    lognormal_sd
  
  analytic_mean_response <- numeric(horizon + 1)
  
  if (horizon >= 1) {
    analytic_mean_response[2:(horizon + 1)] <-
      rho^(0:(horizon - 1)) * (
        beta * delta +
          u_mean * (
            exp(x_baseline + delta) -
              exp(x_baseline)
          )
      )
  }
  
  response <- data.frame(
    horizon = 0:horizon,
    mean_baseline = mean_baseline,
    mean_shocked = mean_shocked,
    mean_response = mean_shocked - mean_baseline,
    mean_response_analytic = analytic_mean_response,
    median_baseline = median_baseline,
    median_shocked = median_shocked,
    median_response = median_shocked - median_baseline,
    mode_baseline = mode_baseline,
    mode_shocked = mode_shocked,
    mode_response = mode_shocked - mode_baseline
  )
  
  list(
    baseline_paths = baseline_paths,
    shocked_paths = shocked_paths,
    response = response,
    parameters = list(
      S = S,
      horizon = horizon,
      y_initial = y_initial,
      x_baseline = x_baseline,
      delta = delta,
      intercept = intercept,
      rho = rho,
      beta = beta,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      sdlog = sdlog
    )
  )
}


# Population mode from a Gaussian KDE

population_mode <- function(x, bw = 1.5 * bw.nrd0(x), n = 8192) {
  
  limits <- quantile(
    x,
    probs = c(0.001, 0.999),
    names = FALSE
  )
  
  density_estimate <- density(
    x,
    kernel = "gaussian",
    bw = bw,
    n = n,
    from = limits[1],
    to = limits[2]
  )
  
  density_estimate$x[
    which.max(density_estimate$y)
  ]
}

