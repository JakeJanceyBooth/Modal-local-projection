# all estimation functions
# LMP/MEM estimator with multi-start, stopping rule, and bandwidth calculation
# mean LP and IRFs
# median/quantile LP and IRFs
# VARs and IRFs

# estimators take arguments y, x, optionally z, and horizons
# estimators return horizons, n_h for each h, coefficients, beta
# then just have one generic LP IRF function. arguments: fitted model, delta

source("0_setup.R")

#construct horizon-h LP data

.prepare_lp_data <- function(y, x, z = NULL, h) {
  
  # input checks
  y <- as.numeric(y)
  x <- as.numeric(x)
  T <- length(y)
  
  if (length(x) != T) {
    stop("y and x must have the same length.")
  }
  
  if (length(h) != 1 || h < 0 || h >= T || h != as.integer(h)) {
    stop("h must be a non-negative integer smaller than T.")
  }
  
  if (is.null(z)) {
    
    z <- matrix(numeric(0), nrow = T, ncol = 0)
    
  } else {
    
    if (is.vector(z)) {
      z <- matrix(z, ncol = 1)
    } else {
      z <- as.matrix(z)
    }
    
    if (nrow(z) != T) {
      stop("z must have T rows.")
    }
    
    if (is.null(colnames(z))) {
      colnames(z) <- paste0("z", seq_len(ncol(z)))
    }
  }
  
  if (anyNA(y) || anyNA(x) || anyNA(z)) {
    stop("Missing values are not currently supported.")
  }
  
  n_h <- T - h
  idx <- seq_len(n_h)
  
  y_h <- y[idx + h]
  x_h <- x[idx]
  z_h <- z[idx, , drop = FALSE]
  
  D_h <- cbind(
    alpha = 1,
    beta = x_h
  )
  
  if (ncol(z_h) > 0) {
    colnames(z_h) <- paste0("gamma_", colnames(z_h))
    D_h <- cbind(D_h, z_h)
  }
  
  list(
    y_h = y_h,
    D_h = D_h,
    n_h = n_h
  )
}

# Mean Local Projection

fit_mean_lp <- function(y, x, z = NULL, horizons) {
  
  horizons <- sort(unique(horizons))
  
  if (length(horizons) == 0) {
    stop("At least one horizon must be supplied.")
  }
  
  fits <- lapply(horizons, function(h) {
    
    dat <- .prepare_lp_data(
      y = y,
      x = x,
      z = z,
      h = h
    )
    
    fit_h <- lm.fit(
      x = dat$D_h,
      y = dat$y_h
    )
    
    if (fit_h$rank < ncol(dat$D_h)) {
      stop(
        paste0(
          "Design matrix is rank deficient at horizon h = ",
          h,
          "."
        )
      )
    }
    
    list(
      coefficients = setNames(
        as.numeric(fit_h$coefficients),
        colnames(dat$D_h)
      ),
      n_h = dat$n_h
    )
  })
  
  # Stack \hat\theta_h across horizons
  coefficients <- do.call(
    rbind,
    lapply(fits, function(fit) fit$coefficients)
  )
  
  rownames(coefficients) <- paste0("h_", horizons)
  
  # Extract \hat\beta_h in particular
  beta <- coefficients[, "beta"]
  
  n_h <- vapply(
    fits,
    function(fit) fit$n_h,
    numeric(1)
  )
  
  out <- list(
    method = "mean_lp",
    horizons = horizons,
    n_h = n_h,
    coefficients = coefficients,
    beta = beta
  )
  
  out
}

# Local projection impulse response

lp_response <- function(fitted_lp, delta = 1) {
  beta <- as.numeric(fitted_lp$beta)
  
  data.frame(
    horizon = fitted_lp$horizons,
    beta = beta,
    delta = delta,
    response = delta * beta
  )
}

# Median local projection

fit_median_lp <- function(y, x, z = NULL, horizons) {
  
  horizons <- sort(unique(horizons))
  
  if (length(horizons) == 0L) {
    stop("At least one horizon must be supplied.")
  }
  
  fits <- lapply(horizons, function(h) {
    
    dat <- .prepare_lp_data(
      y = y,
      x = x,
      z = z,
      h = h
    )
    
    # Check that theta_h is identified
    if (qr(dat$D_h)$rank < ncol(dat$D_h)) {
      stop(
        paste0(
          "Design matrix is rank deficient at horizon h = ",
          h,
          "."
        )
      )
    }
    
    # Horizon-h median regression
    fit_h <- quantreg::rq.fit(
      x = dat$D_h,
      y = dat$y_h,
      tau = 0.5,
      method = "br"
    )
    
    list(
      coefficients = setNames(
        as.numeric(fit_h$coefficients),
        colnames(dat$D_h)
      ),
      n_h = dat$n_h
    )
  })
  
  # Stack \hat\theta_h' across horizons
  coefficients <- do.call(
    rbind,
    lapply(fits, function(fit) fit$coefficients)
  )
  
  rownames(coefficients) <- paste0("h_", horizons)
  
  # Extract \hat\beta_h in particular
  beta <- coefficients[, "beta"]
  
  n_h <- vapply(
    fits,
    function(fit) fit$n_h,
    numeric(1)
  )
  
  list(
    method = "median_lp",
    horizons = horizons,
    n_h = n_h,
    coefficients = coefficients,
    beta = beta
  )
}

# Mean local projection

.modal_objective <- function(theta, y_h, D_h, bandwidth) {
  
  if (bandwidth <= 0) {
    stop("bandwidth must be positive.")
  }
  
  residuals <- y_h - D_h %*% theta
  
  mean(
    dnorm(residuals/bandwidth)
  ) / bandwidth
}

.fit_mem <- function(theta_start, y_h, D_h, bandwidth,
                     tol_theta = 1e-6,
                     tol_objective = 1e-8,
                     max_iter = 1000) {
  
  if (any(!is.finite(theta_start))) {
    stop("theta_start must contain only finite values.")
  }
  
  theta <- theta_start
  objective <- .modal_objective(theta, y_h, D_h, bandwidth)
  
  converged <- FALSE
  failure <- "max_iterations"
  
  for (g in seq_len(max_iter)) {
    
    # E-step
    residuals <- y_h - as.numeric(D_h %*% theta)
    weights <- dnorm(residuals / bandwidth)
    
    # ensure weights are all finite and non-degenerate
    if (any(!is.finite(weights)) ||
        sum(weights) <= .Machine$double.eps) {
      
      return(list(
        coefficients = theta,
        objective = objective,
        iterations = g,
        converged = FALSE,
        failure = "degenerate_weights"
      ))
    }
    
    # M-step
    
    fit <- lm.wfit(
      x = D_h,
      y = y_h,
      w = weights
    )
    
    # ensure non-singular weighted design matrix
    if (fit$rank < ncol(D_h)) {
      return(list(
        coefficients = theta,
        objective = objective,
        iterations = g,
        converged = FALSE,
        failure = "weighted_rank_deficiency"
      ))
    }
    
    theta_new <- fit$coefficients
    
    if (any(!is.finite(theta_new))) {
      return(list(
        coefficients = theta,
        objective = objective,
        iterations = g,
        converged = FALSE,
        failure = "nonfinite_coefficients"
      ))
    }
    
    objective_new <- .modal_objective(
      theta_new,
      y_h,
      D_h,
      bandwidth
    )
    
    # Check stopping criterion
    if (
      sqrt(sum((theta_new - theta)^2)) < tol_theta &&
      abs(objective_new - objective) < tol_objective
    ) {
      theta <- theta_new
      objective <- objective_new
      converged <- TRUE
      failure <- NULL
      break
    }
    
    theta <- theta_new
    objective <- objective_new
  }
  
  list(
    coefficients = theta,
    objective = objective,
    iterations = g,
    converged = converged,
    failure = failure
  )
}

fit_modal_lp <- function(y, x, z = NULL, horizons,
                         bw_constant = 1.6,
                         tol_theta = 1e-6,
                         tol_objective = 1e-8,
                         max_iter = 1000) {
  
  horizons <- sort(unique(horizons))
  
  if (length(horizons) == 0) {
    stop("At least one horizon must be supplied.")
  }
  
  if (!is.finite(bw_constant) || bw_constant <= 0) {
    stop("bw_constant must be positive and finite.")
  }
  
  fits <- lapply(horizons, function(h) {
    
    dat <- .prepare_lp_data(
      y = y,
      x = x,
      z = z,
      h = h
    )
    
    # Ensure theta_h is identified
    if (qr(dat$D_h)$rank < ncol(dat$D_h)) {
      stop(
        paste0(
          "Design matrix is rank deficient at horizon h = ",
          h,
          "."
        )
      )
    }
    
    # OLS starting value + bandwidth selection rule
    
    ols_fit <- lm.fit(
      x = dat$D_h,
      y = dat$y_h
    )
    
    theta_ols <- ols_fit$coefficients
    
    ols_residuals <- dat$y_h -
      as.numeric(dat$D_h %*% theta_ols)
    
    mad_h <- median(
      abs(ols_residuals - median(ols_residuals))
    )
    
    bandwidth <- bw_constant *
      mad_h *
      dat$n_h^(-0.143)
    
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      stop(
        paste0(
          "Invalid bandwidth at horizon h = ",
          h,
          "."
        )
      )
    }

    # Median-regression starting value
    
    median_fit <- quantreg::rq.fit(
      x = dat$D_h,
      y = dat$y_h,
      tau = 0.5,
      method = "br"
    )
    
    theta_median <- median_fit$coefficients
    
    # starting values
    
    starts <- list(
      ols = theta_ols,
      median = theta_median
    )
    
    # run MEM from each starting value
    
    mem_runs <- lapply(
      starts,
      function(theta_start) {
        
        .fit_mem(
          theta_start = theta_start,
          y_h = dat$y_h,
          D_h = dat$D_h,
          bandwidth = bandwidth,
          tol_theta = tol_theta,
          tol_objective = tol_objective,
          max_iter = max_iter
        )
      }
    )
    
    # identify successfully converged runs
    converged <- vapply(
      mem_runs, 
      function(run) {
        isTRUE(run$converged) &&
          is.finite(run$objective)
      },
      logical(1)
    )
    
    #could maybe make it so the function doesn't totally stop just because
    #none of the runs converged at horizon h. Maybe better to recover
    #all the horizon estimates we can and note which ones failed.
    
    if (!any(converged)) {
      stop(
        paste0(
          "MEM failed to converge from all starting values at horizon h = ",
          h,
          "."
        )
      )
    }
    
    # choose converged run with largest objective
    
    converged_indices <- which(converged)
    
    objectives <- vapply(
      mem_runs[converged_indices],
      function(run) run$objective,
      numeric(1)
    )
    
    best_index <- converged_indices[
      which.max(objectives)
    ]
    
    best_fit <- mem_runs[[best_index]]
    
    # Store horizon-specific results
    list(
      coefficients = setNames(
        as.numeric(best_fit$coefficients),
        colnames(dat$D_h)
      ),
      n_h = dat$n_h,
      bandwidth = bandwidth,
      objective = best_fit$objective,
      iterations = best_fit$iterations,
      selected_start = names(starts)[best_index],
      n_converged_starts = sum(converged)
    )
  })
  
  # Combine horizon-specific estimates
  
  coefficients <- do.call(
    rbind,
    lapply(fits, function(fit) fit$coefficients)
  )
  
  rownames(coefficients) <- paste0("h_", horizons)
  
  beta <- coefficients[, "beta"]
  
  n_h <- vapply(
    fits,
    function(fit) fit$n_h,
    numeric(1)
  )
  
  bandwidth <- vapply(
    fits,
    function(fit) fit$bandwidth,
    numeric(1)
  )
  
  objective <- vapply(
    fits,
    function(fit) fit$objective,
    numeric(1)
  )
  
  iterations <- vapply(
    fits,
    function(fit) fit$iterations,
    integer(1)
  )
  
  selected_start <- vapply(
    fits,
    function(fit) fit$selected_start,
    character(1)
  )
  
  n_converged_starts <- vapply(
    fits,
    function(fit) fit$n_converged_starts,
    integer(1)
  )
  
  list(
    method = "modal_lp",
    horizons = horizons,
    n_h = n_h,
    coefficients = coefficients,
    beta = beta,
    bandwidth = bandwidth,
    objective = objective,
    iterations = iterations,
    selected_start = selected_start,
    n_converged_starts = n_converged_starts,
    bw_constant = bw_constant
  )
}

