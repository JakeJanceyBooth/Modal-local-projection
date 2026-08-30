
# Finite-sample blocks-of-blocks bootstrap inference for Modal-LP.
# The resampling unit is the completed, already aligned horizon-specific
# regression row W_{t,h} = (y_{t+h}, x_t, z_t')'.

source("1_estimators.R")
source("2_dgps.R")

# Experiment settings -----------------------------------------------------

dgp_names <- c(
  "gaussian",
  "skewed",
  "disaster"
)

sample_sizes <- c(
  250L,
  1000L
)

horizons <- c(
  1L,
  5L,
  20L
)

n_outer_replications <- 100
n_bootstrap_replications <- 199

confidence_level <- 0.95
delta <- 1
bw_constant <- 2.4
start_quantiles <- 0.5

# Baseline:
block_length_rules <- c("fixed_4")

# Optional sensitivity, approximately doubling the bootstrap work:
# block_length_rules <- c("fixed_4", "cube_root")

modal_tol_theta <- 1e-6
modal_tol_objective <- 1e-8
modal_max_iter <- 1000L

simulation_seed_base <- 812301L
bootstrap_seed_base <- 1700001L

# Optional large-Monte-Carlo comparison from script 4.
# The inference experiment does not depend on this file being present.
response_mc_results_file <-
  file.path("results", "modal_results_R1000.rds")

response_mc_results_file <- NULL

results_directory <- file.path("results", "inference")

dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

if (
  n_outer_replications < 1L ||
  n_bootstrap_replications < 2L
) {
  stop(
    paste(
      "n_outer_replications must be positive and",
      "n_bootstrap_replications must be at least two."
    )
  )
}

if (
  confidence_level <= 0 ||
  confidence_level >= 1
) {
  stop("confidence_level must lie strictly between zero and one.")
}

if (
  any(horizons < 1L) ||
  any(horizons >= min(sample_sizes))
) {
  stop("Each horizon must be positive and smaller than every sample size.")
}

supported_block_length_rules <- c(
  "fixed_4",
  "cube_root"
)

if (
  !all(
    block_length_rules %in%
    supported_block_length_rules
  )
) {
  stop("Unknown block-length rule.")
}

alpha <- 1 - confidence_level
confidence_probabilities <- c(
  alpha / 2,
  1 - alpha / 2
)

run_label <- paste0(
  "blocks_of_blocks_R",
  n_outer_replications,
  "_B",
  n_bootstrap_replications,
  "_h",
  paste(horizons, collapse = "-"),
  "_",
  paste(block_length_rules, collapse = "-")
)

run_started_at <- Sys.time()


# Small bootstrap helpers -------------------------------------------------

block_length_from_rule <- function(
    block_length_rule,
    n_h) {
  if (block_length_rule == "fixed_4") {
    block_length <- 4L
  } else if (block_length_rule == "cube_root") {
    # Simple sensitivity rule only; no optimality claim is attached to it.
    block_length <- as.integer(
      ceiling(n_h^(1 / 3))
    )
  } else {
    stop("Unknown block-length rule.")
  }
  
  if (
    block_length < 1L ||
    block_length > n_h
  ) {
    stop("The selected block length is invalid for this regression sample.")
  }
  
  block_length
}


draw_non_circular_block_indices <- function(
    n_rows,
    block_length) {
  n_candidate_blocks <-
    n_rows - block_length + 1L
  
  n_blocks_to_draw <-
    ceiling(n_rows / block_length)
  
  block_starts <- sample.int(
    n = n_candidate_blocks,
    size = n_blocks_to_draw,
    replace = TRUE
  )
  
  resampled_indices <- unlist(
    lapply(
      block_starts,
      function(block_start) {
        block_start +
          seq.int(0L, block_length - 1L)
      }
    ),
    use.names = FALSE
  )
  
  resampled_indices[seq_len(n_rows)]
}


# This is the horizon-specific part of fit_modal_lp(), applied directly to
# an already aligned outcome and design matrix. It deliberately calls the
# existing MEM routine so that the convergence rule remains unchanged.
fit_aligned_modal_regression <- function(
    y_h,
    D_h,
    bw_constant,
    start_quantiles,
    tol_theta,
    tol_objective,
    max_iter) {
  y_h <- as.numeric(y_h)
  D_h <- as.matrix(D_h)
  
  n_h <- length(y_h)
  
  if (nrow(D_h) != n_h) {
    stop("The aligned outcome and design matrix have different lengths.")
  }
  
  if (qr(D_h)$rank < ncol(D_h)) {
    stop("The aligned Modal-regression design matrix is rank deficient.")
  }
  
  ols_fit <- stats::lm.fit(
    x = D_h,
    y = y_h
  )
  
  theta_ols <- ols_fit$coefficients
  
  if (any(!is.finite(theta_ols))) {
    stop("The OLS starting value is not finite.")
  }
  
  ols_residuals <- ols_fit$residuals
  
  residual_mad <- stats::median(
    abs(
      ols_residuals -
        stats::median(ols_residuals)
    )
  )
  
  bandwidth <-
    bw_constant *
    residual_mad *
    n_h^(-0.143)
  
  if (
    !is.finite(bandwidth) ||
    bandwidth <= 0
  ) {
    stop("The data-dependent Modal-regression bandwidth is not positive.")
  }
  
  start_quantiles <- sort(
    unique(as.numeric(start_quantiles))
  )
  
  if (
    length(start_quantiles) == 0L ||
    any(!is.finite(start_quantiles)) ||
    any(start_quantiles <= 0) ||
    any(start_quantiles >= 1)
  ) {
    stop("All starting quantiles must lie strictly between zero and one.")
  }
  
  quantile_starts <- lapply(
    start_quantiles,
    function(tau) {
      quantreg::rq.fit(
        x = D_h,
        y = y_h,
        tau = tau,
        method = "br"
      )$coefficients
    }
  )
  
  names(quantile_starts) <- paste0(
    "q",
    format(
      start_quantiles,
      trim = TRUE,
      scientific = FALSE
    )
  )
  
  starting_values <- c(
    list(ols = theta_ols),
    quantile_starts
  )
  
  mem_runs <- lapply(
    starting_values,
    function(theta_start) {
      .fit_mem(
        theta_start = theta_start,
        y_h = y_h,
        D_h = D_h,
        bandwidth = bandwidth,
        tol_theta = tol_theta,
        tol_objective = tol_objective,
        max_iter = max_iter
      )
    }
  )
  
  successful_runs <- vapply(
    mem_runs,
    function(mem_run) {
      isTRUE(mem_run$converged) &&
        is.finite(mem_run$objective)
    },
    logical(1)
  )
  
  if (!any(successful_runs)) {
    failure_reasons <- unique(
      vapply(
        mem_runs,
        function(mem_run) {
          failure <- mem_run$failure
          
          if (
            is.null(failure) ||
            length(failure) != 1L ||
            is.na(failure)
          ) {
            "not_converged"
          } else {
            as.character(failure)
          }
        },
        character(1)
      )
    )
    
    stop(
      paste0(
        "No Modal-regression starting value converged: ",
        paste(failure_reasons, collapse = ", ")
      )
    )
  }
  
  successful_objectives <- vapply(
    mem_runs[successful_runs],
    function(mem_run) {
      mem_run$objective
    },
    numeric(1)
  )
  
  selected_successful_run <-
    which.max(successful_objectives)
  
  successful_run_indices <- which(successful_runs)
  selected_run_index <-
    successful_run_indices[selected_successful_run]
  
  selected_run <- mem_runs[[selected_run_index]]
  
  coefficients <- selected_run$coefficients
  names(coefficients) <- colnames(D_h)
  
  list(
    coefficients = coefficients,
    n_h = n_h,
    bandwidth = bandwidth,
    objective = selected_run$objective,
    iterations = selected_run$iterations,
    selected_start = names(mem_runs)[selected_run_index],
    n_converged_starts = sum(successful_runs)
  )
}


# Outer Monte Carlo and bootstrap -----------------------------------------

max_sample_size <- max(sample_sizes)
max_horizon <- max(horizons)

expected_result_rows <-
  n_outer_replications *
  length(dgp_names) *
  length(sample_sizes) *
  length(horizons) *
  length(block_length_rules)

inference_rows <- vector(
  "list",
  expected_result_rows
)

result_row_index <- 0L

dgp_parameter_records <- stats::setNames(
  vector("list", length(dgp_names)),
  dgp_names
)

for (
  outer_replication in
  seq_len(n_outer_replications)
) {
  for (
    dgp_index in
    seq_along(dgp_names)
  ) {
    dgp_name <- dgp_names[dgp_index]
    
    simulation_seed <- as.integer(
      simulation_seed_base +
        100L * (outer_replication - 1L) +
        dgp_index
    )
    
    set.seed(simulation_seed)
    
    if (dgp_name == "gaussian") {
      simulated_dgp <- simulate_gaussian_dgp(
        T = max_sample_size,
        horizon = max_horizon
      )
    } else if (dgp_name == "skewed") {
      simulated_dgp <- simulate_skewed_dgp(
        T = max_sample_size,
        horizon = max_horizon
      )
    } else if (dgp_name == "disaster") {
      simulated_dgp <- simulate_disaster_dgp(
        T = max_sample_size,
        horizon = max_horizon
      )
    } else {
      stop("Unknown DGP.")
    }
    
    if (outer_replication == 1L) {
      dgp_parameter_records[[dgp_name]] <-
        simulated_dgp$parameters
    }
    
    for (
      sample_size_index in
      seq_along(sample_sizes)
    ) {
      sample_size <-
        sample_sizes[sample_size_index]
      
      sample_data <- simulated_dgp$data[
        seq_len(sample_size),
        ,
        drop = FALSE
      ]
      
      conditioning_variables <- as.matrix(
        sample_data[, "y", drop = FALSE]
      )
      
      for (
        horizon_index in
        seq_along(horizons)
      ) {
        horizon <- horizons[horizon_index]
        
        target_index <- which(
          simulated_dgp$true_response$horizon ==
            horizon
        )
        
        if (length(target_index) != 1L) {
          stop("The requested horizon has no unique DGP target.")
        }
        
        target_row <- simulated_dgp$true_response[
          target_index,
          ,
          drop = FALSE
        ]
        
        target_response <-
          delta * as.numeric(target_row$response)
        
        target_type <-
          "formal_linear_modal_response"
        
        coverage_interpretation <-
          "formal_coverage"
        
        aligned_data <- .prepare_lp_data(
          y = sample_data$y,
          x = sample_data$x,
          z = conditioning_variables,
          h = horizon
        )
        
        # The intercept is deterministic and is therefore not included in W.
        completed_regression_rows <- cbind(
          outcome = aligned_data$y_h,
          aligned_data$D_h[
            ,
            -1L,
            drop = FALSE
          ]
        )
        
        original_fit <- tryCatch(
          fit_aligned_modal_regression(
            y_h = aligned_data$y_h,
            D_h = aligned_data$D_h,
            bw_constant = bw_constant,
            start_quantiles = start_quantiles,
            tol_theta = modal_tol_theta,
            tol_objective =
              modal_tol_objective,
            max_iter = modal_max_iter
          ),
          error = function(error) {
            error
          }
        )
        
        if (inherits(original_fit, "error")) {
          original_fit_success <- FALSE
          original_fit_failure <-
            conditionMessage(original_fit)
          
          original_beta <- NA_real_
          original_response_estimate <- NA_real_
          original_bandwidth <- NA_real_
          original_iterations <- NA_integer_
          original_selected_start <- NA_character_
          original_n_converged_starts <- NA_integer_
        } else {
          original_beta <- unname(
            original_fit$coefficients["beta"]
          )
          
          original_fit_success <-
            length(original_beta) == 1L &&
            is.finite(original_beta)
          
          original_fit_failure <- if (
            original_fit_success
          ) {
            NA_character_
          } else {
            "non_finite_shock_coefficient"
          }
          
          original_response_estimate <- if (
            original_fit_success
          ) {
            delta * original_beta
          } else {
            NA_real_
          }
          
          original_bandwidth <-
            original_fit$bandwidth
          
          original_iterations <-
            original_fit$iterations
          
          original_selected_start <-
            original_fit$selected_start
          
          original_n_converged_starts <-
            original_fit$n_converged_starts
        }
        
        for (
          block_rule_index in
          seq_along(block_length_rules)
        ) {
          block_length_rule <-
            block_length_rules[block_rule_index]
          
          block_length <- block_length_from_rule(
            block_length_rule =
              block_length_rule,
            n_h = aligned_data$n_h
          )
          
          bootstrap_seed <- as.integer(
            bootstrap_seed_base +
              100000L *
              (outer_replication - 1L) +
              10000L * dgp_index +
              1000L * sample_size_index +
              100L * horizon_index +
              block_rule_index
          )
          
          if (original_fit_success) {
            n_bootstrap_attempts <-
              n_bootstrap_replications
            
            bootstrap_response_draws <- rep(
              NA_real_,
              n_bootstrap_replications
            )
            
            set.seed(bootstrap_seed)
            
            for (
              bootstrap_replication in
              seq_len(n_bootstrap_replications)
            ) {
              bootstrap_indices <-
                draw_non_circular_block_indices(
                  n_rows = aligned_data$n_h,
                  block_length = block_length
                )
              
              bootstrap_rows <-
                completed_regression_rows[
                  bootstrap_indices,
                  ,
                  drop = FALSE
                ]
              
              bootstrap_y <-
                bootstrap_rows[, "outcome"]
              
              bootstrap_D <- cbind(
                alpha = rep(
                  1,
                  aligned_data$n_h
                ),
                bootstrap_rows[
                  ,
                  -1L,
                  drop = FALSE
                ]
              )
              
              bootstrap_fit <- tryCatch(
                fit_aligned_modal_regression(
                  y_h = bootstrap_y,
                  D_h = bootstrap_D,
                  bw_constant = bw_constant,
                  start_quantiles =
                    start_quantiles,
                  tol_theta =
                    modal_tol_theta,
                  tol_objective =
                    modal_tol_objective,
                  max_iter =
                    modal_max_iter
                ),
                error = function(error) {
                  error
                }
              )
              
              if (!inherits(bootstrap_fit, "error")) {
                bootstrap_beta <- unname(
                  bootstrap_fit$coefficients[
                    "beta"
                  ]
                )
                
                if (
                  length(bootstrap_beta) == 1L &&
                  is.finite(bootstrap_beta)
                ) {
                  bootstrap_response_draws[
                    bootstrap_replication
                  ] <-
                    delta * bootstrap_beta
                }
              }
            }
          } else {
            n_bootstrap_attempts <- 0L
            bootstrap_response_draws <-
              numeric(0)
          }
          
          valid_bootstrap_draws <-
            bootstrap_response_draws[
              is.finite(
                bootstrap_response_draws
              )
            ]
          
          n_bootstrap_success <-
            length(valid_bootstrap_draws)
          
          bootstrap_success_rate <- if (
            n_bootstrap_attempts > 0L
          ) {
            n_bootstrap_success /
              n_bootstrap_attempts
          } else {
            NA_real_
          }
          
          valid_interval <-
            n_bootstrap_success >= 2L
          
          if (valid_interval) {
            bootstrap_se <- stats::sd(
              valid_bootstrap_draws
            )
            
            interval_quantiles <-
              stats::quantile(
                valid_bootstrap_draws,
                probs =
                  confidence_probabilities,
                names = FALSE,
                type = 7
              )
            
            ci_lower <-
              interval_quantiles[1L]
            
            ci_upper <-
              interval_quantiles[2L]
            
            interval_width <-
              ci_upper - ci_lower
          } else {
            bootstrap_se <- NA_real_
            ci_lower <- NA_real_
            ci_upper <- NA_real_
            interval_width <- NA_real_
          }
          
          # Formal coverage is defined only for the correctly specified
          # Gaussian linear Modal-LP.
          coverage_indicator <- if (
            valid_interval
          ) {
            target_response >= ci_lower &&
              target_response <= ci_upper
          } else {
            NA
          }
          # For downside risk this compares a global linear estimator with
          # a pointwise nonlinear response and is therefore diagnostic only.
          pointwise_diagnostic_coverage_indicator <- if (
            valid_interval
          ) {
            target_response >= ci_lower &&
              target_response <= ci_upper
          } else {
            NA
          }
          
          result_row_index <-
            result_row_index + 1L
          
          inference_rows[[result_row_index]] <-
            data.frame(
              dgp = dgp_name,
              sample_size = sample_size,
              outer_replication =
                outer_replication,
              horizon = horizon,
              block_length_rule =
                block_length_rule,
              block_length = block_length,
              n_h = aligned_data$n_h,
              simulation_seed =
                simulation_seed,
              bootstrap_seed =
                bootstrap_seed,
              target_response =
                target_response,
              target_type = target_type,
              coverage_interpretation =
                coverage_interpretation,
              original_beta =
                original_beta,
              original_response_estimate =
                original_response_estimate,
              original_bandwidth =
                original_bandwidth,
              original_iterations =
                original_iterations,
              original_selected_start =
                original_selected_start,
              original_n_converged_starts =
                original_n_converged_starts,
              original_fit_success =
                original_fit_success,
              original_fit_failure =
                original_fit_failure,
              bootstrap_se = bootstrap_se,
              ci_lower = ci_lower,
              ci_upper = ci_upper,
              interval_width =
                interval_width,
              coverage_indicator =
                coverage_indicator,
              pointwise_diagnostic_coverage_indicator =
                pointwise_diagnostic_coverage_indicator,
              n_bootstrap_attempts =
                n_bootstrap_attempts,
              n_bootstrap_success =
                n_bootstrap_success,
              n_bootstrap_failure =
                n_bootstrap_attempts -
                n_bootstrap_success,
              bootstrap_success_rate =
                bootstrap_success_rate,
              valid_interval =
                valid_interval,
              stringsAsFactors = FALSE
            )
        }
      }
    }
  }
  
  if (
    outer_replication == 1L ||
    outer_replication %% 10L == 0L ||
    outer_replication ==
    n_outer_replications
  ) {
    message(
      "Completed outer replication ",
      outer_replication,
      " of ",
      n_outer_replications,
      "."
    )
  }
}

inference_replication_results <-
  dplyr::bind_rows(inference_rows)

if (
  nrow(inference_replication_results) !=
  expected_result_rows
) {
  stop("The inference output has an unexpected number of rows.")
}

run_completed_at <- Sys.time()

inference_run_settings <- list(
  run_label = run_label,
  run_started_at = run_started_at,
  run_completed_at = run_completed_at,
  dgp_names = dgp_names,
  dgp_parameters = dgp_parameter_records,
  sample_sizes = sample_sizes,
  horizons = horizons,
  n_outer_replications =
    n_outer_replications,
  n_bootstrap_replications =
    n_bootstrap_replications,
  confidence_level = confidence_level,
  delta = delta,
  bw_constant = bw_constant,
  bandwidth_rule =
    "bw_constant * OLS_residual_MAD * n_h^(-0.143)",
  start_quantiles = start_quantiles,
  modal_tol_theta = modal_tol_theta,
  modal_tol_objective =
    modal_tol_objective,
  modal_max_iter = modal_max_iter,
  block_length_rules =
    block_length_rules,
  block_length_definitions = list(
    fixed_4 = 4L,
    cube_root =
      "ceiling(n_h^(1/3)); sensitivity only"
  ),
  bootstrap_method =
    "non-circular overlapping moving blocks of completed LP rows",
  bootstrap_horizon_coupling =
    "independent pointwise resampling across horizons",
  interval_method =
    "Efron percentile",
  percentile_quantile_type = 7L,
  simulation_seed_base =
    simulation_seed_base,
  bootstrap_seed_base =
    bootstrap_seed_base,
  response_mc_results_file =
    response_mc_results_file
)


# Checkpoint expensive output before summaries and plots ------------------

replication_results_path <- file.path(
  results_directory,
  paste0(
    "inference_replication_results_",
    run_label,
    ".rds"
  )
)

run_settings_path <- file.path(
  results_directory,
  paste0(
    "inference_run_settings_",
    run_label,
    ".rds"
  )
)

saveRDS(
  inference_replication_results,
  replication_results_path
)

saveRDS(
  inference_run_settings,
  run_settings_path
)


# Across-replication summaries --------------------------------------------

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    NA_real_
  } else {
    mean(x)
  }
}

finite_median <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    NA_real_
  } else {
    stats::median(x)
  }
}

finite_sd <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) <= 1L) {
    NA_real_
  } else {
    stats::sd(x)
  }
}

inference_summary <-
  inference_replication_results |>
  dplyr::group_by(
    dgp,
    sample_size,
    horizon,
    block_length_rule,
    block_length,
    target_type,
    coverage_interpretation
  ) |>
  dplyr::summarise(
    target_response =
      dplyr::first(target_response),
    n_outer_replications =
      dplyr::n(),
    n_original_fit_success =
      sum(original_fit_success),
    outer_fit_success_rate =
      mean(original_fit_success),
    mean_original_response =
      finite_mean(
        original_response_estimate
      ),
    outer_sampling_sd =
      finite_sd(
        original_response_estimate
      ),
    mean_bootstrap_se =
      finite_mean(bootstrap_se),
    median_bootstrap_se =
      finite_median(bootstrap_se),
    total_bootstrap_attempts =
      sum(n_bootstrap_attempts),
    total_bootstrap_success =
      sum(n_bootstrap_success),
    mean_bootstrap_success_rate =
      finite_mean(
        bootstrap_success_rate
      ),
    n_valid_intervals =
      sum(valid_interval),
    mean_interval_width =
      finite_mean(interval_width),
    median_interval_width =
      finite_median(interval_width),
    n_formal_coverage_replications =
      sum(!is.na(coverage_indicator)),
    formal_coverage = if (
      any(!is.na(coverage_indicator))
    ) {
      mean(
        coverage_indicator,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    n_pointwise_diagnostic_replications =
      sum(
        !is.na(
          pointwise_diagnostic_coverage_indicator
        )
      ),
    pointwise_diagnostic_coverage = if (
      any(
        !is.na(
          pointwise_diagnostic_coverage_indicator
        )
      )
    ) {
      mean(
        pointwise_diagnostic_coverage_indicator,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    bootstrap_fit_success_rate = ifelse(
      total_bootstrap_attempts > 0L,
      total_bootstrap_success /
        total_bootstrap_attempts,
      NA_real_
    ),
    bootstrap_se_ratio_to_outer_sd = ifelse(
      is.finite(outer_sampling_sd) &
        outer_sampling_sd > 0,
      mean_bootstrap_se /
        outer_sampling_sd,
      NA_real_
    ),
    formal_coverage_mcse = ifelse(
      n_formal_coverage_replications > 0L,
      sqrt(
        formal_coverage *
          (1 - formal_coverage) /
          n_formal_coverage_replications
      ),
      NA_real_
    ),
    pointwise_diagnostic_coverage_mcse =
      ifelse(
        n_pointwise_diagnostic_replications >
          0L,
        sqrt(
          pointwise_diagnostic_coverage *
            (
              1 -
                pointwise_diagnostic_coverage
            ) /
            n_pointwise_diagnostic_replications
        ),
        NA_real_
      )
  ) |>
  dplyr::arrange(
    match(dgp, dgp_names),
    sample_size,
    horizon,
    match(
      block_length_rule,
      block_length_rules
    )
  )


# Optional comparison with the large script-4 Monte Carlo ----------------

response_mc_sampling_sd_by_cell <- data.frame(
  dgp = character(),
  sample_size = integer(),
  horizon = integer(),
  response_mc_sampling_sd = double(),
  n_response_mc_replications = integer(),
  n_response_mc_success = integer(),
  stringsAsFactors = FALSE
)

if (
  !is.null(response_mc_results_file) &&
  file.exists(response_mc_results_file)
) {
  response_mc_results <- tryCatch(
    readRDS(response_mc_results_file),
    error = function(error) {
      warning(
        "Could not read the optional script-4 results: ",
        conditionMessage(error)
      )
      
      NULL
    }
  )
  
  required_response_mc_columns <- c(
    "dgp",
    "estimator",
    "sample_size",
    "replication",
    "horizon",
    "estimate",
    "failure"
  )
  
  if (
    !is.null(response_mc_results) &&
    all(
      required_response_mc_columns %in%
      names(response_mc_results)
    )
  ) {
    response_mc_sampling_sd_by_cell <-
      response_mc_results |>
      dplyr::filter(
        estimator == "Modal LP",
        dgp %in% dgp_names,
        sample_size %in% sample_sizes,
        horizon %in% horizons
      ) |>
      dplyr::mutate(
        reference_response_estimate =
          delta * estimate,
        reference_fit_success =
          is.na(failure) &
          is.finite(
            reference_response_estimate
          )
      ) |>
      dplyr::group_by(
        dgp,
        sample_size,
        horizon
      ) |>
      dplyr::summarise(
        response_mc_sampling_sd =
          finite_sd(
            reference_response_estimate[
              reference_fit_success
            ]
          ),
        n_response_mc_replications =
          dplyr::n(),
        n_response_mc_success =
          sum(reference_fit_success),
        .groups = "drop"
      )
  } else if (!is.null(response_mc_results)) {
    warning(
      paste(
        "The optional script-4 object lacks the",
        "columns required for the sampling-SD comparison."
      )
    )
  }
}

inference_summary <-
  inference_summary |>
  dplyr::left_join(
    response_mc_sampling_sd_by_cell,
    by = c(
      "dgp",
      "sample_size",
      "horizon"
    )
  ) |>
  dplyr::mutate(
    sampling_sd_for_comparison =
      dplyr::coalesce(
        response_mc_sampling_sd,
        outer_sampling_sd
      ),
    sampling_sd_source =
      dplyr::case_when(
        is.finite(
          response_mc_sampling_sd
        ) ~
          "script_4_response_mc",
        is.finite(outer_sampling_sd) ~
          "current_inference_outer_mc",
        TRUE ~ NA_character_
      ),
    bootstrap_se_ratio_to_sampling_sd =
      ifelse(
        is.finite(
          sampling_sd_for_comparison
        ) &
          sampling_sd_for_comparison > 0,
        mean_bootstrap_se /
          sampling_sd_for_comparison,
        NA_real_
      )
  )

bootstrap_se_calibration <-
  inference_summary |>
  dplyr::select(
    dgp,
    sample_size,
    horizon,
    block_length_rule,
    block_length,
    mean_bootstrap_se,
    median_bootstrap_se,
    outer_sampling_sd,
    response_mc_sampling_sd,
    sampling_sd_for_comparison,
    sampling_sd_source,
    bootstrap_se_ratio_to_outer_sd,
    bootstrap_se_ratio_to_sampling_sd,
    n_original_fit_success,
    n_response_mc_success
  )


# Minimal inferential diagnostic plots ------------------------------------

plot_summary_data <-
  inference_summary |>
  dplyr::mutate(
    dgp_label = dplyr::recode(
      dgp,
      gaussian = "Gaussian",
      skewed = "Skewed",
      disaster = "Disaster"
    ),
    sample_size_label =
      paste0("T = ", sample_size),
    block_rule_label = dplyr::recode(
      block_length_rule,
      fixed_4 = "l = 4",
      cube_root = "ceiling(n_h^(1/3))"
    )
  )

bootstrap_se_plot_data <-
  plot_summary_data |>
  dplyr::select(
    dgp_label,
    sample_size_label,
    horizon,
    block_rule_label,
    mean_bootstrap_se,
    sampling_sd_for_comparison
  ) |>
  tidyr::pivot_longer(
    cols = c(
      mean_bootstrap_se,
      sampling_sd_for_comparison
    ),
    names_to = "uncertainty_measure",
    values_to = "standard_deviation"
  ) |>
  dplyr::mutate(
    uncertainty_measure =
      dplyr::recode(
        uncertainty_measure,
        mean_bootstrap_se =
          "Mean bootstrap SE",
        sampling_sd_for_comparison =
          "Monte Carlo sampling SD"
      )
  )

bootstrap_se_calibration_plot <-
  ggplot2::ggplot(
    bootstrap_se_plot_data,
    ggplot2::aes(
      x = horizon,
      y = standard_deviation,
      color = uncertainty_measure,
      linetype = block_rule_label,
      group = interaction(
        uncertainty_measure,
        block_rule_label
      )
    )
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
    rows = ggplot2::vars(dgp_label),
    cols = ggplot2::vars(
      sample_size_label
    ),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Standard deviation",
    color = NULL,
    linetype = "Block-length rule",
    caption = paste(
      "The Monte Carlo sampling SD comes from script 4",
      "when its saved results are available; otherwise",
      "it is calculated from this experiment's outer replications."
    )
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

formal_coverage_plot_data <-
  plot_summary_data |>
  dplyr::filter(
    n_formal_coverage_replications > 0L
  ) |>
  dplyr::mutate(
    mc_lower = pmax(
      0,
      formal_coverage -
        1.96 * formal_coverage_mcse
    ),
    mc_upper = pmin(
      1,
      formal_coverage +
        1.96 * formal_coverage_mcse
    )
  )

formal_coverage_plot <-
  ggplot2::ggplot(
    formal_coverage_plot_data,
    ggplot2::aes(
      x = horizon,
      y = formal_coverage,
      color = block_rule_label,
      group = block_rule_label
    )
  ) +
  ggplot2::geom_hline(
    yintercept = confidence_level,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mc_lower,
      ymax = mc_upper
    ),
    width = 0.4,
    linewidth = 0.4,
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
    rows = ggplot2::vars(dgp_label),
    cols = ggplot2::vars(sample_size_label)
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Empirical percentile-interval coverage",
    color = "Block-length rule",
    caption = paste(
      "Monte Carlo SE for the estimated coverage probability."
    )
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

interval_width_plot_data <-
  plot_summary_data |>
  dplyr::select(
    dgp_label,
    sample_size_label,
    horizon,
    block_rule_label,
    mean_interval_width,
    median_interval_width
  ) |>
  tidyr::pivot_longer(
    cols = c(
      mean_interval_width,
      median_interval_width
    ),
    names_to = "width_summary",
    values_to = "interval_width"
  ) |>
  dplyr::mutate(
    width_summary = dplyr::recode(
      width_summary,
      mean_interval_width = "Mean",
      median_interval_width = "Median"
    )
  )

interval_width_plot <-
  ggplot2::ggplot(
    interval_width_plot_data,
    ggplot2::aes(
      x = horizon,
      y = interval_width,
      color = width_summary,
      linetype = block_rule_label,
      group = interaction(
        width_summary,
        block_rule_label
      )
    )
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
    rows = ggplot2::vars(dgp_label),
    cols = ggplot2::vars(
      sample_size_label
    ),
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = horizons
  ) +
  ggplot2::labs(
    x = "Horizon",
    y = "Percentile-interval width",
    color = "Across-replication summary",
    linetype = "Block-length rule"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

inference_plots <- list(
  bootstrap_se_calibration =
    bootstrap_se_calibration_plot,
  formal_coverage =
    formal_coverage_plot,
  interval_width =
    interval_width_plot
)


# Save derived output -----------------------------------------------------

summary_path <- file.path(
  results_directory,
  paste0(
    "inference_summary_",
    run_label,
    ".rds"
  )
)

se_calibration_path <- file.path(
  results_directory,
  paste0(
    "bootstrap_se_calibration_",
    run_label,
    ".rds"
  )
)

response_mc_sampling_sd_path <- file.path(
  results_directory,
  paste0(
    "response_mc_sampling_sd_",
    run_label,
    ".rds"
  )
)

plots_path <- file.path(
  results_directory,
  paste0(
    "inference_plots_",
    run_label,
    ".rds"
  )
)

saveRDS(
  inference_summary,
  summary_path
)

saveRDS(
  bootstrap_se_calibration,
  se_calibration_path
)

saveRDS(
  response_mc_sampling_sd_by_cell,
  response_mc_sampling_sd_path
)

saveRDS(
  inference_plots,
  plots_path
)

inference_run_settings$response_mc_comparison_used <-
  any(
    is.finite(
      inference_summary$response_mc_sampling_sd
    )
  )

saveRDS(
  inference_run_settings,
  run_settings_path
)

message(
  "Inference experiment complete. Run label: ",
  run_label
)