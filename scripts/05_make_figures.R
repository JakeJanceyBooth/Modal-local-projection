# Load final saved Monte Carlo objects and produce publication-ready displays.
# This script does not estimate models or rerun simulations.

library(ggplot2)
library(dplyr)

# Paths and final saved results ----
response_results_path <- file.path(
  "results", "responses", "response_plot_data_R1000.rds"
)

forecast_run_label <- paste0(
  "fixed_width_coverage_",
  "R5000_O20_h1-2-5-10"
)

average_forecast_results_path <- file.path(
  "results", "forecasting",
  paste0("average_forecast_summary_", forecast_run_label, ".rds")
)

paired_coverage_results_path <- file.path(
  "results", "forecasting",
  paste0("paired_coverage_summary_", forecast_run_label, ".rds")
)

disaster_coverage_results_path <- file.path(
  "results", "forecasting",
  paste0("disaster_coverage_summary_", forecast_run_label, ".rds")
)

figure_directory <- file.path("results", "figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

response_plot_data <- readRDS(response_results_path)
average_forecast_summary <- readRDS(average_forecast_results_path)
paired_coverage_summary <- readRDS(paired_coverage_results_path)
disaster_coverage_summary <- readRDS(disaster_coverage_results_path)

# Verify the saved summary schemas ----
required_response_columns <- c(
  "dgp",
  "estimator",
  "sample_size",
  "horizon",
  "truth",
  "average_estimate"
)

if (!all(required_response_columns %in% names(response_plot_data))) {
  stop("The final response plotting object has unexpected columns.")
}

required_average_forecast_columns <- c(
  "dgp",
  "sample_size",
  "estimator",
  "horizon",
  "average_forecast"
)

if (!all(required_average_forecast_columns %in%
         names(average_forecast_summary))) {
  stop("The final average-forecast summary has unexpected columns.")
}

required_paired_coverage_columns <- c(
  "dgp",
  "sample_size",
  "comparison",
  "horizon",
  "half_width",
  "mean_coverage_difference"
)

if (!all(required_paired_coverage_columns %in%
         names(paired_coverage_summary))) {
  stop("The final paired-coverage summary has unexpected columns.")
}

required_disaster_coverage_columns <- c(
  "dgp",
  "sample_size",
  "estimator",
  "horizon",
  "disaster_path",
  "half_width",
  "mean_coverage"
)

if (!all(required_disaster_coverage_columns %in%
         names(disaster_coverage_summary))) {
  stop("The final disaster-coverage summary has unexpected columns.")
}

# Reusable figure style ----
theme_modal_paper <- function(base_size = 8.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),

      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.22
      ),

      axis.line = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", size = rel(1)),
      axis.ticks = element_line(colour = "black", linewidth = 0.20),
      axis.ticks.length = grid::unit(2.5, "pt"),

      strip.background = element_blank(),
      strip.text = element_text(colour = "black", face = "plain"),
      strip.text.x = element_text(margin = margin(b = 3)),
      strip.text.y.left = element_text(
        angle = 0,
        hjust = 1,
        margin = margin(r = 4)
      ),
      strip.placement = "outside",
      panel.spacing = grid::unit(6, "pt"),

      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.width = grid::unit(18, "pt"),
      legend.spacing.x = grid::unit(5, "pt"),
      legend.margin = margin(t = 3),

      plot.margin = margin(t = 3, r = 4, b = 2, l = 2)
    )
}

# Stable response-series aesthetics ----
response_series_levels <- c(
  "True response",
  "Mean LP",
  "Median LP",
  "Modal LP",
  "VAR"
)

response_series_colours <- c(
  "True response" = "#000000",
  "Mean LP" = "#3F3F3F",
  "Median LP" = "#2B2B2B",
  "Modal LP" = "#1F5364",
  "VAR" = "#4D4D4D"
)

response_series_linetypes <- c(
  "True response" = "solid",
  "Mean LP" = "longdash",
  "Median LP" = "dotted",
  "Modal LP" = "dotdash",
  "VAR" = "twodash"
)

# ggplot2 linewidths are measured in millimetres.
response_series_linewidths <- c(
  "True response" = 0.70,
  "Mean LP" = 0.50,
  "Median LP" = 0.68,
  "Modal LP" = 0.55,
  "VAR" = 0.50
)

# Common DGP ordering and labels ----------------------------------------------
common_target_dgps <- c("gaussian", "skewed", "disaster", "rich_state")
common_target_dgp_labels <- c(
  "gaussian" = "Gaussian",
  "skewed" = "Skewed",
  "disaster" = "Rare disaster",
  "rich_state" = "Rich state"
)

# Figure 1: common-target response functions ----
figure1_sample_sizes <- c(250, 500, 1000)
common_target_plot_data <- response_plot_data |>
  filter(
    dgp %in% common_target_dgps,
    sample_size %in% figure1_sample_sizes,
    horizon %in% seq_len(20)
  ) |>
  mutate(
    dgp_label = factor(
      unname(common_target_dgp_labels[dgp]),
      levels = unname(common_target_dgp_labels[common_target_dgps])
    ),
    sample_size_label = factor(
      paste0("T = ", sample_size),
      levels = paste0("T = ", figure1_sample_sizes)
    ),
    estimator = factor(estimator, levels = response_series_levels[-1])
  )

# The population response is common across estimators in these four DGPs.
figure1_truth_data <- common_target_plot_data |>
  distinct(
    dgp_label,
    sample_size_label,
    horizon,
    truth
  ) |>
  transmute(
    dgp_label,
    sample_size_label,
    horizon,
    series = factor("True response", levels = response_series_levels),
    response = truth
  )

figure1_estimate_data <- common_target_plot_data |>
  transmute(
    dgp_label,
    sample_size_label,
    horizon,
    series = factor(as.character(estimator), levels = response_series_levels),
    response = average_estimate
  )

figure1_common_target_responses <- ggplot(
  mapping = aes(
    x = horizon,
    y = response,
    group = series,
    colour = series,
    linetype = series,
    linewidth = series
  )
) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.20) +
  geom_line(data = figure1_truth_data, lineend = "round") +
  geom_line(data = figure1_estimate_data, lineend = "round") +
  facet_grid(
    rows = vars(dgp_label),
    cols = vars(sample_size_label),
    switch = "y"
  ) +
  scale_colour_manual(
    values = response_series_colours,
    breaks = response_series_levels,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = response_series_linetypes,
    breaks = response_series_levels,
    drop = FALSE
  ) +
  scale_linewidth_manual(
    values = response_series_linewidths,
    breaks = response_series_levels,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  labs(x = "Horizon", y = "Response", colour = NULL, linetype = NULL) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_modal_paper()

figure1_path <- file.path(figure_directory, "fig_common_target_responses.pdf")

ggsave(
  filename = figure1_path,
  plot = figure1_common_target_responses,
  device = "pdf",
  width = 7.0,
  height = 7.3,
  units = "in",
  bg = "white"
)

# Figure 2: DGP 5 distinct-target response functions ----
# In the saved DGP-5 results, the truth column is estimator specific
figure2_sample_sizes <- c(250, 500, 1000)
figure2_target_levels <- c("Mean target", "Median target", "Modal target")
figure2_target_by_estimator <- c(
  "Mean LP" = "Mean target",
  "VAR" = "Mean target",
  "Median LP" = "Median target",
  "Modal LP" = "Modal target"
)

figure2_plot_data <- response_plot_data |>
  filter(
    dgp == "downside_risk",
    sample_size %in% figure2_sample_sizes,
    horizon %in% seq_len(20)
  ) |>
  mutate(
    target_label = factor(
      unname(figure2_target_by_estimator[as.character(estimator)]),
      levels = figure2_target_levels
    ),
    sample_size_label = factor(
      paste0("T = ", sample_size),
      levels = paste0("T = ", figure2_sample_sizes)
    ),
    estimator = factor(estimator, levels = response_series_levels[-1])
  )

# This also verifies that the Mean-LP and VAR truth entries coincide.
figure2_truth_check <- figure2_plot_data |>
  group_by(
    target_label,
    sample_size_label,
    horizon
  ) |>
  summarise(n_truth_values = n_distinct(truth), .groups = "drop")

if (any(figure2_truth_check$n_truth_values != 1L)) {
  stop(
    paste(
      "The saved DGP-5 truth values do not agree",
      "with the intended estimator-to-target mapping."
    )
  )
}

figure2_truth_data <- figure2_plot_data |>
  distinct(
    target_label,
    sample_size_label,
    horizon,
    truth
  ) |>
  transmute(
    target_label,
    sample_size_label,
    horizon,
    series = factor("True response", levels = response_series_levels),
    response = truth
  )

figure2_estimate_data <- figure2_plot_data |>
  transmute(
    target_label,
    sample_size_label,
    horizon,
    series = factor(as.character(estimator), levels = response_series_levels),
    response = average_estimate
  )

figure2_distinct_target_responses <- ggplot(
  mapping = aes(
    x = horizon,
    y = response,
    group = series,
    colour = series,
    linetype = series,
    linewidth = series
  )
) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.20) +
  geom_line(data = figure2_truth_data, lineend = "round") +
  geom_line(data = figure2_estimate_data, lineend = "round") +
  facet_grid(
    rows = vars(target_label),
    cols = vars(sample_size_label),
    switch = "y"
  ) +
  scale_colour_manual(
    values = response_series_colours,
    breaks = response_series_levels,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = response_series_linetypes,
    breaks = response_series_levels,
    drop = FALSE
  ) +
  scale_linewidth_manual(
    values = response_series_linewidths,
    breaks = response_series_levels,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = c(-0.15, -0.10, -0.05, 0),
    labels = c("-0.15", "-0.10", "-0.05", "0"),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(ylim = c(-0.15, 0.01)) +
  labs(x = "Horizon", y = "Response", colour = NULL, linetype = NULL) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_modal_paper()

figure2_path <- file.path(figure_directory, "fig_distinct_target_responses.pdf")

ggsave(
  filename = figure2_path,
  plot = figure2_distinct_target_responses,
  device = "pdf",
  width = 7.0,
  height = 5.8,
  units = "in",
  bg = "white"
)

# Figure 3: average forecasts ----
# The final saved object contains one average_forecast for each
# DGP / sample size / estimator / forecast-horizon cell.
figure3_sample_size <- 250
figure3_horizons <- c(1, 2, 5, 10)
figure3_plot_data <- average_forecast_summary |>
  filter(
    sample_size == figure3_sample_size,
    dgp %in% common_target_dgps,
    horizon %in% figure3_horizons
  ) |>
  transmute(
    dgp_label = factor(
      unname(common_target_dgp_labels[dgp]),
      levels = unname(common_target_dgp_labels[common_target_dgps])
    ),
    horizon,
    series = factor(estimator, levels = response_series_levels[-1]),
    average_forecast
  )

figure3_average_forecasts <- ggplot(
  data = figure3_plot_data,
  mapping = aes(
    x = horizon,
    y = average_forecast,
    group = series,
    colour = series,
    linetype = series,
    linewidth = series
  )
) +
  geom_line(lineend = "round") +
  facet_wrap(vars(dgp_label), ncol = 2, scales = "free_y") +
  scale_colour_manual(
    values = response_series_colours,
    breaks = response_series_levels[-1],
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = response_series_linetypes,
    breaks = response_series_levels[-1],
    drop = FALSE
  ) +
  scale_linewidth_manual(
    values = response_series_linewidths,
    breaks = response_series_levels[-1],
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = figure3_horizons,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    minor_breaks = NULL,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    x = "Forecast horizon",
    y = "Average forecast",
    colour = NULL,
    linetype = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_modal_paper()

figure3_path <- file.path(figure_directory, "fig_average_forecasts.pdf")

ggsave(
  filename = figure3_path,
  plot = figure3_average_forecasts,
  device = "pdf",
  width = 7.0,
  height = 4.7,
  units = "in",
  bg = "white"
)

# Figure 4: paired fixed-width coverage differences ----
coverage_comparison_levels <- c("Modal - Mean", "Modal - Median")
coverage_comparison_labels <- c(
  "Modal - Mean" = "Modal LP - Mean LP",
  "Modal - Median" = "Modal LP - Median LP"
)

coverage_comparison_colours <- c(
  "Modal - Mean" = "#1F5364",
  "Modal - Median" = "#4D4D4D"
)

coverage_comparison_linetypes <- c(
  "Modal - Mean" = "longdash",
  "Modal - Median" = "dotted"
)

coverage_comparison_linewidths <- c(
  "Modal - Mean" = 0.49,
  "Modal - Median" = 0.45
)

figure4_sample_size <- 250
figure4_horizons <- c(1, 2, 5, 10)
figure4_plot_data <- paired_coverage_summary |>
  filter(
    sample_size == figure4_sample_size,
    dgp %in% common_target_dgps,
    horizon %in% figure4_horizons,
    comparison %in% coverage_comparison_levels
  ) |>
  transmute(
    dgp_label = factor(
      unname(common_target_dgp_labels[dgp]),
      levels = unname(common_target_dgp_labels[common_target_dgps])
    ),
    horizon_label = factor(
      paste0("h = ", horizon),
      levels = paste0("h = ", figure4_horizons)
    ),
    comparison = factor(comparison, levels = coverage_comparison_levels),
    half_width,
    coverage_difference = mean_coverage_difference
  )

fixed_width_axis_breaks <- c(0.5, 1.0, 1.5)
fixed_width_axis_labels <- c("0.5", "1", "1.5")
fixed_width_panel_spacing <- grid::unit(10, "pt")
figure4_paired_coverage_differences <- ggplot(
  data = figure4_plot_data,
  mapping = aes(
    x = half_width,
    y = coverage_difference,
    group = comparison,
    colour = comparison,
    linetype = comparison,
    linewidth = comparison
  )
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.25) +
  geom_line(lineend = "round") +
  facet_grid(rows = vars(dgp_label), cols = vars(horizon_label), switch = "y") +
  scale_colour_manual(
    values = coverage_comparison_colours,
    breaks = coverage_comparison_levels,
    labels = coverage_comparison_labels,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = coverage_comparison_linetypes,
    breaks = coverage_comparison_levels,
    labels = coverage_comparison_labels,
    drop = FALSE
  ) +
  scale_linewidth_manual(
    values = coverage_comparison_linewidths,
    breaks = coverage_comparison_levels,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = fixed_width_axis_breaks,
    labels = fixed_width_axis_labels,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(-0.04, 0, 0.04),
    labels = c("-0.04", "0", "0.04"),
    minor_breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(0.25, 2.00), ylim = c(-0.08, 0.08)) +
  labs(
    x = "Interval half-width",
    y = "Coverage difference",
    colour = NULL,
    linetype = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_modal_paper() +
  theme(
    panel.spacing.x = fixed_width_panel_spacing,
    panel.spacing.y = fixed_width_panel_spacing
  )

figure4_path <- file.path(
  figure_directory,
  "fig_paired_coverage_differences.pdf"
)

ggsave(
  filename = figure4_path,
  plot = figure4_paired_coverage_differences,
  device = "pdf",
  width = 7.0,
  height = 7.2,
  units = "in",
  bg = "white"
)

# Figure 5: rare-disaster conditional coverage decomposition ------------------
figure5_sample_size <- 250
figure5_horizons <- c(1, 2, 5, 10)
figure5_plot_data <- disaster_coverage_summary |>
  filter(
    sample_size == figure5_sample_size,
    dgp == "disaster",
    horizon %in% figure5_horizons
  ) |>
  transmute(
    disaster_path_label = factor(
      disaster_path,
      levels = c("No disaster", "At least one disaster")
    ),
    horizon_label = factor(
      paste0("h = ", horizon),
      levels = paste0("h = ", figure5_horizons)
    ),
    half_width,
    series = factor(estimator, levels = response_series_levels[-1]),
    conditional_coverage = mean_coverage
  )

figure5_disaster_coverage_decomposition <- ggplot(
  data = figure5_plot_data,
  mapping = aes(
    x = half_width,
    y = conditional_coverage,
    group = series,
    colour = series,
    linetype = series,
    linewidth = series
  )
) +
  geom_line(lineend = "round") +
  facet_grid(
    rows = vars(disaster_path_label),
    cols = vars(horizon_label),
    switch = "y"
  ) +
  scale_colour_manual(
    values = response_series_colours,
    breaks = response_series_levels[-1],
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = response_series_linetypes,
    breaks = response_series_levels[-1],
    drop = FALSE
  ) +
  scale_linewidth_manual(
    values = response_series_linewidths,
    breaks = response_series_levels[-1],
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = fixed_width_axis_breaks,
    labels = fixed_width_axis_labels,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    minor_breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(0.25, 2.00), ylim = c(-0.02, 1.02)) +
  labs(
    x = "Interval half-width",
    y = "Conditional coverage",
    colour = NULL,
    linetype = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_modal_paper() +
  theme(panel.spacing.x = fixed_width_panel_spacing)

figure5_path <- file.path(
  figure_directory,
  "fig_disaster_coverage_decomposition.pdf"
)

ggsave(
  filename = figure5_path,
  plot = figure5_disaster_coverage_decomposition,
  device = "pdf",
  width = 7.0,
  height = 4.6,
  units = "in",
  bg = "white"
)
