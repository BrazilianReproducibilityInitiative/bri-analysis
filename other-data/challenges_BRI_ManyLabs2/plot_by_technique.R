### ### ### ### ### ### ### ### ### ### ### ###
#
# plot_by_technique.R
# Builds only the BRI + ManyLabs2 figure with Panel A separated by technique.
#
# Panel A: BRI effect-size dot plot using EPM + MTT + PCR.
# Panel B: ManyLabs2 dot plot restyled to match the BRI panel.
#
# Usage:
#   Rscript "plot_by_technique.R"
#
### ### ### ### ### ### ### ### ### ### ### ###

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

if (!requireNamespace("ggbreak", quietly = TRUE)) {
  stop(
    "Package 'ggbreak' is required to draw the broken x-axis in panel A.",
    call. = FALSE
  )
}


# Configuration -----------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- script_args[str_detect(script_args, "^--file=")]
base_path <- if (length(script_file) > 0) {
  script_file_path <- str_remove(script_file[[1]], "^--file=") |>
    str_replace_all("~\\+~", " ")

  dirname(normalizePath(script_file_path, mustWork = TRUE))
} else {
  getwd()
}

bri_input_path <- file.path(base_path, "inputs", "bri")
raw_data_path <- file.path(base_path, "inputs", "manylabs2")
output_path <- file.path(base_path, "outputs")

output_basename <- "BRI_ManyLabs2_by_technique"


# Shared aesthetics -------------------------------------------------------

bri_color <- c(
  main = "#4292C6",
  dark = "#111111",
  epm = "#2f52a3",
  mtt = "#1cb075",
  pcr = "#d95f02",
  ml2 = "#4292C6",
  grid_major = "#E8E8E8",
  grid_minor = "#F1F1F1",
  zero = "#111111"
)

bri_theme <- list(
  theme_bw(base_size = 9),
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.5, color = "black"),
    panel.grid.major = element_line(color = bri_color[["grid_major"]], linewidth = 0.3),
    panel.grid.minor = element_line(color = bri_color[["grid_minor"]], linewidth = 0.25),
    axis.line = element_blank(),
    axis.text = element_text(color = "#4D4D4D"),
    axis.title = element_text(color = "#222222"),
    strip.background = element_rect(fill = bri_color[["main"]], color = "black", linewidth = 0.5),
    strip.text = element_text(color = "white", face = "plain"),
    plot.title = element_blank()
  )
)


# Helpers -----------------------------------------------------------------

stop_missing_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      label, " not found at:\n  ", normalizePath(path, mustWork = FALSE), "\n",
      "This self-contained folder expects the input files under ./inputs/.",
      call. = FALSE
    )
  }
}

latest_results_path <- function(output_root = file.path(base_path, "output")) {
  if (!dir.exists(output_root)) {
    return(NULL)
  }

  candidates <- list.dirs(output_root, recursive = FALSE, full.names = FALSE)
  candidates <- candidates[str_detect(candidates, "^\\d{4}-\\d{2}-\\d{2}$")]

  if (length(candidates) == 0) {
    return(NULL)
  }

  sort(candidates, decreasing = TRUE)[[1]]
}

citation_only <- function(label) {
  match <- str_match(label, "\\(([^)]+)\\)\\s*$")
  if_else(is.na(match[, 2]), label, match[, 2])
}


# Panel A -----------------------------------------------------------------

technique_levels <- c("EPM", "MTT", "PCR")
legend_levels <- c("Original", "EPM", "MTT", "PCR", "Aggregate")

resolve_bri_analysis_path <- function(filename) {
  file.path(bri_input_path, filename)
}

resolve_bri_replication_path <- function() {
  resolve_bri_analysis_path("Replication Assessment by Replication.tsv")
}

resolve_bri_experiment_path <- function() {
  resolve_bri_analysis_path("Replication Assessment by Experiment.tsv")
}

bri_technique <- function(exp) {
  case_when(
    str_detect(exp, "^EPM") ~ "EPM",
    str_detect(exp, "^MTT") ~ "MTT",
    str_detect(exp, "^PCR") ~ "PCR",
    .default = NA_character_
  )
}

filter_bri_base_techniques <- function(data) {
  data |>
    filter(str_detect(EXP, "^(EPM|MTT|PCR)")) |>
    filter(!str_detect(EXP, "^(ALTMTT|ALTPCR)")) |>
    mutate(
      display_exp = EXP,
      technique = factor(bri_technique(EXP), levels = technique_levels)
    )
}

build_bri_panel_data <- function(replication_path, experiment_path, technique_filter = NULL) {
  stop_missing_file(replication_path, "BRI replication assessment by replication")
  stop_missing_file(experiment_path, "BRI replication assessment by experiment")

  if (!is.null(technique_filter) && !(technique_filter %in% technique_levels)) {
    stop(
      "technique_filter must be NULL, 'EPM', 'MTT', or 'PCR'.",
      call. = FALSE
    )
  }

  rep_raw <- readr::read_tsv(replication_path, show_col_types = FALSE)
  exp_raw <- readr::read_tsv(experiment_path, show_col_types = FALSE)

  required_rep_columns <- c("LAB", "EXP", "original_es", "log_es_ratio_individual")
  missing_rep_columns <- setdiff(required_rep_columns, colnames(rep_raw))

  if (length(missing_rep_columns) > 0) {
    stop(
      "BRI replication assessment by replication is missing required columns for panel A:\n  ",
      paste(missing_rep_columns, collapse = ", "),
      call. = FALSE
    )
  }

  required_exp_columns <- c("EXP", "corrected_sign_original_es", "corrected_sign_replication_es")
  missing_exp_columns <- setdiff(required_exp_columns, colnames(exp_raw))

  if (length(missing_exp_columns) > 0) {
    stop(
      "BRI replication assessment by experiment is missing required columns for panel A:\n  ",
      paste(missing_exp_columns, collapse = ", "),
      call. = FALSE
    )
  }

  rep_data <- rep_raw |>
    filter_bri_base_techniques() |>
    mutate(
      corrected_original_es = abs(original_es),
      corrected_replication_es = corrected_original_es - log_es_ratio_individual
    )

  exp_data <- exp_raw |>
    filter_bri_base_techniques() |>
    transmute(
      EXP,
      display_exp,
      technique,
      original_x = corrected_sign_original_es,
      aggregate_x = corrected_sign_replication_es
    )

  if (!is.null(technique_filter)) {
    rep_data <- rep_data |>
      filter(technique == technique_filter)
    exp_data <- exp_data |>
      filter(technique == technique_filter)
  }

  if (nrow(rep_data) == 0 || nrow(exp_data) == 0) {
    stop(
      "BRI replication assessment exists, but no EPM, MTT, or PCR rows were found.",
      call. = FALSE
    )
  }

  missing_experiments <- rep_data |>
    distinct(EXP) |>
    anti_join(exp_data |> distinct(EXP), by = join_by(EXP))

  if (nrow(missing_experiments) > 0) {
    stop(
      "BRI replication assessment by experiment is missing EXP values present in the by-replication file:\n  ",
      paste(missing_experiments$EXP, collapse = ", "),
      call. = FALSE
    )
  }

  replication_points <- rep_data |>
    mutate(
      point_type = "Replication",
      legend_group = as.character(technique),
      x = corrected_replication_es
    ) |>
    select(EXP, display_exp, technique, point_type, legend_group, x)

  original_points <- exp_data |>
    mutate(
      point_type = "Original",
      legend_group = "Original",
      x = original_x
    ) |>
    select(EXP, display_exp, technique, point_type, legend_group, x)

  aggregate_points <- exp_data |>
    mutate(
      point_type = "Aggregate",
      legend_group = "Aggregate",
      x = aggregate_x
    ) |>
    select(EXP, display_exp, technique, point_type, legend_group, x)

  exp_order <- aggregate_points |>
    arrange(desc(x)) |>
    pull(display_exp)

  bind_rows(original_points, replication_points, aggregate_points) |>
    mutate(
      display_exp = factor(display_exp, levels = rev(unique(exp_order))),
      point_type = factor(point_type, levels = c("Original", "Replication", "Aggregate")),
      technique = factor(technique, levels = technique_levels),
      legend_group = factor(legend_group, levels = legend_levels)
    )
}

order_bri_by_technique <- function(panel_data) {
  top_to_bottom <- panel_data |>
    filter(point_type == "Aggregate") |>
    arrange(technique, desc(x)) |>
    pull(display_exp) |>
    as.character()

  panel_data |>
    mutate(display_exp = factor(as.character(display_exp), levels = rev(unique(top_to_bottom))))
}

bri_technique_separator_yintercepts <- function(panel_data) {
  grouped_counts <- panel_data |>
    filter(point_type == "Aggregate") |>
    count(technique, name = "n") |>
    arrange(technique)

  if (nrow(grouped_counts) < 2) {
    return(NULL)
  }

  total_experiments <- sum(grouped_counts$n)
  separator_after <- cumsum(grouped_counts$n)[-nrow(grouped_counts)]

  total_experiments - separator_after + 0.5
}

add_bri_x_scale <- function(plot, panel_data) {
  if (max(panel_data$x, na.rm = TRUE) > 6) {
    plot +
      scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10)) +
      ggbreak::scale_x_break(c(5, 6), scales = 0.35, space = 0.05)
  } else {
    plot +
      scale_x_continuous(breaks = pretty_breaks(n = 5))
  }
}

build_bri_single_plot <- function(
  panel_data,
  show_legend = TRUE,
  show_x_title = TRUE,
  show_strip = FALSE,
  plot_title = NULL,
  compact = FALSE,
  legend_breaks = NULL,
  separator_yintercepts = NULL
) {
  original_points <- panel_data |> filter(point_type == "Original")
  replication_points <- panel_data |> filter(point_type == "Replication")
  aggregate_points <- panel_data |> filter(point_type == "Aggregate")

  if (is.null(legend_breaks)) {
    legend_breaks <- legend_levels[legend_levels %in% as.character(panel_data$legend_group)]
  }

  x_title <- if (show_x_title) "Effect size (log ROM)" else NULL
  plot_margin <- if (compact) {
    margin(t = 0.2, r = 4, b = 0.2, l = 4)
  } else {
    margin(t = 4, r = 4, b = 4, l = 4)
  }
  plot_title_margin <- if (compact) margin(b = 0.4) else margin(b = 2)
  plot_title_size <- if (compact) 8 else 8.5

  p <- ggplot(panel_data, aes(x = x, y = display_exp)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "#9A9A9A", linewidth = 0.35)

  if (!is.null(separator_yintercepts)) {
    p <- p +
      geom_hline(
        yintercept = separator_yintercepts,
        color = "#999999",
        linewidth = 0.45,
        linetype = "dashed"
      )
  }

  p <- p +
    geom_point(
      data = replication_points,
      aes(color = legend_group, shape = legend_group),
      size = 1.65,
      alpha = 0.5
    ) +
    geom_point(
      data = original_points,
      aes(color = legend_group, shape = legend_group),
      size = 2.2,
      alpha = 0.68
    ) +
    geom_point(
      data = aggregate_points,
      aes(color = legend_group, shape = legend_group),
      size = 3,
      alpha = 0.95,
      stroke = 0.8
    ) +
    scale_color_manual(
      values = c(
        "Original" = bri_color[["dark"]],
        "EPM" = bri_color[["epm"]],
        "MTT" = bri_color[["mtt"]],
        "PCR" = bri_color[["pcr"]],
        "Aggregate" = bri_color[["dark"]]
      ),
      breaks = legend_breaks,
      limits = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Original" = 16,
        "EPM" = 16,
        "MTT" = 16,
        "PCR" = 16,
        "Aggregate" = 124
      ),
      breaks = legend_breaks,
      limits = legend_levels,
      drop = FALSE,
      name = NULL
    )

  if (show_strip) {
    p <- p +
      facet_wrap(~technique, ncol = 1, scales = "free_y")
  }

  p <- add_bri_x_scale(p, panel_data)

  p +
    labs(
      title = plot_title,
      x = x_title,
      y = NULL
    ) +
    bri_theme +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = c("Original" = 16, "EPM" = 16, "MTT" = 16, "PCR" = 16, "Aggregate" = 124)[legend_breaks],
          size = c("Original" = 2.2, "EPM" = 2.0, "MTT" = 2.0, "PCR" = 2.0, "Aggregate" = 3)[legend_breaks],
          alpha = c("Original" = 0.68, "EPM" = 0.5, "MTT" = 0.5, "PCR" = 0.5, "Aggregate" = 0.95)[legend_breaks]
        )
      ),
      shape = "none"
    ) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      legend.text = element_text(size = 7.5),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.spacing.x = unit(0.12, "cm"),
      axis.text.y = element_text(size = 5.8),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 8.5),
      plot.title = element_text(
        size = plot_title_size,
        face = "plain",
        color = "#4D4D4D",
        hjust = 0,
        margin = plot_title_margin
      ),
      plot.margin = plot_margin
    )
}

build_blank_plot <- function() {
  cowplot::ggdraw()
}

build_bri_panel_by_technique_parts <- function(replication_path, experiment_path) {
  all_panel_data <- build_bri_panel_data(replication_path, experiment_path)
  shared_legend <- cowplot::get_legend(
    build_bri_single_plot(
      all_panel_data,
      show_legend = TRUE,
      show_x_title = TRUE,
      show_strip = FALSE,
      legend_breaks = legend_levels
    )
  )

  panel_plots <- lapply(technique_levels, function(technique_name) {
    panel_data <- build_bri_panel_data(
      replication_path,
      experiment_path,
      technique_filter = technique_name
    )

    build_bri_single_plot(
      panel_data,
      show_legend = FALSE,
      show_x_title = technique_name == "PCR",
      show_strip = FALSE,
      plot_title = technique_name,
      compact = TRUE,
      legend_breaks = c("Original", technique_name, "Aggregate")
    )
  })

  stacked_panels <- cowplot::plot_grid(
    panel_plots[[1]],
    build_blank_plot(),
    panel_plots[[2]],
    build_blank_plot(),
    panel_plots[[3]],
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(14, 0.45, 17, 0.45, 15.5)
  )

  list(
    panel = stacked_panels,
    legend = shared_legend
  )
}

build_bri_panel_by_technique <- function(replication_path, experiment_path) {
  parts <- build_bri_panel_by_technique_parts(replication_path, experiment_path)

  bottom_row <- cowplot::plot_grid(
    parts$legend,
    ncol = 1,
    rel_heights = c(1)
  )

  cowplot::plot_grid(
    parts$panel,
    bottom_row,
    ncol = 1,
    rel_heights = c(1, 0.035)
  )
}

build_bri_panel <- function(
  replication_path,
  experiment_path,
  mode = c("combined", "by_technique", "combined_grouped_by_technique"),
  technique_filter = NULL
) {
  mode <- match.arg(mode)

  if (identical(mode, "by_technique")) {
    if (!is.null(technique_filter)) {
      stop(
        "technique_filter must be NULL when mode is 'by_technique'.",
        call. = FALSE
      )
    }

    return(build_bri_panel_by_technique(replication_path, experiment_path))
  }

  if (identical(mode, "combined_grouped_by_technique") && !is.null(technique_filter)) {
    stop(
      "technique_filter must be NULL when mode is 'combined_grouped_by_technique'.",
      call. = FALSE
    )
  }

  panel_data <- build_bri_panel_data(
    replication_path,
    experiment_path,
    technique_filter = technique_filter
  )

  separator_yintercepts <- NULL
  if (identical(mode, "combined_grouped_by_technique")) {
    panel_data <- order_bri_by_technique(panel_data)
    separator_yintercepts <- bri_technique_separator_yintercepts(panel_data)
  }

  legend_breaks <- if (is.null(technique_filter)) {
    legend_levels
  } else {
    c("Original", technique_filter, "Aggregate")
  }

  build_bri_single_plot(
    panel_data,
    show_legend = TRUE,
    show_x_title = TRUE,
    show_strip = FALSE,
    legend_breaks = legend_breaks,
    separator_yintercepts = separator_yintercepts
  )
}


# Panel B -----------------------------------------------------------------

build_ml2_panel_data <- function() {
  sums_path <- file.path(raw_data_path, "Data_Figure_NOweird_sums.csv")
  means_path <- file.path(raw_data_path, "Data_Figure_NOweird_means.csv")
  original_path <- file.path(raw_data_path, "Data_Figure_NOweird_oriEffects.csv")

  stop_missing_file(sums_path, "ManyLabs2 site data")
  stop_missing_file(means_path, "ManyLabs2 mean data")
  stop_missing_file(original_path, "ManyLabs2 original-effects data")

  ml2_sums <- readr::read_csv(sums_path, show_col_types = FALSE)
  ml2_means <- readr::read_csv(means_path, show_col_types = FALSE)
  ml2_original <- readr::read_csv(original_path, show_col_types = FALSE)

  effect_type_map <- ml2_original |>
    distinct(study.labels, esType) |>
    mutate(
      study.labels = str_trim(study.labels),
      esType = str_trim(esType)
    )

  cohens_q_labels <- effect_type_map |>
    filter(esType == "Z.f") |>
    pull(study.labels) |>
    unique()

  original_parsed <- ml2_original |>
    mutate(
      study.labels = str_trim(study.labels),
      usethisES = as.numeric(usethisES),
      oriWEIRD = str_trim(oriWEIRD)
    ) |>
    filter(!is.na(usethisES))

  original_main <- original_parsed |>
    slice_head(n = 1, by = study.labels) |>
    select(study.labels, usethisES)

  original_weird <- original_parsed |>
    filter(oriWEIRD == "WEIRD") |>
    select(study.labels, es_weird = usethisES)

  original_nonweird <- original_parsed |>
    filter(oriWEIRD == "less WEIRD") |>
    select(study.labels, es_nonweird = usethisES)

  split_studies <- intersect(original_weird$study.labels, original_nonweird$study.labels)

  site_data <- ml2_sums |>
    mutate(
      study.labels = str_trim(study.labels),
      loc = as.numeric(loc)
    ) |>
    filter(!is.na(loc))

  mean_data <- ml2_means |>
    mutate(
      study.labels = str_trim(study.labels),
      loc = as.numeric(loc)
    ) |>
    filter(!is.na(loc))

  r_studies <- mean_data |>
    filter(!(study.labels %in% cohens_q_labels)) |>
    arrange(loc) |>
    pull(study.labels)

  q_studies <- mean_data |>
    filter(study.labels %in% cohens_q_labels) |>
    arrange(loc) |>
    pull(study.labels)

  all_studies <- c(r_studies, q_studies)
  cite_labels <- citation_only(all_studies)
  names(cite_labels) <- all_studies

  df_sites <- site_data |>
    filter(study.labels %in% all_studies) |>
    mutate(study_f = factor(study.labels, levels = all_studies)) |>
    select(study_f, x = loc)

  df_means <- mean_data |>
    filter(study.labels %in% all_studies) |>
    mutate(study_f = factor(study.labels, levels = all_studies)) |>
    select(study_f, x = loc)

  df_original_single <- original_main |>
    filter(
      study.labels %in% all_studies,
      !(study.labels %in% split_studies)
    ) |>
    mutate(
      study_f = factor(study.labels, levels = all_studies),
      x = usethisES
    ) |>
    select(study_f, x)

  df_original_weird <- original_weird |>
    filter(study.labels %in% split_studies) |>
    mutate(
      study_f = factor(study.labels, levels = all_studies),
      x = es_weird
    ) |>
    select(study_f, x)

  df_original_nonweird <- original_nonweird |>
    filter(study.labels %in% split_studies) |>
    mutate(
      study_f = factor(study.labels, levels = all_studies),
      x = es_nonweird
    ) |>
    select(study_f, x)

  list(
    sites = df_sites,
    means = df_means,
    original_filled = bind_rows(df_original_single, df_original_weird),
    original_open = df_original_nonweird,
    all_studies = all_studies,
    cite_labels = cite_labels,
    n_r = length(r_studies),
    n_total = length(all_studies)
  )
}

ml2_legend_levels <- c("Original", "Original (non-WEIRD)", "Replication", "Aggregate")

build_ml2_panel <- function(show_x_title = TRUE, show_legend = FALSE) {
  ml2 <- build_ml2_panel_data()
  x_title <- if (show_x_title) "Effect size (r)" else NULL

  # Add legend_group column to each data frame for unified legend
  sites_plot <- ml2$sites |> mutate(legend_group = "Replication")
  means_plot <- ml2$means |> mutate(legend_group = "Aggregate")
  ori_filled_plot <- ml2$original_filled |> mutate(legend_group = "Original")
  ori_open_plot <- ml2$original_open |> mutate(legend_group = "Original (non-WEIRD)")

  ggplot() +
    geom_hline(
      yintercept = seq_len(ml2$n_total),
      color = bri_color[["grid_major"]],
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = 0,
      color = "#777777",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = ml2$n_r + 0.5,
      color = "#999999",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = sites_plot,
      aes(x = x, y = as.numeric(study_f), color = legend_group,
          shape = legend_group, fill = legend_group),
      size = 0.85,
      alpha = 0.38,
      position = position_jitter(width = 0, height = 0.13, seed = 42)
    ) +
    geom_point(
      data = means_plot,
      aes(x = x, y = as.numeric(study_f), color = legend_group,
          shape = legend_group, fill = legend_group),
      size = 3,
      stroke = 0.8
    ) +
    geom_point(
      data = ori_filled_plot,
      aes(x = x, y = as.numeric(study_f), color = legend_group,
          shape = legend_group, fill = legend_group),
      size = 1.8,
      alpha = 0.82
    ) +
    geom_point(
      data = ori_open_plot,
      aes(x = x, y = as.numeric(study_f), color = legend_group,
          shape = legend_group, fill = legend_group),
      size = 1.9,
      stroke = 0.65
    ) +
    scale_color_manual(
      values = c(
        "Original"            = bri_color[["dark"]],
        "Original (non-WEIRD)" = bri_color[["dark"]],
        "Replication"         = bri_color[["ml2"]],
        "Aggregate"           = bri_color[["dark"]]
      ),
      breaks = ml2_legend_levels,
      limits = ml2_legend_levels,
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Original"            = 16,
        "Original (non-WEIRD)" = 21,
        "Replication"         = 16,
        "Aggregate"           = 124
      ),
      breaks = ml2_legend_levels,
      limits = ml2_legend_levels,
      name = NULL
    ) +
    scale_fill_manual(
      values = c(
        "Original"            = bri_color[["dark"]],
        "Original (non-WEIRD)" = "white",
        "Replication"         = bri_color[["ml2"]],
        "Aggregate"           = bri_color[["dark"]]
      ),
      breaks = ml2_legend_levels,
      limits = ml2_legend_levels,
      name = NULL
    ) +
    scale_y_continuous(
      breaks = seq_len(ml2$n_total),
      labels = ml2$cite_labels[ml2$all_studies],
      expand = expansion(add = 0.6)
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      sec.axis = dup_axis(name = "Effect size (Cohen's q)")
    ) +
    labs(x = x_title, y = NULL) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = c(16, 21, 16, 124),
          fill  = c(bri_color[["dark"]], "white", bri_color[["ml2"]], bri_color[["dark"]]),
          color = c(bri_color[["dark"]], bri_color[["dark"]], bri_color[["ml2"]], bri_color[["dark"]]),
          size  = c(2.2, 2.2, 1.6, 3),
          alpha = c(0.82, 0.82, 0.5, 0.95),
          stroke = c(NA, 0.65, NA, 0.8)
        )
      ),
      shape = "none",
      fill  = "none"
    ) +
    bri_theme +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = 7.5),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.spacing.x = unit(0.12, "cm"),
      panel.grid.major.x = element_line(color = bri_color[["grid_major"]], linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 7.4),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 8.5),
      axis.title.x.top = element_text(size = 8.5),
      plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
    )
}


# Output ------------------------------------------------------------------

build_combined_plot <- function(p_bri, p_ml2) {
  cowplot::plot_grid(
    p_bri,
    p_ml2,
    labels = c("A", "B"),
    nrow = 1,
    rel_widths = c(0.95, 1.35),
    label_size = 18,
    label_fontface = "plain",
    align = "h",
    axis = "tb"
  )
}

build_combined_by_technique_plot <- function(bri_parts, p_ml2) {
  rel_widths <- c(bri = 0.95, ml2 = 1.35)
  outer_x <- 0.025
  gap_x <- 0.025
  usable_width <- 1 - (outer_x * 2) - gap_x
  bri_width <- usable_width * rel_widths[["bri"]] / sum(rel_widths)
  ml2_width <- usable_width * rel_widths[["ml2"]] / sum(rel_widths)
  bri_x <- outer_x
  ml2_x <- bri_x + bri_width + gap_x

  main_y <- 0.075
  main_height <- 0.89
  bri_y <- main_y
  legend_y <- 0.028
  legend_height <- 0.040

  # Extract the panel B legend from a version of the ML2 plot with legend enabled
  ml2_legend <- cowplot::get_legend(
    build_ml2_panel(show_x_title = TRUE, show_legend = TRUE)
  )

  cowplot::ggdraw() +
    cowplot::draw_plot(
      bri_parts$panel,
      x = bri_x,
      y = bri_y,
      width = bri_width,
      height = main_height
    ) +
    cowplot::draw_plot(
      p_ml2,
      x = ml2_x,
      y = main_y,
      width = ml2_width,
      height = main_height
    ) +
    cowplot::draw_plot(
      bri_parts$legend,
      x = bri_x,
      y = legend_y,
      width = bri_width,
      height = legend_height
    ) +
    cowplot::draw_plot(
      ml2_legend,
      x = ml2_x,
      y = legend_y,
      width = ml2_width,
      height = legend_height
    ) +
    cowplot::draw_label(
      "A",
      x = 0.006,
      y = 0.99,
      hjust = 0,
      vjust = 1,
      size = 18,
      fontface = "plain"
    ) +
    cowplot::draw_label(
      "B",
      x = ml2_x - 0.018,
      y = 0.99,
      hjust = 0,
      vjust = 1,
      size = 18,
      fontface = "plain"
    )
}

save_combined_figure <- function(plot, basename, width = 14, height = 8.5) {
  png_path <- file.path(output_path, paste0(basename, ".png"))
  pdf_path <- file.path(output_path, paste0(basename, ".pdf"))

  ggsave(
    png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    pdf_path,
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )

  cat("PNG: ", normalizePath(png_path, mustWork = FALSE), "\n", sep = "")
  cat("PDF: ", normalizePath(pdf_path, mustWork = FALSE), "\n", sep = "")
}

cat("-- Resolving BRI replication data --\n")
bri_replication_path <- resolve_bri_replication_path()
bri_experiment_path <- resolve_bri_experiment_path()

cat("-- Building ManyLabs2 panel B --\n")
p_ml2 <- build_ml2_panel()

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

cat("-- Building BRI panel A: by_technique --\n")
bri_parts <- build_bri_panel_by_technique_parts(
  bri_replication_path,
  bri_experiment_path
)

combined_plot <- build_combined_by_technique_plot(
  bri_parts,
  p_ml2
)

cat("-- Saving combined figure: by_technique --\n")
save_combined_figure(
  combined_plot,
  output_basename,
  width = 14,
  height = 11
)

cat("-- Done --\n")
