############################################################
# Quantify associations between spatial predictors and AIV case counts using Negative binomial regression 
#
# Predictors:
# 1. Wild bird-domestic animal overlap
# 2. Wild bird richness
# 3. Domestic animal density
#
# Workflow:
# 1. Load raster predictors and AIV case layers
# 2. Extract spatial values
# 3. Aggregate predictors into equal-width bins
# 4. Calculate q95 case counts within bins
# 5. Fit negative binomial regression models
# 6. Calculate effect sizes and generate prediction curves
#
############################################################

# 1. Load packages ---------------------------------------------------------

library(raster)
library(dplyr)
library(ggplot2)

# 2. Define paths ----------------------------------------------------------

predictor_dir <- file.path("data", "input", "predictors")
case_dir <- file.path("data", "input", "aiv_cases")
output_dir <- file.path("data", "output", "negative_binomial")
figure_dir <- "figures"

invisible(lapply(
  c(output_dir, figure_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

predictor_specs <- list(
  overlap = list(
    label = "Wild bird-domestic animal overlap",
    breeding = file.path(
      predictor_dir,
      "Overlap_Breeding_hostConfirmed_host_cumulative_Mean.tif"
    ),
    nonbreeding = file.path(
      predictor_dir,
      "Overlap_NonBreeding_hostConfirmed_host_cumulative_Mean.tif"
    )
  ),
  wild_bird_richness = list(
    label = "Wild bird richness",
    breeding = file.path(
      predictor_dir,
      paste0(
        "Season_North_Breeding_hostConfirmed_Clipped_fill0_Extended_",
        "QuantileNormalization.tif"
      )
    ),
    nonbreeding = file.path(
      predictor_dir,
      paste0(
        "Season_North_NonBreeding_hostConfirmed_Clipped_fill0_Extended_",
        "QuantileNormalization.tif"
      )
    )
  ),
  domestic_animal_density = list(
    label = "Domestic animal density",
    breeding = file.path(
      predictor_dir,
      "Host_Resampled_QuantileNormalization_Mean.tif"
    ),
    nonbreeding = file.path(
      predictor_dir,
      "Host_Resampled_QuantileNormalization_Mean.tif"
    )
  )
)

case_specs <- list(
  poultry_domestic_mammals = list(
    label = "Poultry and domestic mammals",
    summer = file.path(
      case_dir,
      "AIV_cases_count_PoultryDomesticmammals_Summer.tif"
    ),
    winter = file.path(
      case_dir,
      "AIV_cases_count_PoultryDomesticmammals_Winter.tif"
    )
  )
)

required_paths <- unique(c(
  unlist(lapply(predictor_specs, function(x) c(x$breeding, x$nonbreeding))),
  unlist(lapply(case_specs, function(x) c(x$summer, x$winter)))
))
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_paths, collapse = "\n")
  )
}

# 3. Load input data -------------------------------------------------------

predictor_rasters <- lapply(predictor_specs, function(spec) {
  list(
    breeding = raster(spec$breeding),
    nonbreeding = raster(spec$nonbreeding)
  )
})

case_rasters <- lapply(case_specs, function(spec) {
  list(
    summer = raster(spec$summer),
    winter = raster(spec$winter)
  )
})

analysis_grid <- do.call(
  rbind,
  lapply(names(predictor_specs), function(predictor_id) {
    data.frame(
      predictor_id = predictor_id,
      response_id = names(case_specs),
      stringsAsFactors = FALSE
    )
  })
)
analysis_grid$analysis_id <- paste(
  analysis_grid$predictor_id,
  analysis_grid$response_id,
  sep = "__"
)
analysis_grid$predictor_label <- vapply(
  analysis_grid$predictor_id,
  function(id) predictor_specs[[id]]$label,
  character(1)
)
analysis_grid$response_label <- vapply(
  analysis_grid$response_id,
  function(id) case_specs[[id]]$label,
  character(1)
)

# 4. Extract raster values -------------------------------------------------

assert_aligned <- function(predictor_raster, case_raster, analysis_id, season) {
  aligned <- compareRaster(
    predictor_raster,
    case_raster,
    extent = TRUE,
    rowcol = TRUE,
    crs = TRUE,
    res = TRUE,
    orig = TRUE,
    rotation = TRUE,
    values = FALSE,
    stopiffalse = FALSE
  )

  if (!aligned) {
    stop("Raster geometry mismatch for ", analysis_id, " (", season, ").")
  }
}

extract_raster_values <- function(predictors, cases, analysis_id) {
  assert_aligned(predictors$breeding, cases$summer, analysis_id, "summer")
  assert_aligned(
    predictors$nonbreeding,
    cases$winter,
    analysis_id,
    "winter"
  )

  data.frame(
    breeding_value = getValues(predictors$breeding),
    nonbreeding_value = getValues(predictors$nonbreeding),
    summer_cases = getValues(cases$summer),
    winter_cases = getValues(cases$winter)
  )
}

extracted_values <- setNames(
  lapply(seq_len(nrow(analysis_grid)), function(i) {
    predictor_id <- analysis_grid$predictor_id[i]
    response_id <- analysis_grid$response_id[i]
    extract_raster_values(
      predictor_rasters[[predictor_id]],
      case_rasters[[response_id]],
      analysis_grid$analysis_id[i]
    )
  }),
  analysis_grid$analysis_id
)

# 5. Aggregate predictor bins and calculate q95 cases ---------------------

aggregate_season <- function(
    predictor_values,
    case_values,
    season,
    bin_width = 0.05) {
  season_data <- data.frame(
    predictor_value = predictor_values,
    case_value = case_values
  ) %>%
    filter(!is.na(predictor_value) & !is.na(case_value))

  season_data %>%
    mutate(
      group = cut(
        predictor_value,
        breaks = seq(
          0,
          max(predictor_value, na.rm = TRUE),
          by = bin_width
        ),
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    group_by(group) %>%
    summarise(
      q95_aivcases = quantile(case_value, 0.95, na.rm = TRUE),
      count_aivcases = n(),
      .groups = "drop"
    ) %>%
    mutate(
      group_numeric = as.numeric(
        as.character(gsub("\\[|\\)|,.*", "", group))
      ),
      Season = season
    )
}

aggregate_predictor_bins <- function(values) {
  bind_rows(
    aggregate_season(
      values$breeding_value,
      values$summer_cases,
      "Boreal summer"
    ),
    aggregate_season(
      values$nonbreeding_value,
      values$winter_cases,
      "Boreal winter"
    )
  )
}

binned_data <- lapply(extracted_values, aggregate_predictor_bins)

# 6. Fit negative binomial models ----------------------------------------

fit_negative_binomial_models <- function(data) {
  model_data <- data %>%
    filter(!is.na(group_numeric) & !is.na(q95_aivcases))

  season_data <- split(model_data, model_data$Season)
  lapply(season_data, function(x) {
    MASS::glm.nb(q95_aivcases ~ group_numeric, data = x)
  })
}

extract_model_statistics <- function(models, data) {
  bind_rows(lapply(names(models), function(season) {
    model <- models[[season]]
    season_data <- data %>%
      filter(
        Season == season,
        !is.na(group_numeric),
        !is.na(q95_aivcases)
      )
    model_summary <- summary(model)

    tibble(
      Season = season,
      n = nrow(season_data),
      beta = unname(coef(model)[2]),
      se_beta = model_summary$coefficients[2, 2],
      p_value_nb = model_summary$coefficients[2, 4],
      theta = model$theta,
      AIC = AIC(model)
    )
  }))
}

nb_models <- lapply(binned_data, fit_negative_binomial_models)
model_statistics <- Map(extract_model_statistics, nb_models, binned_data)

# 7. Calculate effect sizes ------------------------------------------------

calculate_effect_sizes <- function(statistics, data) {
  effect_ranges <- bind_rows(lapply(statistics$Season, function(season) {
    season_data <- data %>%
      filter(
        Season == season,
        !is.na(group_numeric),
        !is.na(q95_aivcases)
      )

    tibble(
      Season = season,
      x_p05 = as.numeric(quantile(
        season_data$group_numeric,
        0.05,
        na.rm = TRUE
      )),
      x_p95 = as.numeric(quantile(
        season_data$group_numeric,
        0.95,
        na.rm = TRUE
      ))
    )
  }))

  statistics %>%
    left_join(effect_ranges, by = "Season") %>%
    mutate(
      expected_case_ratio_unit = exp(beta),
      percent_change_unit = (expected_case_ratio_unit - 1) * 100,
      delta_x_p05_p95 = x_p95 - x_p05,
      expected_case_ratio_p05_p95 = exp(beta * delta_x_p05_p95),
      percent_change_p05_p95 =
        (expected_case_ratio_p05_p95 - 1) * 100
    )
}

effect_results <- Map(
  calculate_effect_sizes,
  model_statistics,
  binned_data
)

# 8. Generate predictions and confidence intervals -----------------------

generate_predictions <- function(models, data) {
  bind_rows(lapply(names(models), function(season) {
    model <- models[[season]]
    season_data <- filter(data, Season == season)
    newdata <- data.frame(
      group_numeric = seq(
        min(season_data$group_numeric, na.rm = TRUE),
        max(season_data$group_numeric, na.rm = TRUE),
        length.out = 200
      )
    )
    prediction <- predict(
      model,
      newdata = newdata,
      type = "link",
      se.fit = TRUE
    )

    newdata %>%
      mutate(
        Season = season,
        fit_link = prediction$fit,
        se_link = prediction$se.fit,
        fit = exp(fit_link),
        lower = exp(fit_link - 1.96 * se_link),
        upper = exp(fit_link + 1.96 * se_link)
      )
  }))
}

prediction_data <- Map(generate_predictions, nb_models, binned_data)

# 9. Create figures --------------------------------------------------------

get_significance_code <- function(p_value) {
  case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
}

make_effect_label <- function(value, significance_code) {
  if (significance_code == "") {
    sprintf("exp(beta %%*%% (P[95] - P[5])) == %.2f", value)
  } else {
    sprintf(
      "exp(beta %%*%% (P[95] - P[5])) == %.2f * '%s'",
      value,
      significance_code
    )
  }
}

effect_results <- lapply(effect_results, function(results) {
  results %>%
    mutate(
      significance_code = get_significance_code(p_value_nb),
      linetype = ifelse(p_value_nb < 0.05, "solid", "dashed"),
      label_parse = mapply(
        make_effect_label,
        expected_case_ratio_p05_p95,
        significance_code
      )
    )
})

prediction_data <- Map(function(predictions, results) {
  predictions %>%
    left_join(
      dplyr::select(results, Season, linetype),
      by = "Season"
    )
}, prediction_data, effect_results)

season_colours <- c(
  "Boreal summer" = "#A50F15",
  "Boreal winter" = "#08519C"
)
season_fills <- c(
  "Boreal summer" = "#fc927250",
  "Boreal winter" = "#9ecae150"
)

make_nb_plot <- function(data, predictions, results, x_label) {
  summer_label <- results %>%
    filter(Season == "Boreal summer") %>%
    pull(label_parse)
  winter_label <- results %>%
    filter(Season == "Boreal winter") %>%
    pull(label_parse)

  ggplot(data, aes(x = group_numeric, y = q95_aivcases)) +
    geom_point(
      aes(colour = Season),
      shape = 16,
      alpha = 0.5,
      size = 0.5
    ) +
    geom_ribbon(
      data = predictions,
      aes(
        x = group_numeric,
        ymin = lower,
        ymax = upper,
        fill = Season
      ),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    geom_line(
      data = predictions,
      aes(
        x = group_numeric,
        y = fit,
        colour = Season,
        linetype = linetype
      ),
      linewidth = 0.5
    ) +
    annotate(
      "text",
      x = 0.05,
      y = 57,
      label = summer_label,
      hjust = 0,
      size = 2,
      colour = season_colours[["Boreal summer"]],
      parse = TRUE
    ) +
    annotate(
      "text",
      x = 0.05,
      y = 51,
      label = winter_label,
      hjust = 0,
      size = 2,
      colour = season_colours[["Boreal winter"]],
      parse = TRUE
    ) +
    scale_colour_manual(name = "Season", values = season_colours) +
    scale_fill_manual(name = "Season", values = season_fills) +
    scale_linetype_identity() +
    scale_y_continuous(
      limits = c(0, 60),
      breaks = c(0, 20, 40, 60),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    labs(
      x = x_label,
      y = "AIV cases (domestic animal)"
    ) +
    theme_bw(base_family = "Arial") +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(colour = "black", linewidth = 0.2),
      axis.title = element_text(colour = "black", size = 7),
      axis.text = element_text(colour = "black", size = 7),
      axis.ticks = element_line(colour = "black", linewidth = 0.2),
      plot.margin = margin(t = 5, r = 10, b = 5, l = 0)
    )
}

plot_objects <- setNames(
  lapply(seq_len(nrow(analysis_grid)), function(i) {
    analysis_id <- analysis_grid$analysis_id[i]
    make_nb_plot(
      binned_data[[analysis_id]],
      prediction_data[[analysis_id]],
      effect_results[[analysis_id]],
      analysis_grid$predictor_label[i]
    )
  }),
  analysis_grid$analysis_id
)

combined_plot <- cowplot::plot_grid(
  plotlist = plot_objects,
  ncol = length(case_specs),
  align = "v",
  axis = "lr"
)

# 10. Export results -------------------------------------------------------

add_analysis_metadata <- function(data, analysis_id) {
  metadata <- analysis_grid[
    analysis_grid$analysis_id == analysis_id,
    ,
    drop = FALSE
  ]
  data %>%
    mutate(
      analysis_id = analysis_id,
      predictor_id = metadata$predictor_id,
      response_id = metadata$response_id,
      predictor_label = metadata$predictor_label,
      response_label = metadata$response_label,
      .before = 1
    )
}

binned_table <- bind_rows(Map(
  add_analysis_metadata,
  binned_data,
  names(binned_data)
))
regression_table <- bind_rows(Map(
  add_analysis_metadata,
  effect_results,
  names(effect_results)
))
prediction_table <- bind_rows(Map(
  add_analysis_metadata,
  prediction_data,
  names(prediction_data)
))

figure_stem <- file.path(
  figure_dir,
  "AIVCases_PoultryDomesticmammals_3Predictors_NegativeBinomial"
)
figure_width_cm <- 4.5 * length(case_specs)
figure_height_cm <- 4 * length(predictor_specs)

ggsave(
  filename = paste0(figure_stem, ".png"),
  plot = combined_plot,
  width = figure_width_cm,
  height = figure_height_cm,
  units = "cm",
  dpi = 600,
  bg = "white"
)

svglite::svglite(
  paste0(figure_stem, ".svg"),
  width = figure_width_cm / 2.54,
  height = figure_height_cm / 2.54
)
print(combined_plot)
grDevices::dev.off()

grDevices::cairo_pdf(
  paste0(figure_stem, ".pdf"),
  width = figure_width_cm / 2.54,
  height = figure_height_cm / 2.54,
  family = "Arial"
)
print(combined_plot)
grDevices::dev.off()
