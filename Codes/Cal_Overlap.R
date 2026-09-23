############################################################
# Wild bird-domestic animal overlap
#
# Workflow:
# 1. Load input raster datasets
# 2. Harmonize spatial extent and resolution
# 3. Quantile-normalize raster values
# 4. Calculate seasonal and annual overlap
# 5. Export raster results and figures
#
############################################################

# 1. Load packages ---------------------------------------------------------

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(svglite)

# 2. Define input/output paths --------------------------------------------

domestic_input_dir <- file.path("data", "input", "domestic_animals")
bird_input_dir <- file.path("data", "input", "wild_birds")
resampled_dir <- file.path("data", "output", "resampled")
normalized_dir <- file.path("data", "output", "normalized")
overlap_dir <- file.path("data", "output", "overlap")
figure_dir <- "figures"

output_dirs <- c(resampled_dir, normalized_dir, overlap_dir, figure_dir)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

domestic_paths <- c(
  chicken_2010 = file.path(domestic_input_dir, "5_Ch_2015_Da.tif"),
  duck_2010 = file.path(domestic_input_dir, "5_Dk_2015_Da.tif"),
  chicken_extensive = file.path(domestic_input_dir, "ChExtDn_8k_201507081.tif"),
  chicken_intensive = file.path(domestic_input_dir, "ChIntDn_8k_201507081.tif"),
  cattle_2010 = file.path(domestic_input_dir, "GLW4-2020.D-DA.CTL.tif"),
  pigs_2010 = file.path(domestic_input_dir, "GLW4-2020.D-DA.PGS.tif")
)

bird_paths <- c(
  breeding = file.path(
    bird_input_dir,
    "Season_North_Breeding_hostAll_Clipped_fill0.tif"
  ),
  nonbreeding = file.path(
    bird_input_dir,
    "Season_North_NonBreeding_hostAll_Clipped_fill0.tif"
  ),
  annual = file.path(
    bird_input_dir,
    "Season_North_Global_hostAll_OneTimeCounts_Clipped_fill0.tif"
  )
)

required_paths <- c(domestic_paths, bird_paths)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_paths, collapse = "\n")
  )
}

# 3. Load data -------------------------------------------------------------

domestic_rasters <- lapply(domestic_paths, rast)
bird_rasters <- lapply(bird_paths, rast)

land_medium <- ne_download(
  scale = "medium",
  type = "land",
  category = "physical",
  returnclass = "sf"
)

# 4. Data preprocessing ----------------------------------------------------

fill_land_na <- function(raster_layer, land_mask) {
  raster_layer[is.na(raster_layer) & !is.na(land_mask)] <- 0
  raster_layer
}

align_extent <- function(raster_layer, template) {
  raster_layer <- extend(raster_layer, ext(template))
  crop(raster_layer, ext(template))
}

quantile_normalize <- function(raster_layer) {
  raster_values <- values(raster_layer)
  nonzero_values <- raster_values[
    raster_values > 0 & !is.na(raster_values)
  ]

  if (length(nonzero_values) == 0) {
    return(raster_layer)
  }

  quantiles <- quantile(
    nonzero_values,
    probs = seq(0, 1, length.out = 101),
    na.rm = TRUE,
    names = FALSE
  )
  normalized_values <- approx(
    quantiles,
    seq(0, 1, length.out = 101),
    xout = raster_values,
    rule = 2
  )$y
  normalized_values[raster_values == 0] <- 0

  normalized_raster <- raster_layer
  values(normalized_raster) <- normalized_values
  normalized_raster
}

land_mask_domestic <- rasterize(land_medium, domestic_rasters$chicken_2010, field = 1)
domestic_rasters <- lapply(
  domestic_rasters,
  fill_land_na,
  land_mask = land_mask_domestic
)

bird_rasters <- lapply(
  bird_rasters,
  align_extent,
  template = land_mask_domestic
)

domestic_resampled <- lapply(
  domestic_rasters,
  resample,
  y = bird_rasters$breeding,
  method = "bilinear"
)

resampled_rasters <- c(
  domestic_resampled,
  list(
    season_north_breeding_hostAll_Clipped_fill0_Extended = bird_rasters$breeding,
    season_north_nonbreeding_hostAll_fill0_Extended = bird_rasters$nonbreeding,
    Season_North_Global_hostAll_OneTimeCounts_Clipped_fill0_Extended = bird_rasters$annual
  )
)

normalized_rasters <- lapply(resampled_rasters, quantile_normalize)

host_layer_names <- c(
  "chicken_extensive",
  "chicken_intensive",
  "duck_2010",
  "cattle_2010",
  "pigs_2010"
)
host_stack <- rast(normalized_rasters[host_layer_names])
host_mean <- app(host_stack, function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  x[is.na(x)] <- 0
  sum(x) / length(host_layer_names)
})

# 5. Overlap calculation --------------------------------------------------

breeding_overlap <-
  normalized_rasters$season_north_breeding_hostAll_Clipped_fill0_Extended *
  host_mean
nonbreeding_overlap <-
  normalized_rasters$season_north_nonbreeding_hostAll_fill0_Extended *
  host_mean
annual_overlap <-
  normalized_rasters$Season_North_Global_hostAll_OneTimeCounts_Clipped_fill0_Extended *
  host_mean
delta_overlap <- nonbreeding_overlap - breeding_overlap


# 6. Export results --------------------------------------------------------

resampled_stems <- c(
  chicken_2010 = "chicken_2010_Resampled",
  duck_2010 = "duck_2010_Resampled",
  chicken_extensive = "chicken_extensive_Resampled",
  chicken_intensive = "chicken_intensive_Resampled",
  cattle_2010 = "cattle_2010_Resampled",
  pigs_2010 = "pigs_2010_Resampled",
  season_north_breeding_hostAll_Clipped_fill0_Extended =
    "season_north_breeding_hostAll_Clipped_fill0_Extended",
  season_north_nonbreeding_hostAll_fill0_Extended =
    "season_north_nonbreeding_hostAll_fill0_Extended",
  Season_North_Global_hostAll_OneTimeCounts_Clipped_fill0_Extended =
    "Season_North_Global_hostAll_OneTimeCounts_Clipped_fill0_Extended"
)
resampled_paths <- file.path(
  resampled_dir,
  paste0(unname(resampled_stems[names(resampled_rasters)]), ".tif")
)
normalized_paths <- file.path(
  normalized_dir,
  paste0(
    unname(resampled_stems[names(normalized_rasters)]),
    "_QuantileNormalization.tif"
  )
)

invisible(Map(
  function(raster_layer, output_path) {
    writeRaster(raster_layer, output_path, overwrite = TRUE)
  },
  resampled_rasters,
  resampled_paths
))

invisible(Map(
  function(raster_layer, output_path) {
    writeRaster(raster_layer, output_path, overwrite = TRUE)
  },
  normalized_rasters,
  normalized_paths
))

writeRaster(
  host_mean,
  file.path(
    overlap_dir,
    "Host_Resampled_QuantileNormalization_Mean.tif"
  ),
  overwrite = TRUE
)
writeRaster(
  breeding_overlap,
  file.path(overlap_dir, "Overlap_Breeding_hostAll_host_cumulative_Mean.tif"),
  overwrite = TRUE
)
writeRaster(
  nonbreeding_overlap,
  file.path(overlap_dir, "Overlap_NonBreeding_hostAll_host_cumulative_Mean.tif"),
  overwrite = TRUE
)
writeRaster(
  annual_overlap,
  file.path(overlap_dir, "Overlap_Annual_hostAll_host_cumulative_Mean.tif"),
  overwrite = TRUE
)

# 7. Visualization ---------------------------------------------------------

land_small <- ne_download(
  scale = "small",
  type = "land",
  category = "physical",
  returnclass = "sf"
)

annual_df <- as.data.frame(annual_overlap, xy = TRUE)
names(annual_df)[3] <- "value"

delta_df <- as.data.frame(delta_overlap, xy = TRUE)
names(delta_df)[3] <- "delta"

max_abs_delta <- max(abs(delta_df$delta), na.rm = TRUE)
delta_limits <- c(-max_abs_delta, max_abs_delta)

lon_values <- c(-120, -60, 0, 60, 120)
lat_values <- c(-30, 0, 30, 60)
lon_labels <- c("120°W", "60°W", "0°", "60°E", "120°E")
lat_labels <- c("30°S", "0°", "30°N", "60°N")

graticules_vertical <- data.frame(
  x = lon_values,
  xend = lon_values,
  y = -56,
  yend = 90
)
graticules_horizontal <- data.frame(
  x = -180,
  xend = 180,
  y = lat_values,
  yend = lat_values
)

map_base <- ggplot() +
  geom_segment(
    data = graticules_vertical,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey70",
    linewidth = 0.2,
    linetype = "dotted",
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = graticules_horizontal,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey70",
    linewidth = 0.2,
    linetype = "dotted",
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = lon_values,
    y = -57,
    label = lon_labels,
    size = 3,
    colour = "black",
    vjust = 1
  ) +
  annotate(
    "text",
    x = -200,
    y = lat_values,
    label = lat_labels,
    size = 3,
    colour = "black",
    hjust = 0
  )

map_theme <- theme_bw(base_family = "Arial") +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.background = element_blank(),
    legend.ticks.length = grid::unit(0.03, "cm"),
    legend.key.width = grid::unit(2.0, "cm"),
    legend.key.height = grid::unit(0.1, "cm"),
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0)
  )

map_guide <- guide_colorbar(
  title.position = "top",
  title.theme = element_text(size = 9, angle = 0, hjust = 0.5, vjust = 0.5),
  label.hjust = 0.5,
  ticks = TRUE,
  ticks.colour = "black",
  ticks.linewidth = 0.2,
  frame.colour = "black",
  frame.linewidth = 0.2
)

p_overlap <- map_base +
  geom_raster(data = annual_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    values = scales::rescale(seq(0, 1, by = 0.1), to = c(0, 1)),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),
    na.value = "white"
  ) +
  labs(fill = "AIV wild bird-domestic animal overlap") +
  geom_sf(
    data = land_small,
    colour = "grey70",
    fill = NA,
    linewidth = 0.2,
    alpha = 0.5
  ) +
  coord_sf(xlim = c(-200, 180), ylim = c(-63, 90), expand = FALSE) +
  map_theme +
  guides(fill = map_guide)

p_delta <- map_base +
  geom_raster(data = delta_df, aes(x = x, y = y, fill = delta)) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = delta_limits,
    oob = scales::squish,
    values = scales::rescale(
      seq(delta_limits[1], delta_limits[2], length.out = 11)
    ),
    na.value = "white",
    name = "Seasonal dynamics of overlap"
  ) +
  geom_sf(
    data = land_small,
    colour = "grey70",
    fill = NA,
    linewidth = 0.2,
    alpha = 0.5
  ) +
  coord_sf(xlim = c(-200, 180), ylim = c(-63, 90), expand = FALSE) +
  map_theme +
  guides(fill = map_guide)

land_mask_overlap <- rasterize(land_medium, annual_overlap, field = 1)
annual_overlap_zonal <- fill_land_na(annual_overlap, land_mask_overlap)
delta_overlap_zonal <- fill_land_na(delta_overlap, land_mask_overlap)

calculate_latitude_stats <- function(raster_layer, overlap_type) {
  raster_df <- as.data.frame(raster_layer, xy = TRUE)
  names(raster_df) <- c("lon", "lat", "value")

  raster_df |>
    group_by(lat) |>
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(type = overlap_type)
}

latitude_stats <- bind_rows(
  calculate_latitude_stats(annual_overlap_zonal, "Annual"),
  calculate_latitude_stats(delta_overlap_zonal, "Difference")
)

make_latitude_plot <- function(data, x_limits, x_breaks, x_label) {
  ggplot(data, aes(x = mean_value, y = lat)) +
    geom_ribbon(
      aes(xmin = mean_value - sd_value, xmax = mean_value + sd_value),
      alpha = 0.1,
      fill = "#000000"
    ) +
    geom_path(linewidth = 0.5, colour = "#000000", alpha = 0.5) +
    geom_vline(
      xintercept = 0,
      colour = "grey70",
      linetype = "dotted",
      linewidth = 0.5
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(-60, 90),
      breaks = c(-30, 0, 30, 60),
      labels = c("30°S", "0", "30°N", "60°N"),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = x_label, y = "Latitude") +
    theme_bw(base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 4, 0, 0, unit = "mm"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey70",
        linetype = "dotted",
        linewidth = 0.2
      ),
      panel.grid.minor.y = element_blank(),
      axis.title = element_text(colour = "black", size = 9),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.2),
      axis.ticks = element_line(colour = "black", linewidth = 0.2),
      axis.text = element_text(colour = "black", size = 9)
    )
}

p_lat_delta <- make_latitude_plot(
  filter(latitude_stats, type == "Difference"),
  x_limits = c(-0.3, 0.3),
  x_breaks = c(-0.3, -0.1, 0.1, 0.3),
  x_label = "Seasonal dynamics\nof overlap"
)

p_lat_annual <- make_latitude_plot(
  filter(latitude_stats, type == "Annual"),
  x_limits = c(-0.1, 0.8),
  x_breaks = c(0, 0.3, 0.6),
  x_label = "AIV wild bird-\ndomestic animal overlap"
)

p_annual_delta <- cowplot::plot_grid(
  p_overlap,
  p_delta,
  ncol = 1,
  align = "v"
)

save_publication_figure <- function(plot, filename_stem, width_cm, height_cm) {
  ggsave(
    filename = paste0(filename_stem, ".png"),
    plot = plot,
    height = height_cm,
    width = width_cm,
    units = "cm",
    dpi = 600,
    bg = "white"
  )

  width_in <- width_cm / 2.54
  height_in <- height_cm / 2.54

  svglite::svglite(
    paste0(filename_stem, ".svg"),
    width = width_in,
    height = height_in
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(filename_stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "Arial"
  )
  print(plot)
  grDevices::dev.off()
}

save_publication_figure(
  p_lat_delta,
  file.path(
    figure_dir,
    "Overlap_Annual_hostAll_Delta_Zonal_v1_legend_bottom"
  ),
  width_cm = 4,
  height_cm = 6.75
)
save_publication_figure(
  p_lat_annual,
  file.path(
    figure_dir,
    "Overlap_Annual_hostAll_Annual_Zonal_v1_legend_bottom"
  ),
  width_cm = 4,
  height_cm = 6.75
)
save_publication_figure(
  p_annual_delta,
  file.path(
    figure_dir,
    "Overlap_Annual_hostAll_Annual_Delta_v1_legend_bottom"
  ),
  width_cm = 14,
  height_cm = 15.22
)
