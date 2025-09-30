source("predictions.R")
source("utils.R")

plot_ensemble_vs_SSPs_both <- function(station, method="MSE", from_file=TRUE, start.year=1970, stop.year=2019, show_real=FALSE, ...) {
  #' Plot the predictions from the ensemble method along with the AR6 SSP predictions
  #' @param station station name
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a ggplot object
  #' @examples
  #' plot_grid_ensemble_vs_SSPs(station.names, start.year=1970, stop.year=2019, from_file=TRUE)

  ensemble_no_col <- "ensemble_no"
  ensemble_cop_col <- "ensemble_cop"

  df <- get_predictions_all_method(station, from_file = from_file, start.year = start.year, stop.year = stop.year, ...)
  df <- df %>% dplyr::filter(Year >= 2020)

  ssps <- c("ssp119_medium", "ssp126_medium", "ssp245_medium", "ssp370_medium",
            "ssp585_medium", "ssp126_low", "ssp585_low")
  top_ssps <- ssps

  scenario_labels <- c(
    "ensemble_no" = "Unrestricted",
    "ensemble_cop" = "COP26 Restricted",
    "ssp119_medium" = "SSP1-1.9",
    "ssp126_medium" = "SSP1-2.6",
    "ssp245_medium" = "SSP2-4.5",
    "ssp370_medium" = "SSP3-7.0",
    "ssp585_medium" = "SSP5-8.5",
    "ssp126_low" = "SSP1-2.6 (Low)",
    "ssp585_low" = "SSP5-8.5 (Low)"
  )

  df_melted <- reshape2::melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")
  df_melted <- df_melted[df_melted$Scenario %in% c("ensemble_no", "ensemble_cop", top_ssps), ]

  line_types <- c("solid", "solid", rep("dashed", length(ssps)))
  names(line_types) <- c("ensemble_no", "ensemble_cop", top_ssps)

  line_sizes <- c(1.2, 1.2, rep(0.8, length(ssps)))
  names(line_sizes) <- c("ensemble_no", "ensemble_cop", top_ssps)

  colors <- c(
    "ensemble_no" = "#d62728",
    "ensemble_cop" = "#1f77b4",
    "ssp119_medium" = "#2ecc71",
    "ssp126_low" = "#3498db",
    "ssp126_medium" = "#85c1e9",
    "ssp245_medium" = "#f39c12",
    "ssp370_medium" = "#e74c3c",
    "ssp585_low" = "#8b0000",
    "ssp585_medium" = "#cd5c5c"
  )

  station_name <- get_station_name(station)
  station_id <- get_station_id(station)

  if(show_real){
    real_data <- get_tg_data(station) %>% convert_monthly_to_yearly() %>% dplyr::filter(year >= 2020)
    baseline <- get_baseline(station)
    real_data <- real_data %>% dplyr::mutate(RSL = tg - baseline, Year = year)

  }

  p <- ggplot() +
    geom_line(data = df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors, labels = scenario_labels) +
    scale_linetype_manual(values = line_types, labels = scenario_labels) +
    scale_size_manual(values = line_sizes, labels = scenario_labels) +
    theme_minimal() +
    labs(title = format_station_name(station_id),
         x = "Year",
         y = "RSL (m)",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_text(data = df_melted %>% dplyr::filter(Year == max(Year)),
              aes(x = Year, y = RSL,
                  label = scenario_labels[Scenario],
                  hjust = 1.1),
              size = 3, show.legend = FALSE)

  if(show_real){
    p <- p +
      geom_point(data = real_data,
                aes(x = Year, y = RSL, shape = "observed"),
                color = "black",
                size = 3) +
      scale_shape_manual(values = c(observed = 16),
                        labels = c(observed = "Observed"))
  }

  p
}


plot_compare_global_predictions_ssp <- function(
  #' Plot the global predictions along with the AR6 SSP predictions
  #' @return a ggplot object
  data_path = "data/GMSL_prediction_SSP_perma.csv",
  start.year = 1970, stop.year = 2019, show_real = FALSE, ...) {
  df <- read.csv(data_path)

  df <- df %>% dplyr::filter(Year >= 2020)

  scenarios <- c("mean_no", "mean_cop", "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585")

  scenario_labels <- c(
    "mean_no" = "Unrestricted",
    "mean_cop" = "COP26 Restricted",
    "ssp119" = "SSP1-1.9",
    "ssp126" = "SSP1-2.6",
    "ssp245" = "SSP2-4.5",
    "ssp370" = "SSP3-7.0",
    "ssp460" = "SSP4-6.0",
    "ssp585" = "SSP5-8.5"
  )

  df_melted <- reshape2::melt(df, id.vars = "Year", measure.vars = scenarios, variable.name = "Scenario", value.name = "RSL")
  df_melted <- df_melted %>% dplyr::mutate(RSL = RSL / 1000)

  line_types <- c("solid", "solid", rep("dashed", length(scenarios) - 2))
  names(line_types) <- scenarios

  line_sizes <- c(1.2, 1.2, rep(0.8, length(scenarios) - 2))
  names(line_sizes) <- scenarios

  colors <- c(
    "mean_no" = "#d62728",
    "mean_cop" = "#1f77b4",
    "ssp119" = "#2ecc71",
    "ssp126" = "#85c1e9",
    "ssp245" = "#f39c12",
    "ssp370" = "#e74c3c",
    "ssp460" = "#9b59b6",
    "ssp585" = "#cd5c5c"
  )

  p <- ggplot() +
    geom_line(data = df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors, labels = scenario_labels) +
    scale_linetype_manual(values = line_types, labels = scenario_labels) +
    scale_size_manual(values = line_sizes, labels = scenario_labels) +
    theme_minimal() +
    labs(title = "Global Mean Sea Level Projections: Our Model under Different GWP Scenarios",
         x = "Year",
         y = "RSL (m)",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_text(data = df_melted %>% dplyr::filter(Year == max(Year)),
              aes(x = Year, y = RSL,
                  label = scenario_labels[Scenario],
                  hjust = 1.1),
              size = 3, show.legend = FALSE)
  p
}

plot_compare_global_predictions_ar6 <- function(
  data_path = "data/GMSL_prediction_SSP_perma.csv",
  start.year = 1970, stop.year = 2100, show_real = FALSE, ...) {
  #' Plot the global predictions (mean_no, mean_cop) and AR6 SSP GMSL predictions (columns ending with _50)
  #' @return a ggplot object
  df_my <- read.csv(data_path)
  df_my <- df_my %>% dplyr::filter(Year >= 2020 & Year <= stop.year)
  df_my <- df_my[, c("Year", "mean_no", "mean_cop", "ssp245")]
  df_ar6 <- read.csv("data/ar6/ipcc_ar6_sea_level_projection_psmsl_id_0.csv")
  df_ar6 <- df_ar6 %>% dplyr::filter(Year >= 2020 & Year <= stop.year)
  ar6_cols <- grep("_50$", names(df_ar6), value = TRUE)
  df_ar6 <- df_ar6[, c("Year", ar6_cols)]

  my_2020 <- df_my[df_my$Year == 2020, "ssp245"] / 1000
  ar6_2020 <- df_ar6[df_ar6$Year == 2020, "ssp245_medium_50"]
  baseline_offset <- my_2020 - ar6_2020

  df_my$mean_no <- (df_my$mean_no / 1000) - baseline_offset
  df_my$mean_cop <- (df_my$mean_cop / 1000) - baseline_offset
  df_my_long <- reshape2::melt(df_my[, c("Year", "mean_no", "mean_cop")], id.vars = "Year", variable.name = "Scenario", value.name = "RSL")

  df_ar6_long <- reshape2::melt(df_ar6, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")

  df_all <- dplyr::bind_rows(df_my_long, df_ar6_long)

  scenario_labels <- c(
    "mean_no" = "Unrestricted",
    "mean_cop" = "COP26 Restricted",
    "ssp119_medium_50" = "SSP1-1.9",
    "ssp126_medium_50" = "SSP1-2.6",
    "ssp245_medium_50" = "SSP2-4.5",
    "ssp370_medium_50" = "SSP3-7.0",
    "ssp585_medium_50" = "SSP5-8.5",
    "ssp126_low_50" = "SSP1-2.6 (Low)",
    "ssp585_low_50" = "SSP5-8.5 (Low)"
  )

  scenarios <- intersect(names(scenario_labels), unique(df_all$Scenario))

  line_types <- c("solid", "solid", rep("dashed", length(scenarios) - 2))
  names(line_types) <- scenarios
  line_sizes <- c(1.2, 1.2, rep(0.8, length(scenarios) - 2))
  names(line_sizes) <- scenarios
  colors <- c(
    "mean_no" = "#d62728",
    "mean_cop" = "#1f77b4",
    "ssp119_medium_50" = "#2ecc71",
    "ssp126_medium_50" = "#85c1e9",
    "ssp245_medium_50" = "#f39c12",
    "ssp370_medium_50" = "#e74c3c",
    "ssp585_medium_50" = "#cd5c5c",
    "ssp126_low_50" = "#3498db",
    "ssp585_low_50" = "#8b0000"
  )

  df_all <- df_all[df_all$Scenario %in% scenarios, ]

  p <- ggplot() +
    geom_line(data = df_all, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors, labels = scenario_labels) +
    scale_linetype_manual(values = line_types, labels = scenario_labels) +
    scale_size_manual(values = line_sizes, labels = scenario_labels) +
    theme_minimal() +
    labs(title = "Global Mean Sea Level Projections: Our Model vs AR6",
         x = "Year",
         y = "RSL (m)",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_text(data = df_all %>% dplyr::filter(Year == max(Year)),
              aes(x = Year, y = RSL,
                  label = scenario_labels[Scenario],
                  hjust = 1.1),
              size = 3, show.legend = FALSE)
  p
}

plot_ensemble_vs_SSPs_both_with_zoom <- function(station, method="MSE", from_file=TRUE, start.year=1970, stop.year=2010, ...) {
  #' Plot the predictions from the ensemble method along with the AR6 SSP predictions
  #' @param station station name
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a ggplot object
  #' @examples
  #' plot_grid_ensemble_vs_SSPs(station.names, start.year=1970, stop.year=2019, from_file=TRUE)

  ensemble_no_col <- "ensemble_no"
  ensemble_cop_col <- "ensemble_cop"

  df <- get_predictions_all_method(station, from_file = from_file, start.year = start.year, stop.year = stop.year, ...)
  df <- df %>% dplyr::filter(Year >= 2020)

  ar6_ssps <- c("ssp119_medium", "ssp126_medium", "ssp245_medium", "ssp370_medium",
                "ssp585_medium", "ssp126_low", "ssp585_low")

  ensemble_ssps <- c("ensemble_ssp119", "ensemble_ssp126", "ensemble_ssp245",
                     "ensemble_ssp370", "ensemble_ssp585")

  available_ensemble_ssps <- intersect(ensemble_ssps, names(df))

  all_scenarios <- c("ensemble_no", "ensemble_cop", ar6_ssps, available_ensemble_ssps)

  scenario_labels <- c(
    "ensemble_no" = "Unrestricted",
    "ensemble_cop" = "COP26 Restricted",
    "ssp119_medium" = "AR6 SSP1-1.9",
    "ssp126_medium" = "AR6 SSP1-2.6",
    "ssp245_medium" = "AR6 SSP2-4.5",
    "ssp370_medium" = "AR6 SSP3-7.0",
    "ssp585_medium" = "AR6 SSP5-8.5",
    "ssp126_low" = "AR6 SSP1-2.6 (Low)",
    "ssp585_low" = "AR6 SSP5-8.5 (Low)",
    "ensemble_ssp119" = "Ensemble SSP1-1.9",
    "ensemble_ssp126" = "Ensemble SSP1-2.6",
    "ensemble_ssp245" = "Ensemble SSP2-4.5",
    "ensemble_ssp370" = "Ensemble SSP3-7.0",
    "ensemble_ssp585" = "Ensemble SSP5-8.5"
  )

  df_melted <- reshape2::melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")
  df_melted <- df_melted[df_melted$Scenario %in% all_scenarios, ]

  line_types <- c("solid", "solid", rep("dashed", length(ar6_ssps)), rep("dotted", length(available_ensemble_ssps)))
  names(line_types) <- c("ensemble_no", "ensemble_cop", ar6_ssps, available_ensemble_ssps)

  line_sizes <- c(1.2, 1.2, rep(0.8, length(ar6_ssps)), rep(1.0, length(available_ensemble_ssps)))
  names(line_sizes) <- c("ensemble_no", "ensemble_cop", ar6_ssps, available_ensemble_ssps)

  colors <- c(
    "ensemble_no" = "#d62728",
    "ensemble_cop" = "#1f77b4",
    "ssp119_medium" = "#2ecc71",
    "ssp126_low" = "#3498db",
    "ssp126_medium" = "#85c1e9",
    "ssp245_medium" = "#f39c12",
    "ssp370_medium" = "#e74c3c",
    "ssp585_low" = "#8b0000",
    "ssp585_medium" = "#cd5c5c",
    "ensemble_ssp119" = "#27ae60",
    "ensemble_ssp126" = "#2980b9",
    "ensemble_ssp245" = "#e67e22",
    "ensemble_ssp370" = "#c0392b",
    "ensemble_ssp585" = "#922b21"
  )

  station_name <- get_station_name(station)
  main_plot <- ggplot() +
    geom_line(data = df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors, labels = scenario_labels) +
    scale_linetype_manual(values = line_types, labels = scenario_labels) +
    scale_size_manual(values = line_sizes, labels = scenario_labels) +
    theme_minimal() +
    labs(title = format_station_name(station_name),
         x = "Year",
         y = "RSL (m)",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.margin = margin(5.5, 40, 5.5, 5.5, "pt")) +

    geom_text(data = {
                 final_data <- df_melted %>%
                   dplyr::filter(Year == max(Year)) %>%
                   dplyr::filter(!is.na(RSL)) %>%
                   dplyr::arrange(RSL)

                 if(nrow(final_data) > 1) {
                   rsl_range <- diff(range(final_data$RSL, na.rm = TRUE))
                   offset_size <- rsl_range * 0.02

                   for(i in 2:nrow(final_data)) {
                     if(abs(final_data$RSL[i] - final_data$RSL[i-1]) < offset_size) {
                       final_data$RSL[i] <- final_data$RSL[i-1] + offset_size
                     }
                   }
                 }

                 final_data
               },
              aes(x = Year, y = RSL,
                  label = scenario_labels[Scenario],
                  color = Scenario),
              hjust = -0.05,
              size = 3,
              show.legend = FALSE)

  zoom_data <- df_melted %>% dplyr::filter(Year >= 2090)
  inset_plot <- ggplot() +
    geom_line(data = zoom_data, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors, labels = scenario_labels) +
    scale_linetype_manual(values = line_types, labels = scenario_labels) +
    scale_size_manual(values = line_sizes, labels = scenario_labels) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.background = element_rect(fill = "white", color = "black"),
          plot.title = element_blank()) +
    labs(x = NULL, y = NULL)

  main_plot +
    inset_element(inset_plot, left = 0.01, bottom = 0.4, right = 0.3, top = 1)
}

plot_grid_ensemble_vs_SSPs <- function(station_list=station.names[1:9], method="MSE", from_file=TRUE, start.year=1970, stop.year=2019, ...) {
  #' Create a 3x3 grid of ensemble vs SSP plots for a list of stations
  #' @param station_list List of stations to plot
  #' @return A combined plot with 3x3 layout

  if(length(station_list) != 9) {
    stop("Please provide exactly 9 stations.")
  }

  plot_list <- lapply(station_list, function(station) {
    plot_ensemble_vs_SSPs_both(station, method, from_file, start.year, stop.year, ...)
  })

  extract_legend <- function(a_gplot) {
    g <- ggplotGrob(a_gplot)
    legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
    return(legend)
  }

  legend <- extract_legend(plot_list[[1]])

  plot_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))

  grid_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 3
  )

  final_plot <- cowplot::plot_grid(
    grid_plot,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  final_plot
}

plot_grid_ensemble_vs_SSPs_3_by_2 <- function(station_list=station.names[1:3], method="MSE", from_file=TRUE, start.year=1970, stop.year=2019, ...) {
  #' Create a 3x2 grid of ensemble vs SSP plots for 3 stations, with full range and zoomed views
  #' @param station_list List of 3 stations to plot
  #' @return A combined plot with 3x2 layout

  if(length(station_list) != 3) {
    stop("Please provide exactly 3 stations for the 3x2 grid layout.")
  }

  ssps <- c("ssp119_medium", "ssp126_medium", "ssp245_medium", "ssp370_medium",
            "ssp585_medium", "ssp126_low", "ssp585_low")

  line_types <- c("solid", "solid", rep("dashed", length(ssps)))
  names(line_types) <- c("ensemble_no", "ensemble_cop", ssps)

  line_sizes <- c(1.2, 1.2, rep(0.8, length(ssps)))
  names(line_sizes) <- c("ensemble_no", "ensemble_cop", ssps)

  colors <- c(
    "ensemble_no" = "#d62728",
    "ensemble_cop" = "#1f77b4",
    "ssp119_medium" = "#2ecc71",
    "ssp126_low" = "#3498db",
    "ssp126_medium" = "#85c1e9",
    "ssp245_medium" = "#f39c12",
    "ssp370_medium" = "#e74c3c",
    "ssp585_low" = "#8b0000",
    "ssp585_medium" = "#cd5c5c"
  )

  scenario_labels <- c(
    "ensemble_no" = "Unrestricted",
    "ensemble_cop" = "COP26 Restricted",
    "ssp119_medium" = "SSP1-1.9",
    "ssp126_medium" = "SSP1-2.6",
    "ssp245_medium" = "SSP2-4.5",
    "ssp370_medium" = "SSP3-7.0",
    "ssp585_medium" = "SSP5-8.5",
    "ssp126_low" = "SSP1-2.6 (Low)",
    "ssp585_low" = "SSP5-8.5 (Low)"
  )

  extract_legend <- function(a_gplot) {
    g <- ggplotGrob(a_gplot)
    legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
    return(legend)
  }

  plot_list <- lapply(station_list, function(station) {
    full_plot <- plot_ensemble_vs_SSPs_both(station, method, from_file, start.year, stop.year, ...) +
      theme(legend.position = "none") +
      ggtitle(paste0(format_station_name(get_station_name(station)), "\n(2020-2100)"))

    df <- get_predictions_all_method(station, from_file = from_file, start.year = start.year, stop.year = stop.year, ...)
    df_zoomed <- df %>% dplyr::filter(Year >= 2095 & Year <= 2100)

    df_melted <- reshape2::melt(df_zoomed, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")
    df_melted <- df_melted[df_melted$Scenario %in% c("ensemble_no", "ensemble_cop", ssps), ]

    zoomed_plot <- ggplot() +
      geom_line(data = df_melted,
                aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
      scale_color_manual(values = colors, labels = scenario_labels) +
      scale_linetype_manual(values = line_types, labels = scenario_labels) +
      scale_size_manual(values = line_sizes, labels = scenario_labels) +
      theme_minimal() +
      labs(title = "2095-2100",
           x = "Year") +
      theme(legend.position = "none",
            axis.title.y = element_blank(),
            plot.margin = margin(5, 10, 5, 5),
            plot.title = element_text(hjust = 0.5))

    list(full = full_plot, zoomed = zoomed_plot)
  })

  first_plot <- plot_ensemble_vs_SSPs_both(station_list[1], method, from_file, start.year, stop.year, ...)
  legend <- extract_legend(first_plot)

  plot_grid_list <- list()
  for(i in 1:length(station_list)) {
    plot_grid_list[[2*i-1]] <- plot_list[[i]]$full
    plot_grid_list[[2*i]] <- plot_list[[i]]$zoomed
  }

  grid_plot <- cowplot::plot_grid(
    plotlist = plot_grid_list,
    ncol = 2,
    nrow = 3,
    align = "h",
    rel_widths = c(0.7, 0.3)
  )

  final_plot <- cowplot::plot_grid(
    grid_plot,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  return(final_plot)
}
