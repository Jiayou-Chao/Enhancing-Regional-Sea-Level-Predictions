
library(dplyr)
library(jsonlite)
library(ggplot2)
library(gridExtra)
source("predictions.R")
source("utils.R")
source("generate/compare_prediction_diff.R")
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(reshape2)
library(cowplot)

options(readr.show_col_types = FALSE)

plot_global_station_results <- function(data=get_global_station_results(), from_year=1970, to_year=2019, scenario="no", title="Global Station Predictions", ...) {
  ensemble_no_col <- "ensemble_no"
  ensemble_cop_col <- "ensemble_cop"

  min_y <- min(data %>% dplyr::select(all_of(c(ensemble_no_col, ensemble_cop_col))), na.rm = TRUE)
  max_y <- max(data %>% dplyr::select(all_of(c(ensemble_no_col, ensemble_cop_col))), na.rm = TRUE)

  ensemble_col <- ifelse(scenario == "no", ensemble_no_col, ensemble_cop_col)

  plot_data <- data %>%
    dplyr::select(Year, station, all_of(ensemble_col)) %>%
    dplyr::filter(!is.na(.data[[ensemble_col]]))

  prediction_2100 <- plot_data %>%
    dplyr::filter(Year == 2100) %>%
    dplyr::select(station, prediction_2100 = all_of(ensemble_col))

  plot_data <- plot_data %>%
    dplyr::left_join(prediction_2100, by = "station")

  color_palette <- scales::col_numeric(palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), domain = range(prediction_2100$prediction_2100, na.rm = TRUE))
  plot_data <- plot_data %>%
    dplyr::mutate(color = color_palette(prediction_2100))

  ggplot(plot_data, aes(x = Year, y = get(ensemble_col), color = format_station_name(station), group = station)) +
    geom_line(size = 1) +
    geom_text(data = subset(plot_data, Year == max(Year)),
              aes(label = format_station_name(station)),
              size = 3,
              hjust = 1,
              vjust = -0.1) +
    labs(title = title,
         x = "Year",
         y = "Ensemble Prediction (m)",
         color = "Station") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    scale_y_continuous(limits = c(min_y, max_y)) +
    scale_color_manual(values = setNames(plot_data$color, format_station_name(plot_data$station)))
}

plot_global_scenarios_comparison <- function(data=get_global_station_results(), from_year=1970, to_year=2019, ...) {
  #' Plot Scenario Comparisons Side by Side
  #'
  #' This function generates side-by-side plots for "no" and "cop" scenarios
  #' using the global station results.
  #'
  #' @param data Data frame with global station results.

  no_plot <- plot_global_station_results(data, from_year, to_year, scenario="no", title="", ...) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

  cop_plot <- plot_global_station_results(data, from_year, to_year, scenario="cop", title="", ...) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

  legend <- cowplot::get_legend(no_plot + theme(legend.position = "bottom"))

  combined_plot <- grid.arrange(gridExtra::arrangeGrob(no_plot, cop_plot, ncol = 2), legend, ncol = 1, heights = c(10, 1))
  return(combined_plot)
}

get_valid_global_station_ids <- function(file_path = "global_map_stations.csv",
                                         filters = list(anomalies_detected = FALSE, missing_rate = 0.2)) {
  #' Get Valid Global Station IDs
  #'
  #' This function reads station data from a CSV file and applies customizable filters
  #' to return a list of valid station IDs based on the specified criteria.
  #'
  #' @param file_path Path to the CSV file containing station data. Default is "global_map_stations.csv".
  #' @param filters A named list of filter conditions for each column in the CSV.
  #'        - `anomalies_detected`: Logical filter, e.g., FALSE to exclude anomalies.
  #'        - `missing_rate`: Numeric filter, e.g., < 0.2 to limit missing data.
  #'        - `data_integrity`: Logical filter
  #'        - Additional columns from the global_map_stations.csv can also be used for filtering.
  #' @return A vector of station IDs that meet the specified filtering criteria.
  #' @example
  #' get_valid_global_station_ids(filters = list(anomalies_detected = FALSE, missing_rate = 0.2))
  station_id_data <- read_csv(file_path)

  filtered_data <- station_id_data

  for (filter_name in names(filters)) {
    filter_value <- filters[[filter_name]]
    if (is.logical(filter_value)) {
      filtered_data <- filtered_data %>% dplyr::filter(.data[[filter_name]] == filter_value)
    } else if (is.numeric(filter_value)) {
      filtered_data <- filtered_data %>% dplyr::filter(.data[[filter_name]] < filter_value)
    }
  }

  valid_station_ids <- filtered_data %>% dplyr::pull(station_id)
  valid_station_ids <- c(valid_station_ids, 286, 350, 453, 467, 471, 474, 755, 765, 982, 1109, 1112)
  return(unique(valid_station_ids))
}

get_all_global_station_ids <- function(file_path = "global_map_stations.csv") {
  station_id_data <- read_csv(file_path)
  station_id_data %>% dplyr::pull(station_id)
}

delete_idx_files <- function(directory = "data/era5data") {
  idx_files <- list.files(path = directory, pattern = "\\.idx$", recursive = TRUE, full.names = TRUE)
  file.remove(idx_files)
}

get_all_station_results <- function(save_path="results/global_station_results_2100_full.csv",
                                    from_file=TRUE,
                                    prediction_from_file=TRUE,
                                    add_vlm=FALSE,
                                    log_file="results/log.txt") {
  #' @param save_path path to the results to a CSV file
  #' @param from_file whether to read results directly from save_path
  #' @return a data frame with station_id, lat, and lon, and prediction_2100
  if (from_file && file.exists(save_path)) {
    return(readr::read_csv(save_path))
  } else {
    delete_idx_files()
    all_station_ids <- get_all_global_station_ids()

    location_data <- read_ar6_location_list()

    valid_locations <- location_data %>% dplyr::filter(station_id %in% all_station_ids)

    log_connection <- file(log_file, open = "a")

    predictions <- lapply(valid_locations$station_id, function(station_id) {
      prediction_data <- tryCatch({
        writeLines(paste("Processing station ID:", station_id), log_connection)
        flush(log_connection)
        if (TRUE) {
          get_predictions_all_method(station_id, from_file = prediction_from_file, add_vlm = add_vlm)
        }

      }, error = function(e) {
        writeLines(paste("Error in get_predictions_all_method for station_id:", station_id), log_connection)
        flush(log_connection)
        return(NULL)
      })
      if (!is.null(prediction_data)) {
      prediction_2100 <- prediction_data %>% dplyr::filter(Year == 2100)

      return(cbind(data.frame(
        station_id = station_id,
        prediction_2100_no = prediction_2100 %>% dplyr::mutate(value = ensemble_no) %>% pull(value),
        prediction_2100_cop = prediction_2100 %>% dplyr::mutate(value = ensemble_cop) %>% pull(value),
        ensemble_diff_ipcc_no = prediction_2100 %>% dplyr::mutate(diff_no = ensemble_no - ssp245_medium) %>% pull(diff_no),
        ensemble_diff_ipcc_cop = prediction_2100 %>% dplyr::mutate(diff_cop = ensemble_cop - ssp126_medium) %>% pull(diff_cop)
      ), prediction_2100 %>% dplyr::select(-Year)))
    }
    })

    close(log_connection)

    predictions_df <- do.call(rbind, predictions)

    result <- dplyr::left_join(valid_locations, predictions_df, by = "station_id")

    if (!from_file) {
      write_csv(result, save_path)
    }

    return(result)
  }
}

combine_station_predictions <- function(station_ids) {
  results_list <- list()

  for (station_id in station_ids) {
    prediction_data <- get_AR6_prediction_ipcc_tool_full(station_id, zero_2020 = TRUE)

    filtered_data <- prediction_data %>%
      dplyr::filter(Year == 2100) %>%
      mutate(station_id = station_id) %>%
      dplyr::select(station_id, everything()) %>%
      dplyr::select(-Year)

    results_list[[station_id]] <- filtered_data
  }

  combined_results <- bind_rows(results_list)

  return(combined_results)
}

append_AR6_predictions <- function(df) {

  station_ids <- unique(df$station_id)

  new_df <- combine_station_predictions(station_ids)
  left_join(df, new_df, by = "station_id")

}

get_valid_station_results <- function(from_file=TRUE, filters = list(anomalies_detected = FALSE, missing_rate = 0.2)) {
  #' @param filters A named list of filter conditions for each column in the CSV.
  #'        - `anomalies_detected`: Logical filter, e.g., FALSE to exclude anomalies.
  #'        - `missing_rate`: Numeric filter, e.g., < 0.2 to limit missing data.
  #'        - Additional columns from the CSV can also be used for filtering.
  all_global_results <- get_all_station_results(from_file=from_file)
  valid_station_ids <- get_valid_global_station_ids(filters = filters)
  all_global_results %>%
    dplyr::filter(station_id %in% valid_station_ids) %>%
    tidyr::drop_na(prediction_2100_no) %>%
    dplyr::filter(prediction_2100_no != ensemble_diff_ipcc_no) %>%
    mutate(ipcc_2100_no = prediction_2100_no - ensemble_diff_ipcc_no,
           ipcc_2100_cop = prediction_2100_cop - ensemble_diff_ipcc_cop)
}

plot_global_map <- function(data,
                           scenario="no",
                           color_scale_range=range(-1, 2),
                           show_legend=TRUE) {
  #' This function takes a data frame containing latitude, longitude,
  #' and a value to plot on a global map. The data is visualized with
  #' points colored according to the specified value.
  #'
  #' Args:
  #'   data: A data frame containing the following columns:
  #'     - lat: Numeric values representing the latitude of the points.
  #'     - lon: Numeric values representing the longitude of the points.
  #'     - prediction_2100_no or prediction_2100_cop: Numeric values used for coloring the points on the map based on the scenario.
  #'
  #' Returns:
  #'   A ggplot object displaying the global map with data points.

  predefined_legend_titles <- c(
    "no" = "Sea Level\nChange (m)",
    "cop" = "Sea Level\nChange (m)",
    "diff_no" = "Sea Level\nDifference (m)",
    "diff_cop" = "Sea Level\nDifference (m)",
    "ipcc_no" = "Sea Level\nChange (m)",
    "ipcc_cop" = "Sea Level\nChange (m)"
  )

  predefined_value_columns <- c(
    "no" = "prediction_2100_no",
    "cop" = "prediction_2100_cop",
    "diff_no" = "ensemble_diff_ipcc_no",
    "diff_cop" = "ensemble_diff_ipcc_cop",
    "ipcc_no" = "ipcc_2100_no",
    "ipcc_cop" = "ipcc_2100_cop"
  )

  legend_title <- if (scenario %in% names(predefined_legend_titles)) {
    predefined_legend_titles[[scenario]]
  } else {
    scenario
  }

  value_column <- if (scenario %in% names(predefined_value_columns)) {
    predefined_value_columns[[scenario]]
  } else {
    scenario
  }

  data <- data %>%
    dplyr::rename(latitude = lat, longitude = lon) %>%
    dplyr::mutate(value = .data[[value_column]])
  data_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

  world <- ne_countries(scale = "medium", returnclass = "sf")

  p <- ggplot() +
    geom_sf(data = world, fill = "gray95", color = "gray70", size = 0.2) +
    geom_sf(data = data_sf, aes(color = value), size = 2) +
    scale_color_gradient2(
      low = "darkblue",
      mid = "white",
      high = "darkred",
      midpoint = 0,
      limits = color_scale_range
    ) +
    theme_minimal() +
    labs(color = legend_title) +
    theme(
      legend.title = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = if(show_legend) "right" else "none"
    ) +
    coord_sf(expand = FALSE)

  return(p)
}

create_combined_global_maps <- function(data_no, color_scale_range) {
  scenarios <- list(
    list(scenario = "no", label = "a) Unrestricted Policy"),
    list(scenario = "cop", label = "b) COP26 Policy"),
    list(scenario = "ipcc_no", label = "c) IPCC SSP2-4.5"),
    list(scenario = "ipcc_cop", label = "d) IPCC SSP1-2.6"),
    list(scenario = "diff_no", label = "e) Unrestricted - IPCC SSP2-4.5"),
    list(scenario = "diff_cop", label = "f) COP26 - IPCC SSP1-2.6")
  )

  plots <- lapply(scenarios, function(s) {
    plot_global_map(data = data_no,
                   scenario = s$scenario,
                   color_scale_range = color_scale_range,
                   show_legend = FALSE) +
      labs(title = s$label) +
      theme(
        plot.title = element_text(size = 8, face = "bold"),
        plot.margin = margin(2, 2, 2, 2, "mm")
      )
  })

  legend_plot <- plot_global_map(data = data_no,
                                scenario = "no",
                                color_scale_range = color_scale_range) +
    labs(color = "Sea Level\nChange (m)")

  legend <- cowplot::get_legend(legend_plot)

  plots_grid <- plot_grid(
    plotlist = plots,
    ncol = 2,
    align = 'v',
    axis = 'lr'
  )

  final_plot <- plot_grid(
    plots_grid,
    legend,
    ncol = 2,
    rel_widths = c(0.9, 0.1)
  )

  return(final_plot)
}

plot_predictions_scatter <- function(data_no, x_axis = "prediction_2100_no", y_axis = "ipcc_2100_no",
                                     main = "Scatter Plot with 45-Degree Line") {
  plot(x = data_no[[x_axis]], y = data_no[[y_axis]],
       xlab = x_axis,
       ylab = y_axis,
       main = main)

  abline(a = 0, b = 1, col = "red", lty = 2)
}


plot_ssp_differences <- function(df, x_axis = "prediction_2100_no") {
  #' Plot SSP Differences and Calculate Statistics
  #'
  #' This function calculates differences between predictions and SSP values,
  #' generates a scatter plot of these differences, and prints the mean and
  #' standard deviation of the differences for each SSP.
  #'
  #' @param df Data frame containing prediction and SSP columns.
  #' @param x_axis Column name for the x-axis values. Default is "prediction_2100_no".
  #' @return A ggplot object displaying the scatter plot of differences.

  if (!x_axis %in% colnames(df)) {
    stop(paste("The dataframe must contain a column named", x_axis))
  }

  ssp_columns <- df %>% dplyr::select(starts_with("ssp")) %>%
    dplyr::select(ends_with("medium"), ends_with("low")) %>%
    colnames()

  if (length(ssp_columns) == 0) {
    stop("The dataframe must contain at least one column starting with 'ssp'")
  }

  differences <- data.frame(x_axis = df[[x_axis]])
  for (ssp_col in ssp_columns) {
    differences[[ssp_col]] <- df[[x_axis]] - df[[ssp_col]]
  }

  ssp_stats <- ssp_columns %>%
    purrr::set_names() %>%
    purrr::map_df(~ data.frame(
      SSP = .x,
      MeanDifference = mean(differences[[.x]], na.rm = TRUE),
      SD = sd(differences[[.x]], na.rm = TRUE),
      MeanAbsoluteDifference = mean(abs(differences[[.x]]), na.rm = TRUE)
    ))

  ssp_stats <- ssp_stats %>%
    arrange(SSP)

  print(ssp_stats)
  melted_df <- melt(differences, id.vars = "x_axis", variable.name = "ssp", value.name = "difference")

  plot <- ggplot(melted_df, aes(x = x_axis, y = difference, color = ssp)) +
    geom_point() +
    labs(title = "Scatter Plot of Differences", x = x_axis, y = "Difference (prediction - SSP)") +
    theme_minimal() +
    theme(legend.title = element_blank())

  return(plot)
}

plot_prediction_confidence_vs_ssp <- function(df,
                                              prediction_column = "prediction_2100_no",
                                              ssp_column = "ssp245_medium",
                                              low95_column = "ensemble_no_low95",
                                              up95_column = "ensemble_no_up95") {

  required_columns <- c(prediction_column, low95_column, up95_column, ssp_column)
  if (!all(required_columns %in% colnames(df))) {
    stop(paste("Dataframe must contain the following columns:", paste(required_columns, collapse = ", ")))
  }

  df <- df[complete.cases(df[, required_columns]), ]

  df$color <- ifelse(df[[ssp_column]] >= df[[low95_column]] & df[[ssp_column]] <= df[[up95_column]],
                     "Within 95% Prediction Interval",
                     "Outside 95% Prediction Interval")

  within_count <- sum(df$color == "Within 95% Prediction Interval")
  outside_count <- sum(df$color == "Outside 95% Prediction Interval")
  total_count <- nrow(df)
  within_percentage <- (within_count / total_count) * 100
  outside_percentage <- (outside_count / total_count) * 100

  cat("Count within 95% Prediction Interval:", within_count, "\n")
  cat("Percentage within 95% Prediction Interval:", round(within_percentage, 2), "%\n")
  cat("Count outside 95% Prediction Interval:", outside_count, "\n")
  cat("Percentage outside 95% Prediction Interval:", round(outside_percentage, 2), "%\n")

  ggplot(df, aes(x = .data[[prediction_column]], y = .data[[ssp_column]], color = color)) +
    geom_point() +
    scale_color_manual(values = c("Within 95% Prediction Interval" = "blue", "Outside 95% Prediction Interval" = "red")) +
    labs(title = paste("Scatter Plot of", prediction_column, "vs", ssp_column),
         x = prediction_column,
         y = ssp_column,
         color = "Legend") +
    theme_minimal()
}



if (sys.nframe() == 0) {
  data_no <- get_valid_station_results(from_file=T, filters = list(missing_rate=0.3))
  data_no <- append_AR6_predictions(data_no)
  write_csv(data_no, "results/global_station_results_2100.csv")
  data_no <- read_csv("results/global_station_results_2100.csv")
  color_scale_range <- range(-1, 2)

  combined_plot <- create_combined_global_maps(data_no, color_scale_range)
  combined_plot
  ggsave("results/combined_global_maps.png", combined_plot,
         width = 10, height = 8, dpi = 300)
  ggsave("results/combined_global_maps.jpg", combined_plot,
         width = 10, height = 8, dpi = 300)
}
