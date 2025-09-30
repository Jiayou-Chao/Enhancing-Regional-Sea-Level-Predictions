source("plot_global.R")

data_no <- read_csv("results/global_station_results_2100.csv")
color_scale_range <- range(-1, 2)
Year <- 2050

scenario_pairs <- list(
  c("ensemble_no", "ssp245_medium"),
  c("ensemble_no", "ssp370_medium"),
  c("ensemble_cop", "ssp119_medium"),
  c("ensemble_cop", "ssp126_medium"),
  c("ensemble_cop", "ssp370_medium")
)

for (pair in scenario_pairs) {
  scenario1 <- pair[1]
  scenario2 <- pair[2]
  diff_label <- paste0(scenario1, " - ", scenario2)
  img_path <- paste0("results/prediction_diff/", Year, "_", scenario1, "_", scenario2, ".png")

  data <- data_no$station_id %>%
    generate_prediction_diff(scenario1 = scenario1, scenario2 = scenario2, Year = Year, from_file = F)
  print(paste("The mean of the difference of", scenario1, "and", scenario2, "is", mean(data$diff, na.rm = TRUE)))
  print(paste("The MSE of the difference of", scenario1, "and", scenario2, "is", mean(data$diff^2, na.rm = TRUE)))
  data_no$station_id %>%
    generate_prediction_diff(scenario1 = scenario1, scenario2 = scenario2, Year = Year) %>%
    dplyr::rename(!!diff_label := diff) %>%
    plot_global_map(scenario = diff_label, color_scale_range = color_scale_range)

  ggsave(img_path, width = 10, height = 8)
  print(paste0("Image saved to ", img_path))
}
