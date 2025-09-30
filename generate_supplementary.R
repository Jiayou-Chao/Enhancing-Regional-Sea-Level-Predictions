library(tidyverse)
source("plot_US.R")
source("plot_global.R")

data_no <- read_csv("results/global_station_results_2100.csv")

create_directory("results/supplementary/long_term_prediction")
create_directory("results/supplementary/short_term_backtesting")

for (station_id in unique(data_no$station_id)) {
  long_term_path <- paste0("results/supplementary/long_term_prediction/Long-term_predictions_vs_AR6_station_", station_id, ".jpg")
  if (!file.exists(long_term_path)) {
    plot_ensemble_vs_SSPs_both(station_id, start.year=1970, stop.year=2019, from_file=TRUE, show_real=TRUE)
    ggsave(long_term_path, width = 12, height = 8)
  }

  short_term_path <- paste0("results/supplementary/short_term_backtesting/Backtesting_results_station_", station_id, ".jpg")
  if (!file.exists(short_term_path)) {
    plot_combined_backtesting_results(station_id)
    ggsave(short_term_path, width = 12, height = 8)
  }
}

if(file.exists("results/supplementary/long_term_prediction.zip")){
  file.remove("results/supplementary/long_term_prediction.zip")
}
zip(zipfile = "results/supplementary/long_term_prediction.zip",
    files = list.files("results/supplementary/long_term_prediction", full.names = TRUE, pattern = ".jpg"))

if(file.exists("results/supplementary/short_term_backtesting.zip")){
  file.remove("results/supplementary/short_term_backtesting.zip")
}
zip(zipfile = "results/supplementary/short_term_backtesting.zip",
    files = list.files("results/supplementary/short_term_backtesting", full.names = TRUE, pattern = ".jpg"))
