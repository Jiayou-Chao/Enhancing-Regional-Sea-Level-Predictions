#' If you encounter a bug "Year is not a recognized column name",
#' try change `from_file = TRUE` to `from_file = FALSE` (search `##` in this file)
library(knitr)
library(magrittr)
source("predictions.R")
source("ModelBacktesting.R")
source("data_integrity.R")
options(readr.show_col_types = FALSE)

library(jsonlite)

processed_stations_data <- data.frame(station_id = integer(), data_integrity = logical())

ar6_location_list <- read_ar6_location_list()

for (station_name in ar6_location_list$station_id) {
  integrity_reasons <- check_data_integrity(station_name)
  if (any(unlist(integrity_reasons))) {
    new_record <- data.frame(
      station_id = station_name,
      data_integrity = FALSE,
      data_range_issue = integrity_reasons$data_range_issue,
      missing_values = integrity_reasons$missing_values,
      anomalies_detected = integrity_reasons$anomalies_detected,
      missing_rate = integrity_reasons$missing_rate
    )
  } else {
    new_record <- data.frame(
      station_id = station_name,
      data_integrity = TRUE,
      data_range_issue = FALSE,
      missing_values = FALSE,
      anomalies_detected = FALSE,
      missing_rate = 0
    )
  }
  processed_stations_data <- rbind(processed_stations_data, new_record)
}

write_csv(processed_stations_data, "results/global_map_stations.csv")
