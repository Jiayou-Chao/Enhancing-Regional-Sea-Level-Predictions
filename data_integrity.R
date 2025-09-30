library(tidyverse)
library(anomalize)
library(changepoint)
source("utils/anomaly.R")
source("VLM.R")

.get_sea_level_data <- function(station) {
  get_year_data_raw(station) %>%
    dplyr::rename(Value = tg) %>%
    dplyr::filter(Year >= 1970 & Year < 2020)
}

detect_anomalies_and_shifts <- function(data) {
  z_anomalies <- tryCatch(
    detect_anomalies_zscore(data$Value),
    error = function(e) {
      cat("Error in detecting Z-score anomalies:", e$message, "\n")
      return(numeric(0))
    }
  )

  iqr_anomalies <- tryCatch(
    detect_anomalies_iqr(data$Value),
    error = function(e) {
      cat("Error in detecting IQR anomalies:", e$message, "\n")
      return(numeric(0))
    }
  )
  level_shifts <- tryCatch(
    detect_level_shifts(data$Value),
    error = function(e) {
      cat("Error in detecting level shifts:", e$message, "\n")
      return(numeric(0))
    }
  )
  large_shifts <- tryCatch(
    detect_large_shift(data$Value),
    error = function(e) {
      cat("Error in detecting large shifts:", e$message, "\n")
      return(numeric(0))
    }
  )
  all_shifts <- unique(c(level_shifts, large_shifts))
  list(z_anomalies = z_anomalies, iqr_anomalies = iqr_anomalies, all_shifts = all_shifts)
}

plot_anomaly_year_station <- function(station) {
  sea_level_data <- .get_sea_level_data(station)
  anomalies <- detect_anomalies_and_shifts(sea_level_data)

  cat("Z-score anomalies at indices:", anomalies$z_anomalies, "\n")
  cat("IQR anomalies at indices:", anomalies$iqr_anomalies, "\n")
  cat("Level shifts detected at indices:", anomalies$all_shifts, "\n")

  title <- paste("Sea Level Anomalies for", format_station_name(get_station_name(station)))
  plot <- plot_anomaly_data(sea_level_data, anomalies$z_anomalies, anomalies$iqr_anomalies, anomalies$all_shifts, title=title)
  print(plot)
}

check_data_integrity <- function(station) {
  calculate_missing_rate <- function(data) {
    missing_rate <- sum(data$is_na > 0) / nrow(data)
    return(missing_rate)
  }

  sea_level_data <- .get_sea_level_data(station)
  integrity_reasons <- list(data_range_issue = FALSE, missing_values = FALSE, anomalies_detected = FALSE,
                            missing_rate = 0)

  if (any(is.na(sea_level_data$Value))) {
    integrity_reasons$data_range_issue <- TRUE
    cat("Data range is not 1970-2020 for station", station, "\n")
  }

  sea_level_data$is_na <- ifelse(is.na(sea_level_data$Value), 1, sea_level_data$is_na)
  if (any(sea_level_data$is_na > 0)) {
    integrity_reasons$missing_values <- TRUE
    cat("Missing values detected in station", station, "\n")
  }

  missing_rate <- sum(sea_level_data$is_na > 0) / length(sea_level_data$is_na)
  integrity_reasons$missing_rate <- missing_rate
  cat("Missing rate:", missing_rate, "\n")

  sea_level_data <- sea_level_data %>% drop_na(Value)
  anomalies <- detect_anomalies_and_shifts(sea_level_data)

  if (length(anomalies$z_anomalies) > 0 || length(anomalies$iqr_anomalies) > 0 || length(anomalies$all_shifts) > 0) {
    cat("Anomalies detected in station", station, "\n")
    integrity_reasons$anomalies_detected <- TRUE
  }

  return(integrity_reasons)
}

if (sys.nframe() == 0) {
  for (station in list(1, "SETE")) {
    check_data_integrity(station)
    plot_anomaly_year_station(station)
  }
}
