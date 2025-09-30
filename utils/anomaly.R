library(tidyverse)
library(anomalize)
library(changepoint)
library(changepoint.np)

detect_anomalies_zscore <- function(values, threshold = 3) {
  z_scores <- scale(values)
  anomalies <- which(abs(z_scores) > threshold)
  return(anomalies)
}

detect_anomalies_iqr <- function(values, multiplier = 2) {
  Q1 <- quantile(values, 0.25, na.rm = TRUE)
  Q3 <- quantile(values, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - multiplier * IQR
  upper_bound <- Q3 + multiplier * IQR
  anomalies <- which(values < lower_bound | values > upper_bound)
  return(anomalies)
}

detect_level_shifts <- function(values, minseglen = 20) {
    return(numeric(0))
    cp_pelt <- cpt.np(values, penalty = "MBIC", method = "PELT", test.stat = "empirical_distribution", minseglen = minseglen)
    pelt_shifts <- cpts(cp_pelt)


    return(pelt_shifts)
}
detect_large_shift <- function(values, window_size = 5, threshold_factor = 5) {
  n <- length(values)
  shifts <- numeric()

  med <- median(values, na.rm = TRUE)
  mad <- median(abs(values - med), na.rm = TRUE)

  if (n < 2 * window_size) {
    warning("Data length is too short for the given window size.")
    return(shifts)
  }

  calc_segment_stats <- function(start_idx, end_idx) {
    segment <- values[start_idx:end_idx]
    segment_mean <- mean(segment, na.rm = TRUE)
    segment_mad <- median(abs(segment - segment_mean), na.rm = TRUE)
    return(list(mean = segment_mean, mad = segment_mad))
  }


  for (i in (window_size + 1):(n - window_size)) {
    left_stats <- calc_segment_stats(i - window_size, i - 1)
    right_stats <- calc_segment_stats(i, i + window_size - 1)

    if (abs(right_stats$mean - left_stats$mean) > threshold_factor * max(left_stats$mad, right_stats$mad, mad)) {
      shifts <- c(shifts, i)
    }
  }


  return(unique(shifts))
}




plot_anomaly_data <- function(data, z_anomalies, iqr_anomalies, level_shifts, x_col = "Year", y_col = "Value", title = "Sea Level Records with Anomalies") {
  p <- ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_line(color = "blue") +
    geom_point(aes(color = "Data Points")) +
    geom_point(data = data[iqr_anomalies, ], aes(color = "IQR Anomalies"), size = 3, shape = 16) +
    geom_point(data = data[z_anomalies, ], aes(color = "Z-score Anomalies"), size = 3, alpha = 0.7, shape = 17) +
    labs(title = title, x = x_col, y = y_col) +
    scale_color_manual(values = c("Data Points" = "blue", "IQR Anomalies" = "green", "Z-score Anomalies" = "red", "Level Shifts" = "purple")) +
    theme_minimal() +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = 8)
    )

  if (length(level_shifts) > 0) {
    p <- p + geom_vline(aes(xintercept = data[[x_col]][level_shifts][1], color = "Level Shifts"), linetype = "dashed", size = 1)
  }

  return(p)
}

