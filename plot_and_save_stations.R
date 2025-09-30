library(ggplot2)
library(dplyr)
library(progress)
source("plot_US.R")

mse <- function(x, y) {
  mean((x - y)^2, na.rm = TRUE)
}

plot_and_save_stations <- function(station_list, output_dir = "plots",
                                 method = "MSE", from_file = TRUE,
                                 start.year = 1970, stop.year = 2019) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  total_stations <- length(station_list)
  pb <- progress_bar$new(
    format = "Processing [:bar] :percent | Station :current/:total | Elapsed: :elapsed | ETA: :eta",
    total = total_stations,
    clear = FALSE,
    width = 80
  )

  for (station in station_list) {
    df <- get_predictions_all_method(station, from_file = from_file,
                                   start.year = start.year, stop.year = stop.year,
                                   show_real = TRUE) %>%
            dplyr::select(Year, real_tg, ensemble_no, ensemble_cop, starts_with("ssp"))

    pred_2100 <- df$ensemble_no[df$Year == 2100]

    pred_formatted <- sprintf("%.2f", pred_2100)

    station_name <- get_station_name(station)

    p <- plot_ensemble_vs_SSPs_both(station, method = method,
                                   from_file = TRUE,
                                   start.year = start.year,
                                   stop.year = stop.year,
                                   show_real = T)

    filename <- file.path(output_dir,
                         paste0(pred_formatted, "m_",
                               gsub(" ", "-", station_name),
                               ".png"))
    filename2 <- file.path(output_dir,
                         paste0(pred_formatted, "m_",
                               gsub(" ", "-", station_name),
                               ".jpg"))

    ggsave(filename2, p, width = 10, height = 6, dpi = 300)

    df_future <- df %>%
      dplyr::filter(Year > stop.year & !is.na(real_tg))
    if (nrow(df_future) > 0) {
      mse_cols <- setdiff(colnames(df_future), c("Year", "real_tg"))
      mse_vals <- purrr::map_dbl(mse_cols, ~ mse(df_future[[.x]], df_future$real_tg))
      mse_named <- setNames(as.list(mse_vals), mse_cols)
      ens_vals <- mse_named[c("ensemble_no", "ensemble_cop")]
      smaller_ens <- if (all(!is.na(ens_vals))) {
        names(ens_vals)[which.min(unlist(ens_vals))]
      } else { NA_character_ }
      min_var <- names(mse_named)[which.min(unlist(mse_named))]
      ordered_mse <- paste(
        names(sort(unlist(mse_named))),
        signif(sort(unlist(mse_named)), 4),
        sep = "=", collapse = "; "
      )
      mse_row <- tibble::tibble(
        station = station_name,
        !!!mse_named,
        smaller_ensemble = smaller_ens,
        min_mse_variable = min_var,
        ordered_mse = ordered_mse
      )
      out_csv <- file.path(output_dir, "MSE_summary_all_stations.csv")
      readr::write_csv(
        mse_row,
        out_csv,
        append = file.exists(out_csv),
        col_names = !file.exists(out_csv)
      )
    }

    pb$tick()

    cat(sprintf("Saved plot for %s with 2100 prediction of %sm\n",
                station_name, pred_formatted))
  }
}

data_no <- read_csv("results/global_station_results_2100.csv")
station_list <- data_no$station_id
plot_and_save_stations(station_list, output_dir = "results/prediction_plots", method = "MSE", from_file = T, start.year = 1970, stop.year = 2019)
