source("predictions.R")

generate_prediction_diff <- function(
    station_ids,
    Year = 2100,
    scenario1 = "ensemble_no",
    scenario2 = "ssp245_medium",
    save_path = NA,
    from_file = TRUE
) {
  #' Generate the prediction difference between two scenarios
  #' @param Year the year to compare
  #' @param scenario1 the first scenario, see available scenarios in `names(get_predictions_all_method(station))`
  #' @param scenario2 the second scenario
  #' @param save_path the path to save the results
  #' @param from_file if TRUE, read the results from `save_path` if the file exists. If the file does not exist, the function will call itself with from_file = FALSE
  #' @return a data frame with columns: station_id, diff = scenario1 - scenario2, lat, lon, scenario1, scenario2

    if (is.na(save_path)) {
        create_directory("results/prediction_diff")
        save_path <- paste0("results/prediction_diff/", Year, "_", scenario1, "_", scenario2, ".csv")
    }
    if (from_file && file.exists(save_path)) {
        return(readr::read_csv(save_path, show_col_types = FALSE))
    }

    results <- purrr::map_dfr(
        station_ids,
        function(station) {
            preds <- get_predictions_all_method(station, from_file = TRUE)
            row <- dplyr::filter(preds, Year == !!Year)
            if (nrow(row) == 0 || !(scenario1 %in% names(row)) || !(scenario2 %in% names(row))) {
                coords <- station2lonlat(station)
                return(tibble::tibble(
                    station_id = station,
                    diff = NA_real_,
                    scenario1 = NA_real_,
                    scenario2 = NA_real_,
                    lat = coords$lat,
                    lon = coords$lon
                ))
            }
            val1 <- as.numeric(row[[scenario1]])
            val2 <- as.numeric(row[[scenario2]])
            coords <- station2lonlat(station)
            tibble::tibble(
                station_id = station,
                diff = val1 - val2,
                scenario1 = val1,
                scenario2 = val2,
                lat = coords$lat,
                lon = coords$lon
            )
        }
    )

    readr::write_csv(results, save_path)
    results
}
