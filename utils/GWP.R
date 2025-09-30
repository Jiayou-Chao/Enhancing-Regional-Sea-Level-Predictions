library(timeSeries)
library(tseries)
library(forecast)
library(tidyverse)
library(pracma)
library(astsa)
library(ggpubr)

get_old_GWP_forecast <- function(file_path = 'data/YearData.csv',
                                 end_train_year = 2019){
  #' Function wrapper to get GWP forecast from historical data.
  #' @return
  #' A data frame named `GWP` containing the historical data and the forecasted
  #' values for both scenarios, matching the original script's structure.
  #' Expected columns include: "Year", "GWP_no", "GWP_cop", "GWP_no_80_l", "GWP_no_95_l",
  #' "GWP_no_80_u", "GWP_no_95_u", "GWP_cop_80_l", "GWP_cop_95_l",
  #' "GWP_cop_80_u", "GWP_cop_95_u", "GWP_no_99_u", "GWP_no_99_l",
  #' "GWP_cop_99_u", "GWP_cop_99_l"

  if (!requireNamespace("forecast", quietly = TRUE)) stop("Package 'forecast' is required.")
  climate_full <- read.csv(file_path)
  start_train_year <- min(climate_full$Year)
  forecast_end_year <- 2100
  h <- forecast_end_year - end_train_year

  climate_train <- climate_full[climate_full$Year <= end_train_year, ]

  fit_GWP <- forecast::Arima(climate_train$GWP,
                             order = c(1,1,1),
                             include.drift = TRUE,
                             method = 'ML')
  pred_GWP_no <- forecast::forecast(fit_GWP, h = h)

  pred_GWP_no$x     <- stats::ts(pred_GWP_no$x,
                                 start = start_train_year,
                                 end   = end_train_year)
  pred_GWP_no$mean  <- stats::ts(pred_GWP_no$mean,
                                 start = end_train_year + 1,
                                 end   = forecast_end_year)
  pred_GWP_no$lower <- stats::ts(pred_GWP_no$lower,
                                 start = end_train_year + 1,
                                 end   = forecast_end_year)
  pred_GWP_no$upper <- stats::ts(pred_GWP_no$upper,
                                 start = end_train_year + 1,
                                 end   = forecast_end_year)

  fit_N2O <- forecast::Arima(climate_train$N2O,
                             order = c(1,1,1),
                             include.drift = TRUE,
                             method = 'ML')
  pred_N2O_no <- forecast::forecast(fit_N2O, h = h)$mean

  CO2_2010 <- climate_full$CO2[climate_full$Year == 2010] -
              climate_full$CO2[climate_full$Year == 2009]
  CH4_2020 <- climate_full$CH4[climate_full$Year == 2020] -
              climate_full$CH4[climate_full$Year == 2019]

  last_CO2 <- climate_full$CO2[climate_full$Year == end_train_year]
  last_CH4 <- climate_full$CH4[climate_full$Year == end_train_year]

  years_to_goal <- 2030 - end_train_year
  co2_seq <- last_CO2 + cumsum(
    CO2_2010 - (1:years_to_goal)/years_to_goal*(CO2_2010*0.45)
  )
  ch4_seq <- last_CH4 + cumsum(
    CH4_2020 - (1:years_to_goal)/years_to_goal*(CH4_2020*0.30)
  )
  pred_CO2_cop <- c(co2_seq,
                    rep(tail(co2_seq,1), h - length(co2_seq)))
  pred_CH4_cop <- c(ch4_seq,
                    rep(tail(ch4_seq,1), h - length(ch4_seq)))

  pred_GWP_cop <- pred_GWP_no
  pred_GWP_cop$mean <- pred_GWP_no$mean -
                       pred_GWP_no$mean +
                       stats::ts(pred_CO2_cop +
                                 28/1000*pred_CH4_cop +
                                 265/1000*pred_N2O_no,
                                 start = end_train_year + 1,
                                 frequency = 1)

  ratio <- pred_GWP_cop$mean / pred_GWP_no$mean
  pred_GWP_cop$lower <- pred_GWP_cop$mean -
                        (pred_GWP_no$mean - pred_GWP_no$lower)*ratio
  pred_GWP_cop$upper <- pred_GWP_cop$mean -
                        (pred_GWP_no$mean - pred_GWP_no$upper)*ratio

  n_train <- nrow(climate_train)
  years_all <- seq(start_train_year, forecast_end_year)

  GWP <- data.frame(
    Year      = years_all,
    GWP_no    = c(pred_GWP_no$x, pred_GWP_no$mean),
    GWP_cop   = c(pred_GWP_cop$x, pred_GWP_cop$mean),
    GWP_no_80_l = c(rep(NA, n_train-1),
                    climate_train$GWP[n_train],
                    pred_GWP_no$lower[,1]),
    GWP_no_95_l = c(rep(NA, n_train-1),
                    climate_train$GWP[n_train],
                    pred_GWP_no$lower[,2]),
    GWP_no_80_u = c(rep(NA, n_train-1),
                    climate_train$GWP[n_train],
                    pred_GWP_no$upper[,1]),
    GWP_no_95_u = c(rep(NA, n_train-1),
                    climate_train$GWP[n_train],
                    pred_GWP_no$upper[,2]),
    GWP_cop_80_l = c(rep(NA, n_train-1),
                     climate_train$GWP[n_train],
                     pred_GWP_cop$lower[,1]),
    GWP_cop_95_l = c(rep(NA, n_train-1),
                     climate_train$GWP[n_train],
                     pred_GWP_cop$lower[,2]),
    GWP_cop_80_u = c(rep(NA, n_train-1),
                     climate_train$GWP[n_train],
                     pred_GWP_cop$upper[,1]),
    GWP_cop_95_u = c(rep(NA, n_train-1),
                     climate_train$GWP[n_train],
                     pred_GWP_cop$upper[,2])
  )

  GWP$GWP_no_99_u  <- GWP$GWP_no +
                      (GWP$GWP_no_95_u - GWP$GWP_no)/qnorm(0.975)*qnorm(0.995)
  GWP$GWP_no_99_l  <- GWP$GWP_no -
                      (GWP$GWP_no_95_l - GWP$GWP_no)/qnorm(0.975)*qnorm(0.995)
  GWP$GWP_cop_99_u <- GWP$GWP_cop +
                      (GWP$GWP_cop_95_u - GWP$GWP_cop)/qnorm(0.975)*qnorm(0.995)
  GWP$GWP_cop_99_l <- GWP$GWP_cop -
                      (GWP$GWP_cop_95_l - GWP$GWP_cop)/qnorm(0.975)*qnorm(0.995)

  return(GWP)
}

get_old_GWP_forecast_perma <- function(file='data/GWP_perma_compare.csv') {
  #' returns a data frame with the following columns:
  #' "Year","GWP_no_mean_perma","GWP_no_mean","GWP_cop_mean_perma","GWP_cop_mean"

  read_csv(file)
}

read_co2_monthly <- function(file_path="data/monthly_in_situ_co2_mlo.csv",
                                  from_online=FALSE,
                                  convert_to_yearly=TRUE){
  #' Read and process CO2 data from https://scrippsco2.ucsd.edu/data/atmospheric_co2/mlo.html.
  #' @description
  #' This function is an enhanced version tailored to the Mauna Loa CO2 data
  #' format. It can read from a local file or download the latest version online.
  #' It reads metadata from lines 62 (base headers), 63 (data types),
  #' and 64 (units). It can also aggregate the final data to yearly means.
  #'
  #' @param file_path The path to the input CSV file.s
  #' @param from_online Logical indicating whether to read from the internet.
  #'   If `TRUE`, it downloads the file from the Scripps CO2 website.
  #' @param convert_to_yearly Logical indicating whether to aggregate the monthly
  #'   data into yearly means.
  #'
  #' @return A 'tibble' containing the data, either monthly or yearly.
  #'   Expected columns include: "Yr", "Mn", "Date...3", "Date...4", "CO2", "seasonally (adjusted)",
  #'   "fit", "seasonally (adjusted fit)", "CO2 (filled)", "seasonally (adjusted filled)",
  #'   "X11"
  #'   The units are stored in the 'units' attribute of the tibble.
  #' @examples
  #' \dontrun{
  #'   # Get monthly data from a local file
  #'   monthly_co2_data <- read_co2_monthly(file_path = "monthly_in_situ_co2_mlo.csv", convert_to_yearly = FALSE)
  #'   head(monthly_co2_data)
  #'   attr(monthly_co2_data, "units")
  #' }

  if (from_online || !file.exists(file_path)) {
    online_path <- "https://scrippsco2.ucsd.edu/assets/data/atmospheric/stations/in_situ_co2/monthly/monthly_in_situ_co2_mlo.csv"
    download.file(online_path, file_path, mode = "wb")
    message("Reading data from online source.")
  }

  meta_lines <- read_lines(file_path, skip = 61, n_max = 3)

  parse_meta_line <- function(line) {
    line %>%
      str_split(",") %>%
      unlist() %>%
      str_trim() %>%
      str_replace_all("\"", "")
  }

  base_names <- parse_meta_line(meta_lines[1])
  types      <- parse_meta_line(meta_lines[2])
  units      <- parse_meta_line(meta_lines[3])

  final_col_names <- ifelse(
    !is.na(types) & types != "",
    paste0(base_names, " (", types, ")"),
    base_names
  )

  data <- read_csv(
    file = file_path,
    skip = 64,
    col_names = final_col_names,
    col_types = cols()
  )

  if (length(units) == ncol(data)) {
    names(units) <- final_col_names
  } else {
    warning("Number of units in metadata does not match the number of data columns.")
  }
  attr(data, "units") <- units

  if (convert_to_yearly) {

    year_col <- base_names[1]
    month_col <- base_names[2]

    yearly_data <- data %>%
      group_by(.data[[year_col]]) %>%
      summarise(
        across(
          .cols = where(is.numeric) & !matches(paste0("^", month_col)),
          .fns = ~mean(.x, na.rm = FALSE)
        )
      ) %>%
      ungroup()

    attr(yearly_data, "units") <- attr(data, "units")

    data <- yearly_data
  }

  return(data)
}

read_co2_yearly <- function(file_path="data/yearly_in_situ_co2_mlo.csv",
                                   from_online=FALSE) {

  monthly_data <- read_co2_monthly(file_path, from_online, convert_to_yearly=TRUE)
  return(monthly_data)
}

read_ch4_yearly <- function(file_path = "data/ch4_annmean_gl.csv",
                            from_online = FALSE) {
  #' Read yearly global mean CH4 data from NOAA (https://gml.noaa.gov/ccgg/trends_ch4/)
  #'
  #' @param file_path The path to the input CSV file. This is not used if `from_online` is `TRUE`.
  #' @param from_online Logical indicating whether to read the file from the internet.
  #'
  #' @return A 'tibble' containing the yearly CH4 data. The units ('ppb') are
  #'   stored in the 'units' attribute of the tibble. The expected columns are:
  #'   "year", "mean", "unc"
  #'
  #' @examples
  #' \dontrun{
  #'   ch4_data <- read_ch4_yearly()
  #'   head(ch4_data)
  #'   attr(ch4_data, "units")
  #' }

  if (from_online || !file.exists(file_path)) {
    online_path <- "https://gml.noaa.gov/webdata/ccgg/trends/ch4/ch4_annmean_gl.csv"
    download.file(online_path, file_path, mode = "wb")
    message("Reading CH4 data from online source: ", online_path)
  }

  header_line <- read_lines(file_path, skip = 43, n_max = 1)

  col_names <- header_line %>%
    str_split(",") %>%
    unlist() %>%
    str_trim() %>%
    str_replace_all("\"", "")

  data <- read_csv(
    file = file_path,
    skip = 45,
    col_names = col_names,
    col_types = cols()
  )

  units <- c("NA", "ppb", "ppb")

  if (length(units) == length(col_names)) {
    names(units) <- col_names
  } else {
    warning("The number of hardcoded units does not match the number of columns in the data.")
  }
  attr(data, "units") <- units

  return(data)
}

read_n2o_yearly <- function(file_path = "data/GML_global_N2O.txt",
                            from_online = FALSE,
                            convert_to_yearly = TRUE) {
  #' Read and process N2O data from NOAA/GML.
  #'
  #' @description
  #' Reads Nitrous Oxide (N2O) data from the NOAA Global Monitoring Lab
  #' (https://gml.noaa.gov/hats/combined/N2O.html).
  #' This version uses base R's `read.table` for robust parsing of the
  #' whitespace-delimited source file. The source file contains monthly data,
  #' which this function aggregates into yearly means by default.
  #'
  #' @param file_path The path to the input .txt file. Not used if `from_online` is `TRUE`.
  #' @param from_online Logical indicating whether to read the file from the internet.
  #' @param convert_to_yearly Logical indicating whether to aggregate the monthly
  #'   data into yearly means. Defaults to `TRUE`.
  #'
  #' @return A 'tibble' containing the N2O data. By default, data is yearly.
  #'   The units ('ppb') are stored in the 'units' attribute of the tibble.
  #'
  #' @examples
  #' \dontrun{
  #'   # This should now correctly parse all columns
  #'   n2o_data <- read_n2o_yearly(from_online = TRUE, convert_to_yearly = FALSE)
  #'   print(names(n2o_data))
  #'   print(dim(n2o_data))
  #' }

  if (from_online || !file.exists(file_path)) {
    online_path <- "https://gml.noaa.gov/aftp/data/hats/n2o/combined/GML_global_N2O.txt"
    download.file(online_path, file_path, mode = "wb")
    message("Reading N2O data from online source: ", online_path)
  }

  data <- utils::read.table(
    file = file_path,
    skip = 67,
    header = TRUE,
    sep = "",
    na.strings = "nan"
  ) %>%
    dplyr::as_tibble()

  col_names <- names(data)
  units <- rep("ppb", length(col_names))

  non_measurement_cols <- c("GML_N2O_YYYY", "GML_N2O_MM", "GML_N2O_Programs")
  units[col_names %in% non_measurement_cols] <- "NA"

  names(units) <- col_names
  attr(data, "units") <- units

  if (convert_to_yearly) {
    year_col <- "GML_N2O_YYYY"
    month_col <- "GML_N2O_MM"

    yearly_data <- data %>%
      group_by(.data[[year_col]]) %>%
      summarise(across(
          .cols = where(is.numeric) & !matches(paste0("^", month_col, "$")),
          .fns = ~mean(.x, na.rm = TRUE)
        )) %>%
      ungroup()

    attr(yearly_data, "units") <- attr(data, "units")
    data <- yearly_data
  }

  return(data)
}

read_GWP_yearly <- function() {
  #' GWP is calculated as GWP = CO2 + 28/1000 * CH4 + 265/1000 * N2O. This returns
  #' the real historical data of the GWP values as well as CO2, CH4, and N2O.
  #' @return A tibble containing "Year", "CO2", "CH4", "N2O", and "GWP" columns.

  CO2 <- read_co2_yearly() %>%
    dplyr::select(Yr, CO2) %>%
    dplyr::rename(Year = Yr)
  CH4 <- read_ch4_yearly() %>%
    dplyr::select(year, mean) %>%
    dplyr::rename(Year = year, CH4 = mean)
  N2O <- read_n2o_yearly() %>%
    dplyr::select(GML_N2O_YYYY, GML_Global_N2O) %>%
    dplyr::rename(Year = GML_N2O_YYYY, N2O = GML_Global_N2O)

  GWP <- CO2 %>%
    dplyr::left_join(CH4, by = "Year") %>%
    dplyr::left_join(N2O, by = "Year") %>%
    dplyr::mutate(GWP = CO2 + 28/1000 * CH4 + 265/1000 * N2O)

  return(GWP)
}

plot_gwp_forecast_comparison <- function(old_forecast_file_path = 'data/YearData.csv',
                                         start_year_plot = 2010,
                                         end_year_plot = 2050,
                                         end_train_year = 2019,
                                         ssp_file_path = 'data/Greenhouse/merged_all_ssp.csv') {
  #' Plot GWP forecasts against actual data, including SSP scenarios.
  #'
  #' @description
  #' This function visualizes Global Warming Potential (GWP) forecasts from
  #' `get_old_GWP_forecast` and compares them with more recent actual GWP data
  #' obtained from `read_GWP_yearly`, as well as SSP scenarios from the provided CSV.
  #' It helps assess how well the old predictions align with recent trends and which forecast scenario (No Intervention vs. COP26 vs. SSPs)
  #' appears closer to reality.
  #'
  #' The plot displays:
  #' - Recent actual GWP values.
  #' - The "No Intervention" scenario forecast and its 95% confidence interval.
  #' - The "COP26" scenario forecast and its 95% confidence interval.
  #' - All SSP scenario trajectories as dashed lines with unique colors and shapes.
  #'
  #' @param old_forecast_file_path Path to the CSV file used by `get_old_GWP_forecast`.
  #'                               Defaults to 'data/YearData.csv'.
  #' @param start_year_plot The starting year for the x-axis of the plot.
  #'                        Defaults to 2010.
  #' @param end_year_plot The ending year for the x-axis of the plot.
  #'                      Defaults to 2050.
  #' @param ssp_file_path Path to the SSP scenario CSV file.
  #'                      Defaults to 'data/Greenhouse/merged_all_ssp.csv'.
  #'
  #' @return A ggplot object representing the comparison plot.
  #'
  #' @examples
  #' \dontrun{
  #'   comparison_plot <- plot_gwp_forecast_comparison()
  #'   print(comparison_plot)
  #' }

  actual_gwp_data <- read_GWP_yearly() %>%
    dplyr::select(Year, GWP) %>%
    dplyr::rename(Actual_GWP = GWP) %>%
    dplyr::filter(Year >= start_year_plot, Year <= end_year_plot)

  old_forecasts_data <- get_old_GWP_forecast(file_path = old_forecast_file_path,
                                            end_train_year = end_train_year) %>%
    dplyr::filter(Year >= start_year_plot, Year <= end_year_plot)

  forecast_ci_start_year <- end_train_year + 1

  ssp_data <- readr::read_csv(ssp_file_path, show_col_types = FALSE) %>%
    dplyr::filter(year >= start_year_plot, year <= end_year_plot) %>%
    tidyr::pivot_longer(
      cols = -year,
      names_to = "SSP",
      values_to = "value"
    )

  ssp_scenarios <- unique(ssp_data$SSP)
  ssp_colors <- RColorBrewer::brewer.pal(min(length(ssp_scenarios), 8), "Dark2")
  if (length(ssp_scenarios) > 8) {
    ssp_colors <- c(ssp_colors, RColorBrewer::brewer.pal(length(ssp_scenarios) - 8, "Set2"))
  }
  names(ssp_colors) <- ssp_scenarios
  ssp_shapes <- seq(0, length(ssp_scenarios) - 1) %% 25
  names(ssp_shapes) <- ssp_scenarios

  gwp_plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = Year)) +

    ggplot2::geom_line(data = old_forecasts_data,
                       ggplot2::aes(y = GWP_no, color = "Forecast (No Intervention)"),
                       linewidth = 0.8, linetype = "dashed") +
    ggplot2::geom_point(data = old_forecasts_data,
                        ggplot2::aes(y = GWP_no, color = "Forecast (No Intervention)"),
                        size = 1.5) +
    ggplot2::geom_ribbon(data = old_forecasts_data %>% dplyr::filter(Year >= forecast_ci_start_year),
                         ggplot2::aes(ymin = GWP_no_95_l, ymax = GWP_no_95_u,
                                      fill = "No Intervention (95% CI)"),
                         alpha = 0.2) +

    ggplot2::geom_line(data = old_forecasts_data,
                       ggplot2::aes(y = GWP_cop, color = "Forecast (COP26)"),
                       linewidth = 0.8, linetype = "dotted") +
    ggplot2::geom_point(data = old_forecasts_data,
                        ggplot2::aes(y = GWP_cop, color = "Forecast (COP26)"),
                        size = 1.5) +
    ggplot2::geom_ribbon(data = old_forecasts_data %>% dplyr::filter(Year >= forecast_ci_start_year),
                         ggplot2::aes(ymin = GWP_cop_95_l, ymax = GWP_cop_95_u,
                                      fill = "COP26 (95% CI)"),
                         alpha = 0.2) +

    ggplot2::geom_line(data = actual_gwp_data,
                       ggplot2::aes(y = Actual_GWP, color = "Actual GWP (Recent)"),
                       linewidth = 1.2) +
    ggplot2::geom_point(data = actual_gwp_data,
                        ggplot2::aes(y = Actual_GWP, color = "Actual GWP (Recent)"),
                        size = 2) +

    ggplot2::geom_line(
      data = ssp_data,
      mapping = ggplot2::aes(x = year, y = value, color = SSP, linetype = SSP),
      linewidth = 0.9
    ) +

    ggplot2::scale_color_manual(
      name = "Data Series:",
      values = c(
        "Actual GWP (Recent)" = "black",
        "Forecast (No Intervention)" = "dodgerblue3",
        "Forecast (COP26)" = "forestgreen",
        ssp_colors
      )
    ) +
    ggplot2::scale_shape_manual(
      name = "SSP Scenario:",
      values = ssp_shapes
    ) +
    ggplot2::scale_linetype_manual(
      name = "SSP Scenario:",
      values = c(
        "Forecast (No Intervention)" = "dashed",
        "Forecast (COP26)" = "dotted",
        setNames(rep("longdash", length(ssp_scenarios)), ssp_scenarios)
      ),
      guide = ggplot2::guide_legend(override.aes = list(linewidth = 0.9))
    ) +
    ggplot2::scale_fill_manual(
      name = "Confidence Intervals:",
      values = c(
        "No Intervention (95% CI)" = "lightblue",
        "COP26 (95% CI)" = "lightgreen"
      )
    ) +
    ggplot2::labs(
      title = "GWP Forecast Comparison: Actual, Predicted, and SSP Scenarios",
      subtitle = paste0("Displaying data from ", start_year_plot, " to ", end_year_plot),
      x = "Year",
      y = "Global Warming Potential (GWP)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(xlim = c(start_year_plot, end_year_plot), expand = TRUE)

  return(gwp_plot)
}

plot_gwp_forecast_comparison_perma <- function(old_forecast_file_path = 'data/GWP_perma_compare.csv',
                                               start_year_plot = 2010,
                                               end_year_plot = 2050,
                                               ssp_file_path = 'data/Greenhouse/merged_all_ssp.csv') {
  #' Plot GWP forecasts (including permafrost effects) against actual data and SSPs.
  #'
  #' @description
  #' This function visualizes Global Warming Potential (GWP) forecasts from
  #' `get_old_GWP_forecast_perma` and compares them with recent actual GWP data
  #' from `read_GWP_yearly`, as well as SSP scenarios. It is designed to show
  #' the impact of including permafrost thaw in the forecasts.
  #'
  #' The plot displays:
  #' - Recent actual GWP values.
  #' - "No Intervention" forecast (with and without permafrost).
  #' - "COP26" forecast (with and without permafrost).
  #' - All SSP scenario trajectories.
  #'
  #' @param old_forecast_file_path Path to the CSV file for `get_old_GWP_forecast_perma`.
  #'                               Defaults to 'data/GWP_perma_compare.csv'.
  #' @param start_year_plot The starting year for the plot's x-axis. Defaults to 2010.
  #' @param end_year_plot The ending year for the plot's x-axis. Defaults to 2050.
  #' @param ssp_file_path Path to the SSP scenario CSV file.
  #'                      Defaults to 'data/Greenhouse/merged_all_ssp.csv'.
  #'
  #' @return A ggplot object representing the comparison plot.
  #'
  #' @examples
  #' \dontrun{
  #'   perma_plot <- plot_gwp_forecast_comparison_perma()
  #'   print(perma_plot)
  #' }

  actual_gwp_data <- read_GWP_yearly() %>%
    dplyr::select(Year, GWP) %>%
    dplyr::rename(Actual_GWP = GWP) %>%
    dplyr::filter(Year >= start_year_plot, Year <= end_year_plot)

  old_forecasts_data_long <- get_old_GWP_forecast_perma(file = old_forecast_file_path) %>%
    dplyr::filter(Year >= start_year_plot, Year <= end_year_plot) %>%
    tidyr::pivot_longer(
      cols = -Year,
      names_to = "Forecast_Type",
      values_to = "GWP"
    ) %>%
    dplyr::mutate(
      Scenario = dplyr::case_when(
        stringr::str_detect(Forecast_Type, "_no_") ~ "No Intervention",
        stringr::str_detect(Forecast_Type, "_cop_") ~ "COP26"
      ),
      Model = dplyr::case_when(
        stringr::str_detect(Forecast_Type, "_perma") ~ "Permafrost",
        TRUE ~ "Standard"
      ),
      Legend = paste0(Scenario, " (", Model, ")")
    )

  ssp_data <- readr::read_csv(ssp_file_path, show_col_types = FALSE) %>%
    dplyr::filter(year >= start_year_plot, year <= end_year_plot) %>%
    tidyr::pivot_longer(
      cols = -year,
      names_to = "SSP",
      values_to = "value"
    )

  combined_forecast_data <- old_forecasts_data_long %>%
    dplyr::select(Year, GWP, Legend) %>%
    dplyr::rename(value = GWP, year = Year, series_name = Legend)

  combined_ssp_data <- ssp_data %>%
    dplyr::rename(series_name = SSP)

  combined_actual_data <- actual_gwp_data %>%
    dplyr::mutate(series_name = "Actual GWP (Recent)") %>%
    dplyr::rename(value = Actual_GWP, year = Year)

  all_plot_data <- dplyr::bind_rows(
    combined_forecast_data,
    combined_ssp_data,
    combined_actual_data
  )

  ssp_scenarios <- unique(ssp_data$SSP)
  ssp_colors <- RColorBrewer::brewer.pal(min(length(ssp_scenarios), 8), "Dark2")
  if (length(ssp_scenarios) > 8) {
    ssp_colors <- c(ssp_colors, RColorBrewer::brewer.pal(length(ssp_scenarios) - 8, "Set2"))
  }
  names(ssp_colors) <- ssp_scenarios

  color_values <- c(
    "Actual GWP (Recent)" = "black",
    "No Intervention (Standard)" = "dodgerblue3",
    "No Intervention (Permafrost)" = "deepskyblue",
    "COP26 (Standard)" = "forestgreen",
    "COP26 (Permafrost)" = "limegreen",
    ssp_colors
  )

  linetype_values <- c(
    "Actual GWP (Recent)" = "solid",
    "No Intervention (Standard)" = "solid",
    "No Intervention (Permafrost)" = "solid",
    "COP26 (Standard)" = "solid",
    "COP26 (Permafrost)" = "solid",
    setNames(rep("longdash", length(ssp_scenarios)), ssp_scenarios)
  )

  linewidth_values <- c(
    "Actual GWP (Recent)" = 1.2,
    "No Intervention (Standard)" = 1,
    "No Intervention (Permafrost)" = 1,
    "COP26 (Standard)" = 1,
    "COP26 (Permafrost)" = 1,
    setNames(rep(0.9, length(ssp_scenarios)), ssp_scenarios)
  )

  gwp_plot <- ggplot2::ggplot(
      all_plot_data,
      ggplot2::aes(x = year, y = value, color = series_name, linetype = series_name, size = series_name)
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      data = all_plot_data %>% dplyr::filter(series_name == "Actual GWP (Recent)"),
      size = 2
    ) +

    ggplot2::scale_color_manual(name = "Data Series:", values = color_values) +
    ggplot2::scale_linetype_manual(name = "Data Series:", values = linetype_values) +
    ggplot2::scale_size_manual(name = "Data Series:", values = linewidth_values) +

    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 0.9)),
      linetype = ggplot2::guide_legend(override.aes = list(size = 0.9)),
      size = "none"
    ) +

    ggplot2::labs(
      title = "GWP Forecast Comparison: Permafrost, Actual, and SSP Scenarios",
      subtitle = paste0("Displaying data from ", start_year_plot, " to ", end_year_plot),
      x = "Year",
      y = "Global Warming Potential (GWP)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(xlim = c(start_year_plot, end_year_plot), expand = TRUE)

  return(gwp_plot)
}