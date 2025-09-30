source("ARDL.R")
source("XVARMA.R")
source("utils.R")
source("ModelBacktesting.R")
library(quadprog)

ensemble_prediction <- function(location, measurement.type = "tg", add_VLM = FALSE, VLM_na_interp = TRUE, na.rm = FALSE, equal_time = TRUE, VLM_method = 0, end_date = '2025-01-01', max.p = 2, start.year = -Inf, stop.year = 2019, horizon = 2100, data_dir = DATA_DIR, restriction = "no", selected_vars = c("LSL", "X10m_u_component_of_wind", "X10m_v_component_of_wind", "X2m_temperature", "surface_pressure"), ...) {
  #' Get ARDL predictions
  #' @return A dataframe. Columns: Year, mean, up95, low95, up99, low99
  ardl_pred <- ARDL_prediction(location, measurement.type, add_VLM, VLM_na_interp, na.rm, equal_time, VLM_method, end_date, ...)

  if (restriction == "no") {
    ardl_pred <- ardl_pred %>%
      dplyr::select(Year, mean_no, up95_no, low95_no, up99_no, low99_no)
  } else if (restriction == "cop") {
    ardl_pred <- ardl_pred %>%
      dplyr::select(Year, mean_cop, up95_cop, low95_cop, up99_cop, low99_cop)
  } else {
    stop("Invalid restriction value. Use 'no' or 'cop'.")
  }

  varx_pred <- VARX_prediction(location, measurement.type=measurement.type, max.p=max.p, start.year=start.year, stop.year=stop.year, horizon=horizon, data_dir=data_dir, restriction=restriction, selected_vars=selected_vars, ...)

  merged_pred <- merge(ardl_pred, varx_pred, by = "Year")


  if (restriction == "no") {
    ensemble_pred <- merged_pred %>%
      dplyr::mutate(
        mean = (mean_no + LSL_mean) / 2,
        up95 = (up95_no + up_95) / 2,
        low95 = (low95_no + low_95) / 2,
        up99 = (up99_no + up_99) / 2,
        low99 = (low99_no + low_99) / 2
      ) %>%
      dplyr::select(Year, mean, up95, low95, up99, low99)
  } else if (restriction == "cop") {
    ensemble_pred <- merged_pred %>%
      dplyr::mutate(
        mean = (mean_cop + LSL_mean) / 2,
        up95 = (up95_cop + up_95) / 2,
        low95 = (low95_cop + low_95) / 2,
        up99 = (up99_cop + up_99) / 2,
        low99 = (low99_cop + low_99) / 2
      ) %>%
      dplyr::select(Year, mean, up95, low95, up99, low99)
  }

  return(ensemble_pred)
}


ensemble_prediction_dynamic <- function(location, measurement.type = "tg", add_VLM = FALSE, VLM_na_interp = TRUE, na.rm = FALSE, equal_time = TRUE, VLM_method = 0, end_date = '2025-01-01', max.p = 2, start.year = -Inf, stop.year = 2019, horizon = 2100, data_dir = DATA_DIR, restriction = "no", selected_vars = c("LSL", "X10m_u_component_of_wind", "X10m_v_component_of_wind", "X2m_temperature", "surface_pressure"), ...) {
  #' Get ARDL predictions
  #' @return A dataframe. Columns: Year, mean, up95, low95, up99, low99
  ardl_pred <- ARDL_prediction(location, measurement.type, add_VLM, VLM_na_interp, na.rm, equal_time, VLM_method, end_date, ...)

  if (restriction == "no") {
    ardl_pred <- ardl_pred %>%
      dplyr::select(Year, mean_no, up95_no, low95_no, up99_no, low99_no)
  } else if (restriction == "cop") {
    ardl_pred <- ardl_pred %>%
      dplyr::select(Year, mean_cop, up95_cop, low95_cop, up99_cop, low99_cop)
  } else {
    stop("Invalid restriction value. Use 'no' or 'cop'.")
  }
  ardl_pred[,which(names(ardl_pred)!="Year")] = ardl_pred[,which(names(ardl_pred)!="Year")]-ardl_pred$mean_no[which(ardl_pred$Year==(stop.year+1))]
  varx_pred <- VARX_prediction(location, max.p=max.p, start.year=start.year, stop.year=stop.year+1, horizon=horizon, data_dir=data_dir, restriction=restriction, selected_vars=selected_vars, ...)
  varx_pred[,which(names(varx_pred)!="Year")] = varx_pred[,which(names(varx_pred)!="Year")]-varx_pred$LSL_mean[which(varx_pred$Year==(stop.year+1))]
  merged_pred <- merge(ardl_pred, varx_pred, by = "Year")
  merged_pred <- merged_pred[which(merged_pred$Year>stop.year),]
  varx_weight <- seq(1,by=-0.0125,length.out = length(merged_pred$LSL_mean))
  ardl_weight <- 1-varx_weight

  if (restriction == "no") {
    ensemble_pred <- merged_pred %>%
      dplyr::mutate(
        mean = mean_no*ardl_weight + LSL_mean*varx_weight,
        up95 = (up95_no*ardl_weight + up_95*varx_weight),
        low95 = (low95_no*ardl_weight + low_95*varx_weight),
        up99 = (up99_no*ardl_weight + up_99*varx_weight),
        low99 = (low99_no*ardl_weight + low_99*varx_weight)
      ) %>%
      dplyr::select(Year, mean, up95, low95, up99, low99)
  } else if (restriction == "cop") {
    ensemble_pred <- merged_pred %>%
      dplyr::mutate(
        mean = (mean_cop*ardl_weight + LSL_mean*varx_weight),
        up95 = (up95_cop*ardl_weight + up_95*varx_weight),
        low95 = (low95_cop*ardl_weight + low_95*varx_weight),
        up99 = (up99_cop*ardl_weight + up_99*varx_weight),
        low99 = (low99_cop*ardl_weight + low_99*varx_weight)
      ) %>%
      dplyr::select(Year, mean, up95, low95, up99, low99)
  }

  return(ensemble_pred)
}

ensemble_dynamic_df_1 <- function(df, ardl_col="ARDL_no", varx_col="VARX_no", output_col="ensemble_no", start.year=2020, stop.year=2100){
  #' The hyperparameters are chosen by `utils/dw.py`
  #' @param df A dataframe. Columns: Year, ardl_col, varx_col, ...
  #' @param ardl_col The column name of the ARDL predictions
  #' @param varx_col The column name of the VARX predictions
  #' @param output_col The column name of the ensemble predictions
  #' @return A dataframe. Columns: Year, output_col, ...

  start.year = min(start.year, min(df$Year))
  stop.year = max(stop.year, max(df$Year))
  ardl_starting_weight = 0
  varx_starting_weight = 1
  ardl_ending_weight = 1
  varx_ending_weight = 0
  ardl_weight = seq(ardl_starting_weight, ardl_ending_weight, length.out = stop.year-start.year+1)
  varx_weight = seq(varx_starting_weight, varx_ending_weight, length.out = stop.year-start.year+1)
  df_mod <- df %>% dplyr::filter(Year >= start.year & Year <= stop.year)
  output = df_mod[,ardl_col]*ardl_weight + df_mod[,varx_col]*varx_weight
  df[df$Year >= start.year & df$Year <= stop.year, output_col] = output
  return(df)
}

ensemble_dynamic_df_2 <- function(df, ardl_col="ARDL_no", varx_col="VARX_no", output_col="ensemble_no", start.year=2020, stop.year=2100, alpha=5.021015958254602, beta=0.07012444677753994) {
  #'
  #' @param df A dataframe. Columns: Year, ardl_col, varx_col, ...
  #' @param ardl_col The column name of the ARDL predictions
  #' @param varx_col The column name of the VARX predictions
  #' @param output_col The column name of the ensemble predictions
  #' @param alpha The alpha parameter for the weight calculation
  #' @param beta The beta parameter for the weight calculation
  #' @return A dataframe. Columns: Year, output_col, ...

  start.year = min(start.year, min(df$Year))
  stop.year = max(stop.year, max(df$Year))
  t_0 = 2020

  years = seq(start.year, stop.year)
  w1_t = 1 / (1 + alpha * exp(-beta * (years - t_0)))
  w2_t = 1 - w1_t

  df_mod <- df %>% dplyr::filter(Year >= start.year & Year <= stop.year)

  output = df_mod[,ardl_col] * w1_t + df_mod[,varx_col] * w2_t
  df[df$Year >= start.year & df$Year <= stop.year, output_col] = output

  for (bound in c("low95", "up95", "low99", "up99")) {
    ardl_bound_col = paste0(ardl_col, "_", bound)
    varx_bound_col = paste0(varx_col, "_", bound)
    output_bound_col = paste0(output_col, "_", bound)

    if (ardl_bound_col %in% names(df_mod) && varx_bound_col %in% names(df_mod)) {
      output_bound = df_mod[,ardl_bound_col] * w1_t + df_mod[,varx_bound_col] * w2_t
      df[df$Year >= start.year & df$Year <= stop.year, output_bound_col] = output_bound
    }
  }

  return(df)
}

ensemble_dynamic_df <- function(...) {
  return(ensemble_dynamic_df_2(...))
}

ensemble_lm_weights_location <- function(location, stop.year, ...){
    ardl <- backtesting_ARDL(location, stop.year=stop.year, ...) %>%
        dplyr::select(Year, LSL_true, LSL_mean) %>%
        rename(ardl = LSL_mean)
    sem <- backtesting_SEM(location, stop.year=stop.year, ...) %>%
        dplyr::select(Year, LSL_mean) %>%
        rename(sem = LSL_mean)
    varx <- backtesting_VARX(location, stop.year=stop.year, ...) %>%
        dplyr::select(Year, LSL_mean) %>%
        rename(varx = LSL_mean)
    data <- merge(ardl, sem, by = "Year")
    data <- merge(data, varx, by = "Year")
    data <- data %>% dplyr::filter(Year <= 2023)
    mod <- lm(LSL_true ~ ardl + sem + varx, data = data)
    return(mod)
}

ensemble_lm_weights_locations <- function(locations=c("astoria", "battery_park", "cape_charles", "charleston",
"crescent_city", "elly_oil_platform", "monterey", "newport",
"south_beach", "tofino"), stop.year=2000, ...){
    datasets <- list()
    for (location in locations) {
        ardl <- backtesting_ARDL(location, stop.year=stop.year, ...) %>%
            dplyr::select(Year, LSL_true, LSL_mean) %>%
            rename(ardl = LSL_mean)
        sem <- backtesting_SEM(location, stop.year=stop.year, ...) %>%
            dplyr::select(Year, LSL_mean) %>%
            rename(sem = LSL_mean)
        varx <- backtesting_VARX(location, stop.year=stop.year, ...) %>%
            dplyr::select(Year, LSL_mean) %>%
            rename(varx = LSL_mean)
        data <- merge(ardl, sem, by = "Year")
        data <- merge(data, varx, by = "Year")
        data <- data %>% dplyr::filter(Year <= 2023)
        datasets[[location]] <- data
    }
    data <- Reduce(function(x, y) bind_rows(x, y), datasets)
    data <- data %>% dplyr::filter(Year <= 2023)
    mod <- lm(LSL_true ~ ardl + sem + varx, data = data)
    return(mod)
}

ensemble_qp_weights_locations <- function(locations=c("astoria", "battery_park", "cape_charles", "charleston",
                                             "crescent_city", "elly_oil_platform", "monterey", "newport",
                                             "south_beach", "tofino"), stop.year=2000, ...){
  datasets <- list()
  for (location in locations) {
    ardl <- backtesting_ARDL(location, stop.year=stop.year, ...) %>%
      dplyr::select(Year, LSL_true, LSL_mean) %>%
      rename(ardl = LSL_mean)
    sem <- backtesting_SEM(location, stop.year=stop.year, ...) %>%
      dplyr::select(Year, LSL_mean) %>%
      rename(sem = LSL_mean)
    varx <- backtesting_VARX(location, stop.year=stop.year, ...) %>%
      dplyr::select(Year, LSL_mean) %>%
      rename(varx = LSL_mean)
    data <- merge(ardl, sem, by = "Year")
    data <- merge(data, varx, by = "Year")
    data <- data %>% dplyr::filter(Year <= 2023)
    datasets[[location]] <- data
  }
  data <- Reduce(function(x, y) bind_rows(x, y), datasets)
  data <- data %>% dplyr::filter(Year <= 2023) %>% drop_na()

  x <- as.matrix(data[, c("ardl", "varx")])
  y <- data$LSL_true
  n <- ncol(x)

  Dmat <- t(x) %*% x
  dvec <- t(x) %*% y
  Amat <- cbind(rep(1, n), diag(n))
  bvec <- c(1, rep(0, n))

  fit <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)

  return(fit$solution)
}
