source("data_processing.R")
source("VLM.R")
source("utils.R")
source("models/predict_GMSL.R")
library(broom)

generate_ardl_predictions <- function(fit_NY, model_data, pred_GMSL, use_VLM_as_variable, ar6_VLM, gmsl_col, start_year, end_year) {
    beta1 <- coef(fit_NY)["NY_lag"]
    beta2 <- coef(fit_NY)["GMSL"]
    sigma_NY <- sigma(fit_NY)

    last_obs_year <- max(model_data$Year[!is.na(model_data$NY)])
    last_obs_NY <- model_data$NY[model_data$Year == last_obs_year]

    historical_gmsl_var <- var(model_data$GMSL[!is.na(model_data$GMSL) &
                                              model_data$Year >= start_year &
                                              model_data$Year <= last_obs_year])

    prediction_start <- max(end_year, last_obs_year) + 1
    forecast_years <- prediction_start:2100
    pred_length <- length(forecast_years)

    historical_data <- model_data[!is.na(model_data$NY) & !is.na(model_data$GMSL), ]
    historical_data <- historical_data %>%
        mutate(NY_lag = dplyr::lag(NY)) %>%
        na.omit()
    rho <- cor(historical_data$NY_lag, historical_data$GMSL)

    Y_no <- numeric(pred_length)
    var_Y <- numeric(pred_length)
    sigma_no <- numeric(pred_length)

    initial_idx <- 1
    Y_no[initial_idx] <- last_obs_NY

    new_data <- data.frame(
        NY_lag = Y_no[initial_idx],
        GMSL = pred_GMSL[[gmsl_col]][1]
    )
    if (use_VLM_as_variable) {
        new_data$VLM <- ar6_VLM$VLM[ar6_VLM$Year == forecast_years[1]]
    }

    param_uncertainty <- predict(fit_NY, newdata = new_data, se.fit = TRUE)$se.fit^2
    var_Y[initial_idx] <- param_uncertainty + sigma_NY^2
    sigma_no[initial_idx] <- sqrt(var_Y[initial_idx])

    for (i in (initial_idx+1):pred_length) {
        current_year <- forecast_years[i]
        year_idx <- which(pred_GMSL$Year == current_year)

        if (length(year_idx) == 0) next

        new_data <- data.frame(
            NY_lag = Y_no[i-1],
            GMSL = pred_GMSL[[gmsl_col]][year_idx]
        )
        if (use_VLM_as_variable) {
            new_data$VLM <- ar6_VLM$VLM[ar6_VLM$Year == current_year]
        }

        Y_no[i] <- predict(fit_NY, newdata = new_data)
        param_uncertainty <- predict(fit_NY, newdata = new_data, se.fit = TRUE)$se.fit^2

        if (!is.na(pred_GMSL$up95_no[year_idx]) && !is.na(pred_GMSL$low95_no[year_idx])) {
            var_GMSL <- ((pred_GMSL$up95_no[year_idx] - pred_GMSL$low95_no[year_idx]) / (2 * 1.96))^2
        } else {
            var_GMSL <- historical_gmsl_var
        }

        var_Y[i] <- beta1^2 * var_Y[i-1] +
                    beta2^2 * var_GMSL +
                    2 * rho * beta1 * beta2 * sqrt(var_Y[i-1] * var_GMSL) +
                    param_uncertainty +
                    sigma_NY^2

        sigma_no[i] <- sqrt(var_Y[i])
    }

    predictions <- data.frame(
        Year = forecast_years,
        mean_no = Y_no,
        up95_no = Y_no + qnorm(0.975) * sigma_no,
        low95_no = Y_no - qnorm(0.975) * sigma_no,
        up99_no = Y_no + qnorm(0.995) * sigma_no,
        low99_no = Y_no - qnorm(0.995) * sigma_no
    )

    return(predictions)
}

RSL_predict_location <- function(...) {
  RSL_predict(...)
}

ARDL_prediction <- function(...) {
  ARDL_predictions <- RSL_predict(...)
  return(ARDL_predictions$prediction)
}

RSL_predict_single <- function(location,
                             measurement.type = "tg",
                             add_VLM = FALSE,
                             VLM_na_interp = TRUE,
                             na.rm = FALSE,
                             equal_time = FALSE,
                             VLM_method = 0,
                             start_year = 1970,
                             end_year = 2019,
                             use_VLM_as_variable = FALSE,
                             return_fitted_history_only = FALSE,
                             scenario = "no",
                             ...) {
    #' Function to predict relative sea level (RSL) for a given location and measurement type.
    #' @param scenario Character string indicating the scenario for GMSL prediction. Options are "no", "cop", or a specific SSP scenario.
    #' Possible values include "no", "cop", or and columns in `data/GMSL_prediction_SSP_perma.csv` (ssp119,ssp126,ssp245,ssp370,ssp460,ssp585)
    climate <- get_global_sealevel_raw(return_yearly = TRUE)
    pred_GMSL <- read.csv('data/GMSL_prediction_SSP_perma.csv', header = TRUE)
    if (start_year < 1950) {
        pred_GMSL <- merge_global_sealevel_prediction()
    }

    dat <- get_year_data_raw(location,
                             add_VLM = add_VLM,
                             VLM_na_interp = VLM_na_interp,
                             equal_time = equal_time,
                             VLM_method = VLM_method,
                             use_VLM_as_variable = FALSE,
                             ...)

    annual_data <- dat %>%
        dplyr::filter(Year >= start_year, Year <= end_year) %>%
        mutate(sealevel = case_when(
            measurement.type == "tg" ~ tg,
            measurement.type == "sa" ~ sa,
            TRUE ~ NA_real_
        )) %>%
        group_by(Year) %>%
        summarise(sealevel = mean(sealevel, na.rm = na.rm)) %>%
        ungroup()

    model_data <- climate %>%
        dplyr::select(Year, GMSL) %>%
        right_join(annual_data, by = "Year") %>%
        rename(NY = sealevel)

    if (use_VLM_as_variable) {
        ar6_VLM <- get_ar6_vlm(station = location,
                               start.year = start_year,
                               end.year = 2100)
        model_data <- model_data %>%
            left_join(ar6_VLM, by = "Year")
    }

    if (scenario == "no") {
        gmsl_col <- "mean_no"
    } else if (scenario == "cop") {
        gmsl_col <- "mean_cop"
    } else {
        gmsl_col <- scenario
    }




    model_data_lagged <- model_data %>%
        mutate(
            NY_lag = dplyr::lag(NY, n = 1),
            GMSL_lag = dplyr::lag(GMSL, n = 1)
        ) %>%
        na.omit()

    if (use_VLM_as_variable) {
        fit_NY <- lm(NY ~ NY_lag + GMSL + VLM, data = model_data_lagged)
    } else {
        fit_NY <- lm(NY ~ NY_lag + GMSL, data = model_data_lagged)
    }

    if (return_fitted_history_only) {
        return(augment(fit_NY))
    }

    predictions <- generate_ardl_predictions(
        fit_NY = fit_NY,
        model_data = model_data,
        pred_GMSL = pred_GMSL,
        use_VLM_as_variable = use_VLM_as_variable,
        ar6_VLM = if(use_VLM_as_variable) ar6_VLM else NULL,
        gmsl_col = gmsl_col,
        start_year = start_year,
        end_year = end_year
    )

    return(list(
        prediction = predictions,
        measurement = measurement.type,
        withVLM = add_VLM
    ))
}

RSL_predict <- function(location,
                       measurement.type = "tg",
                       add_VLM = FALSE,
                       VLM_na_interp = TRUE,
                       na.rm = FALSE,
                       equal_time = FALSE,
                       VLM_method = 0,
                       start_year = 1970,
                       end_year = 2019,
                       use_VLM_as_variable = FALSE,
                       return_fitted_history_only = FALSE,
                       ...) {

    pred_no <- RSL_predict_single(
        location = location,
        measurement.type = measurement.type,
        add_VLM = add_VLM,
        VLM_na_interp = VLM_na_interp,
        na.rm = na.rm,
        equal_time = equal_time,
        VLM_method = VLM_method,
        start_year = start_year,
        end_year = end_year,
        use_VLM_as_variable = use_VLM_as_variable,
        return_fitted_history_only = return_fitted_history_only,
        scenario = "no",
        ...
    )

    pred_cop <- RSL_predict_single(
        location = location,
        measurement.type = measurement.type,
        add_VLM = add_VLM,
        VLM_na_interp = VLM_na_interp,
        na.rm = na.rm,
        equal_time = equal_time,
        VLM_method = VLM_method,
        start_year = start_year,
        end_year = end_year,
        use_VLM_as_variable = use_VLM_as_variable,
        return_fitted_history_only = return_fitted_history_only,
        scenario = "cop",
        ...
    )

    combined_pred <- pred_no$prediction %>%
        dplyr::select(Year,
                     mean_no, up95_no, low95_no, up99_no, low99_no) %>%
        dplyr::left_join(
            pred_cop$prediction %>%
                dplyr::rename_with(
                    ~gsub("_no", "_cop", .x),
                    matches("_no")
                ),
            by = "Year"
        )

    return(list(
        prediction = combined_pred,
        measurement = measurement.type,
        withVLM = add_VLM
    ))
}

apply_vlm_adjustment_ardl_backtesting <- function(df, vlm_data) {
  df %>%
    left_join(vlm_data, by = "Year") %>%
        mutate(across(-starts_with("ARDL_"), ~ . + VLM, .names = "{.col}"),
               ) %>%
    dplyr::select(-VLM)
}

backtesting_ARDL_v2 <- function(location,
                            start_year = 1970,
                            end_year = 2000,
                            horizon = 2023,
                            measurement.type = "tg",
                            use_VLM_as_variable = FALSE,
                            ...) {
    historical_fit <- RSL_predict_single(
        location = location,
        measurement.type = measurement.type,
        start_year = start_year,
        end_year = end_year,
        use_VLM_as_variable = use_VLM_as_variable,
        return_fitted_history_only = TRUE,
        ...
    )

    predictions <- RSL_predict_single(
        location = location,
        measurement.type = measurement.type,
        start_year = start_year,
        end_year = end_year,
        use_VLM_as_variable = use_VLM_as_variable,
        scenario = "no",
        ...
    )$prediction

    dat <- get_month_data_raw(location, add_VLM = FALSE, ...)
    years <- as.numeric(substr(dat$month, 1, 4))
    sealev <- if(measurement.type == "tg") dat$tg else if(measurement.type == "sa") dat$sa else NA

    month_sealevel <- data.frame(time = years, sealevel = sealev)
    annual_sealevel <- month_sealevel %>%
        group_by(time) %>%
        summarise(means = mean(sealevel, na.rm = TRUE))
    names(annual_sealevel) <- c('Year', 'NY')


    results <- data.frame(
        Year = start_year:horizon,
        LSL_true = annual_sealevel$NY[match(start_year:horizon, annual_sealevel$Year)]
    )

    historical_years <- start_year:end_year
    fitted_years <- historical_years[-1]
    results$LSL_fitted[results$Year %in% fitted_years] <- historical_fit$.fitted

    results$LSL_fitted[results$Year == start_year] <- NA

    future_years <- (end_year + 1):horizon
    future_idx <- which(results$Year %in% future_years)
    pred_idx <- which(predictions$Year %in% future_years)

    results$LSL_mean[future_idx] <- predictions$mean_no[pred_idx]
    results$low_95[future_idx] <- predictions$low95_no[pred_idx]
    results$up_95[future_idx] <- predictions$up95_no[pred_idx]
    results$low_99[future_idx] <- predictions$low99_no[pred_idx]
    results$up_99[future_idx] <- predictions$up99_no[pred_idx]

    return(results)
}