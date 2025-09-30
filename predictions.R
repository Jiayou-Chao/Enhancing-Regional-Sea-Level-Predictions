source("SEM_Fitting.R")
source("XVARMA.R")
source("VLM.R")
library(reshape2)
source("ensemble.R")
source("data_integrity.R")

ger_ar6_vlm_prediction <- function(station){
  loc_id <- find_closest_ar6_station(station)
  df <- read_csv(paste0("data/ar6/ipcc_ar6_vlm_projection_psmsl_id_", loc_id, ".csv"))
  return(df)
}

.get_ARDL_prediction <- function(station, measurement.type = "sa", add_VLM=TRUE, add_VLM_afterwards=FALSE, start.year=1993, stop.year=2019, ...){
    end_date <- paste0(stop.year, "-12-31")
    vlm_no_2050 <- 0
    vlm_no_2100 <- 0
    vlm_cop_2050 <- 0
    vlm_cop_2100 <- 0
    if (add_VLM_afterwards) {
      add_VLM <- FALSE
      ar6_vlm_prediction <- ger_ar6_vlm_prediction(station)
      vlm_no_2050 <- ar6_vlm_prediction %>% dplyr::filter(Year == 2050) %>% dplyr::select(ssp245_medium_50) %>% pull()
      vlm_no_2100 <- ar6_vlm_prediction %>% dplyr::filter(Year == 2100) %>% dplyr::select(ssp245_medium_50) %>% pull()
      vlm_cop_2050 <- ar6_vlm_prediction %>% dplyr::filter(Year == 2050) %>% dplyr::select(ssp245_medium_50) %>% pull()
      vlm_cop_2100 <- ar6_vlm_prediction %>% dplyr::filter(Year == 2100) %>% dplyr::select(ssp245_medium_50) %>% pull()
    }

    results  <- RSL_predict(station, measurement.type = measurement.type, add_VLM=add_VLM, end_date = end_date, equal_time = TRUE, ...)$prediction
    rsl <- get_baseline(station)
    mean_no_2050 <- results %>% dplyr::filter(Year == 2050) %>% dplyr::select(mean_no) %>% pull() - rsl_2020
    mean_no_2100 <- results %>% dplyr::filter(Year == 2100) %>% dplyr::select(mean_no) %>% pull() - rsl_2020
    mean_cop_2050 <- results %>% dplyr::filter(Year == 2050) %>% dplyr::select(mean_cop) %>% pull() - rsl_2020
    mean_cop_2100 <- results %>% dplyr::filter(Year == 2100) %>% dplyr::select(mean_cop) %>% pull() - rsl_2020
    up95_no_2050 <- results %>% dplyr::filter(Year == 2050) %>% dplyr::select(up95_no) %>% pull() - rsl_2020
    up95_no_2100 <- results %>% dplyr::filter(Year == 2100) %>% dplyr::select(up95_no) %>% pull() - rsl_2020
    low95_no_2050 <- results %>% dplyr::filter(Year == 2050) %>% dplyr::select(low95_no) %>% pull() - rsl_2020
    low95_no_2100 <- results %>% dplyr::filter(Year == 2100) %>% dplyr::select(low95_no) %>% pull() - rsl_2020
    up95_cop_2050 <- results %>% dplyr::filter(Year == 2050) %>% dplyr::select(up95_cop) %>% pull() - rsl_2020
    up95_cop_2100 <- results %>% dplyr::filter(Year == 2100) %>% dplyr::select(up95_cop) %>% pull() - rsl_2020
    low95_cop_2050 <- results %>% dplyr::filter(Year == 2050) %>% dplyr::select(low95_cop) %>% pull() - rsl_2020
    low95_cop_2100 <- results %>% dplyr::filter(Year == 2100) %>% dplyr::select(low95_cop) %>% pull() - rsl_2020


    notes <- "no VLM"
    if (add_VLM) {
      notes <- "VLM corrected"
    }
    if (add_VLM_afterwards) {
      notes <- "VLM after"
    }

    result = data.frame(station = station,
                        model = "ARDL",
                        rsl_2020 = rsl_2020,
                        source = measurement.type,
                        notes = notes,
                        mean_no_2050 = mean_no_2050 - vlm_no_2050,
                        mean_no_2100 = mean_no_2100 - vlm_no_2100,
                        mean_cop_2050 = mean_cop_2050 - vlm_cop_2050,
                        mean_cop_2100 = mean_cop_2100 - vlm_cop_2100,
                        up95_no_2050 = up95_no_2050 - vlm_no_2050,
                        up95_no_2100 = up95_no_2100 - vlm_no_2100,
                        low95_no_2050 = low95_no_2050 - vlm_no_2050,
                        low95_no_2100 = low95_no_2100 - vlm_no_2100,
                        up95_cop_2050 = up95_cop_2050 - vlm_cop_2050,
                        up95_cop_2100 = up95_cop_2100 - vlm_cop_2100,
                        low95_cop_2050 = low95_cop_2050 - vlm_cop_2050,
                        low95_cop_2100 = low95_cop_2100 - vlm_cop_2100)
    return(result)
}

get_ARDL_prediction <- function(station, start.year=1993, stop.year=2019){
    results_tg <- .get_ARDL_prediction(station, add_VLM = F, measurement.type = "tg", start.year=start.year, stop.year=stop.year)
    results_sa <- .get_ARDL_prediction(station, add_VLM = F, measurement.type = "sa", start.year=start.year, stop.year=stop.year)
    results_sa_vlm <- .get_ARDL_prediction(station, add_VLM = T, measurement.type = "sa", VLM_na_interp = T, VLM_retreive_method = "circlemean", choose_which = 5, max_distance = 300, VLM_method = c('ULR', 'NGL', 'JPL', 'GFZ'), start.year=start.year, stop.year=stop.year)
    results_sa_vlm_afterwards <- .get_ARDL_prediction(station, add_VLM = F, add_VLM_afterwards = T, start.year=start.year, stop.year=stop.year)
    results <- rbind(results_tg, results_sa, results_sa_vlm, results_sa_vlm_afterwards)
    return(results)
}

.get_SEM_prediction <- function(station, rsl_2020=0, measurement.type = "sa", add_VLM=TRUE, start.year=1993, stop.year=2019, ...){
    results_no <- SEM_prediction(station, restriction = "no", measurement = measurement.type, add_VLM = add_VLM, start.year=start.year, stop.year=stop.year, ...)
    results_cop <- SEM_prediction(station, restriction = "cop", measurement = measurement.type, add_VLM = add_VLM, start.year=start.year, stop.year=stop.year, ...)
    results_no <- data.frame(results_no)
    results_cop <- data.frame(results_cop)
    rsl_2020 <- get_baseline(station)
    mean_no_2050 <- results_no %>% dplyr::filter(Year == 2050) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
    mean_cop_2050 <- results_cop %>% dplyr::filter(Year == 2050) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
    mean_no_2100 <- results_no %>% dplyr::filter(Year == 2100) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
    mean_cop_2100 <- results_cop %>% dplyr::filter(Year == 2100) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020

    result = data.frame(station = station,
                        model = "SEM",
                        rsl_2020 = rsl_2020,
                        source = measurement.type,
                        notes = if (add_VLM) "VLM corrected" else NA,
                        mean_no_2050 = mean_no_2050,
                        mean_no_2100 = mean_no_2100,
                        mean_cop_2050 = mean_cop_2050,
                        mean_cop_2100 = mean_cop_2100,
                        up95_no_2050 = NA,
                        up95_no_2100 = NA,
                        low95_no_2050 = NA,
                        low95_no_2100 = NA,
                        up95_cop_2050 = NA,
                        up95_cop_2100 = NA,
                        low95_cop_2050 = NA,
                        low95_cop_2100 = NA)
    return(result)

}

get_SEM_prediction <- function(station, rsl_2020, start.year=1993, stop.year=2019){
    results_tg <- .get_SEM_prediction(station, rsl_2020[1], add_VLM = F, measurement.type = "tg", start.year=start.year, stop.year=stop.year)
    results_sa <- .get_SEM_prediction(station, rsl_2020[2], add_VLM = F, measurement.type = "sa", start.year=start.year, stop.year=stop.year)
    results_sa_vlm <- .get_SEM_prediction(station, rsl_2020[3], add_VLM = T, measurement.type = "sa", VLM_na_interp = T, VLM_retreive_method = "circlemean", choose_which = 5, max_distance = 300, VLM_method = c('ULR', 'NGL', 'JPL', 'GFZ'), start.year=start.year, stop.year=stop.year)
    results <- rbind(results_tg, results_sa, results_sa_vlm)
    return(results)
}

.get_VARX_prediction <- function(station, rsl_2020=0, measurement.type = "sa", add_VLM=TRUE, ...){
  results_no <- VARX_prediction(station, restriction = "no", measurement = measurement.type, add_VLM = add_VLM, ...)
  results_cop <- VARX_prediction(station, restriction = "cop", measurement = measurement.type, add_VLM = add_VLM, ...)
  results_no <- data.frame(results_no)
  results_cop <- data.frame(results_cop)
  rsl_2020 <- get_baseline(station)
  mean_no_2050 <- results_no %>% dplyr::filter(Year == 2050) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
  mean_cop_2050 <- results_cop %>% dplyr::filter(Year == 2050) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
  mean_no_2100 <- results_no %>% dplyr::filter(Year == 2100) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
  mean_cop_2100 <- results_cop %>% dplyr::filter(Year == 2100) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020

  result = data.frame(station = station,
                      model = "VARX",
                      rsl_2020 = rsl_2020,
                      source = measurement.type,
                      notes = if (add_VLM) "VLM corrected" else NA,
                      mean_no_2050 = mean_no_2050,
                      mean_no_2100 = mean_no_2100,
                      mean_cop_2050 = mean_cop_2050,
                      mean_cop_2100 = mean_cop_2100,
                      up95_no_2050 = NA,
                      up95_no_2100 = NA,
                      low95_no_2050 = NA,
                      low95_no_2100 = NA,
                      up95_cop_2050 = NA,
                      up95_cop_2100 = NA,
                      low95_cop_2050 = NA,
                      low95_cop_2100 = NA)
  return(result)

}

get_VARX_prediction <- function(station, rsl_2020=0, start.year=1950, stop.year=2019,...){
  results_tg <- .get_VARX_prediction(station, rsl_2020[1], add_VLM = F, measurement.type = "tg", start.year=start.year, stop.year=stop.year, ...)
  results <- rbind(results_tg)
  return(results)
}

.get_XVARMA_prediction <- function(station, rsl_2020, measurement.type = "sa", add_VLM=TRUE, ...){
  results_no <- XVARMA_prediction(station, restriction = "no", measurement = measurement.type, add_VLM = add_VLM, ...)
  results_cop <- XVARMA_prediction(station, restriction = "cop", measurement = measurement.type, add_VLM = add_VLM, ...)
  results_no <- data.frame(results_no)
  results_cop <- data.frame(results_cop)
  mean_no_2050 <- results_no %>% dplyr::filter(Year == 2050) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
  mean_cop_2050 <- results_cop %>% dplyr::filter(Year == 2050) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
  mean_no_2100 <- results_no %>% dplyr::filter(Year == 2100) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020
  mean_cop_2100 <- results_cop %>% dplyr::filter(Year == 2100) %>% dplyr::select(LSL_mean) %>% pull() - rsl_2020

  result = data.frame(station = station,
                      model = "XVARMA",
                      rsl_2020 = rsl_2020,
                      source = measurement.type,
                      notes = if (add_VLM) "VLM corrected" else NA,
                      mean_no_2050 = mean_no_2050,
                      mean_no_2100 = mean_no_2100,
                      mean_cop_2050 = mean_cop_2050,
                      mean_cop_2100 = mean_cop_2100,
                      up95_no_2050 = NA,
                      up95_no_2100 = NA,
                      low95_no_2050 = NA,
                      low95_no_2100 = NA,
                      up95_cop_2050 = NA,
                      up95_cop_2100 = NA,
                      low95_cop_2050 = NA,
                      low95_cop_2100 = NA)
  return(result)

}

get_XVARMA_prediction <- function(station, rsl_2020){
  results_tg <- .get_XVARMA_prediction(station, rsl_2020[1], add_VLM = F, measurement.type = "tg")
  results_sa <- .get_XVARMA_prediction(station, rsl_2020[2], add_VLM = F, measurement.type = "sa")
  results_sa_vlm <- .get_XVARMA_prediction(station, rsl_2020[3], add_VLM = T, measurement.type = "sa", VLM_na_interp = T, VLM_retreive_method = "circlemean", choose_which = 5, max_distance = 300, VLM_method = c('ULR', 'NGL', 'JPL', 'GFZ'))
  results <- rbind(results_tg, results_sa, results_sa_vlm)
  return(results)
}


get_AR6_prediction_github <- function(station){
    #' Get AR6 projections from https://github.com/Rutgers-ESSP/IPCC-AR6-Sea-Level-Projections
    #' @param station station name, which will be mapped using station_loc_id_map
    #' @return data frame
    #' @examples
    #' get_AR6_prediction_github('battery_park')
    #'

    loc_id <- find_closest_ar6_station(station)

    df_no <- read_csv(paste0("data/ar6/sea_level_change_loc", loc_id, "_ssp245.csv"))
    df_cop <- read_csv(paste0("data/ar6/sea_level_change_loc", loc_id, "_ssp126.csv"))

    results_no <- df_no %>% dplyr::filter(years %in% c(2050, 2100)) %>%
        dplyr::group_by(workflow_id) %>%
        dplyr::summarize(
            mean_no_2050 = dplyr::first(sea_level_change[years == 2050]),
            mean_no_2100 = dplyr::first(sea_level_change[years == 2100])
        )

    results_cop <- df_cop %>% dplyr::filter(years %in% c(2050, 2100)) %>%
        dplyr::group_by(workflow_id) %>%
        dplyr::summarize(
            mean_cop_2050 = dplyr::first(sea_level_change[years == 2050]),
            mean_cop_2100 = dplyr::first(sea_level_change[years == 2100])
        )

    results <- dplyr::full_join(results_no, results_cop, by = "workflow_id") %>%
        dplyr::mutate(
            station = station,
            model = "AR6",
            rsl_2020 = 0,
            source = "github",
            add_VLM = F,
            notes = workflow_id
        ) %>%
        dplyr::select(station, model, source, notes, mean_no_2050, mean_no_2100, mean_cop_2050, mean_cop_2100)

    return(results)
}

get_AR6_prediction_ipcc_tool_full <- function(station, zero_2020=FALSE){
    #' @example
    #' get_AR6_prediction_ipcc_tool_full("battery_park", zero_2020=TRUE)
    #'
    #'
    if (is.numeric(station)){
      loc_id <- as.integer(station)
    } else {
      if (station != "global") {
        loc_id <- find_closest_ar6_station(station)
      } else {
        loc_id <- station
      }
    }
    if (station == "global") {
      df <- read_csv(paste0("data/ar6/ipcc_ar6_sea_level_projection_psmsl_id_", loc_id, ".csv"))
    } else {
      df <- read_csv(paste0("data/ar6/ipcc_ar6_sea_level_projection_psmsl_id_", loc_id, ".csv"))
    }

    if (FALSE) {
        rsl_2020 <- get_baseline(station)
        df <- df %>% dplyr::mutate(across(starts_with("ssp"), ~ . - rsl_2020))
    }
    return(df)
}

get_AR6_prediction_ipcc_tool <- function(station){
    #' Get AR6 projections from IPCC AR6 projection tools (https://sealevel.nasa.gov/ipcc-ar6-sea-level-projection-tool)
    #' @param station station name, which will be mapped using station_loc_id_map
    #' @return data frame
    #' @examples
    #' get_AR6_prediction_ipcc_tool('battery_park')
    #'

    loc_id <- find_closest_ar6_station(station)

    df <- read_csv(paste0("data/ar6/ipcc_ar6_sea_level_projection_psmsl_id_", loc_id, ".csv"))

    rsl_2020 <- get_baseline(station)

get_results <- function(df, confidence_level, include_no) {
        df %>% dplyr::filter(Year %in% c(2050, 2100)) %>%
            dplyr::summarize(
                mean_no_2050 = if (include_no) dplyr::first(df[[paste0("ssp245_", confidence_level, "_50")]][df$Year == 2050]) - rsl_2020 else NA,
                mean_no_2100 = if (include_no) dplyr::first(df[[paste0("ssp245_", confidence_level, "_50")]][df$Year == 2100]) - rsl_2020 else NA,
                mean_cop_2050 = dplyr::first(df[[paste0("ssp126_", confidence_level, "_50")]][df$Year == 2050]) - rsl_2020,
                mean_cop_2100 = dplyr::first(df[[paste0("ssp126_", confidence_level, "_50")]][df$Year == 2100]) - rsl_2020
            ) %>%
            dplyr::mutate(
                station = station,
                model = "AR6",
                rsl_2020 = rsl_2020,
                source = "ipcc_tool",
                add_VLM = F,
                notes = confidence_level
            ) %>%
            dplyr::select(station, model, source, notes, mean_no_2050, mean_no_2100, mean_cop_2050, mean_cop_2100)
    }

    results_medium <- get_results(df, "medium", TRUE)
    results_low <- get_results(df, "low", FALSE)

    results <- bind_rows(results_medium, results_low)

    return(results)
}

get_AR6_prediction <- function(station){
    #' Merge AR6 projections from both GitHub and IPCC tool
    #' @param station station name, which will be mapped using station_loc_id_map
    #' @return data frame
    #' @examples
    #' get_AR6_prediction('battery_park')
    #'

    results_ipcc_tool <- get_AR6_prediction_ipcc_tool(station)

    results <- results_ipcc_tool
    return(results)
}

get_predictions_comparison <- function(station, from_file = FALSE, start.year=1993, stop.year=2019, save=TRUE,...) {
  #' Get the predictions comparison for a station
  #' @param station station name
  #' @param from_file if TRUE, read the results from a file. If the file does not exist, the function will call itself with from_file = FALSE
  #' @return data frame
  #' @examples
  #' get_predictions_comparison("battery_park")
  file_path <- paste0("results/", station, "/prediction_results.csv")
  if (from_file) {

    if (file.exists(file_path)) {
      results <- read_csv(file_path)
      return(results)
    }
  } else {
    results_AR6 <- get_AR6_prediction(station)
    results_ARDL <- get_ARDL_prediction(station, start.year=start.year, stop.year=stop.year)
    results_SEM <- get_SEM_prediction(station, results_ARDL$rsl_2020, start.year=start.year, stop.year=stop.year)
    results_VARX <- get_VARX_prediction(station, results_ARDL$rsl_2020, start.year=start.year, stop.year=stop.year,...)
    results_ensemble <- data.frame(station = station, model="ensemble", source="tg", notes="",
                                  mean_no_2050=(results_ARDL$mean_no_2050[1] + results_VARX$mean_no_2050[1]) / 2,
                                  mean_no_2100=(results_ARDL$mean_no_2100[1] + results_VARX$mean_no_2100[1]) / 2,
                                  mean_cop_2050=(results_ARDL$mean_cop_2050[1] + results_VARX$mean_cop_2050[1]) / 2,
                                  mean_cop_2100=(results_ARDL$mean_cop_2100[1] + results_VARX$mean_cop_2100[1]) / 2)
    results <- bind_rows(results_ARDL, results_SEM, results_VARX, results_ensemble, results_AR6) %>% dplyr::select(-rsl_2020)
    if (save) {
      write_csv(results, file_path, na = "")
    }
    return(results)
  }
}

prediction_statistics <- function(station, notes="VLM corrected") {
  predictions <- get_predictions_comparison(station, from_file=TRUE)
  station_name <- predictions$station[1]
  vlm_notes <- notes

  baseline <- predictions %>% dplyr::filter(station == station_name & model == "AR6" & source == "ipcc_tool" & notes == "medium")

  ardl_sa <- predictions %>% dplyr::filter(station == station_name & model == "ARDL" & source == "sa" & notes == "no VLM")
  ardl_sa_vlm <- predictions %>% dplyr::filter(station == station_name & model == "ARDL" & source == "sa" & notes == vlm_notes)

  results <- list(
    mean_no_2050_improved = 0,
    mean_no_2100_improved = 0,
    mean_cop_2050_improved = 0,
    mean_cop_2100_improved = 0,
    ipcc_in_ardl_sa_ci_no_2050 = 0,
    ipcc_in_ardl_sa_vlm_ci_no_2050 = 0,
    ipcc_in_ardl_sa_ci_no_2100 = 0,
    ipcc_in_ardl_sa_vlm_ci_no_2100 = 0,
    ipcc_in_ardl_sa_ci_cop_2050 = 0,
    ipcc_in_ardl_sa_vlm_ci_cop_2050 = 0,
    ipcc_in_ardl_sa_ci_cop_2100 = 0,
    ipcc_in_ardl_sa_vlm_ci_cop_2100 = 0
  )

  if (abs(ardl_sa$mean_no_2050 - baseline$mean_no_2050) > abs(ardl_sa_vlm$mean_no_2050 - baseline$mean_no_2050)) {
    results$mean_no_2050_improved <- 1
  }

  if (abs(ardl_sa$mean_no_2100 - baseline$mean_no_2100) > abs(ardl_sa_vlm$mean_no_2100 - baseline$mean_no_2100)) {
    results$mean_no_2100_improved <- 1
  }

  if (abs(ardl_sa$mean_cop_2050 - baseline$mean_cop_2050) > abs(ardl_sa_vlm$mean_cop_2050 - baseline$mean_cop_2050)) {
    results$mean_cop_2050_improved <- 1
  }

  if (abs(ardl_sa$mean_cop_2100 - baseline$mean_cop_2100) > abs(ardl_sa_vlm$mean_cop_2100 - baseline$mean_cop_2100)) {
    results$mean_cop_2100_improved <- 1
  }

  if (baseline$mean_no_2050 <= ardl_sa$up95_no_2050 & baseline$mean_no_2050 >= ardl_sa$low95_no_2050) {
    results$ipcc_in_ardl_sa_ci_no_2050 <- 1
  }

  if (baseline$mean_no_2050 <= ardl_sa_vlm$up95_no_2050 & baseline$mean_no_2050 >= ardl_sa_vlm$low95_no_2050) {
    results$ipcc_in_ardl_sa_vlm_ci_no_2050 <- 1
  }

  if (baseline$mean_no_2100 <= ardl_sa$up95_no_2100 & baseline$mean_no_2100 >= ardl_sa$low95_no_2100) {
    results$ipcc_in_ardl_sa_ci_no_2100 <- 1
  }

  if (baseline$mean_no_2100 <= ardl_sa_vlm$up95_no_2100 & baseline$mean_no_2100 >= ardl_sa_vlm$low95_no_2100) {
    results$ipcc_in_ardl_sa_vlm_ci_no_2100 <- 1
  }

  if (baseline$mean_cop_2050 <= ardl_sa$up95_cop_2050 & baseline$mean_cop_2050 >= ardl_sa$low95_cop_2050) {
    results$ipcc_in_ardl_sa_ci_cop_2050 <- 1
  }

  if (baseline$mean_cop_2050 <= ardl_sa_vlm$up95_cop_2050 & baseline$mean_cop_2050 >= ardl_sa_vlm$low95_cop_2050) {
    results$ipcc_in_ardl_sa_vlm_ci_cop_2050 <- 1
  }

  if (baseline$mean_cop_2100 <= ardl_sa$up95_cop_2100 & baseline$mean_cop_2100 >= ardl_sa$low95_cop_2100) {
    results$ipcc_in_ardl_sa_ci_cop_2100 <- 1
  }

  if (baseline$mean_cop_2100 <= ardl_sa_vlm$up95_cop_2100 & baseline$mean_cop_2100 >= ardl_sa_vlm$low95_cop_2100) {
    results$ipcc_in_ardl_sa_vlm_ci_cop_2100 <- 1
  }

  return(results)

}



prediction_statistics_stations <- function(stations=station.names, ...) {
    #' Calculate the prediction statistics for a list of stations
    #' @param stations list of station names
    #' @return list of prediction statistics
    #' @examples
    #' prediction_statistics_stations()

    first_station_results <- prediction_statistics(stations[1], ...)
    summed_results <- lapply(first_station_results, function(x) 0)

    for (station in stations) {
        results <- prediction_statistics(station, ...)

        if (length(results) != length(summed_results)) {
            stop("Results structure mismatch between stations")
        }

        summed_results <- mapply(function(sum, res) sum + res, summed_results, results, SIMPLIFY = FALSE)
    }
    summed_results$total_stations <- length(stations)
    return(summed_results)
}



get_GWP_our <- function(save_path="data/Greenhouse/our_GWP.csv", from_file=FALSE) {
    #' Get the global GWP.
    #' auto.arima to predict the GWP until 2100.
    #' @return a data frame with columns: Year, etc
    #'
    #'

    if (from_file) {
        return(read.csv(save_path))
    }

    climate <- read.csv('data/YearData.csv')
    climate.l <- cbind(climate[-1,], climate[-72,])
    names(climate.l)[12:22] <- paste(names(climate.l)[12:22],'.l', sep = '')

    fit_GWP <- Arima(climate$GWP, order = c(1,1,1), include.drift = TRUE, method = 'ML')
    pred_GWP_no <- forecast::forecast(fit_GWP, h = 79)
    pred_GWP_no$x <- ts(pred_GWP_no$x, start = 1950, end = 2021)
    pred_GWP_no$mean <- ts(pred_GWP_no$mean, start = 2022, end = 2100)
    pred_GWP_no$lower <- ts(pred_GWP_no$lower, start = 2022, end = 2100)
    pred_GWP_no$upper <- ts(pred_GWP_no$upper, start = 2022, end = 2100)

    fit_N2O <- Arima(climate$N2O, order = c(1,1,1), include.drift = TRUE, method = 'ML')
    pred_N2O_no <- forecast(fit_N2O, h = 79)$mean

    CO2_2010 <- climate$CO2[61] - climate$CO2[60]
    pred_CO2_cop <- c(climate$CO2[72] + cumsum(CO2_2010 - (1:9)/9*(CO2_2010*0.45)))
    pred_CO2_cop <- c(pred_CO2_cop, rep(pred_CO2_cop[9], 70))

    CH4_2020 <- climate$CH4[71] - climate$CH4[70]
    pred_CH4_cop <- c(climate$CH4[72] + cumsum(CH4_2020 - (1:9)/9*(CH4_2020*0.3)))
    pred_CH4_cop <- c(pred_CH4_cop, rep(pred_CH4_cop[9], 70))

    pred_GWP_cop_mean <- pred_CO2_cop + 28/1000 * pred_CH4_cop + 265/1000 * pred_N2O_no

    pred_GWP_cop <- pred_GWP_no
    pred_GWP_cop$mean <- pred_GWP_cop$mean - pred_GWP_no$mean + ts(pred_GWP_cop_mean, start = 2022)

    pred_GWP_cop$lower[,1] <- pred_GWP_cop$mean - (pred_GWP_no$mean - pred_GWP_no$lower[,1]) * pred_GWP_cop$mean / pred_GWP_no$mean
    pred_GWP_cop$lower[,2] <- pred_GWP_cop$mean - (pred_GWP_no$mean - pred_GWP_no$lower[,2]) * pred_GWP_cop$mean / pred_GWP_no$mean
    pred_GWP_cop$upper[,1] <- pred_GWP_cop$mean - (pred_GWP_no$mean - pred_GWP_no$upper[,1]) * pred_GWP_cop$mean / pred_GWP_no$mean
    pred_GWP_cop$upper[,2] <- pred_GWP_cop$mean - (pred_GWP_no$mean - pred_GWP_no$upper[,2]) * pred_GWP_cop$mean / pred_GWP_no$mean

    GWP <- data.frame(
        Year = 1950:2100,
        GWP_no = c(pred_GWP_no$x, pred_GWP_no$mean),
        GWP_cop = c(pred_GWP_cop$x, pred_GWP_cop$mean),
        GWP_no_80_l = c(rep(NA, 71), climate$GWP[72], pred_GWP_no$lower[,1]),
        GWP_no_95_l = c(rep(NA, 71), climate$GWP[72], pred_GWP_no$lower[,2]),
        GWP_no_80_u = c(rep(NA, 71), climate$GWP[72], pred_GWP_no$upper[,1]),
        GWP_no_95_u = c(rep(NA, 71), climate$GWP[72], pred_GWP_no$upper[,2]),
        GWP_cop_80_l = c(rep(NA, 71), climate$GWP[72], pred_GWP_cop$lower[,1]),
        GWP_cop_95_l = c(rep(NA, 71), climate$GWP[72], pred_GWP_cop$lower[,2]),
        GWP_cop_80_u = c(rep(NA, 71), climate$GWP[72], pred_GWP_cop$upper[,1]),
        GWP_cop_95_u = c(rep(NA, 71), climate$GWP[72], pred_GWP_cop$upper[,2])
    )

    GWP$GWP_no_99_u <- GWP$GWP_no + (GWP$GWP_no_95_u - GWP$GWP_no) / qnorm(0.975) * qnorm(0.995)
    GWP$GWP_no_99_l <- GWP$GWP_no - (GWP$GWP_no_95_u - GWP$GWP_no) / qnorm(0.975) * qnorm(0.995)
    GWP$GWP_cop_99_u <- GWP$GWP_cop + (GWP$GWP_cop_95_u - GWP$GWP_cop) / qnorm(0.975) * qnorm(0.995)
    GWP$GWP_cop_99_l <- GWP$GWP_cop - (GWP$GWP_cop_95_u - GWP$GWP_cop) / qnorm(0.975) * qnorm(0.995)

    if (!is.na(save_path)) {
        write.csv(GWP, save_path, row.names = FALSE)
    }

    return(GWP)
}

get_GWP_SSPs <- function() {
    #' Get the GWP SSPs
    #' @return a data frame with columns: Year, etc
    #' @examples
    #' get_GWP_SSPs()
    df <- read_csv("data/Greenhouse/merged_all_ssp.csv")

    df %>% rename(Year = year)
}

get_GWP <- function(save_path="data/Greenhouse/full_GWP.csv", from_file=FALSE,
                    cut_off=2100) {
    #' Get the GWP, merged by get_GWP_our and get_GWP_SSPs
    #' @return a data frame with columns: Year, etc
    #' @examples
    #' get_GWP()

    if (from_file) {
        return(read.csv(save_path))
    }

    df1 <- get_GWP_our()
    df2 <- get_GWP_SSPs()
    df <- dplyr::full_join(df1, df2, by = "Year")
    df <- df %>% dplyr::filter(Year <= cut_off)
    if (!is.na(save_path)) {
        write.csv(df, save_path, row.names = FALSE)
    }
    return(df)
}

compare_GWP <- function(scenario) {
    #' Compare the mean squared error of GWP of the scenario with the GWP of the SSPs
    #' @param scenario scenario name, "no" or "cop"
    #' @return a data frame. The columns are the SSPs and MSE, ordered by MSE ascending
    #' @examples
    #' compare_GWP("no")

    df <- get_GWP()

    gwp_col <- ifelse(scenario == "no", "GWP_no", "GWP_cop")

    mse_list <- sapply(names(df)[grepl("ssp", names(df))], function(ssp) {
        mean((df[[gwp_col]] - df[[ssp]])^2, na.rm = TRUE)
    })

    mse_df <- data.frame(SSP = names(mse_list), MSE = mse_list)

    mse_df <- mse_df[order(mse_df$MSE), ]

    return(mse_df)
}

plot_compare_GWP <- function(scenario, top_n = 3) {
    #' Visualize the time series data of GWP under different scenarios. Highlight method:
    #' Line Thickness: Increase the line thickness for the GWP_no or GWP_cop and the closest SSP.
    #' Line Type: Use different line types (e.g., solid for GWP_no or GWP_cop and the closest SSP, dashed for other SSPs).
    #' Color: Use distinct colors for GWP_no or GWP_cop and the closest SSP.
    #' @param scenario scenario name, "no" or "cop"
    #' @param top_n number of closest SSPs to show
    #' @examples
    #' plot_compare_GWP("no", top_n = 5)

    df <- get_GWP()

    gwp_col <- ifelse(scenario == "no", "GWP_no", "GWP_cop")

    top_ssps <- compare_GWP(scenario)$SSP[1:top_n]

    df_melted <- melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "GWP")

    df_melted <- df_melted[df_melted$Scenario %in% c(gwp_col, top_ssps), ]

    line_types <- c("solid", "solid", rep("dashed", top_n - 1))
    line_sizes <- c(1.2, 1.2, rep(0.8, top_n - 1))
    colors <- c("red", "blue", rep("grey", top_n - 1))

    names(line_types) <- c(gwp_col, top_ssps[1], top_ssps[-1])
    names(line_sizes) <- c(gwp_col, top_ssps[1], top_ssps[-1])
    names(colors) <- c(gwp_col, top_ssps[1], top_ssps[-1])

    ggplot(df_melted, aes(x = Year, y = GWP, color = Scenario, linetype = Scenario, size = Scenario)) +
        geom_line() +
        scale_color_manual(values = colors) +
        scale_linetype_manual(values = line_types) +
        scale_size_manual(values = line_sizes) +
        theme_minimal() +
        labs(title = paste("Global Warming Potential (GWP) under", scenario, "Scenario."),
             x = "Year",
             y = "GWP",
             color = "Scenario") +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom",
              legend.title = element_blank()) +
        geom_text(data = df_melted[df_melted$Year == max(df_melted$Year), ],
                  aes(label = Scenario, hjust = 1.1), size = 3, show.legend = FALSE)
}

read_ar6_tg <- function(id, ...) {
  #' The function has been moved to `data_processing.R` and renamed to `get_tg_data()`
  #' @examples
  #' read_ar6_tg(1, convert_to_monthly = TRUE, from_file = FALSE)
  #'
  #'
  get_tg_data(id, ...)

}

get_ARDL_AR6_prediction <- function(station, end_date="2019-12-31") {
  #' Get the ARDL prediction and the AR6 SSP prediction
  #' @param station station name
  #' @return a data frame
  #' @examples
  #' get_ARDL_AR6_prediction("battery_park")
  #' get_ARDL_AR6_prediction(336)

  ARDL_tg_prediction <- RSL_predict(station, measurement.type = "tg", add_VLM=FALSE, end_date = end_date)$prediction
  ARDL_tg_prediction  <- ARDL_tg_prediction %>% dplyr::select(Year, mean_no, mean_cop)
  rsl_2020 <- get_baseline(station)
  ARDL_tg_prediction <- ARDL_tg_prediction %>% dplyr::mutate(across(starts_with("mean"), ~ . - rsl_2020))

  AR6_prediction <- get_AR6_prediction_ipcc_tool_full(station, zero_2020 = TRUE)
  AR6_prediction <- AR6_prediction %>% dplyr::select(Year, ends_with("_50"))
  AR6_prediction <- AR6_prediction %>% dplyr::rename_with(~ gsub("_50", "", .), ends_with("_50"))

  df <- dplyr::left_join(ARDL_tg_prediction, AR6_prediction, by = "Year")
  return(df)
}

apply_vlm_adjustment <- function(df, vlm_data) {
  df %>%
    left_join(vlm_data, by = "Year") %>%
        mutate(across(starts_with("ARDL_"), ~ . + VLM, .names = "{.col}"),
               across(starts_with("SEM_"), ~ . + VLM, .names = "{.col}"),
               across(starts_with("VARX_"), ~ . + VLM, .names = "{.col}"),
               across(starts_with("ensemble_"), ~ . + VLM, .names = "{.col}")) %>%
    dplyr::select(-VLM)
}

rename_and_adjust <- function(df, prefix, rsl_2020) {
  names(df) <- gsub("mean_", paste0(prefix, "_"), names(df))
  df <- df %>% dplyr::rename(
    !!paste0(prefix, "_no_up95") := up95_no,
    !!paste0(prefix, "_no_low95") := low95_no,
    !!paste0(prefix, "_no_up99") := up99_no,
    !!paste0(prefix, "_no_low99") := low99_no,
    !!paste0(prefix, "_cop_up95") := up95_cop,
    !!paste0(prefix, "_cop_low95") := low95_cop,
    !!paste0(prefix, "_cop_up99") := up99_cop,
    !!paste0(prefix, "_cop_low99") := low99_cop
  )
  df %>% dplyr::mutate(across(starts_with(prefix), ~ . - rsl_2020))
}

rename_and_adjust_ssp <- function(df, prefix, ssp, rsl_2020) {
  df <- df %>% dplyr::select(Year, mean_no)
  colnames(df)[colnames(df) == "mean_no"] <- paste0(prefix, "_", ssp)
  df[[paste0(prefix, "_", ssp)]] <- df[[paste0(prefix, "_", ssp)]] - rsl_2020
  return(df)
}

get_model_col <- function(model, scenario) {
  if (scenario %in% c("no", "cop")) {
    return(paste0(model, "_", scenario))
  } else {
    return(paste0(model, "_", scenario))
  }
}

get_ssp_ensemble_predictions <- function(station, start.year, stop.year, rsl_2020, df, ...) {
  #' @return a data frame with "Year", "ensemble_ssp119", "ensemble_ssp126", "ensemble_ssp245",
  #' "ensemble_ssp370", "ensemble_ssp460", "ensemble_ssp585"
  ssp_scenarios <- c("ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585")
  all_ssp_ensemble_predictions <- list()

  VARX_col <- "VARX_no"

  for (ssp in ssp_scenarios) {
    ardl_ssp_prediction <- RSL_predict_single(station, measurement.type = "tg", add_VLM=FALSE, scenario = ssp, start_year = start.year, end_year = stop.year, use_VLM_as_variable = FALSE, ...)$prediction
    ardl_ssp_increase_2100_no <- ardl_ssp_prediction %>% dplyr::filter(Year == 2100) %>% dplyr::select(mean_no) %>% pull() - rsl_2020
    ardl_ssp_prediction <- rename_and_adjust_ssp(ardl_ssp_prediction, "ARDL", ssp, rsl_2020)

    ardl_tg_vlm_prediction <- RSL_predict_single(station, measurement.type = "tg", add_VLM=FALSE, use_VLM_as_variable = TRUE, scenario = ssp, start_year = start.year, end_year = stop.year, ...)$prediction
    ardl_tg_vlm_increase_2100_no <- ardl_tg_vlm_prediction %>% dplyr::filter(Year == 2100) %>% dplyr::select(mean_no) %>% pull() - rsl_2020
    ardl_tg_vlm_prediction <- rename_and_adjust_ssp(ardl_tg_vlm_prediction, "ARDL_VLM", ssp, rsl_2020)

    ensemble_models <- select_ensemble_models(
      ardl_ssp_increase_2100_no,
      ardl_tg_vlm_increase_2100_no,
      df %>% dplyr::filter(Year == 2100) %>% dplyr::select(VARX_no) %>% pull() - rsl_2020
    )

    ardl_col <- get_model_col("ARDL", ssp)
    ardl_vlm_col <- get_model_col("ARDL_VLM", ssp)

    temp_df <- df %>% dplyr::select(Year, !!VARX_col)
    temp_df <- dplyr::left_join(temp_df, ardl_ssp_prediction, by = "Year")
    temp_df <- dplyr::left_join(temp_df, ardl_tg_vlm_prediction, by = "Year")

    if (length(ensemble_models$no) == 2) {
      col1 <- ifelse(ensemble_models$no[1] == "ARDL_no", ardl_col,
                     ifelse(ensemble_models$no[1] == "ARDL_VLM_no", ardl_vlm_col, VARX_col))
      col2 <- ifelse(ensemble_models$no[2] == "ARDL_no", ardl_col,
                     ifelse(ensemble_models$no[2] == "ARDL_VLM_no", ardl_vlm_col, VARX_col))
      temp_df <- temp_df %>%
        ensemble_dynamic_df(col1, col2, paste0("ensemble_", ssp), start.year=2020, stop.year=2100)
    } else {
      ensemble_col_name <- paste0("ensemble_", ssp)
      col1 <- ifelse(ensemble_models$no == "ARDL_no", ardl_col,
                     ifelse(ensemble_models$no == "ARDL_VLM_no", ardl_vlm_col, VARX_col))
      temp_df[, ensemble_col_name] <- temp_df[, col1]
    }

    all_ssp_ensemble_predictions[[ssp]] <- temp_df %>% dplyr::select(Year, starts_with(paste0("ensemble_", ssp)))
  }

  final_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Year"), all_ssp_ensemble_predictions)
  return(final_df)
}

get_predictions_all_method <- function(station, start.year=1970, stop.year=2019, from_file=FALSE,
                                        save_dir='results/predictions_all',
                                        save=TRUE,
                                        add_vlm=FALSE,
                                        show_real=FALSE,
                                        ...) {
  #' @return a data frame with predictions from all methods ("Year", "ARDL_no", "ARDL_cop", "ARDL_no_up95", "ARDL_no_low95",
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  save_path <- file.path(save_dir, paste0(station, ".csv"))
  if (from_file && file.exists(save_path)) {
    df <- read_csv(save_path)
    if (add_vlm) {
      vlm_data <- get_ar6_vlm(station, zero_year = stop.year)
      df <- apply_vlm_adjustment(df, vlm_data)
    }
    if (show_real) {
      real_data <- get_year_data_raw(station, adjust_baseline=TRUE, use_VLM_as_variable=FALSE) %>%
        dplyr::filter(Year > stop.year) %>%
        dplyr::select(Year, tg) %>%
        dplyr::rename(real_tg = tg)
      df <- dplyr::left_join(df, real_data, by = "Year")
    }
    return(df)
  }

  calculate_missing_rate <- function(data) {
    missing_rate <- sum(data$is_na > 0) / nrow(data)
    return(missing_rate)
  }
  .get_sea_level_data <- function(station) {
    get_year_data_raw(station) %>%
      dplyr::rename(Value = tg) %>%
      dplyr::filter(Year >= 1970 & Year < 2020)
  }
  sea_level_data <- .get_sea_level_data(station)
  missing_rate <- calculate_missing_rate(sea_level_data)

  if (missing_rate > 0.99) {
    stop("Data integrity check failed: all data are missing.")
  }

  ARDL_tg_prediction <- RSL_predict(station, measurement.type = "tg", add_VLM=FALSE, start_year = start.year, end_year = stop.year)$prediction
  rsl_2020 <- get_baseline(station)
  ARDL_tg_increase_2100_no <- ARDL_tg_prediction %>% dplyr::filter(Year == 2100) %>% dplyr::select(mean_no) %>% pull() - rsl_2020
  ARDL_tg_prediction <- ARDL_tg_prediction %>% dplyr::select(Year, mean_no, mean_cop, up95_no, low95_no, up99_no, low99_no, up95_cop, low95_cop, up99_cop, low99_cop)
  ARDL_tg_prediction <- rename_and_adjust(ARDL_tg_prediction, "ARDL", rsl_2020)

  ARDL_tg_vlm_prediction <- RSL_predict(station, measurement.type = "tg", add_VLM=FALSE, use_VLM_as_variable = T, equal_time = F, start_year = start.year, end_year = stop.year)$prediction
  ARDL_tg_vlm_increase_2100_no <- ARDL_tg_vlm_prediction %>% dplyr::filter(Year == 2100) %>% dplyr::select(mean_no) %>% pull() - rsl_2020
  ARDL_tg_vlm_prediction <- ARDL_tg_vlm_prediction %>% dplyr::select(Year, mean_no, mean_cop, up95_no, low95_no, up99_no, low99_no, up95_cop, low95_cop, up99_cop, low99_cop)
  ARDL_tg_vlm_prediction <- rename_and_adjust(ARDL_tg_vlm_prediction, "ARDL_VLM", rsl_2020)

  SEM_tg_prediction_no <- tryCatch({
    SEM_prediction(station, measurement = "tg", start.year = start.year, stop.year = 2019, restriction = "no")
  }, error = function(e) {
    data.frame(Year = integer(),
               SEM_no = numeric(),
               SEM_no_low95 = numeric(),
               SEM_no_up95 = numeric(),
               SEM_no_low99 = numeric(),
               SEM_no_up99 = numeric())
  })

  if(nrow(SEM_tg_prediction_no) > 0) {
    SEM_tg_prediction_no <- SEM_tg_prediction_no %>% dplyr::select(Year, LSL_mean, low_95, up_95, low_99, up_99) %>%
      dplyr::rename(
        SEM_no = LSL_mean,
        SEM_no_low95 = low_95,
        SEM_no_up95 = up_95,
        SEM_no_low99 = low_99,
        SEM_no_up99 = up_99
      ) %>%
      dplyr::mutate(across(starts_with("SEM_"), ~ . - rsl_2020))
  }

  SEM_tg_prediction_cop <- tryCatch({
    SEM_prediction(station, measurement = "tg", start.year = start.year, stop.year = 2019, restriction = "cop")
  }, error = function(e) {
    data.frame(Year = integer(),
               SEM_cop = numeric(),
               SEM_cop_low95 = numeric(),
               SEM_cop_up95 = numeric(),
               SEM_cop_low99 = numeric(),
               SEM_cop_up99 = numeric())
  })

  if(nrow(SEM_tg_prediction_cop) > 0) {
    SEM_tg_prediction_cop <- SEM_tg_prediction_cop %>% dplyr::select(Year, LSL_mean, low_95, up_95, low_99, up_99) %>%
      dplyr::rename(
        SEM_cop = LSL_mean,
        SEM_cop_low95 = low_95,
        SEM_cop_up95 = up_95,
        SEM_cop_low99 = low_99,
        SEM_cop_up99 = up_99
      ) %>%
      dplyr::mutate(across(starts_with("SEM_"), ~ . - rsl_2020))
  }

  SEM_tg_prediction <- dplyr::left_join(SEM_tg_prediction_no, SEM_tg_prediction_cop, by = "Year", suffix = c(".no", ".cop"))

  VARX_tg_prediction_no <- VARX_prediction(station, measurement = "tg", start.year = start.year, stop.year = stop.year, restriction = "no", ...)
  VARX_tg_prediction_no <- VARX_tg_prediction_no %>% dplyr::select(Year, LSL_mean, low_95, up_95, low_99, up_99) %>%
    dplyr::rename(
      VARX_no = LSL_mean,
      VARX_no_low95 = low_95,
      VARX_no_up95 = up_95,
      VARX_no_low99 = low_99,
      VARX_no_up99 = up_99
    ) %>%
    dplyr::mutate(across(starts_with("VARX_"), ~ . - rsl_2020))
  VARX_tg_increase_2100_no <- VARX_tg_prediction_no %>% dplyr::filter(Year == 2100) %>% dplyr::select(VARX_no) %>% pull() - rsl_2020

  VARX_tg_prediction_cop <- VARX_prediction(station, measurement = "tg", start.year = start.year, stop.year = stop.year, restriction = "cop", ...)
  VARX_tg_prediction_cop <- VARX_tg_prediction_cop %>% dplyr::select(Year, LSL_mean, low_95, up_95, low_99, up_99) %>%
    dplyr::rename(
      VARX_cop = LSL_mean,
      VARX_cop_low95 = low_95,
      VARX_cop_up95 = up_95,
      VARX_cop_low99 = low_99,
      VARX_cop_up99 = up_99
    ) %>%
    dplyr::mutate(across(starts_with("VARX_"), ~ . - rsl_2020))

  VARX_tg_prediction <- dplyr::left_join(VARX_tg_prediction_no, VARX_tg_prediction_cop, by = "Year")

  AR6_prediction <- get_AR6_prediction_ipcc_tool_full(station, zero_2020 = FALSE)
  AR6_prediction <- AR6_prediction %>% dplyr::select(Year, ends_with("_50"))
  AR6_prediction <- AR6_prediction %>% dplyr::rename_with(~ gsub("_50", "", .), ends_with("_50"))

  df <- dplyr::left_join(ARDL_tg_prediction, AR6_prediction, by = "Year")
  df <- dplyr::left_join(df, ARDL_tg_vlm_prediction, by = "Year")
  df <- dplyr::left_join(df, SEM_tg_prediction, by = "Year")
  df <- dplyr::left_join(df, VARX_tg_prediction, by = "Year")

  ensemble_models <- select_ensemble_models(
    ARDL_tg_increase_2100_no,
    ARDL_tg_vlm_increase_2100_no,
    VARX_tg_increase_2100_no
  )

  if (length(ensemble_models$no) == 2) {
    df <- df %>% ensemble_dynamic_df(ensemble_models$no[1], ensemble_models$no[2], "ensemble_no", start.year=2020, stop.year=2100)
  } else {
    df[, "ensemble_no"] <- df[, ensemble_models$no]
    for(b in c("low95", "up95", "low99", "up99")) {
      src_col <- paste0(ensemble_models$no, "_", b)
      if (src_col %in% names(df)) {
        df[, paste0("ensemble_no_", b)] <- df[, src_col]
      }
    }
  }

  if (length(ensemble_models$cop) == 2) {
    df <- df %>% ensemble_dynamic_df(ensemble_models$cop[1], ensemble_models$cop[2], "ensemble_cop", start.year=2020, stop.year=2100)
  } else {
    df[, "ensemble_cop"] <- df[, ensemble_models$cop]
    for(b in c("low95", "up95", "low99", "up99")) {
      src_col <- paste0(ensemble_models$cop, "_", b)
      if (src_col %in% names(df)) {
        df[, paste0("ensemble_cop_", b)] <- df[, src_col]
      }
    }
  }

  ssp_ensemble_preds <- get_ssp_ensemble_predictions(station, start.year, stop.year, rsl_2020, df, rename_and_adjust, ...)
  df <- dplyr::left_join(df, ssp_ensemble_preds, by = "Year")

  if (save) {
    write_csv(df, save_path)
  }
  if (show_real) {
    real_data <- get_year_data_raw(station, adjust_baseline=TRUE, use_VLM_as_variable=FALSE) %>%
      dplyr::filter(Year > stop.year) %>%
      dplyr::select(Year, tg) %>%
      dplyr::rename(real_tg = tg)
    df <- dplyr::left_join(df, real_data, by = "Year")
  }
  return(df)
}

generate_ensemble_summary_table <- function(stations=c("astoria", "battery_park",
        "cape_charles", "charleston", "crescent_city",
        "monterey", "newport", "south_beach", "tofino"), start.year=-Inf, stop.year=2100, from_file=TRUE, save_dir='results/predictions_all', save=TRUE, ...) {
  #' Generate a summary table with lower bound, mean, and upper bound values for the ensemble model in 2050 and 2100
  #' @param stations a vector of station names
  #' @return a data frame with columns: Location, Lower bound in 2050, Mean in 2050, Upper bound in 2050, Lower bound in 2100, Mean in 2100, Upper bound in 2100
  #' @examples
  #' generate_ensemble_summary_table(c("battery_park", "astoria"))

  summary_list <- lapply(stations, function(station) {
    df <- get_predictions_all_method(station, start.year=start.year, stop.year=stop.year, from_file=from_file, save=TRUE, ...)

    summary_2050 <- df %>% dplyr::filter(Year == 2050) %>% dplyr::select(
      `Lower bound in 2050` = ensemble_no_low95,
      `Mean in 2050` = ensemble_no,
      `Upper bound in 2050` = ensemble_no_up95
    )

    summary_2100 <- df %>% dplyr::filter(Year == 2100) %>% dplyr::select(
      `Lower bound in 2100` = ensemble_no_low95,
      `Mean in 2100` = ensemble_no,
      `Upper bound in 2100` = ensemble_no_up95
    )

    summary_table <- dplyr::bind_cols(summary_2050, summary_2100)
    summary_table <- summary_table %>% dplyr::mutate(Location = station)

    summary_table <- summary_table %>% dplyr::select(
      Location,
      `Lower bound in 2050`,
      `Mean in 2050`,
      `Upper bound in 2050`,
      `Lower bound in 2100`,
      `Mean in 2100`,
      `Upper bound in 2100`
    )

    return(summary_table)
  })

  final_summary_table <- dplyr::bind_rows(summary_list)

  if (save) {
    write_csv(final_summary_table, file.path(save_dir, "ensemble_summary_table.csv"))
  }

  return(final_summary_table)
}

compare_ARDL_prediction_with_SSP <- function(station, scenario="no", method="MSE") {
  #' Compare the ARDL prediction with the SSPs
  #' @param station station name
  #' @param scenario scenario name, "no" or "cop"
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a data frame
  #' @examples
  #' compare_ARDL_prediction_with_SSP("battery_park", method="correlation")

  ARDL_col <- ifelse(scenario == "no", "mean_no", "mean_cop")
  df <- get_ARDL_AR6_prediction(station)

  similarity_metric <- function(ARDL_series, SSP_series, method) {
    if (method == "MSE") {
      return(mean((ARDL_series - SSP_series)^2, na.rm = TRUE))
    } else if (method == "correlation") {
      return(-cor(ARDL_series, SSP_series, use = "complete.obs"))
    } else if (method == "cosine") {
      return(-sum(ARDL_series * SSP_series) / (sqrt(sum(ARDL_series^2)) * sqrt(sum(SSP_series^2))))
    } else if (method == "KL_divergence") {
      ARDL_series <- ARDL_series / sum(ARDL_series, na.rm = TRUE)
      SSP_series <- SSP_series / sum(SSP_series, na.rm = TRUE)
      return(sum(ARDL_series * log(ARDL_series / SSP_series), na.rm = TRUE))
    } else if (method == "diff_2050") {
      return(abs(ARDL_series[df$Year == 2050] - SSP_series[df$Year == 2050]))
    } else if (method == "diff_2100") {
      return(abs(ARDL_series[df$Year == 2100] - SSP_series[df$Year == 2100]))
    } else {
      stop("Invalid method")
    }
  }

  similarity_list <- sapply(names(df)[grepl("ssp", names(df))], function(ssp) {
    similarity_metric(df[[ARDL_col]], df[[ssp]], method)
  })

  similarity_df <- data.frame(SSP = names(similarity_list), Similarity = similarity_list)

  similarity_df <- similarity_df[order(similarity_df$Similarity), ]

  return(similarity_df)
}

plot_compare_ARDL_SSP <- function(station, scenario="no", top_n = 3, method="MSE"){
  #' Plot the ARDL prediction and the AR6 SSP prediction
  #' @param station station name
  #' @param scenario scenario name, "no" or "cop"
  #' @param top_n number of closest SSPs to show
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a ggplot object
  #' @examples
  #' plot_compare_ARDL_SSP("battery_park", "no", top_n = 3, method="correlation")

  df <- get_ARDL_AR6_prediction(station)
  df <- df %>% dplyr::filter(Year >= 2020)
  ARDL_col <- ifelse(scenario == "no", "mean_no", "mean_cop")

  top_ssps <- compare_ARDL_prediction_with_SSP(station, scenario, method)$SSP[1:top_n]

  highlight_ssp <- ifelse(scenario == "no", "ssp245_medium", "ssp126_medium")
  if (!(highlight_ssp %in% top_ssps)) {
    top_ssps <- c(top_ssps, highlight_ssp)
  }

  df_melted <- melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")

  df_melted <- df_melted[df_melted$Scenario %in% c(ARDL_col, top_ssps), ]

  line_types <- rep("dashed", length(top_ssps) + 1)
  line_sizes <- rep(0.8, length(top_ssps) + 1)
  colors <- rep("grey", length(top_ssps) + 1)

  line_types[1] <- "solid"
  line_sizes[1] <- 1.2
  colors[1] <- "red"

  closest_ssp <- top_ssps[1]
  closest_index <- which(top_ssps == closest_ssp) + 1
  line_types[closest_index] <- "solid"
  line_sizes[closest_index] <- 1.2
  colors[closest_index] <- "blue"

  if (highlight_ssp %in% top_ssps) {
    highlight_index <- which(top_ssps == highlight_ssp) + 1
    line_types[highlight_index] <- "dashed"
    line_sizes[highlight_index] <- 1.2
    colors[highlight_index] <- "green"
  }

  names(line_types) <- c(ARDL_col, top_ssps)
  names(line_sizes) <- c(ARDL_col, top_ssps)
  names(colors) <- c(ARDL_col, top_ssps)

  ggplot(df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    geom_line() +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = line_types) +
    scale_size_manual(values = line_sizes) +
    theme_minimal() +
    labs(title = paste("RSL under", scenario, "Scenario for", station, "using", method, "method."),
         x = "Year",
         y = "RSL",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_text(data = df_melted[df_melted$Year == max(df_melted$Year), ],
              aes(label = Scenario, hjust = 1.1), size = 3, show.legend = FALSE)
}

plot_compare_VARX_SSP <- function(station, scenario="no", top_n = 3, method="MSE"){
  #' Plot the ARDL prediction and the AR6 SSP prediction
  #' @param station station name
  #' @param scenario scenario name, "no" or "cop"
  #' @param top_n number of closest SSPs to show
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a ggplot object
  #' @examples
  #' plot_compare_ARDL_SSP("battery_park", "no", top_n = 3, method="correlation")

  pred <- VARX_prediction(station,stop.year = 2020)

  df <- get_ARDL_AR6_prediction(station)
  df <- df %>% dplyr::filter(Year >= 2020)
  df$mean_no = pred$LSL_mean - pred$LSL_mean[1]
  ARDL_col <- "mean_no"

  top_ssps <- compare_ARDL_prediction_with_SSP(station, scenario, method)$SSP[1:top_n]

  highlight_ssp <- ifelse(scenario == "no", "ssp245_medium", "ssp126_medium")
  if (!(highlight_ssp %in% top_ssps)) {
    top_ssps <- c(top_ssps, highlight_ssp)
  }

  df_melted <- melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")

  df_melted <- df_melted[df_melted$Scenario %in% c(ARDL_col, top_ssps), ]

  line_types <- rep("dashed", length(top_ssps) + 1)
  line_sizes <- rep(0.8, length(top_ssps) + 1)
  colors <- rep("grey", length(top_ssps) + 1)

  line_types[1] <- "solid"
  line_sizes[1] <- 1.2
  colors[1] <- "red"

  closest_ssp <- top_ssps[1]
  closest_index <- which(top_ssps == closest_ssp) + 1
  line_types[closest_index] <- "solid"
  line_sizes[closest_index] <- 1.2
  colors[closest_index] <- "blue"

  if (highlight_ssp %in% top_ssps) {
    highlight_index <- which(top_ssps == highlight_ssp) + 1
    line_types[highlight_index] <- "dashed"
    line_sizes[highlight_index] <- 1.2
    colors[highlight_index] <- "green"
  }

  names(line_types) <- c(ARDL_col, top_ssps)
  names(line_sizes) <- c(ARDL_col, top_ssps)
  names(colors) <- c(ARDL_col, top_ssps)

  ggplot(df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    geom_line() +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = line_types) +
    scale_size_manual(values = line_sizes) +
    theme_minimal() +
    labs(title = paste("RSL under", scenario, "Scenario for", station, "using", method, "method."),
         x = "Year",
         y = "RSL",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_text(data = df_melted[df_melted$Year == max(df_melted$Year), ],
              aes(label = Scenario, hjust = 1.1), size = 3, show.legend = FALSE)
}


plot_compare_all_methods <- function(station, scenario="no", top_n = 3, method="MSE", from_file=TRUE, start.year=1993, stop.year=2019, ...) {
  #' Plot the predictions from ARDL, SEM, VARX, and ensemble methods along with the AR6 SSP predictions
  #' @param station station name
  #' @param scenario scenario name, "no" or "cop"
  #' @param top_n number of closest SSPs to show
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a ggplot object
  #' @examples
  #' plot_compare_all_methods("battery_park", "no", top_n = 3, method="correlation")

  df <- get_predictions_all_method(station, from_file = from_file, start.year = start.year, stop.year = stop.year, ...)
  df <- df %>% dplyr::filter(Year >= 2020)
  ARDL_col <- ifelse(scenario == "no", "ARDL_no", "ARDL_cop")
  ARDL_VLM_col <- ifelse(scenario == "no", "ARDL_VLM_no", "ARDL_VLM_cop")
  SEM_col <- ifelse(scenario == "no", "SEM_no", "SEM_cop")
  VARX_col <- ifelse(scenario == "no", "VARX_no", "VARX_cop")
  ensemble_col <- ifelse(scenario == "no", "ensemble_no", "ensemble_cop")

  top_ssps <- compare_ARDL_prediction_with_SSP(station, scenario, method)$SSP[1:top_n]

  highlight_ssp <- ifelse(scenario == "no", "ssp245_medium", "ssp126_medium")
  if (!(highlight_ssp %in% top_ssps)) {
    top_ssps <- c(top_ssps, highlight_ssp)
  }

  ci_data <- df %>%
    dplyr::select(Year,
                  !!ensemble_col,
                  !!paste0(ensemble_col, "_low95"),
                  !!paste0(ensemble_col, "_up95"),
                  !!paste0(ensemble_col, "_low99"),
                  !!paste0(ensemble_col, "_up99"))

  df_melted <- reshape2::melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")

  df_melted <- df_melted[df_melted$Scenario %in% c(ARDL_col, ARDL_VLM_col, SEM_col, VARX_col, ensemble_col, top_ssps), ]

  line_types <- rep("dashed", length(top_ssps) + 5)
  line_sizes <- rep(0.8, length(top_ssps) + 5)
  colors <- rep("grey", length(top_ssps) + 5)

  line_types[1:5] <- "solid"
  line_sizes[1:5] <- 1.2
  colors[1:5] <- c("red", "black", "blue", "green", "brown")

  closest_ssp <- top_ssps[1]
  closest_index <- which(top_ssps == closest_ssp) + 5
  line_types[closest_index] <- "dashed"
  line_sizes[closest_index] <- 1.2
  colors[closest_index] <- "purple"

  if (highlight_ssp %in% top_ssps) {
    highlight_index <- which(top_ssps == highlight_ssp) + 5
    line_types[highlight_index] <- "dashed"
    line_sizes[highlight_index] <- 1.2
    colors[highlight_index] <- "orange"
  }

  names(line_types) <- c(ARDL_col, ARDL_VLM_col, SEM_col, VARX_col, ensemble_col, top_ssps)
  names(line_sizes) <- c(ARDL_col, ARDL_VLM_col, SEM_col, VARX_col, ensemble_col, top_ssps)
  names(colors) <- c(ARDL_col, ARDL_VLM_col, SEM_col, VARX_col, ensemble_col, top_ssps)

  ggplot() +
    geom_ribbon(data = ci_data, aes(x = Year, ymin = !!sym(paste0(ensemble_col, "_low99")), ymax = !!sym(paste0(ensemble_col, "_up99"))), fill = "grey", alpha = 0.3) +
    geom_ribbon(data = ci_data, aes(x = Year, ymin = !!sym(paste0(ensemble_col, "_low95")), ymax = !!sym(paste0(ensemble_col, "_up95"))), fill = "grey", alpha = 0.5) +
    geom_line(data = df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = line_types) +
    scale_size_manual(values = line_sizes) +
    theme_minimal() +
    labs(title = paste("RSL under", scenario, "Scenario for", station, "using", method, "method."),
         x = "Year",
         y = "RSL",
         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_text(data = df_melted %>% dplyr::filter(Year == max(Year)),
              aes(x = Year, y = RSL, label = Scenario, hjust = 1.1), size = 3, show.legend = FALSE)
}

plot_compare_ensemble_and_SSPs <- function(station, scenario="no", top_n = 3, method="MSE", from_file=TRUE, start.year=1993, stop.year=2019, ...) {
  #' Plot the predictions from the ensemble method along with the AR6 SSP predictions
  #' @param station station name
  #' @param scenario scenario name, "no" or "cop"
  #' @param top_n number of closest SSPs to show
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a ggplot object
  #' @examples
  #' plot_compare_ensemble_and_SSPs("battery_park", "no", top_n = 3, method="correlation")

  df <- get_predictions_all_method(station, from_file = from_file, start.year = start.year, stop.year = stop.year, ...)
  df <- df %>% dplyr::filter(Year >= 2020)
  ensemble_col <- ifelse(scenario == "no", "ensemble_no", "ensemble_cop")

  top_ssps <- compare_ARDL_prediction_with_SSP(station, scenario, method)$SSP[1:top_n]

  highlight_ssp <- ifelse(scenario == "no", "ssp245_medium", "ssp126_medium")
  if (!(highlight_ssp %in% top_ssps)) {
    top_ssps <- c(top_ssps, highlight_ssp)
  }

  ci_data <- df %>%
    dplyr::select(Year,
                  !!ensemble_col,
                  !!paste0(ensemble_col, "_low95"),
                  !!paste0(ensemble_col, "_up95"),
                  !!paste0(ensemble_col, "_low99"),
                  !!paste0(ensemble_col, "_up99"))

  df_melted <- reshape2::melt(df, id.vars = "Year", variable.name = "Scenario", value.name = "RSL")

  df_melted <- df_melted[df_melted$Scenario %in% c(ensemble_col, top_ssps), ]

  line_types <- rep("dashed", length(top_ssps) + 1)
  line_sizes <- rep(0.8, length(top_ssps) + 1)
  colors <- c("#0048ff", "green", "orange", "yellow", "brown", "pink", "gray", "coral")
  ssp_list <- c("ssp119_medium", "ssp126_medium", "ssp126_low", "ssp245_medium", "ssp370_medium",
  "ssp585_medium", "ssp585_low")
  top_ssps <- ssp_list

  line_types[1] <- "solid"
  line_sizes[1] <- 1.2
  colors[1] <- "#0048ff"

  closest_ssp <- top_ssps[1]
  closest_index <- which(top_ssps == closest_ssp) + 1
  line_types[closest_index] <- "dashed"
  line_sizes[closest_index] <- 1.2


  names(line_types) <- c(ensemble_col, top_ssps)
  names(line_sizes) <- c(ensemble_col, top_ssps)
  names(colors) <- c(ensemble_col, top_ssps)

  ggplot() +
    geom_line(data = df_melted, aes(x = Year, y = RSL, color = Scenario, linetype = Scenario, size = Scenario)) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = line_types) +
    scale_size_manual(values = line_sizes) +
    theme_minimal() +
    labs(title = paste("RSL under", scenario, "Scenario for", station, "using", method, "method."),
         x = "Year",
         y = "RSL",

         color = "Scenario") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",

          legend.title = element_blank()) +

    geom_text(data = df_melted %>% dplyr::filter(Year == max(Year)),

              aes(x = Year, y = RSL, label = Scenario, hjust = 1.1), size = 3, show.legend = FALSE)

}


compare_all_predictions_with_SSP <- function(station, scenario="no", method="MSE", from_file=FALSE, ...) {
  #' Compare the predictions (ARDL) with the SSPs
  #' @param station station name
  #' @param scenario scenario name, "no" or "cop"
  #' @param method similarity metric, one of "MSE", "correlation", "cosine", "KL_divergence", "diff_2050", "diff_2100"
  #' @return a data frame
  #' @examples
  #' compare_all_predictions_with_SSP("battery_park", method="correlation")

  ARDL_col <- ifelse(scenario == "no", "ARDL_no", "ARDL_cop")

  df <- get_predictions_all_method(station, from_file = from_file, ...)

  similarity_metric <- function(pred_series, SSP_series, method) {
    if (method == "MSE") {
      return(mean((pred_series - SSP_series)^2, na.rm = TRUE))
    } else if (method == "correlation") {
      return(-cor(pred_series, SSP_series, use = "complete.obs"))
    } else if (method == "cosine") {
      return(-sum(pred_series * SSP_series) / (sqrt(sum(pred_series^2)) * sqrt(sum(SSP_series^2))))
    } else if (method == "KL_divergence") {
      pred_series <- pred_series / sum(pred_series, na.rm = TRUE)
      SSP_series <- SSP_series / sum(SSP_series, na.rm = TRUE)
      return(sum(pred_series * log(pred_series / SSP_series), na.rm = TRUE))
    } else if (method == "diff_2050") {
      return(abs(pred_series[df$Year == 2050] - SSP_series[df$Year == 2050]))
    } else if (method == "diff_2100") {
      return(abs(pred_series[df$Year == 2100] - SSP_series[df$Year == 2100]))
    } else {
      stop("Invalid method")
    }
  }

  similarity_list <- sapply(names(df)[grepl("ssp", names(df))], function(ssp) {
    similarity_metric(df[[ARDL_col]], df[[ssp]], method)
  })

  similarity_df <- data.frame(SSP = names(similarity_list), Similarity = similarity_list)

  similarity_df <- similarity_df[order(similarity_df$Similarity), ]

  return(similarity_df)
}

