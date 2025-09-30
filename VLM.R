library(maps)
library(ggplot2)
library(ggrepel)
library(lubridate)
library(zoo)
library(tidyverse)
library(docstring)
source("utils.R")
source("data_processing.R")
library(geosphere)
library(leaflet)
library(readxl)
library(dplyr)

vlm_rename_dict <- function(type = "method") {
  #' Get the rename dictionary for VLM data
  #' @param type character, "method" or "rename". The names of the methods are "1", "2", "3", "4", which correspond to "ULR", "NGL", "JPL", "GFZ".
  #' @return a dictionary
  #' @export
  #' @examples
  #' vlm_rename_dict("method")
  method_dict <- c("1" = "ULR", "2" = "NGL", "3" = "JPL", "4" = "GFZ")
  rename_dict <- c("ULR" = "1", "NGL" = "2", "JPL" = "3", "GFZ" = "4")
  if (type == "method") {
    return(method_dict)
  }
  else if (type == "rename") {
    return(rename_dict)
  }
  else {
    stop("Invalid type. Please use 'method' or 'rename'.")
  }
}

get_VLM_old <- function(lan_input,lon_input,verbose=FALSE){
  location<-read.table('data/VLM/JPL_location.txt',head=TRUE)
  location<-location[location$X2=="POS",]
  distance=(location$N-lan_input)^2+(location$E-lon_input)^2
  VLM_station_name=location[which.min(distance),]$Name
  if(verbose){
    print('The info of the nearest VLM station: ')
    print(location[which.min(distance),])
  }

  VLM_ts=read.table(paste('https://sideshow.jpl.nasa.gov/pub/JPL_GPS_Timeseries/repro2018a/post/point/',VLM_station_name,'.series',sep='')
                    ,head=FALSE)

  VLM_ts=data.frame(as.Date(date_decimal(VLM_ts$V1)),VLM_ts$V4)
  names(VLM_ts)=c('Date','VLM')
  VLM_ts$month <- format(VLM_ts$Date, '%Y-%m')
  VLM_ts_monthly<- data.frame(VLM = tapply(VLM_ts$VLM,VLM_ts$month, mean))


  all_month=paste(rownames(VLM_ts_monthly),'-01',sep='')

  VLM_ts_monthly=data.frame(month=all_month,
                            VLM=VLM_ts_monthly$VLM)
  VLM_ts_monthly <- VLM_ts_monthly %>% mutate(month = as.Date(month))

  return(VLM_ts_monthly)

}

read_neu <- function(neu_path, convert_decimal_year = T) {
  neu <- read.table(neu_path, header = FALSE, comment.char = "#")
  names(neu) <- c("Year", "DN", "DE", "DU", "SDN", "SDE", "SDU")
  if (convert_decimal_year) {
    neu$Year <- date_decimal(neu$Year)
  }
  return(neu)
}

read_station_table <- function(table_path) {
  station_table <- read.table(table_path, header = FALSE, skip = 13)
  names(station_table) <- c("Site", "DOMES", "Lon", "Lat", "T_GPS", "Data", "V_GPS", "S_GPS", "MODEL")
  return(station_table)
}

read_ULR_table <- function(ULR_table_path = 'data/VLM/ulr7_vertical_velocities.txt') {
  ULR_table <- read_station_table(ULR_table_path)
  return(ULR_table)
}

read_NGL_table <- function(NGL_table_path = 'data/VLM/vertical_velocities_table_ngl14.txt') {
  NGL_table <- read_station_table(NGL_table_path)
  return(NGL_table)
}

read_JPL_table <- function(JPL_table_path = 'data/VLM/vertical_velocities_table_jpl14.txt') {
  JPL_table <- read_station_table(JPL_table_path)
  return(JPL_table)
}

read_GFZ_table <- function(GFZ_table_path = 'data/VLM/vertical_velocities_table_gt3.txt') {
  GFZ_table <- read_station_table(GFZ_table_path)
  return(GFZ_table)
}

.get_VLM <- function(lat, lon, method = 'ULR', choose_which = 1,
                     max_distance = 100000,
                     ULR_data_dir = 'data/VLM/ULR7a_neu/',
                     NGL_data_dir = 'data/VLM/NGL14/',
                     JPL_data_dir = 'data/VLM/JPL14/',
                     GFZ_data_dir = 'data/VLM/GT3_neu/',
                     return_VLM_only = T,
                     convert_daily_to_monthly = T,
                     choose_multiple = FALSE,
                     verbose = F,
                     return_VLM_station_info_only = F,
                     return_yearly = FALSE) {
  #' Get VLM data from the nearest station
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character or int. 'default', 'ULR', 'NGL', 'JPL', 'GFZ'. Alternatively, you can use 0, 1, 2, 3, 4.
  #' The `default` method uses the code from Guanchao's original code to keep the backward compatibility.
  #' @param choose_which integer, choose the nth nearest station
  #' @param max_distance numeric, the maximum distance to search for the nearest station.
  #' If there is no station within the max_distance, an empty data frame will be returned.
  #' @param return_VLM_only logical, if TRUE, only return VLM data, if FALSE, return more details including distance and speed in 3 directions.
  #' @param convert_daily_to_monthly logical, if TRUE, convert daily VLM data to monthly VLM data.
  #' @param choose_multiple logical, if TRUE, choose multiple stations from 1 to choose_which. Otherwise, choose the choose_whichth nearest station.
  #' @param verbose logical, if TRUE, print the info of the station.
  #' @param return_VLM_station_info_only logical, if TRUE, only return the info of the station.
  #' @param return_yearly logical, if TRUE, return yearly data instead of monthly data.
  #' @return VLM data. If `return_yearly` is TRUE, the data will be converted to yearly data.
  #' @details Distance is calculated using the distHaversine function from the geosphere package. The distHaversine function uses the Haversine formula, which is a good general-purpose formula for distances on the Earth's surface.
  #' The alternative distVincentySphere function is a bit more accurate for longer distances.
  #' @export
  #' @examples
  #' get_VLM(lat = 40, lon = 40, method = 'ULR')
  #' get_VLM_by_station(station.names[3], method = 'ULR', choose_which = 5, choose_multiple = T, return_VLM_station_info_only = T, max_distance = 50)
  #'
  if (method == 'default' | method == 0) {
    return(get_VLM_old(lat, lon, verbose))
  }
  else if (method == 'ULR' | method == 1) {
    station_table <- read_ULR_table()
    file_suffix <- paste0('_ULR7.neu')
    data_dir <- ULR_data_dir
  }
  else if (method == 'NGL' | method == 2) {
    station_table <- read_NGL_table()
    file_suffix <- paste0('_NGL14.neu')
    data_dir <- NGL_data_dir
  }
  else if (method == 'JPL' | method == 3) {
    station_table <- read_JPL_table()
    file_suffix <- paste0('_JPL14.neu')
    data_dir <- JPL_data_dir
  }
  else if (method == 'GFZ' | method == 4) {
    station_table <- read_GFZ_table()
    file_suffix <- paste0('_GT3.neu')
    data_dir <- GFZ_data_dir
  }
  else {
    stop(paste0('Invalid method. Please use "default", "ULR", "NGL", "JPL", "GFZ", or 0, 1, 2, 3, 4. You used: ', method, '.'))
  }

  station_table$distance <- distHaversine(c(lon, lat), station_table[, c('Lon', 'Lat')]) / 1000
  station_table <- station_table[station_table$distance < max_distance, ]
  empty_data_frame <- data.frame(month = character(), VLM = numeric())
  if (nrow(station_table) == 0) {
    return(empty_data_frame)
  }
  if (choose_which > nrow(station_table)) {
    warning(paste0('There are only ', nrow(station_table), ' stations within ', max_distance, ' km.'))
    if (choose_multiple) {
      choose_which <- nrow(station_table)
    }
    else {
      return(empty_data_frame)
    }
  }
  station_table <- station_table[order(station_table$distance), ]
  get_station_path <- function(choose_which) {
    station_name <- paste0('d', station_table$Site[choose_which], '_', station_table$DOMES[choose_which])
    station_path <- paste0(data_dir, station_name, file_suffix)
    return(station_path)
  }

  check_file_exist <- function(choose_which) {
    station_path <- get_station_path(choose_which)
    return(file.exists(station_path))
  }

  while((!check_file_exist(choose_which)) && (choose_which > 1)) {
    choose_which <- choose_which - 1
  }

  station_path <- get_station_path(choose_which)
  is_file_exist <- file.exists(station_path)
  while ((!is_file_exist) && (choose_which < nrow(station_table))) {
    choose_which <- choose_which + 1
    is_file_exist <- check_file_exist(choose_which)
    if (verbose) {
      print(paste0('The ', choose_which, 'th nearest station is not available. Try the ', choose_which + 1, 'th nearest station.'))
    }
  }


  station_path <- get_station_path(choose_which)

  if (verbose) {
    print(paste0('The requested location: ', lat, ', ', lon))
    print(paste0('The info of the nearest station: '))
    print(station_table[choose_which, ])
  }
  if (return_VLM_station_info_only) {
    if (choose_multiple) {
      station_table <- station_table[1:choose_which, ]
      station_table <- station_table %>%
        mutate(site = paste(method, 1:choose_which, sep = '_'))
      return(station_table)
    }
    station_table <- station_table[choose_which, ]
    station_table$site <- paste(method, choose_which, sep = '_')
    return(station_table)
  }

  neu <- read_neu(station_path)

  if (return_VLM_only) {
    neu  <- neu %>%
      rename(VLM = DU, month = Year) %>%
      dplyr::select(month, VLM)

    if (convert_daily_to_monthly | choose_multiple | return_yearly) {
      neu <- convert_daily_to_monthly(neu)
    }
  }

  if (choose_multiple) {

    for (i in 1:(choose_which - 1)) {
      if (check_file_exist(i)) {
        station_path <- get_station_path(i)
        neu1 <- read_neu(station_path)
        neu1  <- neu1 %>%
          rename(!!(paste("VLM", i, sep = "_")) := DU, month = Year) %>%
          dplyr::select(month, !!paste("VLM", i, sep = "_"))
        if (convert_daily_to_monthly) {
          neu1 <- convert_daily_to_monthly(neu1)
        }
        neu <- merge(neu, neu1, by = 'month', all = TRUE)
      }
    }
  }

  if (return_yearly) {
    neu <- convert_monthly_to_yearly(neu)
  }

  return(neu)
}

get_VLM_longest <- function(lat, lon,
                method = 'default',
                choose_which = 3,
                convert_daily_to_monthly = TRUE,
                return_VLM_station_info_only = FALSE,
                return_yearly = FALSE,
                ...) {
  #' Get VLM data with the longest time series from the nearby stations
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character or int, 'default', 'ULR', 'NGL', 'JPL', 'GFZ', or 0, 1, 2, 3, 4.
  #' @param choose_which integer, the maximum number of stations to choose from. When choose_which = 1, it is the same as get_VLM.
  #' @param convert_daily_to_monthly logical, it has to be TRUE for get_VLM_longest.
  #' @param ... other parameters passed to .get_VLM. Run `docstring(.get_VLM)` or `?.get_VLM` for details.
  #' @return VLM data
  #' @export
  #' @examples
  #' get_VLM_longest(lat = 40, lon = 40, method = 'ULR', choose_which = 3, convert_daily_to_monthly = F)
  #'
  #'




  longest_station <- 0
  if (choose_which == 1) {
    return(.get_VLM(lat=lat, lon=lon, method = method, choose_which = choose_which, ...))
  }
  for (.choose_which in 1:choose_which) {
    vlm_data <- .get_VLM(lat, lon, method = method, choose_which = .choose_which,
                        return_VLM_station_info_only = FALSE,
                        ...)
    if (nrow(vlm_data) > longest_station) {
      longest_station <- nrow(vlm_data)
      longest_station_data <- vlm_data
      if (return_VLM_station_info_only) {
        longest_station_data <- .get_VLM(lat=lat, lon=lon, method = method, choose_which = .choose_which,
                                        return_VLM_station_info_only = TRUE,
                                        ...)
        longest_station_data$site <- paste(method, .choose_which, sep = '_')
        if (.choose_which > 1) {

          closest_station <- .get_VLM(lat=lat, lon=lon, method = method, choose_which = 1,
                                    return_VLM_station_info_only = TRUE,
                                    ...)
          closest_station$site <- paste(method, 1, sep = '_')
          longest_station_data <- rbind(closest_station, longest_station_data)
        }

      }
    }
  }
  if (return_yearly) {
    longest_station_data <- convert_monthly_to_yearly(longest_station_data)
  }
  return(longest_station_data)
}

get_VLM_mean <- function(lat, lon,
            methods = c('ULR', 'NGL', 'JPL', 'GFZ'),
            choose_which = 1,
            return_VLM_only = TRUE,
            convert_daily_to_monthly = TRUE,
            return_VLM_station_info_only = FALSE,
            return_yearly = FALSE,
            ...) {
  #' Get VLM data Using the four methods and take the mean.
  #' If one method is not available, impute it with the mean of the other methods.
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param methods character vector, the methods to use.
  #' @param choose_which integer, the maximum number of stations to choose from
  #' @param convert_daily_to_monthly logical, it has to be TRUE for get_VLM_mean. It doesn't really have any effect.
  #' @param ... other parameters passed to get_VLM
  #' @return VLM data
  #' @examples
  #' get_VLM_mean(lat = 40, lon = 40, methods = c('ULR', 'NGL', 'JPL', 'GFZ'))
  #' get_VLM_mean(lat = 30, lon = 30, methods = c('ULR', 'NGL', 'JPL', 'GFZ'), return_VLM_station_info_only = T)
  #'

  for (method in methods) {
    vlm_data <- .get_VLM(lat=lat, lon=lon, method = method, choose_which = choose_which, convert_daily_to_monthly = T, return_VLM_only = TRUE, return_VLM_station_info_only = return_VLM_station_info_only, ...)
    if (method == methods[1]) {
      if (return_VLM_station_info_only) {
        vlm_data_all <- vlm_data
        vlm_data_all$site <- paste(method, choose_which, sep = '_')
      }
      else {
        vlm_data_all <- vlm_data  %>%
          rename(!!(method) := VLM)
      }
    }
    else {
      if (return_VLM_station_info_only) {
        vlm_data$site <- paste(method, 1, sep = '_')
        vlm_data_all <- rbind(vlm_data_all, vlm_data)
      }
      else {
        vlm_data_all <- merge(vlm_data_all, vlm_data, by = 'month', all = TRUE) %>%
          rename(!!(method) := VLM)
      }
    }
  }

  if (return_VLM_station_info_only) {
    return(vlm_data_all)
  }

  vlm_data_all <- vlm_data_all %>%
    mutate_at(vars(methods), ~ifelse(is.na(.), mean(., na.rm = T), .))
  vlm_data_all <- vlm_data_all %>%
    mutate(VLM = rowMeans(.[methods], na.rm = T)) %>%
    dplyr::select(month, VLM)

  if (return_yearly) {
    vlm_data_all <- convert_monthly_to_yearly(vlm_data_all)
  }
  return(vlm_data_all)

}

get_VLM_circlemean <- function(lat, lon,
                               method = c('ULR', 'NGL', 'JPL', 'GFZ'),
                               max_distance = 300,
                               choose_multiple = TRUE,
                               convert_daily_to_monthly = TRUE,
                               na.rm = TRUE,
                               choose_which = 5,
                               return_VLM_station_info_only = FALSE,
                               return_yearly = FALSE,
                               verbose = FALSE,
                               ...){
  #' Get VLM data Using the mean of the stations within a circle.
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method character or character vector, 'default', 'ULR', 'NGL', 'JPL', 'GFZ', or 0, 1, 2, 3, 4.
  #' @param max_distance numeric, the maximum distance to search for the nearest station. The unit is km.
  #' @param ... other parameters passed to .get_VLM
  #' @return data frame with VLM data
  #' @examples
  #' get_VLM_circlemean(lat = 30, lon = 30, method = 'ULR', max_distance = 50)
  #' get_VLM_circlemean(lat = 30, lon = 30, method = c('ULR', 'NGL', 'JPL', 'GFZ'), max_distance = 50)
  #'
  #'
  datalist <- list()
  for (.method in method) {
    data <- .get_VLM(lat=lat, lon=lon, method = .method, max_distance = max_distance,
                     choose_which = choose_which, choose_multiple = TRUE,
                     return_VLM_station_info_only = return_VLM_station_info_only, ...)
    if (return_VLM_station_info_only) {
      data <- data %>%
        mutate(site = paste0(.method, '_', row_number()))
    }
    else {
      data <- data %>%
        rename_with(~ paste(.method, .x, sep = '_'), -month)
    }
    datalist[[length(datalist) + 1]] <- data
  }
  if (return_VLM_station_info_only) {
    data <- Reduce(function(x, y) rbind(x, y), datalist)
    return(data)
  }
  data <- Reduce(function(x, y) merge(x, y, by = 'month', all = TRUE), datalist)
  data <- data %>%
    rowwise() %>%
    mutate(VLM = mean(c_across(ends_with("VLM")), na.rm = na.rm)) %>%
    dplyr::select(month, VLM)
  if (return_yearly) {
    data <- convert_monthly_to_yearly(data)
  }
  return(data)
}

IDW_weight <- function(d, d_max, p = 2) {
  return(1 / (d^p))
}

get_weighted_circle_mean <- function(lat, lon,
                               method = c('ULR', 'NGL', 'JPL', 'GFZ'),
                               max_distance = 300,
                               choose_multiple = TRUE,
                               convert_daily_to_monthly = TRUE,
                               na.rm = TRUE,
                               choose_which = 100,
                               return_VLM_station_info_only = FALSE,
                               verbose = FALSE,
                               ...){

}

get_VLM <- function(lat, lon,
                    method = 'ULR',
                    VLM_retreive_method = 'deterministic',
                    return_yearly = FALSE,
                    ...) {
  #' Get VLM data from the nearest station
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method Character, integer, or vector. Accepts 'default', 'ULR', 'NGL', 'JPL', 'GFZ' or the integers 0 through 4. The argument name aligns with VLM_method used elsewhere.
  #' When `VLM_retreive_method` is `mean` or `circlemean`, `method` can be a vector of methods that are used to calculate the mean.
  #' The `default` method uses the code from Guanchao's original code to keep the backward compatibility.
  #' @param VLM_retreive_method: character, 'deterministic', 'longest', 'mean' or 'circlemean'. 'deterministic' uses the exact specified station.
  #' 'longest' uses the station with the longest time series constrained by `choose_which`. Alternatively, you can `get_VLM_longest()`
  #' 'mean' uses the mean of the VLM data from different stations constrained by `choose_which`. Alternatively, you can `get_VLM_mean()`.
  #' 'circlemean' uses the mean of the VLM data from different stations within a circle constrained by `max_distance`. Alternatively, you can `get_VLM_circlemean()`.
  #' @param choose_which integer, choose the nth nearest station. When `VLM_retreive_method` is `longest`, it is the maximum number of stations to choose from.
  #' @param return_VLM_only logical, if TRUE, only return VLM data, if FALSE, return more details including distance and speed in 3 directions.
  #' @param convert_daily_to_monthly logical, if TRUE, convert daily VLM data to monthly VLM data.
  #' @param verbose logical, if TRUE, print the info of the station.
  #' @param return_VLM_station_info_only logical, if TRUE, only return the info of the station.
  #' @param return_yearly logical, if TRUE, return yearly data instead of monthly data.
  #' @param ... other parameters passed to `.get_VLM`. Run `docstring(.get_VLM)` or `?.get_VLM` for details.
  #' @return VLM data
  #' @export
  #' @examples
  #' get_VLM(lat = 40, lon = 40, method = 'ULR')
  #' get_VLM(lat = 40, lon = 40, method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = 'mean')
  #' get_VLM_mean(lat = 40, lon = 40, methods = c('ULR', 'NGL', 'JPL', 'GFZ'))
  #' get_VLM(lat = 40, lon = 40, method = 'NGL', VLM_retreive_method = 'longest', choose_which = 3, verbose = T)
  #' get_VLM_longest(lat = 40, lon = 40, method = 'NGL', choose_which = 3, verbose = T)
  #' get_VLM(lat = 30, lon = 30, method = 'ULR', VLM_retreive_method = 'circlemean', max_distance = 500)
  #' get_VLM(lat = 30, lon = 30, method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = 'circlemean', max_distance = 500)
  #' get_VLM(lat = 40, lon = 40, method = 'ULR', return_yearly = TRUE)
  result <- if (VLM_retreive_method == 'deterministic') {
    .get_VLM(lat=lat, lon=lon, method = method, return_yearly = return_yearly, ...)
  }
  else if (VLM_retreive_method == 'longest') {
    get_VLM_longest(lat=lat, lon=lon, method = method, return_yearly = return_yearly, ...)
  }
  else if (VLM_retreive_method == 'mean') {
    get_VLM_mean(lat=lat, lon=lon, methods = method, return_yearly = return_yearly, ...)
  }
  else if (VLM_retreive_method == 'circlemean') {
    get_VLM_circlemean(lat=lat, lon=lon, method = method, return_yearly = return_yearly, ...)
  }
  else {
    stop(paste0('Invalid VLM_retreive_method. Please use "deterministic", "longest", "mean", or "circlemean". You used: ', VLM_retreive_method, '.'))
  }

  return(result)
}


get_VLM_by_station <- function(station_name, data_dir=DATA_DIR, return_yearly = FALSE, ...){
  #' get VLM data by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory containing the data.
  #' @param return_yearly logical, if TRUE, return yearly data instead of monthly data.
  #' @param ... other parameters passed to get_VLM. See `docstring(get_VLM)` or `?get_VLM` for details.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' get_VLM_by_station("astoria")
  #' get_VLM_by_station("astoria", method="JPL", choose_which=2)
  #' get_VLM_by_station("astoria", return_yearly = TRUE)
  station_id <- get_station_id(station_name)
  station_info <- station2lonlat(station_id, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  return(get_VLM(lat=lat, lon=lon, return_yearly = return_yearly, ...))
}

vlm_plot <- function(lat, lon, method = 1, choose_which = 1, ...) {
  #' plot VLM data based on the location
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character or int, 'default', 'ULR', 'NGL', 'JPL', 'GFZ', or 0, 1, 2, 3, 4.
  #' @param choose_which integer, choose the nth nearest station
  #' @return a ggplot object
  vlm_data <- get_VLM(lat=lat, lon=lon, method = method, choose_which = choose_which, ...)
  p <- ggplot(vlm_data, aes(x = month, y = VLM)) +
    geom_line() +
    theme_bw() +
    labs(x = "Month", y = "VLM (m)")
  return(p)
}

vlm_plot_full <- function(lat, lon, methods = c("ULR", "NGL", "JPL", "GFZ"), ...) {
  vlm_data_list <- list()

  for (method in methods) {
    vlm_data_list[[method]] <- get_VLM(lat = lat, lon = lon, method = method)
  }

  circlemean_data <- get_VLM_circlemean(lat = lat, lon = lon, ...)
  vlm_data_list[["circlemean"]] <- circlemean_data

  combined_data <- do.call(rbind, lapply(names(vlm_data_list), function(method) {
    data.frame(month = vlm_data_list[[method]]$month, VLM = vlm_data_list[[method]]$VLM, method = method)
  }))

  y_min <- min(combined_data$VLM, na.rm = TRUE)
  y_max <- max(combined_data$VLM, na.rm = TRUE)

  y_range <- y_max - y_min
  y_min <- y_min - 0.1 * y_range
  y_max <- y_max + 0.1 * y_range

  break_size <- ifelse(y_range > 1, 0.2, ifelse(y_range > 0.5, 0.1, 0.05))

  custom_colors <- c("ULR" = "#1f77b4", "NGL" = "#ff7f0e", "JPL" = "#2ca02c",
                     "GFZ" = "#d62728", "circlemean" = "#000000")

  p <- ggplot(combined_data, aes(x = month, y = VLM, color = method, group = method)) +
    geom_line(data = subset(combined_data, method != "circlemean"), size = 1, alpha = 0.7) +
    geom_line(data = subset(combined_data, method == "circlemean"), size = 1.5) +
    scale_color_manual(values = custom_colors) +
    scale_x_date(limits = as.Date(c("1995-01-01", "2025-01-01")),
                date_breaks = "5 years", date_labels = "%Y") +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = seq(floor(y_min/break_size)*break_size,
                  ceiling(y_max/break_size)*break_size,
                  by = break_size)
    ) +
    theme_minimal(base_size = 14) +
    labs(x = "Year", y = "VLM (m)",
        title = "Vertical Land Motion (VLM) Data for Different Methods",
        subtitle = paste("Location:", round(lat, 2), "°N,", round(lon, 2), "°E")) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "darkgrey"),
      axis.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))

  return(p)
}



vlm_plot_station <- function(station_name, data_dir=DATA_DIR, method = "ULR", ...) {
  #' get VLM data by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @param ... other parameters passed to get_VLM_by_station and vlm_source_plot.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' vlm_source_plot_by_station("astoria")
  #' vlm_source_plot_by_station("astoria",method="JPL", choose_which=2)
  station_info <- station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  p <- vlm_plot(lat, lon, method = method, ...)
  method_dict <- vlm_rename_dict("method")
  if (method %in% names(method_dict)) {
    method <- method_dict[[method]]
  }
  p <- p + labs(title = paste(station_name, 'VLM (m) from', method))
  return(p)
}

vlm_plot_full_station <- function(station_name, data_dir=DATA_DIR, ...) {
  #' get VLM data by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @param ... other parameters passed to get_VLM_by_station and vlm_source_plot.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' vlm_source_plot_by_station("astoria")
  #' vlm_source_plot_by_station("astoria",method="JPL", choose_which=2)
  station_name = get_station_name(station_name)
  station_info <- get_ar6_station_latlon(station_name)
  lon <- station_info$lon
  lat <- station_info$lat
  p <- vlm_plot_full(lat=lat, lon=lon, ...)
  return(p)
}

vlm_sources_plots <- function(lat, lon, methods = c(1, 2, 3, 4), choose_which = 1) {
  vlm_data <- get_VLM(lat=lat, lon=lon, method = methods[1], choose_which = choose_which)
  vlm_data <- vlm_data %>%
    rename(!!as.character(methods[1]) := VLM)
  for (method in methods[-1]) {
    vlm_data <- merge(vlm_data, get_VLM(lat, lon, method = method, choose_which = choose_which), by = "month", all = TRUE)
    vlm_data <- vlm_data %>%
      rename(!!as.character(method) := VLM)
  }

  method_dict <- c("1" = "ULR", "2" = "NGL", "3" = "JPL", "4" = "GFZ")
  rename_dict <- c("ULR" = "1", "NGL" = "2", "JPL" = "3", "GFZ" = "4")

  for (m in methods) {
    if (as.character(m) %in% names(method_dict)) {
      vlm_data <- vlm_data %>%
        rename(!!method_dict[[as.character(m)]] := !!as.character(m))
    }
  }

  plot_data <- vlm_data %>%
    dplyr::select(month, any_of(c('ULR', 'NGL', 'JPL', 'GFZ'))) %>%
    pivot_longer(cols = -month, names_to = "source", values_to = "vlm")

  custom_colors <- c("ULR" = "#1f77b4", "NGL" = "#ff7f0e", "JPL" = "#2ca02c",
                    "GFZ" = "#d62728", "circlemean" = "#000000")

  y_min <- min(plot_data$vlm, na.rm = TRUE)
  y_max <- max(plot_data$vlm, na.rm = TRUE)
  y_range <- y_max - y_min
  y_min <- y_min - 0.1 * y_range
  y_max <- y_max + 0.1 * y_range
  break_size <- ifelse(y_range > 1, 0.2, ifelse(y_range > 0.5, 0.1, 0.05))

  p <- ggplot(plot_data, aes(x = month, y = vlm, color = source, group = source)) +
    geom_line(size = 1, alpha = 0.7) +
    scale_color_manual(values = custom_colors) +
    scale_x_date(limits = as.Date(c("1995-01-01", "2025-01-01")),
                date_breaks = "5 years", date_labels = "%Y") +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = seq(floor(y_min/break_size)*break_size,
                  ceiling(y_max/break_size)*break_size,
                  by = break_size)
    ) +
    theme_minimal(base_size = 14) +
    labs(x = "Year", y = "VLM (m)",
        title = "Vertical Land Motion (VLM) Data from Different Sources",
        subtitle = paste("Location:", round(lat, 2), "°N,", round(lon, 2), "°E")) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "darkgrey"),
      axis.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))

  return(p)
}

vlm_sources_plots_by_station <- function(station_name, data_dir=DATA_DIR, include_ar6_vlm = TRUE, ...) {
  #' Compare VLM data from different sources by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @param include_ar6_vlm logical, if TRUE, include AR6 VLM data
  #' @param ... other parameters passed to get_VLM_by_station and vlm_sources_plots.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' vlm_sources_plots_by_station("astoria")
  #' vlm_sources_plots_by_station("astoria", method="JPL", choose_which=2)
  station_info <- station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  p <- vlm_sources_plots(lat=lat, lon=lon, ...)

  if (include_ar6_vlm) {
    ar6_vlm_data <- get_ar6_vlm(station_name) %>%
      dplyr::filter(Year < 2025) %>%
      mutate(month = as.Date(paste0(Year, "-01-01"))) %>%
      dplyr::select(month, VLM) %>%
      dplyr::mutate(VLM = -VLM, source = "AR6 VLM")

    p <- p + geom_line(data = ar6_vlm_data,
                       aes(x = month, y = VLM, color = source),
                       size = 1)

    custom_colors <- c("ULR" = "#1f77b4", "NGL" = "#ff7f0e", "JPL" = "#2ca02c",
                       "GFZ" = "#d62728", "circlemean" = "#000000", "AR6 VLM" = "#9467bd")

    p <- p + scale_color_manual(values = custom_colors)
  }

  return(p)
}

.station2color <- function(station) {
  if (startsWith(station, 'ULR')) {
    return("#1b9e77")
  }
  else if (startsWith(station, 'NGL')) {
    return("#d95f02")
  }
  else if (startsWith(station, 'JPL')) {
    return("#7570b3")
  }
  else if (startsWith(station, 'GFZ')) {
    return("#e7298a")
  }
  else {
    return("#000000")
  }
}

station2color <- function(station) {
  #' get the color of the station
  #' @param station vector of characters, the name of the station
  #' @return a vector of colors
  #'
  #'
  if (is.character(station)) {
    return(sapply(station, .station2color))
  }
  else {
    return(.station2color(station))
  }
}

get_VLM_coords <- function(lat,lon, methods = c('ULR','NGL','JPL','GFZ'), ...){
  coords_df = data.frame(site = NULL,lat = NULL,lon=NULL,distance=NULL)
  i = 1
  for (meth in methods){
    station_table = get_VLM(lat=lat,lon=lon,method=meth,return_VLM_station_info_only=TRUE, ...)
    if ('site' %in% names(station_table)){
      tempdf = data.frame('site'=station_table$site,lat=station_table$Lat,lon=station_table$Lon,distance = station_table$distance)
    } else {

    }
    coords_df = rbind(coords_df,tempdf)
    i=i+1
  }
  coords_df <- coords_df  %>%
    mutate(color = station2color(site))

  return(coords_df)
}

get_loc_info <- function(station_name,data_dir=DATA_DIR,save=FALSE, ...){
  #' get VLM station info by station.
  #' @param station_name character, the name of the station
  #' @param save logical, if TRUE, save the data to a csv file
  #' @param ... other parameters passed to get_VLM_by_station and get_VLM_coords.
  #' @return a data frame with the VLM station information for the station. An example:
  #' @examples
  #' > get_loc_info("astoria")
  #'      site      lat       lon    distance   color
  #' 1 astoria 46.20667 -123.7683  0.00000000 #000000
  #' 2   ULR_1 46.20740 -123.7684  0.08176028 #1b9e77
  #' 3   NGL_3 46.11820 -122.8961 67.96752197 #d95f02
  #' 4   JPL_2 46.20490 -123.9561 14.46707541 #7570b3
  #' get_VLM_coords_by_station(station.names[6], method=1,choose_which = 11, VLM_retreive_method = "longest")
  station_info = station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  df = data.frame(site = station_name,lat = lat,lon=lon,distance=0, color="#000000")
  df = rbind(df,get_VLM_coords(lat=lat,lon=lon, ...))

  if (save){
    write.csv(df,file.path(data_dir,station_name,"VLM_distances.csv"))
  }

  return(df)
}

get_VLM_coords_by_station <- function(...) get_loc_info(...)

get_snapshot<-function(df,centlat,centlon){
  #' get a snapshot of the VLM stations around a location
  #' @param df data frame, the output of `get_VLM_coords` or `get_loc_info`.
  #' @param centlat numeric, the latitude of the center of the snapshot
  #' @param centlon numeric, the longitude of the center of the snapshot
  #' @return a data frame with the boundaries of the snapshot
  #' @examples
  #' df = get_VLM_coords(30,30)
  #' get_snapshot(df,30,30)
  radius_lat = 1.1*max(abs(df$lat-centlat))
  radius_lon = 1.1*max(abs(df$lon-centlon))
  r = max(radius_lon,radius_lat)

  df = data.frame(uplat = centlat+r,lowlat = centlat-r,uplon = centlon+r,lowlon=centlon-r)
  df
}

get_snapshot_station <- function(station_name, data_dir=DATA_DIR){
  #' get VLM station info by station.
  #' @param station_name character, the name of the station
  df <- get_loc_info(station_name, data_dir=data_dir)
  centlat <- df$lat[1]
  centlon <- df$lon[1]
  get_snapshot(df,centlat,centlon)
}

vlm_map_df <- function(coords){
  #' plots a map of the VLM stations around a location.
  #' @param coords data frame, usually the output of `get_VLM_coords` or `get_loc_info`.
  #' @return a ggplot object
  #' @examples
  #' get_VLM_coords_by_station(station.names[6], method=1,choose_which = 11, VLM_retreive_method = "longest")
  #'    %>% vlm_map_df()
  #' get_VLM_coords_by_station(station.names[6],choose_which = 11, VLM_retreive_method = "mean") %>% vlm_map_df()
  #' get_loc_info("astoria") %>% vlm_map_df()
  us_map = map_data("world")
  snap = get_snapshot(coords,coords$lat[1],coords$lon[1])
  print(coords$color)
  base_map <- ggplot() + geom_polygon(data=us_map, aes(x=long, y=lat, group=group),
                                      color="black", fill="white") +
    coord_cartesian(xlim = c(snap$lowlon,snap$uplon),ylim = c(snap$lowlat,snap$uplat))
  print(unique(coords$color))
  map_w_points <- base_map + geom_point(data=coords,aes(x=lon, y=lat, color = color),show.legend = FALSE) + scale_color_manual(values = unique(coords$color))
  map_w_points+ geom_text_repel(data = coords, max.overlaps = Inf, aes(x = lon, y = lat, label = site),
                                                                                                 box.padding = 0.25, point.padding = 0.25)
}

vlm_map <- function(station_name, data_dir = DATA_DIR,savedistance=FALSE,saveplot=FALSE){
  coords = get_loc_info(station_name,data_dir = data_dir,save=savedistance)

  map_w_points <- vlm_map_df(coords)

  if(saveplot){
    ggsave(file.path(data_dir,station_name,"VLM_Locations.jpg"))
  }

  map_w_points

}

vlm_map_df_live <- function(coords) {
  pal <- colorFactor(palette = "Set1", domain = coords$site)

  map <- leaflet(coords) %>%
    addTiles() %>%
    addProviderTiles("CartoDB.Positron") %>%

    addCircleMarkers(
      ~lon, ~lat,
      color = ~pal(site),
      radius = 8,
      fillOpacity = 0.8,
      stroke = FALSE,
      label = ~site,
      labelOptions = labelOptions(noHide = TRUE, direction = "auto", textOnly = TRUE)
    ) %>%

    addLegend("bottomright",
              pal = pal,
              values = ~site,
              title = "VLM Stations",
              opacity = 1
    ) %>%

    fitBounds(
      min(coords$lon) - 0.5,
      min(coords$lat) - 0.5,
      max(coords$lon) + 0.5,
      max(coords$lat) + 0.5
    )

  return(map)
}

vlm_map_station <- function(...) vlm_map(...)

vlm_map_live_station <- function(station_name, data_dir=DATA_DIR, ...) {
  coords = get_loc_info(station_name, data_dir=data_dir, ...)
  vlm_map_df_live(coords)
}

merge_VLM_df <- function(datasets, col_names = c('ULR', 'NGL', 'JPL', 'GFZ')) {
  #' Merge monthly VLM data
  #' @param datasets a list of data frames with monthly VLM data. The datasets should come from `get_VLM` or `get_VLM_by_station`.
  #' The column names should be `month` and `VLM`.
  #' @param col_names character vector, the col_names of the column names. It should have the same length as `datasets`.
  #' @return a data frame with merged monthly VLM data
  #' @export
  #' @examples
  #' monthly_VLM_1 <- get_VLM(lat = 40, lon = 40, method = 'ULR')
  #' monthly_VLM_2 <- get_VLM(lat = 40, lon = 40, method = 'NGL')
  #' merge_VLM(list(monthly_VLM_1, monthly_VLM_2), col_names = c('ULR', 'NGL', 'JPL', 'GFZ'))
  #'
  #'
  for (i in 1:length(datasets)) {
    datasets[[i]] <- datasets[[i]] %>%
      rename(!!(col_names[i]) := VLM) %>%
      mutate(month = as.Date(month))
  }
  Reduce(function(x, y) merge(x, y, by = 'month', all = TRUE), datasets)
}

merge_VLM_station <- function(station, data_dir = DATA_DIR, methods = c('ULR', 'NGL', 'JPL', 'GFZ')) {
  #' Merge monthly VLM data by station
  #' @param station character, the name of the station
  #' @param data_dir character, the directory containing the data.
  #' @param methods character vector, the methods to use. Default is c('ULR', 'NGL', 'JPL', 'GFZ').
  #' @return a data frame with merged monthly VLM data
  #' @export
  #' @examples
  #' merge_VLM_station(station.names[6])
  #'
  #'
  datasets <- list()
  for (method in methods) {
    datasets[[method]] <- get_VLM_by_station(station, data_dir = data_dir, method = method)
  }
  merge_VLM_df(datasets, col_names = methods)
}

if (sys.nframe() == 0){
  test_data <- get_VLM(lan_input=40,lon_input=40)
  test_data
}


get_VLM_ar6_prediction <- function(station, confidence="medium", scenario="ssp245", quantiles=c(50), interpolate_to_yearly=TRUE) {
  station_id <- get_station_id(station)
  file_path <- paste0("data/ar6/projection_tool_raw/ipcc_ar6_sea_level_projection_psmsl_id_", station_id, ".xlsx")

  if (!file.exists(file_path)) {
    if (stringr::str_detect(getwd(), "chaoj")) {
      system(paste0('C:/Users/chaoj/miniforge3/envs/RSL/python.exe "c:/Users/chaoj/OneDrive/R/Regional Sea Level/GNSSR-Plotting/utils/ar6.py" ', station_id, " -u"))
    }
    else if (stringr::str_detect(getwd(), "jichao")) {
      system(paste0('/gpfs/projects/ZhuGroup/miniconda3/bin/python utils/sa_data.py" ', station_id, " -u"))
    }
    else {
      system(paste0("python utils/ar6.py ", station_id, " -u"))
    }
  }
  vlm_data <- read_excel(file_path, sheet = "VerticalLandMotion")

  filtered_data <- vlm_data %>%
    dplyr::filter(confidence == !!confidence, scenario == !!scenario) %>%
    dplyr::filter(quantile %in% !!quantiles)

  long_data_list <- lapply(quantiles, function(q) {
    df <- filtered_data %>%
      dplyr::filter(quantile == q) %>%
      pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "VLM") %>%
      mutate(Year = as.numeric(Year),
             quantile = q)

    if (interpolate_to_yearly) {
      df <- df %>%
        ungroup() %>%
        complete(Year = full_seq(Year, 1), nesting(psmsl_id, process, confidence, scenario, quantile)) %>%
        arrange(Year) %>%
        mutate(VLM = zoo::na.approx(VLM, na.rm = FALSE))
    }

    return(df)
  })

  final_data <- bind_rows(long_data_list) %>%
    mutate(station = station)

  return(final_data %>% dplyr::select(Year, VLM, station, quantile))
}

plot_ar6_VLM_predictions <- function(station, confidence = "medium", scenario = "ssp245", quantiles = c(5, 50, 95), interpolate_to_yearly = TRUE) {
  #' Plot AR6 VLM predictions with mean and confidence intervals
  #' @param station character, name of the station.
  #' @param confidence character, confidence level (e.g., "medium")
  #' @param scenario character, scenario name (e.g., "ssp245")
  #' @param quantiles numeric vector, quantiles to use for the confidence intervals and mean; default is c(5, 50, 95)
  #' @param interpolate_to_yearly logical, if TRUE, interpolate data to yearly.
  #' @return a ggplot object with the VLM predictions
  #' @export
  #' @examples
  #' plot_ar6_VLM_predictions("astoria")

  prediction_data <- get_VLM_ar6_prediction(station, confidence = confidence, scenario = scenario, quantiles = quantiles, interpolate_to_yearly = interpolate_to_yearly)

  quantile_low <- min(quantiles)
  quantile_high <- max(quantiles)
  quantile_mean <- 50

  plot_data <- prediction_data %>%
    filter(quantile %in% c(quantile_low, quantile_mean, quantile_high)) %>%
    mutate(quantile = factor(quantile, levels = c(quantile_low, quantile_mean, quantile_high), labels = c("Low", "Mean", "High")))

  p <- ggplot(plot_data, aes(x = Year, y = VLM, color = quantile, group = quantile)) +
    geom_line(size = 1) +
    geom_ribbon(data = plot_data %>% filter(quantile %in% c("Low", "High")) %>% pivot_wider(names_from = quantile, values_from = VLM),
                aes(x = Year, ymin = Low, ymax = High), fill = "grey70", alpha = 0.5, inherit.aes = FALSE) +
    geom_line(data = plot_data %>% filter(quantile == "Mean"), size = 1, color = "blue") +
    labs(
      title = paste0("Vertical Land Motion Predictions for Station: ", station),
      subtitle = paste0("Confidence: ", confidence, ", Scenario: ", scenario),
      x = "Year",
      y = "Vertical Land Motion (m)",
      caption = "Source: IPCC AR6"
    ) +
    scale_x_continuous(breaks = seq(2020, 2100, by = 5), limits = c(2020, 2100)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key = element_rect(fill = "white", color = "grey80")
    ) +
    scale_color_manual(values = c("Low" = "red", "Mean" = "blue", "High" = "green")) +
    theme(legend.position = "bottom")

  return(p)
}

get_VLM_full <- function(station, ...) {
  #' Get full VLM data by station, including both historical and future data
  #' @param station character, the name of the station
  #' @return a data frame with full VLM data (always yearly data)
  #' @export
  #' @examples
  #' get_VLM_full(station.names[6])
  #'

  historical_VLM <- get_VLM_by_station(station, return_yearly = TRUE, ...)
  future_VLM <- get_VLM_ar6_prediction(station, ...)

  full_VLM <- bind_rows(historical_VLM, future_VLM)

  return(full_VLM)
}

get_ar6_vlm_rate <- function(station, confidence = "medium", scenario = "ssp245", quantile = 50) {
  #' Get the VLM rate of change
  #' @param station character, name of the station
  #' @param confidence character, confidence level of the prediction
  #' @param scenario character, climate scenario
  #' @param quantiles numeric vector, quantiles of the prediction
  #' @return float, the constant VLM rate of change

  vlm_predictions <- get_VLM_ar6_prediction(station, quantiles = c(quantile), interpolate_to_yearly = FALSE)

  vlm_predictions <- vlm_predictions %>%
    mutate(lagged_diff = VLM - dplyr::lag(VLM))

  lagged_diffs <- na.omit(vlm_predictions$lagged_diff)

  if(length(unique(lagged_diffs)) == 1) {
    return(unique(lagged_diffs) / 10)
  } else {
    warning("The VLM rate of change is not constant across the years.")
    return(mean(lagged_diffs, na.rm = TRUE) / 10)
  }
}

get_ar6_vlm <- function(station, start.year = 1993, end.year = 2100, confidence = "medium", scenario = "ssp245", quantile = 50, zero_year = -1) {
  #' Get the VLM based on AR6 predictions. It uses the value of 2020 and the VLM rate
  #' @param station character, name of the station
  #' @param start.year integer, the start year for the VLM calculation
  #' @param end.year integer, the end year for the VLM calculation
  #' @param confidence character, confidence level of the prediction
  #' @param scenario character, climate scenario
  #' @param quantiles numeric vector, quantiles of the prediction
  #' @param zero_year integer, the base year to make VLM zero. If -1, no adjustment is made.
  #' @return data.frame, VLM values from start year to end year. Yearly data. Column name is "VLM" and "Year"

  vlm_rate <- get_ar6_vlm_rate(station, confidence=confidence, scenario=scenario, quantile=quantile)

  years <- seq(start.year, end.year)

  vlm_2020 <- get_VLM_ar6_prediction(station, confidence=confidence, scenario=scenario, quantiles=c(quantile)) %>% dplyr::filter(Year == 2020) %>% pull(VLM)

  vlm_values <- rep(NA, length(years))
  vlm_values[years == 2020] <- vlm_2020

  for (i in 1:length(years)) {
    if (years[i] < 2020) {
      vlm_values[i] <- vlm_2020 - (2020 - years[i]) * vlm_rate
    } else if (years[i] > 2020) {
      vlm_values[i] <- vlm_2020 + (years[i] - 2020) * vlm_rate
    }
  }

  if (zero_year > -1 && zero_year %in% years) {
    zero_year_value <- vlm_values[years == zero_year]
    vlm_values <- vlm_values - zero_year_value
  }

  vlm_df <- data.frame(Year = years, VLM = vlm_values)

  return(vlm_df)
}

get_month_data_raw <- function(station,
                               data_dir = DATA_DIR,
                               adjust_by_mean=FALSE,
                               add_VLM = FALSE,
                               VLM_na_interp = TRUE,
                               equal_time=NULL,
                               drop_na_columns = c(),
                               VLM_method = 1,
                               from_file=TRUE,
                               ...) {
  #' Compared with `get_month_data_raw_old()`, there are 2 changes:
  #' 1. GNSSR data is no longer used in the analysis.
  #' 2. The data sources of TG and SA has been changed to support more stations globally and align with AR6.
  #' 3. Use `from_file=FALSE` to forcefully update the data from the Internet.
  month_sa <- get_sa_data(station, from_file = from_file)
  if ("sla" %in% colnames(month_sa)) {
    month_sa %<>% rename(sa = sla)
  } else {
    month_sa$sa <- NA
  }
  month_tg <- get_tg_data(station, from_file = from_file)
  month_tg$month <- as.Date(month_tg$month)

  sealev.df <- merge(month_sa, month_tg, all=TRUE)


  if(add_VLM){
    VLM_data <- get_VLM_by_station(station, data_dir=data_dir, method = VLM_method, convert_daily_to_monthly = T, ...)

    sealev.df %<>% left_join(VLM_data,by="month")
    if(VLM_na_interp){
      VLM_start <- min(which(!is.na(sealev.df$VLM)))
      VLM_end <- length(sealev.df$VLM)
      start_year <- as.numeric(sealev.df$month[VLM_start],0,4)
      start_month <- as.numeric(sealev.df$month[VLM_start],6,7)
      sealev.df$VLM[VLM_start:VLM_end] <- na.StructTS(ts(sealev.df$VLM[VLM_start:VLM_end],frequency = 12,start(start_year,start_month)))
    }
    sealev.df %<>% rename(sa_original=sa)
    sealev.df %<>% mutate(sa=sa_original-VLM)
  }

  if (!is.null(equal_time)) {
    if (equal_time) {
      drop_na_columns = c("sa", "tg")
    }
    else {
      drop_na_columns = c()
    }
    if (adjust_by_mean){
      sealev.df$sa=adjust_baseline(sealev.df$sa,sealev.df$tg,method='mean')$data1
      sealev.df$tg=adjust_baseline(sealev.df$tg,sealev.df$tg,method='mean')$data2
    }
  }



  sealev.df %>% drop_na({{drop_na_columns}})
}

get_year_data_raw <- function(station, data_dir = DATA_DIR, use_VLM_as_variable=TRUE, adjust_baseline=FALSE, ...) {
  #' Get yearly data by station
  #' @param station character, name of the station
  #' @param data_dir character, directory of the data
  #' @param use_VLM_as_variable logical, if TRUE, use VLM as a variable
  #' @param adjust_baseline logical. Must be set to FALSE. Currently only works on TG data.
  #' @param ... additional arguments to pass to `get_month_data_raw()` (get_VLM_station)
  #' @return a data frame with yearly data
  #' @export
  #' @examples
  #' get_year_data_raw(station.names[6])
  month_data <- get_month_data_raw(station, data_dir = data_dir, ...)
  year_data <- convert_monthly_to_yearly(month_data) %>% dplyr::rename(Year = year)
  if (use_VLM_as_variable) {
    start.year <- min(year_data$Year)
    end.year <- max(year_data$Year)
    VLM_ar6 <- get_ar6_vlm(station, start.year = start.year, end.year = end.year)
    VLM_ar6 <- VLM_ar6 %>% dplyr::rename(VLM_ar6 = VLM)
    year_data <- year_data %>% left_join(VLM_ar6, by = "Year")
  }
  if (adjust_baseline) {
    baseline <- get_baseline(station, method = "ar6")
    year_data <- year_data %>%
      mutate(tg = tg - baseline)
  }

  year_data
}

get_baseline <- function(station, year=2019, method="ar6") {
  #' Get the baseline sea level for a station
  if (method == "ar6") {
    get_year_data_raw(station, use_VLM_as_variable = FALSE) %>%
      dplyr::filter(Year >= 1995 & Year <= 2014) %>%
      summarise(mean_tg = mean(tg, na.rm = TRUE)) %>%
      pull(mean_tg) %>%
      return()
  } else{
    get_year_data_raw(station, use_VLM_as_variable = FALSE) %>%
      dplyr::filter(Year == {{year}}) %>%
      pull(tg) %>%
      return()
  }

}
