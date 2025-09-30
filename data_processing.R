source("utils.R")
library(lubridate)
library(readr)
library(dplyr)

station_loc_id_map <- list(
    astoria = 256,
    battery_park = 12,
    cape_charles = 636,
    charleston = 234,
    crescent_city = 378,
    elly_oil_platform = 245,
    south_beach = 1196,
    monterey = 1352,
    newport = 351,
    tofino = 165
)

read_ar6_location_list <- function(file_path="data/ar6/location_list.txt", full=FALSE) {
  #' Read the location list from a file
  #' @param file_path character, the path to the location list file
  #' @param full logical, when FALSE, return only the location with names (station_id < 1000000000).
  #' These stations have data from PSMSL.
  #' @return data frame with columns: station_name, station_id, lat, lon

  location_list <- readr::read_delim(file_path, delim = "\t", col_names = c("station_name", "station_id", "lat", "lon"), show_col_types = FALSE)
  if (!full) {
    return(location_list %>% dplyr::filter(station_id < 1000000000))
  } else {
    return(location_list)
  }
}

get_ar6_station_id <- function(station_name) {
  if (station_name == "global") {
    return(0)
  }
  location_list <- read_ar6_location_list()
  station_info <- dplyr::filter(location_list, station_name == !!station_name)
  if (nrow(station_info) == 0) {
    stop(paste0("Station not found in the location list: ", station_name))
  }
  station_info$station_id
}

get_ar6_station_name <- function(station_id) {
  location_list <- read_ar6_location_list()
  station_info <- dplyr::filter(location_list, station_id == !!{{station_id}})
  if (nrow(station_info) == 0) {
    stop(paste0("Station not found in the location list: ", station_id))
  }
  station_info$station_name
}

get_ar6_station_info <- function(station_id) {
  location_list <- read_ar6_location_list()
  if (is.character(station_id)) {
    station_info <- dplyr::filter(location_list, station_name == !!{{station_id}})
  } else {
    station_info <- dplyr::filter(location_list, station_id == !!{{station_id}})
  }
  if (nrow(station_info) == 0) {
    stop(paste0("Station not found in the location list: ", station_id))
  }
  station_info
}

get_ar6_station_latlon <- function(station_id) {
  station_info <- get_ar6_station_info(station_id)
  list(lat = station_info$lat, lon = station_info$lon)
}

read_location_list <- function(...) {
  read_ar6_location_list(...)
}

find_closest_ar6_station_latlon <- function(lat, lon, location_list=read_location_list(), max_distance=300, top_k=1) {
  #' Find the closest station based on latitude and longitude
  #' @param lat numeric, the latitude of the point of interest
  #' @param lon numeric, the longitude of the point of interest
  #' @param location_list data frame, the location list with columns: station_name, station_id, lat, lon
  #' @param max_distance numeric, maximum allowed distance
  #' @param top_k integer, number of closest stations to return

  input_lat <- lat
  input_lon <- lon

  closest_stations <- location_list %>%
    dplyr::mutate(distance = sqrt((input_lat - lat)^2 + (input_lon - lon)^2)) %>%
    dplyr::arrange(distance) %>%
    dplyr::slice(1:top_k)

  if (any(closest_stations$distance > max_distance)) {
    stop("No station found within the max distance")
  }

  return(closest_stations %>% dplyr::pull(station_id))
}

find_closest_ar6_station <- function(station, location_list=read_location_list(), max_distance=2) {
  loc_id <- ifelse(station %in% names(station_loc_id_map), station_loc_id_map[[station]], -1)
  if (loc_id>0) {
    return(loc_id)
  }
  station_coords <- station2lonlat(station)
  return(find_closest_ar6_station_latlon(station_coords$lat, station_coords$lon, location_list, max_distance))
}

get_station_name <- function(station_id) {
  if (is.numeric(station_id)) {
    station_id = as.integer(station_id)
    return(get_ar6_station_name(station_id))
  }
  if (station_id %in% names(station_loc_id_map)) {
    return(station_id)
  }
  return(station_id)
}

get_station_id <- function(station_name) {
  if (is.numeric(station_name)) {
    return(as.integer(station_name))
  }
  if (station_name %in% names(station_loc_id_map)) {
    return(station_loc_id_map[[station_name]])
  }
  get_ar6_station_id(station_name)
}

get_tg_data <- function(id, convert_to_monthly=TRUE, from_file = TRUE, save_dir='data/tidalgauge/raw', min_end_date="2019-12-31") {
  #' Get TG data from PSMSL (https://psmsl.org/data/obtaining/).
  #' @param id character or numeric, the station id or name. Use `read_ar6_location_list()` to get the list of station names.
  #' @param convert_to_monthly logical, whether to convert the data to monthly. Default is TRUE.
  #' @param from_file logical, whether to read the data from file. Default is TRUE to save time.
  #' It will automatically download the data if the file does not exist even set to TRUE.
  #' Use `from_file = FALSE` to force update the data.
  #' @param save_dir character, the directory to save the data. Default is 'data/tidalgauge/raw'.
  #' @return data frame with columns: month, tg
  #' @examples
  #' get_tg_data(1, convert_to_monthly = TRUE, from_file = FALSE)
  if (is.character(id)) {
    if (id %in% names(station_loc_id_map)) {
      id <- station_loc_id_map[[id]]
    } else {
      id <- get_ar6_station_id(id)
    }
  }
  filename <- paste0(id, ".metdata")

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  if (!file.exists(file.path(save_dir, filename)) | from_file == FALSE) {
    url <- paste0("https://psmsl.org/data/obtaining/rlr.monthly.data/", id, ".rlrdata")
    download.file(url, destfile = file.path(save_dir, filename), mode = "wb")
  }

  if (file.info(file.path(save_dir, filename))$size == 0) {
    warning("Failed to download data or file is empty. Please check the URL or network connection.")
    return(data.frame(month=character(), tg=numeric()))
  }

  df <- read.table(file.path(save_dir, filename), sep=";", header=FALSE, strip.white=TRUE)
  names(df) <- c("Year", "tg", "tag1", "tag2")

  df$Year <- lubridate::date_decimal(df$Year)

  if (convert_to_monthly) {
    df$Year <- floor_date(df$Year, "month")
  }

  df$tg <- ifelse(df$tg == -99999, NA, df$tg)

  df <- df %>% dplyr::select(Year, tg) %>% rename(month = Year)

  end_date <- as.Date(min_end_date)

  df$month <- as.Date(df$month)

  if (max(df$month, na.rm = TRUE) < end_date) {
    df <- df %>%
      complete(month = seq.Date(from = min(df$month, na.rm = TRUE),
                                to = end_date,
                                by = "month"))
  }

  VLM_start <- min(which(!is.na(df$tg)))
  VLM_end <- length(df$tg)
  start_year <- as.numeric(format(df$month[VLM_start], "%Y"))
  start_month <- as.numeric(format(df$month[VLM_start], "%m"))

  df$is_na <- is.na(df$tg)
  df$tg[VLM_start:VLM_end] <- na.StructTS(ts(df$tg[VLM_start:VLM_end], frequency = 12,
                                  start = c(start_year, start_month)))
  df$tg <- df$tg / 1000


  df
}

get_sa_data <- function(station_id, from_file=TRUE, data_dir='data/satellite/processed', convert_to_monthly=TRUE, end_date="2025-01-01") {
  #' This function is no longer used.
  return(data.frame("month"=c(1993-01-01), "sla"=c(NA)))
}

get_era5_yearly_data <- function(station_id, from_file=TRUE, data_dir='data/era5data_processed', convert_to_monthly=TRUE, end_date="2025-01-01") {
  #' Get ERA5 data from the AR6 database
  #' @param station_id character or numeric, the station id or name. Use `read_ar6_location_list()` to get the list of station names.
  #' @param from_file logical, whether to read the data from file. Default is TRUE to save time.
  #' @param data_dir character, the directory to save the data. Default is 'data/era5/processed'.
  #' @return data frame with columns: Year, and other variables. The unit is dependent on the variable.
  #' @examples
  #' get_era5_yearly_data(1, from_file = FALSE)
  #'
  if (is.character(station_id)) {
    if (station_id %in% names(station_loc_id_map)) {
      station_id <- station_loc_id_map[[station_id]]
    } else {
      station_id <- get_ar6_station_id(station_id)
    }
  }

  file_path <- file.path(data_dir, paste0(station_id, "_era5_data.csv"))
  if (file.exists(file_path) && from_file) {
    df <- readr::read_csv(file_path, show_col_types = FALSE)
  } else {
    system(paste0("python utils/Reading_ERA5_data.py ", station_id, " -u"))
    df <- readr::read_csv(file_path, show_col_types = FALSE)
  }
  df <- df %>% dplyr::filter(Year < as.numeric(format(as.Date(end_date), "%Y")))
  df
}

station2lonlat <- function(station_name, data_dir=DATA_DIR){
  #' returns a list with the longitude and latitude of the station
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @export
  #' @return a list with the longitude and latitude of the station. The list has two elements: lon and lat.
  station_id <- get_station_id(station_name)
  station_name <- get_station_name(station_name)
  possible_station_names <- c("astoria", "battery_park", "cape_charles", "charleston",
    "crescent_city", "elly_oil_platform", "monterey", "newport",
    "south_beach", "tofino")
  if(!(station_name %in% possible_station_names)){
    station_info <- get_ar6_station_latlon(station_id)
  } else {
    station_info <- read_csv(file.path(data_dir,station_name,"GMR.csv"),show_col_types = FALSE)
  }
  return(list("lon"=station_info$lon,"lat"=station_info$lat))
}

global_station_latlon <- list(
  "london" = c(51.5074, -0.1278),
  "tokyo" = c(35.6895, 139.6917),
  "beijing" = c(39.9042, 116.4074),
  "new_york" = c(40.7128, -74.0060),
  "san_francisco" = c(37.7749, -122.4194),
  "sydney" = c(-33.8688, 151.2153),
  "pretoria" = c(-25.7479, 28.1881),
  "sao_paulo" = c(-23.5505, -46.6333)
)

convert_dates <- function(dates, row_numbers) {
  mapply(function(date_str, row_num) {
    if (!grepl("/", date_str)) return(as.character(ymd(date_str)))

    parts <- as.numeric(strsplit(date_str, "/")[[1]])
    year <- parts[3]

    if (year < 100) {
      century <- if(row_num >= 80 && row_num < 480) "19" else "20"
      date_str <- sprintf("%d/%d/%s%02d", parts[1], parts[2], century, year)
    }

    as.character(mdy(date_str))
  }, dates, row_numbers) %>% ymd()
}

get_global_sealevel_raw <- function(file_path = "data/global_sea_level_raw.txt",
                                     return_yearly = TRUE) {
  df <- read_tsv(file_path, comment = "#") %>%
    dplyr::select(-`...4`)

  df <- df %>%
    mutate(
      row_num = row_number(),
      date = convert_dates(Date, row_num),
      CW_2011_filled = coalesce(CW_2011, UHSLC_FD),
      UHSLC_FD_filled = coalesce(UHSLC_FD, CW_2011),
      GMSL = (CW_2011_filled + UHSLC_FD_filled) / 2
    ) %>%
    dplyr::select(-CW_2011_filled, -UHSLC_FD_filled, -row_num)

  if (return_yearly) {
    df <- df %>%
      group_by(Year = year(date)) %>%
      summarise(
        date = first(date),
        GMSL = mean(GMSL, na.rm = TRUE)
      ) %>%
      dplyr::select(-date) %>%
      mutate(GMSL = GMSL + 10.76)
  }

  return(df)
}

merge_global_sealevel_prediction <- function(save_csv = FALSE) {
  #' merge prediction data and historical data

  historical_data <- get_global_sealevel_raw(return_yearly = TRUE)
  prediction_data <- read_csv("data/GMSL_prediction_SSP_perma.csv", show_col_types = FALSE)
  merged_data <- historical_data %>%
    dplyr::select(Year, GMSL) %>%
    dplyr::rename(mean_no = GMSL) %>%
    dplyr::filter(Year < 1950) %>%
    bind_rows(prediction_data %>% dplyr::filter(Year >= 1950)) %>%
    dplyr::arrange(Year)
  if (save_csv) {
    write_csv(merged_data, "data/GMSL_historical_prediction_merged.csv")
  }
  return(merged_data)
}

select_ensemble_models <- function(ARDL_tg_increase_2100_no,
                                    ARDL_tg_vlm_increase_2100_no,
                                    VARX_tg_increase_2100_no) {
  model_mapping <- list(
    default = list(
      no = c("ARDL_no", "VARX_no"),
      cop = c("ARDL_cop", "VARX_cop")
    ),
    vlm_positive = list(
      no = c("ARDL_VLM_no", "VARX_no"),
      cop = c("ARDL_VLM_cop", "VARX_cop")
    )
  )

  if (ARDL_tg_increase_2100_no >= 0) {
    return(model_mapping$default)
  }

  if (ARDL_tg_vlm_increase_2100_no > 0) {
    return(model_mapping$vlm_positive)
  }

  increases <- c(
    ARDL_tg_increase_2100_no = ARDL_tg_increase_2100_no,
    ARDL_tg_vlm_increase_2100_no = ARDL_tg_vlm_increase_2100_no,
    VARX_tg_increase_2100_no = VARX_tg_increase_2100_no
  )

  max_increase_model <- names(which.max(increases))

  model_selection <- list(
    ARDL_tg_increase_2100_no = list(no = "ARDL_no", cop = "ARDL_cop"),
    ARDL_tg_vlm_increase_2100_no = list(no = "ARDL_VLM_no", cop = "ARDL_VLM_cop"),
    VARX_tg_increase_2100_no = list(no = "VARX_no", cop = "VARX_cop")
  )

  return(model_selection[[max_increase_model]])
}


create_map_station <- function(station) {
  latlon <- get_ar6_station_latlon(station)
  create_map_with_point(lat=latlon$lat, lon=latlon$lon, label=station)
}

create_map_static_station <- function(station) {
  latlon <- get_ar6_station_latlon(station)
  create_map_with_point_static(lat=latlon$lat, lon=latlon$lon, label=station)
}