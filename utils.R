library(magrittr)
library(tidyverse)
library(corrplot)
library(leaflet)
library(ggplot2)
library(maps)
source("utils/country_name.R")

DATA_DIR <- "data/plots"
station.names <- c("astoria", "battery_park", "cape_charles", "charleston", "crescent_city",
"elly_oil_platform", "monterey", "newport", "south_beach",
"tofino")

convert_daily_VLM_to_monthly <- function(data, date_col = 'month', value_col = 'VLM', na.rm = F) {
  #' Converts daily VLM data to monthly VLM data. The output is YYYY-MM-01 format.
  #' Deprecated. Use convert_daily_to_monthly() instead.
  warning("Deprecated. Use convert_daily_to_monthly() instead.")
  if (value_col == 'VLM') {
    data1 <- data %>%
              mutate(month1 = month(get(date_col)), year = year(get(date_col))) %>%
              group_by(year, month1) %>%
              summarise_all(mean)
    data1$month <- as.Date(paste0(data1$year,'-',data1$month1,'-01'))
    data1 %<>%
      ungroup() %>%
      select(month, VLM)
    return(data1)

  }
  data <- data %>%
    mutate(month1 = month(get(date_col)), year = year(get(date_col))) %>%
    group_by(year, month1) %>%
    summarise(value = mean(get(value_col), na.rm = na.rm))
  data %>% select(-month1)
}

convert_daily_to_monthly <- function(data, month_col = 'month', na.rm = F) {
  #' Converts daily VLM data to monthly VLM data. The output is YYYY-MM-01 format.
  #' Columns with numeric values will be averaged. Columns with non-numeric values will be dropped.

  data_num <- data %>%
    dplyr::select(where(is.numeric), all_of(month_col))

  data_num <- data_num %>%
    mutate(month1 = lubridate::month(get(month_col)), year = lubridate::year(get(month_col))) %>%
    group_by(year, month1) %>%
    summarise_all(mean, na.rm = na.rm)

  data_num$month <- as.Date(NA_character_)

  valid_rows <- !is.na(data_num$year) & !is.na(data_num$month1)

  if(any(valid_rows)) {
    date_strings_valid <- paste(data_num$year[valid_rows],
                                data_num$month1[valid_rows],
                                "01", sep = "-")
    data_num$month[valid_rows] <- as.Date(date_strings_valid, format = "%Y-%m-%d")
  }

  data_num %>%
    ungroup() %>%
    dplyr::select(-month1, -year) %>%
    dplyr::select(month, everything())
}

merge_monthly_data <- function(datasets, by = 'month') {
  #' Merge monthly VLM data
  #' @param datasets a list of data frames with monthly VLM data
  #' @param by the column name to merge by
  #' @return a data frame with merged monthly VLM data
  #' @export
  #' @examples
  #' monthly_VLM_1 <- get_VLM_by_station(station.names[1], method = 1)
  #' monthly_VLM_2 <- get_VLM_by_station(station.names[1], method = 2)
  #' merge_monthly_data(c(monthly_VLM_1, monthly_VLM_2), by = 'month')
  for (dataset in datasets) {
    dataset[by] <- as.Date(dataset[[by]])
  }
  Reduce(function(x, y) merge(x, y, by = by, all = TRUE), datasets)
}

beautiful_corplot <- function(cor_data) {
  #' Make a beautiful correlation plot
  #' @param cor_data a data frame with correlation data, usually the output of cor().
  #' @return a correlation plot
  #' @export
  #' @examples
  #' data <- merge_VLM_station("south_beach")
  #' cor_data <- data %>% select(-month) %>% cor(use = "complete.obs")
  #' beautiful_corplot(cor_data)
  corrplot(cor_data,
         method = 'circle',
         type = 'lower',
         sig.level = c(.001, .01, .05),
         tl.pos="lt",
         tl.col="black", tl.cex=1.3,
         tl.offset=0.2,
         cl.pos="r",
         insig = "label_sig",
         pch.cex = 1.3,
         pch.col="red",
         cl.cex = 1.3)
  corrplot(cor_data,  type="upper", method="number",
          col="coral4",  tl.pos="n", cl.pos="n", number.cex = 1.2, add=T,diag=F)
  recordPlot()
}


create_map_with_point <- function(lat, lon, label = "Point") {
  map <- leaflet() %>%
    addTiles() %>%
    setView(lng = lon, lat = lat, zoom = 10) %>%
    addMarkers(lng = lon, lat = lat, label = label)

  return(map)
}

create_map_with_point_static <- function(lat, lon, label = "Point") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("maps", quietly = TRUE)) {
    install.packages("maps")
  }

  world_map <- map_data("world")

  static_map <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
    geom_point(aes(x = lon, y = lat), color = "red", size = 3) +
    geom_text(aes(x = lon, y = lat, label = label), vjust = -1, hjust = 0.5, size = 5) +
    coord_fixed(1.3) +
    labs(title = "Static Map with Point", x = "Longitude", y = "Latitude") +
    theme_minimal()

  return(static_map)
}

create_map_with_points <- function(points) {

  map <- leaflet() %>%
    addTiles()

  for (i in 1:nrow(points)) {
    map <- map %>%
      addMarkers(lng = points$lon[i], lat = points$lat[i], label = points$label[i])
  }

  return(map)
}

simple_monthly_plot <- function(data){
  data %>% pivot_longer(!month) %>% ggplot(aes(x=month,y=value)) + geom_line(aes(color=name),size=0.75)
}

simple_yearly_plot <- function(data){
  data %>% pivot_longer(!Year) %>% ggplot(aes(x=Year,y=value)) + geom_line(aes(color=name),size=0.75)
}

monthly_plot <- function(data, stat, method){
  simple_monthly_plot(data) + labs(title=paste(str_to_title(stat),"Monthly Mean Sea Level Measurement,",method,"adjustment"),x="Time",y="Height (m)")+
    theme(plot.title = element_text(hjust=0.5,size=32))
}

convert_monthly_to_yearly <- function(data, year_col = 'Year', na.rm = FALSE) {
  #' Converts monthly data to yearly data.
  #'
  #' @param data A data frame containing monthly data.
  #' @param year_col The name of the column to use for grouping by year. Default is 'year'.
  #' @param na.rm A logical value indicating whether to remove NA values when calculating means. Default is FALSE.
  #' @return A data frame with yearly data.
  #'
  #' @examples
  #' monthly_data <- data.frame(
  #'   month = seq(as.Date("2020-01-01"), as.Date("2021-12-31"), by = "month"),
  #'   value = runif(24)
  #' )
  #' yearly_data <- convert_monthly_to_yearly(monthly_data)

  if (!year_col %in% names(data)) {
    data[[year_col]] <- lubridate::year(data$month)
  }

  data_num <- data %>%
    dplyr::select(where(is.numeric), all_of(year_col))

  yearly_data <- data_num %>%
    group_by(!!sym(year_col)) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = na.rm)))

  return(yearly_data)
}

convert_monthly_to_yearly <- function(data, date_col = 'month', value_col='', na.rm=F) {
  if (value_col == '') {
    return(data %>%
             mutate(year = year(get(date_col))) %>%
             group_by(year) %>%
             summarise_all(mean))
  }
  data <- data %>%
    mutate(year = year(month)) %>%
    group_by(year) %>%
    summarise(value = mean(get(value_col), na.rm=na.rm))
  data
}

format_station_name <- function(station) {
  #' Convert station name to proper format with country name
  #' south_beach -> South Beach, United States
  #' elly_oil_platform -> Elly Oil Platform, United States
  #' HANASAKI II -> Hanasaki II, Japan
  if (station == "global") {
    return("Global")
  }
  station_id <- get_station_id(station)
  station <- get_station_name(station)
  station_coords <- get_ar6_station_latlon(station_id)
  country_name <- get_country_name(lat = station_coords$lat, lon = station_coords$lon)

  station <- gsub("_", " ", station)

  station <- tools::toTitleCase(tolower(station))

  roman_numerals <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")
  for (numeral in roman_numerals) {
    pattern <- paste0("\\b", tools::toTitleCase(tolower(numeral)), "\\b")
    station <- gsub(pattern, numeral, station)
  }

  common_abbrevs <- c("DC", "USA", "UK")
  for (abbrev in common_abbrevs) {
    pattern <- paste0("\\b", tools::toTitleCase(tolower(abbrev)), "\\b")
    station <- gsub(pattern, abbrev, station)
  }

  station <- paste0(station, ", ", country_name)
  return(station)
}

distance_to_interval <- function(lower, upper, value) {
  ifelse(value < lower, lower - value,
         ifelse(value > upper, value - upper, 0))
}

scatter_plot_df <- function(df, x_col, y_col) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)

  plot <- ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point() +
    theme_minimal() +
    labs(
      x = x_col,
      y = y_col,
      title = paste("Scatter Plot of", y_col, "vs", x_col)
    )

  print(plot)
}

create_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}
