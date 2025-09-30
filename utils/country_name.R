library(sf)
library(rnaturalearth)
library(bdc)
library(nngeo)
library(dplyr)

sf_use_s2_original <- sf::sf_use_s2()

get_country_name_1 <- function(lat, lon) {
  #' Get country name from latitude and longitude using the bdc package
  #'
  #' @param lat Latitude of the location
  #' @param lon Longitude of the location
  #' @return A character string containing the country name

  data <- data.frame(
    decimalLatitude = lat,
    decimalLongitude = lon,
    country = NA_character_
  )

  result <- bdc::bdc_country_from_coordinates(
    data = data,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    country = "country"
  )

  return(result$country)
}

get_country_name_2 <- function(lat, lon) {
  #' Get country name from latitude and longitude using spatial join
  #'
  #' @param lat Latitude of the location
  #' @param lon Longitude of the location
  #' @return A character vector containing the country abbreviation

  original_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)

  result <- tryCatch({
    countries <- ne_countries(returnclass = "sf")

    countries <- st_make_valid(countries)

    point <- st_as_sf(
      data.frame(lon = lon, lat = lat),
      coords = c("lon", "lat"),
      crs = 4326
    )

    st_join(point, countries)
  },
  error = function(e) {
    message("Error in spatial join: ", e$message)
    return(NULL)
  },
  finally = {
    sf::sf_use_s2(original_s2)
  })

  return(result)
}

get_country_abbreviation <- function(lat, lon) {
  result <- get_country_name_2(lat, lon)
  if(is.null(result)) return(NA_character_)
  return(result$abbrev)
}

get_country_with_eez <- function(lon, lat, buffer_km = 10) {
  #' Get country information considering marine Exclusive Economic Zones

  original_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)

  result <- tryCatch({
    countries <- ne_countries(
      scale = 50,
      type = "countries",
      returnclass = "sf"
    )

    countries <- st_make_valid(countries)

    point <- st_as_sf(
      data.frame(lon = lon, lat = lat),
      coords = c("lon", "lat"),
      crs = 4326
    )

    country_buffers <- st_buffer(
      st_transform(countries, 3857),
      dist = buffer_km * 1000
    ) %>% st_transform(4326)

    result <- st_join(point, country_buffers)

    if (is.na(result$name[1])) {
      nearest <- st_nn(point, countries, k = 1, progress = FALSE)
      result <- countries[unlist(nearest), ]
    }

    return(result)
  },
  error = function(e) {
    message("Error in EEZ search: ", e$message)
    return(NULL)
  },
  finally = {
    sf::sf_use_s2(original_s2)
  })

  return(result)
}

get_country_name <- function(lat, lon, buffer_km = 100) {
  #' Primary function to get country name from coordinates
  #' Uses multiple methods with fallbacks

  result <- get_country_name_2(lat=lat, lon=lon)

  if(!is.null(result) && !is.na(result$name_long[1])) {
    return(result$name_long[1])
  }

  result2 <- get_country_with_eez(lon=lon, lat=lat, buffer_km=buffer_km)
  if(!is.null(result2) && !is.na(result2$name_long[1])) {
    return(result2$name_long[1])
  }

  result3 <- get_country_name_1(lat=lat, lon=lon)
  if(!is.null(result3) && !is.na(result3[1])) {
    return(result3[1])
  }

  return(NA_character_)
}
