r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

required_packages <- c(
  "tidyverse", "maps", "ggplot2", "ggrepel", "lubridate", "zoo", "tidyverse", "docstring", "geosphere",
  "corrplot", "lmodel2", "tseries", "forecast", "magrittr", "testthat",
  'RgoogleMaps', "reticulate", "leaflet", "timeSeries", "pracma", "astsa",
  "ggpubr", "lavaan", "regsem", "MTS", "MASS", "reshape2",
  "rnaturalearth", "rnaturalearthdata", "anomalize", "changepoint", "changepoint.np",
  "sf", "bdc", "nngeo"
)

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    renv::install(package)
  }
}

