#!/bin/bash

conda install -c conda-forge cmake -y
conda install -c conda-forge \
r-base r-essentials \
r-tidyverse r-maps r-ggplot2 r-ggrepel r-lubridate r-zoo  \
r-geosphere r-corrplot r-lmodel2 r-tseries r-forecast r-magrittr \
r-testthat r-rgooglemaps r-reticulate r-leaflet r-timeseries r-pracma \
r-astsa r-ggpubr r-lavaan r-regsem r-mass r-reshape2 \
r-rnaturalearth r-rnaturalearthdata r-changepoint \
r-changepoint.np r-sf -y

conda install conda-forge::r-devtools -y
conda install conda-forge::r-rnaturalearth -y
