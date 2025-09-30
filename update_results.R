library(knitr)
library(magrittr)
source("plot_US.R")
source("plot_global.R")

station_list <- c(832, 448, 1099, 12, 10, 167, 62, 982, 47)
plot_grid_backtesting_results(station_list = station_list)
ggsave("results/grid_backtesting_results.png", width = 12, height = 8)
plot_grid_ensemble_vs_SSPs(station_list, start.year=1970, stop.year=2019, from_file=TRUE, show_real=TRUE)
ggsave("results/grid_ensemble_vs_SSPs.png", width = 12, height = 8)
c(9, 10, 12, 13, 58, 111, 118, 127, 135) %>% plot_grid_long_backtesting_results()
ggsave("results/grid_long_backtesting_results.png", width = 12, height = 8)
ggsave("results/grid_long_backtesting_results.jpeg", width = 12, height = 8)
source("generate_global_data_integrity_results.R")
data_no <- get_valid_station_results(from_file=T, filters = list(missing_rate=0.3))
data_no <- append_AR6_predictions(data_no)
write_csv(data_no, "results/global_station_results_2100.csv")
data_no <- read_csv("results/global_station_results_2100.csv")
color_scale_range <- range(-1, 2)

combined_plot <- create_combined_global_maps(data_no, color_scale_range)
combined_plot
ggsave("results/combined_global_maps.png", combined_plot,
        width = 10, height = 8, dpi = 300)
ggsave("results/combined_global_maps.jpg", combined_plot,
        width = 10, height = 8, dpi = 300)
