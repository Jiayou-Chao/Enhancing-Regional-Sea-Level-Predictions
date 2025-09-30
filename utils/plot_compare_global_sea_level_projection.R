library(docstring)
library(tidyverse)

plot_compare_global_sea_level_projection <- function(
  file_no_perma = "data/GMSL_prediction_SSP.csv",
  file_with_perma = "data/GMSL_prediction_SSP_perma.csv",
  title = "Global Mean Sea Level Projection Comparison",
  ylab = "GMSL anomaly (mm)",
  xlab = "Year"
) {
  "Plots and compares global mean sea level projections, showing the impact of permafrost feedback across two different scenarios (e.g., mitigation vs. no mitigation).

  @param file_no_perma Path to CSV file with GMSL predictions without permafrost feedback.
  @param file_with_perma Path to CSV file with GMSL predictions with permafrost feedback.
  @param title Plot title.
  @param ylab Y-axis label.
  @param xlab X-axis label.
  @return A ggplot object comparing the four projection series.

  Data frame columns used from each file:
    - Year: integer, year of projection
    - mean_no, up95_no, low95_no: numeric, projection for one scenario (e.g., no mitigation)
    - mean_cop, up95_cop, low95_cop: numeric, projection for another scenario (e.g., with mitigation)
  "

  process_file <- function(filepath, permafrost_label) {
    readr::read_csv(filepath, show_col_types = FALSE) %>%
      dplyr::select(Year, starts_with("mean_"), starts_with("low95_"), starts_with("up95_")) %>%
      tidyr::pivot_longer(
        cols = -Year,
        names_to = c(".value", "scenario"),
        names_pattern = "([a-z0-9]+)_(no|cop)"
      ) %>%
      dplyr::mutate(
        comparison_label = paste(permafrost_label, scenario, sep = " - ")
      )
  }

  df_no_perma <- process_file(file_no_perma, "Without Permafrost")
  df_with_perma <- process_file(file_with_perma, "With Permafrost")

  df_plot <- dplyr::bind_rows(df_no_perma, df_with_perma) %>%
    dplyr::mutate(
        comparison_label = case_when(
            comparison_label == "Without Permafrost - no" ~ "Without Permafrost - No Scenario",
            comparison_label == "Without Permafrost - cop" ~ "Without Permafrost - COP26 Scenario",
            comparison_label == "With Permafrost - no" ~ "With Permafrost - No Scenario",
            comparison_label == "With Permafrost - cop" ~ "With Permafrost - COP26 Scenario",
            TRUE ~ comparison_label
        )
    )


  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Year, y = mean, color = comparison_label)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = low95, ymax = up95, fill = comparison_label),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::labs(
      title = title,
      y = ylab,
      x = xlab,
      color = "Scenario",
      fill = "Scenario"
    ) +
    ggplot2::theme_minimal()

  return(p)
}