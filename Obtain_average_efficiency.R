
library(lmodel2)


Obtain_average_efficiency <- function(x, y) {
  OLS_y = lm(y ~ x)
  y_hat = OLS_y$coefficients[1] + OLS_y$coefficients[2] * x
  e1 = 1
  OLS_x = lm(x ~ y)

  x_hat = OLS_x$coefficients[1] + OLS_x$coefficients[2] * y
  x_wave_OLS_y = (y - OLS_y$coefficients[1]) / OLS_y$coefficients[2]
  e2 = sum((x - x_hat) ^ 2, na.rm = TRUE) / sum((x - x_wave_OLS_y) ^ 2, na.rm =
                                                  TRUE)

  average_efficiency_OLS = (e1 + e2) / 2


  dataused = data.frame(a = y, b = x)
  model2 <- lmodel2(a ~ b, data = dataused)
  intercept_value = model2$regression.results[3, 'Intercept']
  slope_value = model2$regression.results[3, 'Slope']
  y_wave_GMR <- intercept_value + slope_value * x
  x_wave_GMR <- (y - intercept_value) / slope_value
  e1 = sum((y - y_hat) ^ 2, na.rm = TRUE) / sum((y - y_wave_GMR) ^ 2, na.rm =
                                                  TRUE)
  e2 = sum((x - x_hat) ^ 2, na.rm = TRUE) / sum((x - x_wave_GMR) ^ 2, na.rm =
                                                  TRUE)


  average_efficiency_GMR = (e1 + e2) / 2


  return(data.frame(average_efficiency_OLS, average_efficiency_GMR))
}

