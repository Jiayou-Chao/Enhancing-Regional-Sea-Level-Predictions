library(MASS)
library(timeSeries)
library(tseries)
library(forecast)
library(ggplot2)
library(pracma)
library(astsa)
library(ggpubr)
library(dplyr)


predict_GMSL <- function(end_year) {
  pred_steps <- 2100 - end_year
  climate <- read.csv('data/YearData.csv')
  climate.l <- cbind(climate[-1, ], climate[-72, ])
  names(climate.l)[12:22] <-
    paste(names(climate.l)[12:22], '.l', sep = '')

  climate_train <-  climate.l %>%
    dplyr::filter(Year <= end_year)
  n <- nrow(climate_train)
  climate_test <-
    climate.l[(n + 1):nrow(climate.l), ]

  fit_TEMP_train <-
    lm(TEMP ~ TEMP.l + humidity + humidity.l + GWP.l, data = climate_train)
  fit_humidity_train <-
    lm(humidity ~ TEMP.l + humidity.l, data = climate_train)
  fit_SeaIce_train <-
    lm(SeaIce ~ SeaIce.l + Mass.l + TEMP, data = climate_train)
  fit_Mass_train <- lm(Mass ~ Mass.l + TEMP, data = climate_train)
  fit_GMSL_train <-
    lm(GMSL ~ GMSL.l + Mass.l + SeaIce.l + TEMP, data = climate_train)
  fit_TEMP_all <-
    lm(TEMP ~ TEMP.l + humidity + humidity.l + GWP.l, data = climate.l)
  fit_humidity_all <-
    lm(humidity ~ TEMP.l + humidity.l, data = climate.l)
  fit_SeaIce_all <-
    lm(SeaIce ~ SeaIce.l + Mass.l + TEMP, data = climate.l)
  fit_Mass_all <- lm(Mass ~ Mass.l + TEMP, data = climate.l)
  fit_GMSL_all <-
    lm(GMSL ~ GMSL.l + Mass.l + SeaIce.l + TEMP, data = climate.l)
  fit_GWP_train <- Arima(climate$GWP[1:(n+1)], order = c(1,1,1), include.drift = T, method = 'ML')


  pred_GWP_no <- forecast::forecast(fit_GWP_train, h = pred_steps)$mean

  newdata <- data.frame('GWP' = c(climate$GWP[n + 1], pred_GWP_no))
  pred_TEMP <- c(climate$TEMP[n + 1])
  pred_humidity <- c(climate$humidity[n + 1])
  pred_SeaIce <- c(climate$SeaIce[n + 1])
  pred_Mass <- c(climate$Mass[n + 1])
  pred_GMSL <- c(climate$GMSL[n + 1])
  for (i in 1:pred_steps) {
    pred_humidity <-
      c(pred_humidity,
        predict(
          fit_humidity_train,
          newdata = data.frame(TEMP.l = pred_TEMP[i], humidity.l = pred_humidity[i])
        ))
    pred_TEMP <-
      c(pred_TEMP, predict(
        fit_TEMP_train,
        newdata = data.frame(
          GWP.l = newdata[i, 'GWP'],
          TEMP.l = pred_TEMP[i],
          humidity = pred_humidity[i + 1],
          humidity.l = pred_humidity[i]
        )
      ))
    pred_Mass <-
      c(pred_Mass, predict(
        fit_Mass_train,
        newdata = data.frame(
          TEMP = pred_TEMP[i + 1],
          TEMP.l = pred_TEMP[i],
          SeaIce.l = pred_SeaIce[i],
          Mass.l = pred_Mass[i]
        )
      ))
    pred_SeaIce <-
      c(pred_SeaIce,
        predict(
          fit_SeaIce_train,
          newdata = data.frame(
            TEMP = pred_TEMP[i + 1],
            TEMP.l = pred_TEMP[i],
            SeaIce.l = pred_SeaIce[i],
            Mass.l = pred_Mass[i],
            Mass = pred_Mass[i + 1]
          )
        ))
    pred_GMSL <-
      c(pred_GMSL, predict(
        fit_GMSL_train,
        newdata = data.frame(
          TEMP = pred_TEMP[i + 1],
          TEMP.l = pred_TEMP[i],
          Mass = pred_Mass[i + 1],
          Mass.l = pred_Mass[i],
          GMSL.l = pred_GMSL[i],
          SeaIce = pred_SeaIce[i + 1],
          SeaIce.l = pred_SeaIce[i]
        )
      ))
  }
  vcv <-
    cov(
      cbind(
        fit_GWP_train$residuals[-1],
        fit_humidity_train$residuals,
        fit_TEMP_train$residuals,
        fit_Mass_train$residuals,
        fit_SeaIce_train$residuals,
        fit_GMSL_train$residuals
      )
    )
  vcv[1, 1] <- fit_GWP_train$sigma2
  vcv <- rbind(matrix(rep(0, 6 * (pred_steps - n)), ncol = 6), vcv)
  vcv <- cbind(matrix(rep(0, (pred_steps - n) * (pred_steps + 6 - n)), ncol = (pred_steps - n)), vcv)
  A <- eye(pred_steps + 6 - n)
  A[pred_steps + 3 - n, pred_steps + 2 - n] <- -fit_TEMP_train$coefficients['humidity']
  A[pred_steps + 4 - n, pred_steps + 3 - n] <- -fit_Mass_train$coefficients['TEMP']
  A[pred_steps + 5 - n, pred_steps + 3 - n] <- -fit_SeaIce_train$coefficients['TEMP']
  A[pred_steps + 6 - n, pred_steps + 3 - n] <- -fit_GMSL_train$coefficients['TEMP']



  B <- zeros(pred_steps + 6 - n)
  for (i in 1:(pred_steps - n))
    B[i, i + 1] <- 1
  B[(pred_steps + 1 - n), 1:(pred_steps + 1 - n)] <-
    rev(-ARMAtoAR(
      ar = c(1 + fit_GWP_train$coef[1], -fit_GWP_train$coef[1]),
      ma = fit_GWP_train$coef[2],
      lag.max = (pred_steps + 1 - n)
    ))


  B[(pred_steps + 2 - n), (pred_steps + 2 - n)] <- fit_humidity_train$coefficients['humidity.l']
  B[(pred_steps + 2 - n), (pred_steps + 3 - n)] <- fit_humidity_train$coefficients['TEMP.l']
  B[(pred_steps + 3 - n), (pred_steps + 1 - n)] <- fit_TEMP_train$coefficients['GWP.l']
  B[(pred_steps + 3 - n), (pred_steps + 2 - n)] <- fit_TEMP_train$coefficients['humidity.l']
  B[(pred_steps + 3 - n), (pred_steps + 3 - n)] <- fit_TEMP_train$coefficients['TEMP.l']
  B[(pred_steps + 4 - n), (pred_steps + 4 - n)] <- fit_Mass_train$coefficients['Mass.l']
  B[(pred_steps + 5 - n), (pred_steps + 4 - n)] <- fit_SeaIce_train$coefficients['Mass.l']
  B[(pred_steps + 5 - n), (pred_steps + 5 - n)] <- fit_SeaIce_train$coefficients['SeaIce.l']
  B[(pred_steps + 6 - n), (pred_steps + 4 - n)] <- fit_GMSL_train$coefficients['Mass.l']
  B[(pred_steps + 6 - n), (pred_steps + 5 - n)] <- fit_GMSL_train$coefficients['SeaIce.l']
  B[(pred_steps + 6 - n), (pred_steps + 6 - n)] <- fit_GMSL_train$coefficients['GMSL.l']

  sigma <- zeros(pred_steps + 6 - n)


  tmp <- inv(A) %*% vcv %*% t(inv(A))

  sigma_TEMP <- c()
  sigma_GMSL <- c()
  for(i in 1:(pred_steps)){
    sigma <- inv(A) %*% B %*% sigma %*% t(B) %*% t(inv(A)) + tmp
    sigma_TEMP <- c(sigma_TEMP, sigma[pred_steps + 3 - n, pred_steps + 3 - n])
    sigma_GMSL <- c(sigma_GMSL, sigma[pred_steps + 6 - n, pred_steps + 6 - n])
  }

    sigma_TEMP <- sqrt(sigma_TEMP)
    sigma_GMSL <- sqrt(sigma_GMSL)

    TEMP_multi <-
      data.frame(
        Year = end_year:2100,
        up95 = pred_TEMP + c(0, qnorm(0.975) *
                                                  sigma_TEMP),
        low95 = pred_TEMP - c(0, qnorm(0.975) *
                                                   sigma_TEMP),
        up99 =  pred_TEMP + c(0, qnorm(0.995) *
                                                  sigma_TEMP),
        low99 = pred_TEMP - c(0, qnorm(0.995) *
                                                   sigma_TEMP)
      )

    GMSL_multi <-
      data.frame(
        Year = end_year:2100,
        up95 = pred_GMSL + c(0, qnorm(0.975) *
                                                  sigma_GMSL),
        low95 = pred_GMSL - c(0, qnorm(0.975) *
                                                   sigma_GMSL),
        up99 = pred_GMSL + c(0, qnorm(0.995) *
                                                  sigma_GMSL),
        low99 = pred_GMSL - c(0, qnorm(0.995) *
                                                   sigma_GMSL),
        sigma = c(0,sigma_GMSL),
        mean = pred_GMSL
      )

    GMSL_multi

}

.unfinished <- function(end_year) {
  climate <- read.csv('data/YearData.csv')
  climate <- climate %>%
    filter(Year <= end_year)
  climate.l <- cbind(climate[-1,], climate[-nrow(climate),])
  names(climate.l)[12:22] <- paste(names(climate.l)[12:22],'.l', sep = '')


  fit_TEMP <- lm(TEMP ~ TEMP.l + humidity + humidity.l + GWP.l, data = climate.l)
  fit_humidity <- lm(humidity ~ TEMP.l + humidity.l, data = climate.l)
  fit_SeaIce <- lm(SeaIce ~ SeaIce.l + Mass.l + TEMP, data = climate.l)
  fit_Mass <- lm(Mass ~ Mass.l + TEMP, data = climate.l)
  fit_GMSL <- lm(GMSL ~ GMSL.l + Mass.l + SeaIce.l + TEMP, data = climate.l)

  fit_GWP <- Arima(climate$GWP,order = c(1,1,1), include.drift =  T, method = 'ML')
  pred_GWP_no <- forecast::forecast(fit_GWP, h = 2100 - end_year)
  pred_GWP_no$x <- ts(pred_GWP_no$x, start = 1950, end = end_year)
  pred_GWP_no$mean <- ts(pred_GWP_no$mean, start = end_year + 1, end = 2100)
  pred_GWP_no$lower <- ts(pred_GWP_no$lower, start = end_year + 1, end = 2100)
  pred_GWP_no$upper <- ts(pred_GWP_no$upper, start = end_year + 1, end = 2100)

  fit_N2O <- Arima(climate$N2O, order = c(1,1,1), include.drift = T, method = 'ML')
  pred_N2O_no <- forecast(fit_N2O, h = 2100 - end_year)$mean

  CO2_2010 <- climate %>% filter(Year == 2010) %>% select(CO2) %>% as.numeric()
  CO2_2009 <- climate %>% filter(Year == 2009) %>% select(CO2) %>% as.numeric()

  CO2_2010 <- CO2_2010 - CO2_2009
  pred_CO2_cop <- c(climate$CO2[72] + cumsum(CO2_2010 - (1:9)/9*(CO2_2010*0.45)))
  pred_CO2_cop <- c(pred_CO2_cop, rep(pred_CO2_cop[9],70))

  CH4_2020 <- climate %>% filter(Year == 2020) %>% select(CH4) %>% as.numeric()
  CH4_2019 <- climate %>% filter(Year == 2019) %>% select(CH4) %>% as.numeric()

  CH4_2020 <- climate$CH4[71] - climate$CH4[70]
  pred_CH4_cop <- c(climate$CH4[72] + cumsum(CH4_2020 - (1:9)/9*(CH4_2020*0.3)))
  pred_CH4_cop <- c(pred_CH4_cop, rep(pred_CH4_cop[9],70))

  pred_GWP_cop_mean <- pred_CO2_cop + 28/1000* pred_CH4_cop + 265/1000 * pred_N2O_no

}

forecast_multi_step <- function(end_year) {
  climate <- read.csv('data/YearData.csv')

  climate.l <- cbind(climate[-1,], climate[-72,])
  names(climate.l)[12:22] <- paste(names(climate.l)[12:22], '.l', sep = '')
  climate.l <- climate[-1, ]
  names(climate.l) <- paste(names(climate.l), 'l', sep = '')

  fit_TEMP <- lm(TEMP ~ TEMP.l + humidity + humidity.l + GWP.l, data = climate.l)
  fit_humidity <- lm(humidity ~ TEMP.l + humidity.l, data = climate.l)
  fit_SeaIce <- lm(SeaIce ~ SeaIce.l + Mass.l + TEMP, data = climate.l)
  fit_Mass <- lm(Mass ~ Mass.l + TEMP, data = climate.l)
  fit_GMSL <- lm(GMSL ~ GMSL.l + Mass.l + SeaIce.l + TEMP, data = climate.l)


}

forecast_multi_step <- function(end_year) {
  climate <- read.csv('data/YearData.csv')

  climate <- climate[climate$Year <= end_year, ]

  climate.l <- climate[-1, ]
  names(climate.l)[2:ncol(climate.l)] <- paste(names(climate.l)[2:ncol(climate.l)], 'l', sep = '')

  fit_TEMP <- lm(TEMP ~ TEMP.l + humidity + humidity.l + GWP.l, data = climate.l)
  fit_humidity <- lm(humidity ~ TEMP.l + humidity.l, data = climate.l)

  fit_GWP <- Arima(climate$GWP, order = c(1,1,1), include.drift = TRUE, method = 'ML')

  forecast_horizon <- end_year - 1950

  pred_GWP_no <- forecast::forecast(fit_GWP, h = forecast_horizon)

  pred_GWP_no$x <- ts(pred_GWP_no$x, start = 1950, end = end_year)
  pred_GWP_no$mean <- ts(pred_GWP_no$mean, start = end_year + 1, end = 2100)
  pred_GWP_no$lower <- ts(pred_GWP_no$lower, start = end_year + 1, end = 2100)
  pred_GWP_no$upper <- ts(pred_GWP_no$upper, start = end_year + 1, end = 2100)



  fit_N2O <- Arima(climate$N2O, order = c(1,1,1), include.drift = TRUE, method = 'ML')
  h_forecast <- end_year - max(climate$Year) + (2100 - end_year)
  pred_N2O_no <- forecast(fit_N2O, h = h_forecast)$mean

  CO2_2010_diff <- climate$CO2[61] - climate$CO2[60]
  pred_CO2_cop <- c(climate$CO2[72] + cumsum(CO2_2010_diff - (1:9)/9 * (CO2_2010_diff * 0.45)))
  pred_CO2_cop <- c(pred_CO2_cop, rep(pred_CO2_cop[9], 2100 - end_year))

  CH4_2020_diff <- climate$CH4[71] - climate$CH4[70]
  pred_CH4_cop <- c(climate$CH4[72] + cumsum(CH4_2020_diff - (1:9)/9 * (CH4_2020_diff * 0.3)))
  pred_CH4_cop <- c(pred_CH4_cop, rep(pred_CH4_cop[9], 2100 - end_year))

  pred_GWP_cop_mean <- pred_CO2_cop + (28/1000) * pred_CH4_cop + (265/1000) * pred_N2O_no

  pred_GWP_cop <- pred_GWP_no
  pred_GWP_cop$mean <- pred_GWP_cop$mean - pred_GWP_no$mean + ts(pred_GWP_cop_mean, start = 2022)

  pred_GWP_cop$lower[,1] <- pred_GWP_cop$mean - (pred_GWP_no$mean - pred_GWP_no$lower[,1]) * pred_GWP_cop$mean / pred_GWP_no$mean
  pred_GWP_cop$lower[,2] <- pred_GWP_cop$mean - (pred_GWP_no$mean - pred_GWP_no$lower[,2]) * pred_GWP_cop$mean / pred_GWP_no$mean
  pred_GWP_cop$upper[,1] <- pred_GWP_cop$mean + (pred_GWP_no$upper[,1] - pred_GWP_no$mean) * pred_GWP_cop$mean / pred_GWP_no$mean
  pred_GWP_cop$upper[,2] <- pred_GWP_cop$mean + (pred_GWP_no$upper[,2] - pred_GWP_no$mean) * pred_GWP_cop$mean / pred_GWP_no$mean



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

  estimate_no <- matrix(c(1, climate$GWP[71], climate$GWP[72], climate$humidity[72], climate$TEMP[72],
                          climate$Mass[72], climate$SeaIce[72], climate$GMSL[72]), nrow = 8)

  for (i in 1:79) {
    estimate_no <- cbind(estimate_no, c(1, inv(A) %*% B %*% estimate_no[, ncol(estimate_no)]))
    if (i == 1) estimate_no[3, ncol(estimate_no)] <- pred_GWP_no$mean[1]
  }

  pred_TEMP_no <- estimate_no[5,] + 0.306
  pred_GMSL_no <- estimate_no[8,] - mean(climate$GMSL[37:56])


}
