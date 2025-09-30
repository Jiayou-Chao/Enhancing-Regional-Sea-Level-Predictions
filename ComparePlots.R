library(ggplot2)
library(tidyverse)
library(lmodel2)
library(tseries)
library(zoo)
source("VLM.R")
source("Obtain_average_efficiency.R")
source("utils.R")
library(magrittr)
source("models/predict_GMSL.R")
source("data_processing.R")
source("ARDL.R")

library(httr)
library(lubridate)




meths = c("gmr","mean","zscore")
measures = c("tg","sa","gnssr")

adjust_baseline <- function(data1, data2, method="gmr"){
  nas = is.na(data1)|is.na(data2)
  data1.filt = data1[!nas]
  data2.filt = data2[!nas]
  m=NA
  ml=NA
  mu=NA

  if(method=="zscore"){
    transf1 = (data1.filt-mean(data1.filt))/sd(data1.filt)
    transf2 = (data2.filt-mean(data2.filt))/sd(data2.filt)
  }
  else if(method=="gmr"){
    fit.gmr <- lmodel2(data2.filt~data1.filt)
    m = fit.gmr$regression.results[3,'Slope']
    b = fit.gmr$regression.results[3,'Intercept']

    ml = fit.gmr$confidence.interval[3,4]
    mu = fit.gmr$confidence.interval[3,5]

    transf1 = m*data1.filt+b
    transf2 = data2.filt
  }
  else if(method=="mean"){
    transf1 = data1.filt-mean(data1.filt)
    transf2 = data2.filt-mean(data2.filt)
  }
  else{
    stop("unknown method")
  }
  data1.fin <- data1
  data1.fin[nas]<-NA
  data1.fin[!nas]<-transf1
  data2.fin <- data2
  data2.fin[nas]<-NA
  data2.fin[!nas]<-transf2
  return(list("data1"=data1.fin,"data2"=data2.fin,"slope"=m,"slope.l"=ml,"slope.u"=mu))
}


month_plot_by_station <- function(station, data_dir=DATA_DIR, add_VLM=F, VLM_na_interp=T,equal_time=F){
  get_month_data_raw(station, adjust_by_mean = T, add_VLM = add_VLM, VLM_na_interp=VLM_na_interp, data_dir = data_dir, equal_time=equal_time) %>%
    dplyr::select(-is_na) %>%
    monthly_plot(station, 'mean')
}

plot_tides <- function(station,method, data_dir=DATA_DIR,save_plot=T,only_data=F,add_VLM=F,VLM_na_interp=T, equal_time=F){
  sealev.df <- get_month_data_raw(station, data_dir=data_dir, adjust_by_mean=T,add_VLM=add_VLM, VLM_na_interp = VLM_na_interp,equal_time=equal_time)
  if (method != 'mean'){
    sealev.df$sa=adjust_baseline(sealev.df$sa,sealev.df$tg,method=method)$data1
    sealev.df$gnssr=adjust_baseline(sealev.df$gnssr,sealev.df$tg,method=method)$data1
  }

  if(only_data) return(sealev.df)


  r_tgsa = signif(cor(sealev.df$sa,sealev.df$tg,use="complete.obs"),4)
  r_tggnssr = signif(cor(sealev.df$tg,sealev.df$gnssr,use="complete.obs"),4)
  r_sagnssr = signif(cor(sealev.df$sa,sealev.df$gnssr,use="complete.obs"),4)


  lab_point1 = max(sealev.df$tg[!is.na(sealev.df$tg)])
  lab_point2 = lab_point1*.95+.05*min(sealev.df$tg[!is.na(sealev.df$tg)])

  station.plot=monthly_plot(sealev.df,station,method=method) + annotate("text",x=as.Date("1995-01-01"),y=lab_point1,label=paste("TG v SA cor=",r_tgsa))+
    annotate("text",x=as.Date("1995-01-01"),y=lab_point2,label=paste("TG v GNSSR cor=",r_tggnssr))

  if(save_plot){
    ggsave(filename = paste0(station,"_",method,"_plot.png"),plot = station.plot,path = file.path(data_dir,station),width=16,height = 8 )
  }


  return(station.plot)
}



get_data <- function(station, data_dir=DATA_DIR, add_VLM=FALSE, VLM_na_interp=TRUE,equal_time=FALSE){
  station_data <- get_month_data_raw(station, data_dir, adjust_by_mean = T, add_VLM=add_VLM, VLM_na_interp = VLM_na_interp,equal_time=equal_time)
  station_data
}



get_slope_and_cor <- function(station, add_VLM = F, adjust_by_mean=FALSE) {
  return(regression_analysis_by_station(station, add_VLM = add_VLM, adjust_by_mean=adjust_by_mean))
}

regression_analysis_by_station <- function(station,
                                           data_dir=DATA_DIR,
                                           add_VLM=FALSE,VLM_na_interp=TRUE,
                                           adjust_by_mean=TRUE,
                                           equal_time = FALSE,
                                           time_unit = "month"){

  if(time_unit=="year"){
    station_data <- get_month_data_raw(station, data_dir, adjust_by_mean = adjust_by_mean, add_VLM=add_VLM, VLM_na_interp=VLM_na_interp,equal_time=equal_time)
    sa_year=convert_monthly_to_yearly(station_data,value_col='sa',na.rm=T)
    tg_year=convert_monthly_to_yearly(station_data,value_col='tg',na.rm=T)
    gnssr_year = convert_monthly_to_yearly(station_data,value_col='gnssr',na.rm=T)
    station_data <- merge(sa_year,tg_year,by="year")
    station_data <- merge(station_data,gnssr_year, by="year")
    names(station_data)<- c("year","sa","tg","gnssr")

  }else if(time_unit=="month"){
    station_data <- get_month_data_raw(station, data_dir, adjust_by_mean = adjust_by_mean, add_VLM=add_VLM, VLM_na_interp=VLM_na_interp,equal_time=equal_time)
  }else{
    stop("Unrecognized time units.")
  }

  sa.gmr <- lmodel2(station_data$tg~station_data$sa)
  gnssr.gmr <- lmodel2(station_data$tg~station_data$gnssr)
  sa.gmr.intercept <- sa.gmr$regression.results[3,'Intercept']
  sa.gmr.slope <- sa.gmr$regression.results[3,'Slope']
  sa.ols.intercept <- sa.gmr$regression.results[1,'Intercept']
  sa.ols.slope <- sa.gmr$regression.results[1,'Slope']
  gnssr.gmr.intercept <- gnssr.gmr$regression.results[3,'Intercept']
  gnssr.gmr.slope <- gnssr.gmr$regression.results[3,'Slope']
  sa.gmr.slope.l <- sa.gmr$confidence.interval[3,4]
  sa.gmr.slope.u <- sa.gmr$confidence.interval[3,5]

  average_efficiency_sa <- Obtain_average_efficiency(station_data$tg, station_data$sa)
  average_efficiency_gnssr <- Obtain_average_efficiency(station_data$tg, station_data$gnssr)
  sa.ols.average_efficiency <- average_efficiency_sa$average_efficiency_OLS
  sa.gmr.average_efficiency <- average_efficiency_sa$average_efficiency_GMR
  gnssr.ols.average_efficiency <- average_efficiency_gnssr$average_efficiency_OLS
  gnssr.gmr.average_efficiency <- average_efficiency_gnssr$average_efficiency_GMR

  data.frame(
    "station" = station,
    "VLM.Corrected" = add_VLM,
    "tg.variance" = var(station_data$tg,na.rm=T),
    "sa.variance" = var(station_data$sa,na.rm=T),
    "sa.gmr.intercept" = sa.gmr.intercept,
    "sa.gmr.slope" = sa.gmr.slope,
    "sa.gmr.slope.l" = sa.gmr.slope.l,
    "sa.gmr.slope.u" = sa.gmr.slope.u,
    "sa.ols.intercept" = sa.ols.intercept,
    "sa.ols.slope" = sa.ols.slope,
    "gnssr.gmr.intercept" = gnssr.gmr.intercept,
    "gnssr.gmr.slope" = gnssr.gmr.slope,
    "sa.ols.average_efficiency" = sa.ols.average_efficiency,
    "sa.gmr.average_efficiency" = sa.gmr.average_efficiency,
    "gnssr.ols.average_efficiency" = gnssr.ols.average_efficiency,
    "gnssr.gmr.average_efficiency" = gnssr.gmr.average_efficiency
  )

  }


backtesting_prediction <- function(location, measurement.type = "tg", endyear = 2009, horizon = 2022,na.rm=F, add_VLM=T,...){
  #' Returns a data frame. measurement.type can be 'tg', 'sa' or 'gnssr'.
  #' @param location The name of the station. It must be one of c("astoria", "battery_park", "cape_charles", "charleston", "crescent_city",
  #' "elly_oil_platform", "monterey", "newport", "south_beach", "tofino"
  #' )
  #' @param measurement.type The type of measurement to use. It can be one of "tg", "sa", or "gnssr".
  #' @param add_VLM Whether to add VLM to the data.
  #' @param endyear When to "stop" recording real values, and begin predicting sealevel values based on data up to that point
  #' @param horizon How far out to predict. Default 2022
  #' @param na_rm Whether to remove NA's in the data. Usually set to FALSE to avoid bias in the mean calculation.
  #' Default is FALSE.
  #' @param equal_time Whether to use the same time range for all data.
  #' @param VLM_method The method to use for VLM calculation. See `docstring(get_VLM)` or `?get_VLM` for details.
  #' @param ... Other parameters to pass to `get_month_data_raw`.
  dat = get_month_data_raw(location,add_VLM = add_VLM, ...)
  climate <- read.csv('data/YearData.csv')[,c("Year","GMSL")]

  pred_GMSL <- predict_GMSL(endyear)


  times = dat$month

  if(measurement.type=="tg"){
    sealev = dat$tg
  }else if(measurement.type=="sa"){
    sealev = dat$sa
  }else if(measurement.type=="gnssr"){
    sealev = dat$gnssr
  }else{
    stop("Unknown measurement type")
  }

  years <- as.numeric(substr(times,1,4))
  month_sealevel <- data.frame(time=years,sealevel=sealev)
  annual_sealevel <- month_sealevel %>% group_by(time) %>% summarise(means = mean(sealevel,na.rm=na.rm))
  names(annual_sealevel)=c('Year','NY')
  Year=data.frame(Year=1950:endyear)
  Year_station=merge(Year,annual_sealevel,by='Year')
  Year_station[Year_station<(-500)]=NA
  Year_station$NY[is.nan(Year_station$NY)]=NA
  GMSL <- climate$GMSL

  NYSL <- merge(Year_station, climate, by="Year")



  baseline <- min(which(!is.na(NYSL$NY)))

  NYSL$GMSL <- NYSL$GMSL - NYSL$GMSL[baseline]
  annual_sealevel$NY <- annual_sealevel$NY - NYSL$NY[baseline]
  NYSL$NY <- NYSL$NY - NYSL$NY[baseline]
  NYSL.l <- cbind(NYSL[-1,],NYSL[-nrow(NYSL),])
  names(NYSL.l)[4:6] <- paste0(names(NYSL.l)[4:6], '.l')


  fit_NY <- lm(NY ~ NY.l + GMSL, data = NYSL.l)


  return(fit_NY)

  endindex = which(NYSL.l$Year == endyear)
  timespan = horizon-NYSL.l$Year[endindex]



  sigma_NY <- summary(fit_NY)$sigma
  sigma_GMSL <- pred_GMSL$sigma[1:timespan+1]
  sigma_both <- c(0)

  for(i in 1:timespan){
    sigma_both <- c(sigma_both,
                    sqrt(fit_NY$coefficients['NY.l']^2 * sigma_both[i]^2 + fit_NY$coefficients['GMSL']^2 * sigma_GMSL[i]^2 + sigma_NY^2))

  }
  sigma_both <- sigma_both[-1]




  pred_NY_no <- c(NYSL.l$NY[endindex])
  for(i in 1:timespan){
    pred_NY_no <- c(pred_NY_no, predict(fit_NY, newdata = data.frame(NY.l = pred_NY_no[i], GMSL = pred_GMSL$mean[which(pred_GMSL$Year==endyear)+i])))
  }
  pred_NY_no <- pred_NY_no[-1]

  pred_NY <- data.frame(Year = endyear:horizon,
                        mean = c(NYSL$NY[endindex+1], pred_NY_no),
                        up95 = c(NYSL$NY[endindex+1], pred_NY_no + qnorm(0.975) * sigma_both),
                        low95 = c(NYSL$NY[endindex+1], pred_NY_no - qnorm(0.975) * sigma_both),
                        up99 = c(NYSL$NY[endindex+1], pred_NY_no + qnorm(0.995) * sigma_both),
                        low99 = c(NYSL$NY[endindex+1], pred_NY_no - qnorm(0.995) * sigma_both)
  )

  pred_NY <- merge(pred_NY,annual_sealevel,by="Year",all.x = TRUE)


  return(pred_NY)



}


plot_backtesting <- function(location, measurement.type = "tg", endyear = 2009, horizon = 2022, na.rm=F, add_VLM=T, ...){
  pred_NY <- backtesting_prediction(location, measurement.type = measurement.type, endyear = endyear, horizon = horizon, na.rm=na.rm, add_VLM=add_VLM, ...) %>%
    rename(Observations = NY)
  ggplot(data=pred_NY,aes(x=Year))+
    geom_line(aes(y=Observations, color="Observations"))+
    geom_line(aes(y=mean, color="Mean"))+
    geom_ribbon(aes(ymin=low95,ymax=up95, fill="95% CI"),alpha=0.2)+
    geom_ribbon(aes(ymin=low99,ymax=up99, fill="99% CI"),alpha=0.1)+
    scale_color_manual(name = "Lines", values = c("Observations" = "black", "Mean" = "blue")) +
    scale_fill_manual(name = "Ribbons", values = c("95% CI" = "blue", "99% CI" = "blue")) +
    labs(title=paste(location, "Backtesting Result"),x="Year",y="Relative Sea Level (m)")+
    theme_minimal()
}

predict_RSL_data <- function(data, measurement.type = "tg", steps=1) {
  #' Reimplementation of RSL_predict_data. Returns a data frame. measurement.type can be 'tg', 'sa' or 'gnssr'.
  #' @param data A data frame with the raw data. The data frame must have a column named 'Year' and another column with the name of the measurement type.
  #' Usually, the data is obtained from `get_month_data_raw`.
  #' @param measurement.type The type of measurement to use. It can be one of "tg", "sa", or "gnssr".
  #' @return A data frame with the prediction.
  #'
  #'
  data <- data %>%
    rename(Year = year) %>%
    left_join(climate, by = "Year")
  data <- data %>% filter(!is.na(get(measurement.type)))
  data <- data %>%
    rename(RSL = !! {measurement.type})
  fit_ARDL <- lm(RSL ~ lag(RSL,1) + GMSL, data = data)
  fit_NY <- fit_ARDL
  sigma_NY <- summary(fit_NY)$sigma
  sigma_GMSL_no <- pred_GMSL$sigma_no[(startyear+1):150]
  sigma_GMSL_no[is.na(sigma_GMSL_no)] <- 0
  sigma_GMSL_cop <- (pred_GMSL$up95_cop-pred_GMSL$mean_cop)[(startyear+1):150]/qnorm(0.975)
  sigma_GMSL_cop[is.na(sigma_GMSL_cop)] <- 0
  sigma_both_no <- c(0)
  sigma_both_cop <- c(0)
}

RSL_predict_backtest_data <- function(data, measurement.type, splitratio = 0.8) {
  data <- left_join(data, climate, by = "Year") %>%
    rename(RSL = !! {measurement.type})
  data <- data %>% mutate(lag_RSL = lag(RSL,1))
  total_rows <- nrow(data)
  train_row <- round(total_rows * splitratio)
  train_data <- data[1:train_row,]
  test_data <- data[(train_row):total_rows,]
  Year <- data$Year[(train_row+1):total_rows]
  fit_ARDL <- lm(RSL ~ lag_RSL + GMSL, data = train_data)
  one_step_ahead <- predict(fit_ARDL, newdata = test_data)
  multi_step_ahead <- c(one_step_ahead[1])
  for (row in (train_row+2):(total_rows)) {
    fit_new <- lm(RSL ~ lag_RSL + GMSL, data = data[1:(row-1),])
    multi_step_ahead <- c(multi_step_ahead, predict(fit_new, newdata = data[row,]))
  }
  result <- data.frame(Year = Year,
  RSL = test_data$RSL[1:nrow(test_data)-1],
  one_step_ahead = one_step_ahead[1:nrow(test_data)-1],
  multi_step_ahead = multi_step_ahead
  )
  return(result)


}

RSL_predict_backtest_station <- function(station, measurement.type = "tg", splitratio = 0.8, ...) {
  "Returns a data frame. measurement.type can be 'tg', 'sa' or 'gnssr'.
  Time range is the maximum possible."
  dat = get_month_data_raw(station, ...)
  dat <- dat %>% convert_monthly_to_yearly %>% rename(Year = year)
  RSL_predict_backtest_data(dat, measurement.type = measurement.type, splitratio=splitratio)
}

plot_backtesting_data <- function(data, measurement.type = "tg", savePlot = T, save_location = './results/', save_name = NULL) {
  "Plots the backtesting results."
  p <- ggplot(data, aes(x = Year)) +
    geom_line(aes(y = RSL, color = "Observed")) +
    geom_line(aes(y = one_step_ahead, color = "One step ahead")) +
    geom_line(aes(y = multi_step_ahead, color = "Multi step ahead")) +
    ylab("Sea Level (m)") +
    scale_color_manual(values = c("Observed" = "black", "One step ahead" = "red", "Multi step ahead" = "blue")) +
    theme(legend.position = c(0.2,0.7), legend.title = element_blank(), legend.text = element_text(size = 5), legend.key.size = unit(0.5, "lines"), legend.spacing.y = unit(0, "pt"))
  if (savePlot) {
    if (is.null(save_name)) {
      save_name <- paste0(measurement.type, "_backtesting.png")
    }
    ggsave(filename = save_name, plot = p, path = save_location, width = 2.36 * 3, height = 1.965 * 3, device = 'pdf')
  }
  return(p)
}

plot_backtesting_station <- function(station, measurement.type = "tg", splitratio = 0.8, savePlot = F, save_location = './', save_name = NULL, ...) {
  "Plots the backtesting results."
  dat <- RSL_predict_backtest_station(station, measurement.type = measurement.type, splitratio = splitratio, ...)
  plot_backtesting_data(dat, measurement.type = measurement.type, savePlot = savePlot, save_location = save_location, save_name = save_name)
}

round2 <- function(x){
  tmp <- round(x,2)

  if(tmp-round(x,0) == 0){
    return(paste0(as.character(tmp),'.00'))
  }
  if(tmp-round(x,1) == 0){
    return(paste0(as.character(tmp),'0'))
  }
  return(as.character(tmp))
}

plot_pred <- function(prediction,title,data_dir=DATA_DIR,savePlot=T,
                      width=2.36 * 3,height=1.965*3, device = 'pdf',
                      legend.position = c(0.2,0.70),
                      ci_type = 'no',
                      save_location = './',
                      save_name = NULL){
  "Set ci_type to NA to not include confidence intervals"
  pred=prediction$prediction
  method = prediction$measurement
  addVLM = prediction$withVLM

  colors <- c(SSP585 = '#F28522FF',SSP370 = '#91D1C2FF', SSP460 = '#7E6148FF',SSP245 = '#00A087FF', Without= "#DC0000FF", With = "#3C5488FF", SSP126 = '#8491B4FF',SSP119 = '#4DBBD5FF')
  labels<-c('SSP5-8.5','SSP3-7.0','SSP4-6.0','SSP2-4.5','Unrestricted scenario','COP26 scenario','SSP1-2.6','SSP1-1.9')
  picture_labels<- c('Unrestricted', 'COP26', 'SSP1-1.9', 'SSP1-2.6',  'SSP2-4.5', 'SSP3-7.0', 'SSP4-6.0', 'SSP5-8.5')
  picture_colors<- c(Without= "#DC0000FF", With = "#3C5488FF", SSP119 = '#4DBBD5FF', SSP126 = '#8491B4FF', SSP245 = '#00A087FF', SSP370 = '#91D1C2FF', SSP460 = '#7E6148FF',  SSP585 = '#F28522FF')
  numbersize = 2.5

  plt <- ggplot(pred)+
    geom_line(aes(x=Year,y=mean_no, color = 'Without')) +
    geom_line(aes(x=Year,y=mean_cop, color = 'With')) +
    geom_line(aes(x = Year, y = ssp119, col = 'SSP119')) +
    geom_line(aes(x = Year, y = ssp126, col = 'SSP126')) +
    geom_line(aes(x = Year, y = ssp245, col = 'SSP245')) +
    geom_line(aes(x = Year, y = ssp370, col = 'SSP370')) +
    geom_line(aes(x = Year, y = ssp460, col = 'SSP460')) +
    geom_line(aes(x = Year, y = ssp585, col = 'SSP585'))
  if (is.na(ci_type)) {
  } else if (ci_type == 'no'){
    plt <- plt +
    geom_ribbon(aes(x=Year,ymin = low95_no, ymax = up95_no),  fill="#DC0000FF", alpha = 0.2) +
    geom_ribbon(aes(x=Year,ymin = low99_no, ymax = up99_no),  fill="#DC0000FF", alpha = 0.1)
  } else if (ci_type == 'cop'){
    plt <- plt +
    geom_ribbon(aes(x=Year,ymin = low95_cop, ymax = up95_cop),  fill="#3C5488FF", alpha = 0.2) +
    geom_ribbon(aes(x=Year,ymin = low99_cop, ymax = up99_cop),  fill="#3C5488FF", alpha = 0.1)
  }
  plt <- plt +
    geom_line(aes(x= c(Year[1:12],rep(NA,79)), y = c(mean_no[1:12],rep(NA,79)))) +
    ylab(expression('Sea Level (m)')) +
    scale_color_manual(values = colors,
                       labels = labels) +
    theme(legend.position=legend.position,legend.title=element_blank(),legend.text = element_text(size=5),legend.key.size = unit(0.5, "lines"),legend.spacing.y = unit(0, "pt")) +
    annotate(geom="text", x= 2100, y=pred$mean_no[91], label = round2(pred$mean_no[91]),size = numbersize,  col = "#DC0000FF") +
    annotate(geom="text", x= 2100, y=pred$up95_no[91], label = round2(pred$up95_no[91]),size = numbersize, col = "#E64B35FF")+
    annotate(geom="text", x= 2100, y=pred$low95_no[91], label = round2(pred$low95_no[91]),size = numbersize, col = "#E64B35FF")+
    annotate(geom="text", x= 2100, y=pred$up99_no[91], label = round2(pred$up99_no[91]),size = numbersize, col = "#E64B35FF", alpha = 0.8)+
    annotate(geom="text", x= 2100, y=pred$low99_no[91], label = round2(pred$low99_no[91]),size = numbersize, col = "#E64B35FF", alpha = 0.8)+
    annotate(geom="text", x= 2100, y=pred$mean_cop[91], label = round2(pred$mean_cop[91]),size = numbersize,col = "#3C5488FF")+
    annotate(geom="text", x = 2050, y=pred$mean_no[41], label = round2(pred$mean_no[41]),size = numbersize, col = "#DC0000FF") +
    annotate(geom="text", x = 2050, y=pred$mean_cop[41], label = round2(pred$mean_cop[41]),size = numbersize,col = "#3C5488FF")+
    annotate(geom="text", x = 2050, y=pred$up99_no[41], label = round2(pred$up99_no[41]),size = numbersize, col = "#E64B35FF", alpha = 0.8) +
    annotate(geom="text", x = 2050, y=pred$low99_no[41], label = round2(pred$low99_no[41]),size = numbersize, col = "#E64B35FF", alpha = 0.8) +
    annotate(geom="text", x = 2021, y=pred$mean_cop[12]*1.1, label = round2(pred$mean_cop[12]),size = numbersize) +
    annotate(geom="text", x = 2110, y=c(pred[91,2],as.numeric(pred[91, c(6,12:17)])), label = picture_labels,  col = picture_colors, alpha = 0.8, size = 1.5)

  if (is.null(save_name)){
    save_name = paste0(title,'_',method,'_prediction_plot.',device)
  }
  if(savePlot){
    if(addVLM){
      if (is.null(save_location)){
        save_location = paste0(title,'_',method,'_with_VLM','_prediction_plot.',device)
      }
      ggsave(save_name,plot=plt,path=file.path(data_dir,title),
             width = width, height = height)
    }
    else{

      ggsave(save_name,plot=plt,path=save_location,
             width = width, height = height)
    }
  }

  return(plt)
}

plot_pred_by_station <- function(station, method = "tg", data_dir=DATA_DIR,add_VLM=F, VLM_na_interp=T, na.rm=T,
                                  savePlot=F, device='pdf', height=5, width=8, legend.position = 'right', ci_type = 'no',
                                  save_location = NULL){
  prediction<-RSL_predict(station,method,add_VLM = add_VLM, VLM_na_interp=T, na.rm=na.rm)
  if((!is.nan(prediction$prediction$up95_no[15]))&(!is.na(prediction$prediction$up95_no[15]))){
    pred_plot <- plot_pred(prediction,station, savePlot = T, device=device,
                           height=height, width = width,
                           legend.position = legend.position,
                           ci_type = ci_type,
                           save_location = save_location)
    pred_plot

  }
}

plot_multi_prediction_1_fig <- function(..., station, method = "tg", data_dir=DATA_DIR,add_VLM=F, na.rm=T,
                                        savePlot=F, device='pdf', height=5, width=8, legend.position = 'right', ci_type = 'no',
                                        save_location = NULL){
  "predictions should be a list of predictions from the RSL_predict function"
  predictions <- list(...)
  combined_predictions <- data.frame(Year = 2023)
  for (i in 1:length(predictions)){
    prediction <- predictions[[i]]
    measurement_type <- prediction$measurement
    if (measurement_type == 'tg') {
      if (prediction$withVLM) with_VLM <- "_with_VLM"
      else with_VLM <- "_without_VLM"
    } else {
      with_VLM <- ""
    }
    prediction <- prediction$prediction %>%
      select(Year | contains('mean_no') | contains('mean_cop')) %>%
      rename_with(~paste0(., "_", measurement_type, with_VLM), -Year)
    combined_predictions <- merge(combined_predictions, prediction, by = "Year", all = T)
  }
  combined_plot <- combined_predictions %>%
    pivot_longer(!Year) %>%
    ggplot(aes(x = Year, y = value, color = name)) +
    geom_line(aes(color=name),size=0.75)
  if (savePlot) {
    if (is.null(save_location)){
      save_location = paste0(station,'_',method,'_prediction_combined_plot.',device)
    }
    ggsave(save_location,plot=combined_plot,path=file.path(data_dir,station),
            width = width, height = height)
  }

  return(combined_plot)
}

compare_pred <- function(station, add_VLM = FALSE, VLM_na_interp = TRUE, save_csv = FALSE, data_dir=DATA_DIR){
  meths <- c("tg","sa","gnssr")
  df <- NULL
  for (measure in meths){
    res <- RSL_predict(station,measurement.type=measure, add_VLM = add_VLM)
    prediction <- res$prediction
    prediction$method <- rep(measure, length(prediction[,1]))
    prediction <- cbind(prediction[,length(prediction)],prediction[,-length(prediction)])
    names(prediction)[1] <- "measure"
    df <- rbind(df,prediction[which(prediction$Year %in% c(2050,2075,2100)),1:12])
  }

  df = df %>% add_column(width95_no=(df$up95_no-df$low95_no),.after = "low95_no")
  df = df %>% add_column(width99_no=(df$up99_no-df$low99_no),.after = "low99_no")
  df = df %>% add_column(width95_cop=(df$up95_cop-df$low95_cop),.after = "low95_cop")
  df = df %>% add_column(width99_cop=(df$up99_cop-df$low99_cop),.after = "low99_cop")

  if(save_csv){
    path <- file.path(data_dir,station,"compare_pred.csv")
    write.csv(df,path)
  }

  return(df)
}


group_data_type <- function(input_str){
  if (endsWith(input_str, "gnssr")) {
    return("gnssr")
  } else if (endsWith(input_str, "sa")) {
    return("sa")
  } else if (endsWith(input_str, "tg")) {
    return("tg")
  } else {
    return("no")
  }
}

compare_pred_plot <- function(data, station, title=NULL) {
  if(is.null(title)) {
    title = paste("99% Confidence Bands for Sea Level Predictions \nStation: ", str_to_title(station))
  }
  data %>%
    dplyr::select(Year | contains("99_no") | contains("99_cop")) %>%
    pivot_longer(!Year) %>%
    mutate(type = purrr::map_chr(name, group_data_type)) %>%
    ggplot(aes(x = Year, y = value, color = type)) +


    geom_line(aes(color = type, group = name), size = 0.9) +



    geom_line(aes(), size=0.3, linetype = 6) +

    labs(title=title, x="Time",y="Height (m)",
         color = "Measurement\nMethod")+
    theme(plot.title = element_text(hjust=0.5, size=14, face = "bold"),
          axis.title = element_text(face = "bold"))
}

compare_pred_plot_by_data <- function(pred_tg,
                                      pred_sa,
                                      pred_gnssr,
                                      station = "",
                                      adjust_by_mean = T,
                                      data_dir=DATA_DIR,
                                      add_VLM = T,
                                      VLM_na_interp = TRUE,
                                      na.rm = T,
                                      data_type = "no",
                                      return_data_only = F,
                                      title=NULL) {
  if (is.null(pred_gnssr)) {
    pred_gnssr <- NULL
  } else {
    pred_gnssr <- pred_gnssr$prediction %>%
      dplyr::select(Year | contains(data_type)) %>%
      rename_with(~paste0(., "_gnssr"), -Year)
  }
  pred_sa <- pred_sa$prediction %>%
    dplyr::select(Year | contains(data_type)) %>%
    rename_with(~paste0(., "_sa"), -Year)
  pred_tg <- pred_tg$prediction %>%
    dplyr::select(Year | contains(data_type)) %>%
    rename_with(~paste0(., "_tg"), -Year)

  if (is.null(pred_gnssr)) {
    pred <- pred_sa %>% merge(pred_tg, by = "Year", all = T)
  } else {
    pred <- merge(pred_gnssr, pred_sa, by = "Year", all = T) %>% merge(pred_tg, by = "Year", all = T)
  }
  if (return_data_only) {
    return(pred)
  }

  compare_pred_plot(pred, station, title=title)
}


compare_pred_plot_by_station <- function(station,
                                         adjust_by_mean = T,
                                         data_dir=DATA_DIR,
                                         add_VLM = T,
                                         VLM_na_interp = TRUE,
                                         na.rm = T,
                                         data_type = "no",
                                         return_data_only = F) {
  pred_gnssr <- RSL_predict(station,measurement.type="gnssr", add_VLM = add_VLM, VLM_na_interp = VLM_na_interp, na.rm = na.rm)
  Sys.sleep(1)
  pred_sa <- RSL_predict(station,measurement.type="sa", add_VLM = add_VLM, VLM_na_interp = VLM_na_interp, na.rm = na.rm)
  Sys.sleep(1)
  pred_tg <- RSL_predict(station,measurement.type="tg", add_VLM = add_VLM, VLM_na_interp = VLM_na_interp, na.rm = na.rm)

  compare_pred_plot_by_data(pred_tg, pred_sa, pred_gnssr, station = station, adjust_by_mean = adjust_by_mean, data_dir = data_dir, add_VLM = add_VLM, VLM_na_interp = VLM_na_interp, na.rm = na.rm, data_type = data_type, return_data_only = return_data_only)
}


station_variable_selection <- function(station, data_dir=DATA_DIR){

  station_data = get_month_data_raw(station,data_dir=data_dir)

  New = convert_monthly_to_yearly(station_data, value_col="tg")


  names(New)[2] = "sealevel"


  New$sealevel=New$sealevel*1000
  up_rows=data.frame(X=1950:1992,x=rep(NA,length(1950:1992)))
  names(New)=c('X','x')
  NYMSL=rbind(up_rows,New[which(New$X>1992 & New$X<=2021),])
  GMSL <- climate$GMSL


  Mass <- climate$Mass
  SeaIce <- climate$SeaIce
  Temp <- climate$TEMP
  Humidity <- climate$humidity
  SSN <- climate$SSN
  CH4 <- climate$CH4
  N2O <- climate$N2O
  CO2 <- climate$CO2
  GWP <- climate$GWP


  NYSL <- data.frame(Year = 1950:2021, NY = NYMSL$x, GMSL = GMSL, Mass = Mass, SeaIce = SeaIce,
                     Temp = Temp, Humidity = Humidity, SSN = SSN, CH4 = CH4, N2O = N2O, CO2 = CO2, GWP = GWP)


  NYSL$GMSL <- NYSL$GMSL - mean(NYSL$GMSL[37:56],na.rm = TRUE)
  NYSL$NY <- NYSL$NY - mean(NYSL$NY[37:56],na.rm = TRUE)
  NYSL.l <- cbind(NYSL[-1,],NYSL[-72,])
  names(NYSL.l)[13:15] <- paste0(names(NYSL.l)[13:15], '.l')
  NYSL.l <- NYSL.l[,-(16:24)]

  NYSL.l <- NYSL.l[,-1]
  NYSL.l <- NYSL.l[,-12]


  fit_NY <- lm(NY ~ ., data = na.omit(NYSL.l))
  summary(fit_NY)
  stepAIC(fit_NY, direction = 'both', trace = F)

}

plots_compare_measuremtns <- function(station, title = waiver(), measures = c("tg", "sa", "gnssr"), ...) {
  #' Plots the predictions for the different measurement types
  #' @param station The name of the station.
  #' @param measures The measurement types to be compared. Default is c("tg", "sa", "gnssr").
  #' @param ... Additional arguments to be passed to `RSL_predict`.
  arg_list = list(...)
  predictions = c()
  for (meth in measures) {
    tryCatch(
      {
        prediction <- RSL_predict(station, meth, ...)
        predictions <- c(predictions, list(prediction))
      }, error = function(e) {
        print(e)
        predictions <- c(predictions, list(NULL))
      }
    )
  }
  if (length(predictions) == 2) {
    compare_pred_plot_by_data(predictions[[1]], predictions[[2]], NULL, station = station, title = title)
  } else if (length(predictions) == 3) {
    compare_pred_plot_by_data(predictions[[1]], predictions[[2]], predictions[[3]], station = station, title = title)
  }
}











regression_analysis_by_station_for_fit_plot <- function(station,
                                                        data_dir=DATA_DIR,
                                                        add_VLM=FALSE, VLM_na_interp=TRUE,
                                                        adjust_by_mean=TRUE,
                                                        equal_time = TRUE,
                                                        time_unit = "month",
                                                        VLM_method = 0,
                                                        VLM_retreive_method = 'deterministic',
                                                        ...){
  #' Returns a data frame.
  #' @examples
  #' regression_analysis_by_station_for_fit_plot("tofino", add_VLM = TRUE, VLM_na_interp = FALSE, equal_time = NULL, VLM_method = 1)
  #' regression_analysis_by_station_for_fit_plot("tofino", add_VLM = TRUE, VLM_na_interp = FALSE, equal_time = NULL, VLM_method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = "mean")
  #' regression_analysis_by_station_for_fit_plot("tofino", add_VLM = TRUE, VLM_na_interp = FALSE, equal_time = NULL, VLM_method = 'ULR', VLM_retreive_method = "longest", choose_which = 5)
  station_data <- get_month_data_raw(station, data_dir, adjust_by_mean = adjust_by_mean, add_VLM=add_VLM, VLM_na_interp=VLM_na_interp,equal_time=equal_time, VLM_method = VLM_method, VLM_retreive_method = VLM_retreive_method, ...)
  if(time_unit=="year"){
    sa_year=convert_monthly_to_yearly(station_data,value_col='sa',na.rm=T)
    tg_year=convert_monthly_to_yearly(station_data,value_col='tg',na.rm=T)
    gnssr_year = convert_monthly_to_yearly(station_data,value_col='gnssr',na.rm=T)
    station_data <- merge(sa_year,tg_year,by="year")
    station_data <- merge(station_data,gnssr_year, by="year")
    names(station_data)<- c("year","sa","tg","gnssr")
    year_start = station_data$year[1]
    year_end = station_data$year[nrow(station_data)]
  }else if(time_unit=="month"){
    year_start = year(station_data$month[1])
    year_end = year(station_data$month[nrow(station_data)])
  }else{
    stop("Unrecognized time units.")
  }

  VLM_source <- do.call(paste, as.list(VLM_method))
  if (VLM_retreive_method == 'longest') {
    VLM_source <- get_VLM_by_station(station, method = VLM_method, VLM_retreive_method = VLM_retreive_method, return_VLM_station_info_only = TRUE, ...)
    VLM_source <- tail(VLM_source$site, n=1)
  }

  sa.gmr <- lmodel2(station_data$tg~station_data$sa)
  sa.gmr.intercept <- sa.gmr$regression.results[3,'Intercept']
  sa.gmr.slope <- sa.gmr$regression.results[3,'Slope']
  sa.ols.intercept <- sa.gmr$regression.results[1,'Intercept']
  sa.ols.slope <- sa.gmr$regression.results[1,'Slope']
  sa.ols.slope.l <- sa.gmr$confidence.interval[1,4]
  sa.ols.slope.u <- sa.gmr$confidence.interval[1,5]
  sa.gmr.slope.l <- sa.gmr$confidence.interval[3,4]
  sa.gmr.slope.u <- sa.gmr$confidence.interval[3,5]

  tg.gmr <- lmodel2(station_data$sa~station_data$tg)
  tg.ols.intercept <- tg.gmr$regression.results[1,'Intercept']
  tg.ols.slope <- tg.gmr$regression.results[1,'Slope']
  tg.ols.slope.l <- tg.gmr$confidence.interval[1,4]
  tg.ols.slope.u.tmp <- tg.gmr$confidence.interval[1,5]
  tg.ols.intercept <- -1*(tg.ols.intercept/tg.ols.slope)
  tg.ols.slope <- 1/tg.ols.slope
  tg.ols.slope.u <- 1/tg.ols.slope.l
  tg.ols.slope.l <- 1/tg.ols.slope.u.tmp

  sa.or.intercept <- sa.gmr$regression.results[2, 'Intercept']
  sa.or.slope <- sa.gmr$regression.results[2, 'Slope']
  sa.or.slope.l <- sa.gmr$confidence.interval[2,4]
  sa.or.slope.u <- sa.gmr$confidence.interval[2,5]


  average_efficiency_sa <- Obtain_average_efficiency(station_data$sa, station_data$tg)
  sa.gmr.average_efficiency <- average_efficiency_sa$average_efficiency_GMR

  data.frame(
    "station" = station,
    "tg.variance" = var(station_data$tg,na.rm=T),
    "sa.variance" = var(station_data$sa,na.rm=T),
    "sa.gmr.intercept" = sa.gmr.intercept,
    "sa.gmr.slope" = sa.gmr.slope,
    "sa.gmr.slope.l" = sa.gmr.slope.l,
    "sa.gmr.slope.u" = sa.gmr.slope.u,
    "sa.ols.intercept" = sa.ols.intercept,
    "sa.ols.slope" = sa.ols.slope,
    "sa.ols.slope.l" = sa.ols.slope.l,
    "sa.ols.slope.u" = sa.ols.slope.u,
    "sa.gmr.average_efficiency" = sa.gmr.average_efficiency,

    "tg.ols.intercept" = tg.ols.intercept,
    "tg.ols.slope" = tg.ols.slope,
    "tg.ols.slope.l" = tg.ols.slope.l,
    "tg.ols.slope.u" = tg.ols.slope.u,


    "sa.or.intercept" = sa.or.intercept,
    "sa.or.slope" = sa.or.slope,
    "sa.or.slope.l" = sa.or.slope.l,
    "sa.or.slope.u" = sa.or.slope.u,

    "year_start" = year_start,
    "year_end" = year_end,

    "VLM_source" = VLM_source,
    "VLM_retreive_method" = VLM_retreive_method
  )
}


plot_calibration_fitted_lines_df <- function(res1, add_VLM = TRUE, VLM_method = 0, station = '') {

  year_start <- res1$year_start
  year_end <- res1$year_end

  VLM_source = res1$VLM_source

  slopes <- c(rep(res1$sa.gmr.slope,5),rep(res1$sa.ols.slope,5),rep(res1$sa.or.slope,5),rep(res1$tg.ols.slope,5))
  intercepts <- c(rep(res1$sa.gmr.intercept,5),rep(res1$sa.ols.intercept,5),rep(res1$sa.or.intercept,5),rep(res1$tg.ols.intercept,5))
  fit.method <- c(rep('GMR',5), rep('OLS: TG on SA', 5), rep('OR', 5), rep('OLS: SA on TG', 5))
  x <- c(seq(1:5),seq(1:5),seq(1:5),seq(1:5))
  y <- c(seq(1:20))
  df1 <- data.frame(cbind(slopes,intercepts,fit.method,x,y))
  df1$x <- as.numeric(df1$x)
  df1$y <- as.numeric(df1$y)
  df1$slopes <- as.numeric(df1$slopes)
  df1$intercepts <- as.numeric(df1$intercepts)

  ggplot(data = df1, aes(x = x, y = y)) +
    geom_blank() +
    scale_y_continuous(limits=c(-0.05,0.05)) +
    scale_x_continuous(limits=c(-0.05,0.05)) +
    geom_abline(data = df1, aes(slope = slopes, intercept = intercepts, color = fit.method), lwd = 1.3) +
    geom_abline(color = 'black', lwd = 1, linetype = "longdash") +
    labs(title = paste("Station: ", str_to_title(station)),
         subtitle = paste("VLM data source: ", VLM_source,
                          '\nData range: ', year_start, '-', year_end),
         caption = "Dashed black line is reference line Y = X\nSA has been corrected for VLM",
         color = "Fitting Method") +
    theme(
      plot.title = element_text(size=14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size=12, hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, face = "bold")
    )

}

plot_calibration_fitted_lines <- function(station, add_VLM = TRUE, VLM_method = 0, time_unit = "year", ...) {
  res1 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = VLM_method, ...)
  plot_calibration_fitted_lines_df(res1, add_VLM = add_VLM, VLM_method = VLM_method, station = station)
}

regression_analysis_with_plots_station <- function(station, VLM_na_interp = FALSE,
                                                   save_dir = NULL, save_plots = TRUE,
                                                   width = 8, height = 8, time_unit = "year") {
  if (is.null(save_dir)) {
    save_dir = file.path('results', station)
  }

  res0 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = FALSE, time_unit = time_unit, VLM_method = 'No', VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res0)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_no_VLM.png')), plot = p, width = width, height = height)
  }

  res1 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'ULR', VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res1)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_ULR.png')), plot = p, width = width, height = height)
  }

  res2 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'NGL', VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res2)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_NGL.png')), plot = p, width = width, height = height)
  }

  res3 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'JPL', VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res3)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_JPL.png')), plot = p, width = width, height = height)
  }

  res4 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'GFZ', VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res4)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_GFZ.png')), plot = p, width = width, height = height)
  }

  res5 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'ULR', VLM_retreive_method = 'longest', choose_which = 5, VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res5)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_ULR_longest.png')), plot = p, width = width, height = height)
  }

  res6 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'NGL', VLM_retreive_method = 'longest', choose_which = 5, VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res6)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_NGL_longest.png')), plot = p, width = width, height = height)
  }

  res7 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'JPL', VLM_retreive_method = 'longest', choose_which = 5, VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res7)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_JPL_longest.png')), plot = p, width = width, height = height)
  }

  res8 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = 'GFZ', VLM_retreive_method = 'longest', choose_which = 5, VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res8)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_GFZ_longest.png')), plot = p, width = width, height = height)
  }

  res9 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = 'mean', VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res9)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_mean.png')), plot = p, width = width, height = height)
  }

  res10 <- regression_analysis_by_station_for_fit_plot(station, add_VLM = T, time_unit = time_unit, VLM_method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = 'circlemean', choose_which = 5, max_distance = 500, VLM_na_interp = VLM_na_interp)
  p <- plot_calibration_fitted_lines_df(res10)
  if (save_plots) {
    ggsave(file.path(save_dir, paste0(station, '_calibration_fitted_lines_circlemean.png')), plot = p, width = width, height = height)
  }

  resutls_df <- rbind(res0, res1, res2, res3, res4, res5, res6, res7, res8, res9, res10)

  resutls_df

}


