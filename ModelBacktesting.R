source("ComparePlots.R")
source("models/predict_GMSL.R")
source("XVARMA.R")
source("SEM_data_preprocessing.R")
source("SEM_Fitting.R")

backtesting_ARDL <- function(location, measurement.type = "tg", start.year = 1950, stop.year = 2019, horizon = 2100,na.rm=T, add_VLM=F,...){
  #' Returns a data frame. measurement.type can be 'tg', 'sa' or 'gnssr'.
  #' @param location The name of the station. It must be one of c("astoria", "battery_park", "cape_charles", "charleston", "crescent_city",
  #' "elly_oil_platform", "monterey", "newport", "south_beach", "tofino"
  #' )
  #' @param measurement.type The type of measurement to use. It can be one of "tg", "sa", or "gnssr".
  #' @param add_VLM Whether to add VLM to the data.
  #' @param stop.year When to "stop" recording real values, and begin predicting sealevel values based on data up to that point
  #' @param horizon How far out to predict. Default 2022
  #' @param na_rm Whether to remove NA's in the data. Usually set to FALSE to avoid bias in the mean calculation.
  #' Default is FALSE.
  #' @param equal_time Whether to use the same time range for all data.
  #' @param VLM_method The method to use for VLM calculation. See `docstring(get_VLM)` or `?get_VLM` for details.
  #' @param ... Other parameters to pass to `get_month_data_raw`.
  dat = get_month_data_raw(location,add_VLM = add_VLM, ...)
  climate <- read.csv('data/YearData.csv')[,c("Year","GMSL")]



  pred_GMSL <- predict_GMSL(stop.year)


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
  Year=data.frame(Year=start.year:stop.year)
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



  endindex = which(NYSL.l$Year == stop.year)
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
    pred_NY_no <- c(pred_NY_no, predict(fit_NY, newdata = data.frame(NY.l = pred_NY_no[i], GMSL = pred_GMSL$mean[which(pred_GMSL$Year==stop.year)+i])))
  }
  pred_NY_no <- pred_NY_no[-1]

  pred_NY <- data.frame(Year = stop.year:horizon,
                        LSL_mean = c(NYSL$NY[endindex+1], pred_NY_no),

                        low_95 = c(NYSL$NY[endindex+1], pred_NY_no - qnorm(0.975) * sigma_both),
                        up_95 = c(NYSL$NY[endindex+1], pred_NY_no + qnorm(0.975) * sigma_both),
                        low_99 = c(NYSL$NY[endindex+1], pred_NY_no - qnorm(0.995) * sigma_both),
                        up_99 = c(NYSL$NY[endindex+1], pred_NY_no + qnorm(0.995) * sigma_both)

  )

  pred_NY <- merge(annual_sealevel,pred_NY,by="Year",all.y = TRUE)

  names(pred_NY)[2] = "LSL_true"
  baseline = pred_NY$LSL_true[1]
  pred_NY = pred_NY - baseline
  pred_NY$Year = pred_NY$Year + baseline

  return(pred_NY)



}

backtesting_SEM <- function(location,measurement = "tg", max.p = 2, start.year = -Inf, stop.year = 2019, horizon = 2100,data_dir = DATA_DIR, restriction = "no",...){
  local_data <- preprocess_SEM_data(location,start.year = start.year,stop.year = Inf, measurement.type = measurement,...)
  recovery_data = read.csv(file.path(data_dir,location,"era5data/recovery_data.csv"),header = T)
  true_LSL <- local_data[,c("Year","LSL")]
  true_LSL$LSL <- true_LSL$LSL*recovery_data$LSL_SD+recovery_data$LSL_Mean
  prediction <- SEM_prediction(location, measurement = measurement, start.year = start.year, stop.year = stop.year,restriction = restriction, horizon = horizon,...)

  out <- merge(true_LSL,prediction,by = "Year", all.y = T)
  baseline <- out$LSL[1]
  out <- out[,which(names(out)!="LSL_var")]-baseline
  out$Year <- out$Year + baseline

  names(out)[2] <- "LSL_true"

  return(out)
}

backtesting_VARX <- function(location,measurement = "tg", max.p = 1, start.year = -Inf, stop.year = 2019, horizon = 2100,data_dir = DATA_DIR, restriction = "no",...){
  local_data <- preprocess_SEM_data(location,start.year = start.year,stop.year = Inf, measurement.type = measurement,...)
  recovery_data = read.csv(file.path(data_dir,location,"era5data/recovery_data.csv"),header = T)
  true_LSL <- local_data[,c("Year","LSL")]
  true_LSL$LSL <- true_LSL$LSL*recovery_data$LSL_SD+recovery_data$LSL_Mean
  prediction <- VARX_prediction(location, measurement = measurement, start.year = start.year, stop.year = stop.year,restriction = restriction, horizon = horizon,...)


  out <- merge(true_LSL,prediction,by = "Year", all.y = T)
  baseline <- out$LSL[1]
  out <- out[,which(names(out)!="LSL_var")]-baseline
  out$Year <- out$Year + baseline

  names(out)[2] <- "LSL_true"


  return(out)
}

get_prediction_accuracy <- function(location,stop.year = 2010, horizon = 2021,...){
  bt_ARDL = backtesting_ARDL(location,stop.year = stop.year, horizon = horizon)
  bt_ARDL = bt_ARDL[2:length(bt_ARDL[,1]),]
  n_succ_ARDL95 = sum((bt_ARDL$low_95 <= bt_ARDL$LSL_true & bt_ARDL$LSL_true <= bt_ARDL$up_95))
  ARDL_95acc = n_succ_ARDL95/length(bt_ARDL[,1])
  ARDL_95pval = pbinom(n_succ_ARDL95,length(bt_ARDL[,1]),0.95)
  ARDL_mse = sum((bt_ARDL$LSL_mean-bt_ARDL$LSL_true)^2)/length(bt_ARDL[,1])
  ARDL_results = data.frame("Meth"="ARDL","mse"=ARDL_mse,"acc.95"=ARDL_95acc,"p.value.95"=ARDL_95pval)
  ARDL_misses = bt_ARDL$Year[which(!(bt_ARDL$low_95 <= bt_ARDL$LSL_true & bt_ARDL$LSL_true <= bt_ARDL$up_95))]


  bt_SEM = backtesting_SEM(location,stop.year = stop.year, horizon = horizon)
  bt_SEM = bt_SEM[2:length(bt_SEM[,1]),]
  n_succ_SEM95 = sum((bt_SEM$low_95 <= bt_SEM$LSL_true & bt_SEM$LSL_true <= bt_SEM$up_95))
  SEM_95acc = n_succ_SEM95/length(bt_SEM[,1])
  SEM_95pval = pbinom(n_succ_SEM95,length(bt_SEM[,1]),0.95)
  SEM_mse = sum((bt_SEM$LSL_mean-bt_SEM$LSL_true)^2)/length(bt_SEM[,1])
  SEM_results = data.frame("Meth"="SEM","mse"=SEM_mse,"acc.95"=SEM_95acc,"p.value.95"=SEM_95pval)
  SEM_misses = bt_SEM$Year[which(!(bt_SEM$low_95 <= bt_SEM$LSL_true & bt_SEM$LSL_true <= bt_SEM$up_95))]



  bt_VARX = backtesting_VARX(location,stop.year = stop.year, horizon = horizon,...)
  bt_VARX = bt_VARX[2:length(bt_VARX[,1]),]
  n_succ_VARX95 = sum((bt_VARX$low_95 <= bt_VARX$LSL_true & bt_VARX$LSL_true <= bt_VARX$up_95))
  VARX_95acc = n_succ_VARX95/length(bt_VARX[,1])
  VARX_95pval = pbinom(n_succ_VARX95,length(bt_VARX[,1]),0.95)
  VARX_mse = sum((bt_VARX$LSL_mean-bt_VARX$LSL_true)^2)/length(bt_VARX[,1])
  VARX_results = data.frame("Meth"="VARX","mse"=VARX_mse,"acc.95"=VARX_95acc,"p.value.95"=VARX_95pval)
  VARX_misses = bt_VARX$Year[which(!(bt_VARX$low_95 <= bt_VARX$LSL_true & bt_VARX$LSL_true <= bt_VARX$up_95))]

  bt_ensam = data.frame(Year = bt_VARX$Year, low_95 = (bt_VARX$low_95+bt_ARDL$low_95)/2,up_95 = (bt_VARX$up_95+bt_ARDL$up_95)/2,LSL_mean = (bt_VARX$LSL_mean+bt_ARDL$LSL_mean)/2)
  n_succ_ensam95 = sum((bt_ensam$low_95 <= bt_VARX$LSL_true & bt_VARX$LSL_true <= bt_ensam$up_95))
  ensam_95acc = n_succ_ensam95/length(bt_ensam[,1])
  ensam_95pval = pbinom(n_succ_ensam95,length(bt_ensam[,1]),0.95)
  ensam_mse = sum((bt_ensam$LSL_mean-bt_VARX$LSL_true)^2)/length(bt_ensam[,1])
  ensam_results = data.frame("Meth"="Ensam","mse"=ensam_mse,"acc.95"=ensam_95acc,"p.value.95"=ensam_95pval)
  ensam_misses = bt_ensam$Year[which(!(bt_ensam$low_95 <= bt_VARX$LSL_true & bt_VARX$LSL_true <= bt_ensam$up_95))]


  results = rbind(ARDL_results,SEM_results, VARX_results,ensam_results)
  misses = list("ARDL"=ARDL_misses, "SEM"=SEM_misses, "VARX"=VARX_misses,"Ensam"=ensam_misses)

  return(list("acc"=results,"misses"=misses,"first.pred.year"=stop.year+1,"last.pred.year"=horizon))

}


var_selection_backtesting <- function(location_list, selected_vars = c("X10m_u_component_of_wind","X10m_v_component_of_wind","X2m_temperature","surface_pressure","X2m_dewpoint_temperature","total_precipitation"),
                                      stop.year = 2000,horizon = 2020,...){
  n_succ_VARX95 = 0
  total = 0
  accs = c()


  for(location in location_list){
    bt_VARX = backtesting_VARX(location,stop.year = stop.year, horizon = horizon,selected_vars=selected_vars,...)
    bt_VARX = bt_VARX[2:length(bt_VARX[,1]),]
    this_succ = sum((bt_VARX$low_95 <= bt_VARX$LSL_true & bt_VARX$LSL_true <= bt_VARX$up_95))
    n_succ_VARX95 = n_succ_VARX95 + this_succ
    this_total = length(bt_VARX[,1])
    total = total + this_total
    accs = c(accs,this_succ)

  }
  return(list("vars"= selected_vars,"acc" = n_succ_VARX95/total, "p.value"=pbinom(n_succ_VARX95,total,0.95)))

}

plot_backtesting_results <- function(data, title = NULL) {
  #' Plot the backtesting results with confidence intervals
  #' @param data data frame with columns: Year, LSL_true, LSL_mean, low_95, up_95, low_99, up_99, LSL_fitted (optional)
  #' @param title Optional title for the plot
  #' @return a ggplot object

  required_columns <- c("Year", "LSL_true", "LSL_mean", "low_95", "up_95", "low_99", "up_99")
  if (!all(required_columns %in% colnames(data))) {
    stop("Data must contain the following columns: Year, LSL_true, LSL_mean, low_95, up_95, low_99, up_99")
  }

  ci_data <- data %>%
    dplyr::select(Year, LSL_mean, low_95, up_95, low_99, up_99)

  p <- ggplot() +
    geom_ribbon(data = ci_data, aes(x = Year, ymin = low_99, ymax = up_99), fill = "grey", alpha = 0.3) +
    geom_ribbon(data = ci_data, aes(x = Year, ymin = low_95, ymax = up_95), fill = "grey", alpha = 0.5) +
    geom_line(data = data, aes(x = Year, y = LSL_true, color = "Observation"), size = 1.2, linetype = "solid") +
    geom_line(data = data, aes(x = Year, y = LSL_mean, color = "Prediction"), size = 1.2, linetype = "solid") +
    scale_color_manual(values = c("Observation" = "black", "Prediction" = "blue")) +
    theme_minimal() +
    labs(x = "Year",
         y = "RSL (m)",
         color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank())

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  if ("LSL_fitted" %in% colnames(data)) {
    p <- p + geom_line(data = data, aes(x = Year, y = LSL_fitted, color = "Fitted"), size = 1.2, linetype = "solid") +
      scale_color_manual(values = c("Observation" = "black", "Prediction" = "blue", "Fitted" = "green"))
  }

  return(p)
}

plot_combined_backtesting_results <- function(station, stop.year = 2010, horizon = 2023, start.year = 1970, measurement.type = "tg", use_VLM_as_variable = TRUE, equal_time = FALSE) {
  ardl_results <- backtesting_ARDL(
    station,
    measurement.type = measurement.type,
    use_VLM_as_variable = use_VLM_as_variable,
    stop.year = stop.year,
    horizon = horizon,
    start.year = start.year,
    equal_time = equal_time
  )

  varx_results <- backtesting_VARX(
    station,
    stop.year = stop.year,
    horizon = horizon,
    start.year = start.year
  )

  required_columns <- c("Year", "LSL_true", "LSL_mean", "low_95", "up_95", "low_99", "up_99")
  if (!all(required_columns %in% colnames(ardl_results)) || !all(required_columns %in% colnames(varx_results))) {
    stop("Data must contain the following columns: Year, LSL_true, LSL_mean, low_95, up_95, low_99, up_99")
  }
  station_name <- get_station_name(station)
  ggplot() +
    geom_ribbon(data = ardl_results, aes(x = Year, ymin = low_99, ymax = up_99), fill = "lightblue", alpha = 0.2, color = "skyblue", linetype = "dotted") +
    geom_ribbon(data = ardl_results, aes(x = Year, ymin = low_95, ymax = up_95), fill = "lightblue", alpha = 0.3) +
    geom_ribbon(data = varx_results, aes(x = Year, ymin = low_99, ymax = up_99), fill = "lightgreen", alpha = 0.2, color = "lightgreen", linetype = "dotted") +
    geom_ribbon(data = varx_results, aes(x = Year, ymin = low_95, ymax = up_95), fill = "lightgreen", alpha = 0.3) +
    geom_line(data = varx_results, aes(x = Year, y = LSL_mean, color = "dSEM"), size = 1, linetype = "solid") +
    geom_line(data = ardl_results, aes(x = Year, y = LSL_mean, color = "ARDL"), size = 1, linetype = "solid") +
    geom_line(data = ardl_results, aes(x = Year, y = LSL_true, color = "Observation"), size = 1.2, linetype = "dashed") +
    scale_color_manual(values = c("Observation" = "black", "ARDL" = "blue", "dSEM" = "green"),
                       breaks = c("Observation", "ARDL", "dSEM")) +
    theme_minimal() +
    labs(title = paste(format_station_name(station)),
         x = "Year",
         y = "RSL (m)",
         color = "Legend") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "longdash"))))
}

plot_grid_backtesting_results <- function(station_list=station.names[1:9],
                                          stop.year = 2010, horizon = 2023, start.year = 1970,
                                          measurement.type = "tg", use_VLM_as_variable = TRUE, equal_time = FALSE) {
  #' Create a 3x3 grid of backtesting plots for a list of stations
  #' @param station_list List of stations to plot
  #' @return A combined plot with 3x3 layout

  if(length(station_list) != 9) {
    stop("Please provide exactly 9 stations.")
  }

  plot_list <- lapply(station_list, function(station) {
    plot_combined_backtesting_results(station, stop.year, horizon, start.year,
                                      measurement.type, use_VLM_as_variable, equal_time)
  })

  extract_legend <- function(a_gplot) {
    g <- ggplotGrob(a_gplot)
    legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
    return(legend)
  }

  legend <- extract_legend(plot_list[[1]])

  plot_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))

  grid_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 3
  )

  final_plot <- cowplot::plot_grid(
    grid_plot,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  final_plot
}

plot_grid_long_backtesting_results <- function(station_list = station.names[1:9],
                                             start_year = 1900,
                                             end_year = 1950,
                                             horizon = 2023,
                                             measurement.type = "tg",
                                             use_VLM_as_variable = TRUE) {
  #' Create a 3x3 grid of long-term backtesting plots for a list of stations
  #' @param station_list List of stations to plot
  #' @param start_year Start year for the backtesting period
  #' @param end_year End year for the model training period
  #' @param horizon End year for predictions
  #' @return A combined plot with 3x3 layout
  #' @examples
  #' For single plot, use `backtesting_ARDL_v2(station, start_year=1900, end_year=1950, horizon=2023) %>% plot_backtesting_results`

  if(length(station_list) != 9) {
    stop("Please provide exactly 9 stations.")
  }

  plot_list <- lapply(station_list, function(station) {
    results <- backtesting_ARDL_v2(
      station,
      start_year = start_year,
      end_year = end_year,
      horizon = horizon,
      measurement.type = measurement.type,
      use_VLM_as_variable = use_VLM_as_variable
    )

    station_name <- get_station_name(station)
    plot_backtesting_results(results, title = format_station_name(station_name))
  })

  extract_legend <- function(a_gplot) {
    g <- ggplotGrob(a_gplot)
    legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
    return(legend)
  }

  shared_legend <- extract_legend(plot_list[[1]])

  plot_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))

  grid_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 3
  )

  final_plot <- cowplot::plot_grid(
    grid_plot,
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  return(final_plot)
}