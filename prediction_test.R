RSL_predict_test<-function(location, measurement.type = "tg"){
  dat = get_month_data_raw(location)
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
  annual_sealevel <- month_sealevel %>% group_by(time) %>% summarise(means = mean(sealevel,na.rm=T))
  names(annual_sealevel)=c('Year','sealevel')
  Year=data.frame(Year=1950:2021)
  Year_station=merge(Year,annual_sealevel,by='Year',all.x=TRUE)
  Year_station[Year_station<(-500)]=NA
  Year_station$sealevel[is.nan(Year_station$sealevel)]=NA
  GMSL <- climate$GMSL
  NYSL <- data.frame(Year = 1950:2021, NY = Year_station$sealevel, GMSL = GMSL)

  baseline <- min(which(!is.na(NYSL$NY)))

  NYSL$GMSL <- NYSL$GMSL - NYSL$GMSL[baseline]
  NYSL$NY <- NYSL$NY - NYSL$NY[baseline]
  NYSL.l <- cbind(NYSL[-1,],NYSL[-72,])
  names(NYSL.l)[4:6] <- paste0(names(NYSL.l)[4:6], '.l')



  fit_NY <- lm(NY ~ NY.l + GMSL, data = NYSL.l)
  adf.test(fit_NY$residuals, k = 2)

  startyear = max(which(!is.na(NYSL.l$NY)))
  timespan = 2100-NYSL.l$Year[startyear]


  sigma_NY <- summary(fit_NY)$sigma
  sigma_GMSL_no <- pred_GMSL$sigma_no[(startyear+1):150]
  sigma_GMSL_cop <- (pred_GMSL$up95_cop-pred_GMSL$mean_cop)[(startyear+1):150]/qnorm(0.975)
  sigma_both_no <- c(0)
  sigma_both_cop <- c(0)


  for(i in 1:timespan){
    sigma_both_no <- c(sigma_both_no,
                       sqrt(fit_NY$coefficients['NY.l']^2 * sigma_both_no[i]^2 + fit_NY$coefficients['GMSL']^2 * sigma_GMSL_no[i]^2 + sigma_NY^2))
    sigma_both_cop <- c(sigma_both_cop,
                        sqrt(fit_NY$coefficients['NY.l']^2 * sigma_both_cop[i]^2 + fit_NY$coefficients['GMSL']^2 * sigma_GMSL_cop[i]^2 + sigma_NY^2))
  }
  sigma_both_no <- sigma_both_no[-1]
  sigma_both_cop <- sigma_both_cop[-1]




  pred_NY_no <- c(NYSL.l$NY[startyear])
  pred_NY_cop <- c(NYSL.l$NY[startyear])


  ssp119<- c(NYSL.l$NY[startyear])
  ssp126<- c(NYSL.l$NY[startyear])
  ssp245<- c(NYSL.l$NY[startyear])
  ssp370<- c(NYSL.l$NY[startyear])
  ssp460<- c(NYSL.l$NY[startyear])
  ssp585<- c(NYSL.l$NY[startyear])



  for(i in 1:timespan){
    pred_NY_no <- c(pred_NY_no, predict(fit_NY, newdata = data.frame(NY.l = pred_NY_no[i], GMSL = pred_GMSL$mean_no[startyear-1+i])))
    pred_NY_cop <- c(pred_NY_cop, predict(fit_NY, newdata = data.frame(NY.l = pred_NY_cop[i], GMSL = pred_GMSL$mean_cop[startyear-1+i])))
    ssp119 <- c(ssp119, predict(fit_NY, newdata = data.frame(NY.l =ssp119[i], GMSL = pred_GMSL$ssp119[startyear-1+i])))
    ssp126 <- c(ssp126, predict(fit_NY, newdata = data.frame(NY.l =ssp126[i], GMSL = pred_GMSL$ssp126[startyear-1+i])))
    ssp245 <- c(ssp245, predict(fit_NY, newdata = data.frame(NY.l =ssp245[i], GMSL = pred_GMSL$ssp245[startyear-1+i])))
    ssp370 <- c(ssp370, predict(fit_NY, newdata = data.frame(NY.l =ssp370[i], GMSL = pred_GMSL$ssp370[startyear-1+i])))
    ssp460 <- c(ssp460, predict(fit_NY, newdata = data.frame(NY.l =ssp460[i], GMSL = pred_GMSL$ssp460[startyear-1+i])))
    ssp585 <- c(ssp585, predict(fit_NY, newdata = data.frame(NY.l =ssp585[i], GMSL = pred_GMSL$ssp585[startyear-1+i])))

  }
  pred_NY_no <- pred_NY_no[-1]
  pred_NY_cop <- pred_NY_cop[-1]
  ssp119 <- ssp119[-1]
  ssp126 <- ssp126[-1]
  ssp245 <- ssp245[-1]
  ssp370 <- ssp370[-1]
  ssp460 <- ssp460[-1]
  ssp585 <-ssp585[-1]


  index2010=which(NYSL$Year==2010)



  startyear=startyear+1

  num.nas = length(NYSL$NY[index2010:startyear])-1
  if(num.nas<0){
    num.nas=0
  }

  pred_NY <- data.frame(Year = 2010:2100,
                        mean_no = c(NYSL$NY[index2010:startyear], pred_NY_no),
                        up95_no = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_no + qnorm(0.975) * sigma_both_no),
                        low95_no = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_no - qnorm(0.975) * sigma_both_no),
                        up99_no = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_no + qnorm(0.995) * sigma_both_no),
                        low99_no = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_no - qnorm(0.995) * sigma_both_no),
                        mean_cop = c(NYSL$NY[index2010:startyear], pred_NY_cop),
                        up95_cop = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_cop + qnorm(0.975) * sigma_both_cop),
                        low95_cop = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_cop - qnorm(0.975) * sigma_both_cop),
                        up99_cop = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_cop + qnorm(0.995) * sigma_both_cop),
                        low99_cop = c(rep(NA,num.nas),NYSL$NY[startyear], pred_NY_cop - qnorm(0.995) * sigma_both_cop),
                        ssp119 = c(NYSL$NY[index2010:startyear], ssp119),
                        ssp126 = c(NYSL$NY[index2010:startyear], ssp126),
                        ssp245 = c(NYSL$NY[index2010:startyear], ssp245),
                        ssp370 = c(NYSL$NY[index2010:startyear], ssp370),
                        ssp460 = c(NYSL$NY[index2010:startyear], ssp460),
                        ssp585 = c(NYSL$NY[index2010:startyear], ssp585))


  return(pred_NY)


}
