source("ComparePlots.R")
source("data_processing.R")
source("VLM.R")

DATA_DIR = "data/plots"

normalize <- function(x) return((x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE))



preprocess_SEM_data <- function(location, start.year = -Inf, stop.year = Inf, measurement.type = "tg",na.rm=T, add_VLM=F,...){
  #' @param start.year Inclusive start year.
  #' @param stop.year Inclusive stop year.

  dat = get_month_data_raw(location,add_VLM = add_VLM, equal_time=F, ...)
  climate <- read.csv('data/YearData.csv')[,c("Year","GMSL")]



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
  names(annual_sealevel)=c('Year','LSL')
  Year_station=annual_sealevel
  Year_station[Year_station<(-500)]=NA
  Year_station$LSL[is.nan(Year_station$LSL)]=NA
  GMSL <- climate$GMSL

  NYSL <- merge(Year_station, climate, by="Year")



  NYSL.l <- cbind(NYSL[-1,],NYSL[-nrow(NYSL),])
  names(NYSL.l)[4:6] <- paste0(names(NYSL.l)[4:6], '.l')


  all_data = NYSL.l

  datapath = file.path(DATA_DIR, location,"era5data")

  era_data = get_era5_yearly_data(location)

  era_data = era_data[which(names(era_data)!="X")]

  era_lag = era_data
  era_lag$Year = era_lag$Year+1
  names(era_lag)=paste0(names(era_lag),".l")
  names(era_lag)[1]="Year"


  era_data = merge(era_data,era_lag,by = "Year", all = FALSE)


  all_data = merge(all_data,era_data, all=FALSE)

  VLM_data = get_ar6_vlm(location,start.year = min(all_data$Year),end.year = max(all_data$Year))


  all_data = merge(all_data,VLM_data,all=FALSE, by = "Year")

  all_data = all_data[which(all_data$Year <= stop.year & all_data$Year >= start.year),]

  year = all_data$Year

  all_data = subset(all_data, select = -c(Year,Year.l))

  sea_level_recovery = data.frame("GMSL_Mean" = mean(all_data$GMSL, na.rm=T),"GMSL_SD" = sd(all_data$GMSL,na.rm=T),"LSL_Mean" = mean(all_data$LSL, na.rm=T),"LSL_SD" = sd(all_data$LSL,na.rm=T),
                                  "VLM_Mean"=mean(all_data$VLM, na.rm = T),"VLM_SD"=sd(all_data$VLM,na.rm=T))

  all_data = as.data.frame(sapply(all_data,normalize))

  all_data$Year = year

  names(all_data)[grepl("^[0-9]", names(all_data))] = paste0("X",names(all_data)[grepl("^[0-9]", names(all_data))])

  if(!dir.exists(datapath)){
    dir.create(datapath, recursive = TRUE)
  }

  write.csv(all_data, file = file.path(datapath,"SEM_vars.csv"))
  write.csv(sea_level_recovery, file = file.path(datapath,"recovery_data.csv"))

  return(all_data)
}
