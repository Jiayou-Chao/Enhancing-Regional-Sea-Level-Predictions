source("C:/Users/andre/OneDrive/Desktop/Thesis/scratch.R")

DATA_DIR_SSE = "C:\\Users\\andre\\OneDrive\\Desktop\\Thesis\\result\\"

SSE <- function(x,y){
  return(sum((x-y)^2,na.rm = T))
}

callibration_SSE<-function(station, station_list, name_list, data_dir=DATA_DIR_SSE){
  my_data = NULL
  them_data = data.frame("month"=NULL,"sa"=NULL,"tg"=NULL)
  for (s in station_list){
    data = get_month_data_raw2(s,drop_na_columns = c(sa,tg),adjust_by_mean = T,data_dir = data_dir,add_VLM=T)
    if(length(data$month)>50){
      new_dat = data[sample(length(data$month),50),]
      if(min(new_dat$sa)>-0.5 & s != station){
        them_data = rbind(them_data,new_dat)

      }
      if(station==s){
        my_data = data
      }
    }
  }

  result = data.frame("NoCali_SSE"=SSE(my_data$tg,my_data$sa))

  mods = lmodel2(tg ~ sa, them_data)

  my_GMR_preds = my_data$sa*mods$regression.results[3,"Slope"]+mods$regression.results[3,"Intercept"]

  result$GMR_SSE = SSE(my_GMR_preds,my_data$tg)

  my_OLS_preds = my_data$sa*mods$regression.results[1,"Slope"]+mods$regression.results[1,"Intercept"]

  result$OLS_SSE = SSE(my_OLS_preds,my_data$tg)

  return(result)

}


sse_stuff = data.frame("NoCali_SSE"=NULL,"GMR_SSE"=NULL,"OLS_SSE"=NULL)
for(i in 1:length(west_coast_stations)){
  sse_stuff = rbind(sse_stuff,callibration_SSE(west_coast_stations[i],west_coast_stations,west_coast_names))

}