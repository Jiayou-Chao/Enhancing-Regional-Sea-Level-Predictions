library(MTS)
library(dplyr)
source("SEM_data_preprocessing.R")
source("SEM_Fitting.R")
source("VLM.R")
library(MASS)

calculate_varma_aic <- function(p, q, data) {
  fit <- refVARMA(VARMA(data, p, q))
  if (is.nan(fit$aic)){
    return(Inf)
  }
  return(fit$aic)
}

calculate_varx_criteria <- function(p, y, x, p.value, method = "AIC", criteria_index = 1) {
  mod = VARX(y, p, x)
  df = length(x[,1]) - length(mod$se.coef)
  t = abs(qt(p.value/2, df))
  if(length(y)==1){
    fit <- capture.output(VARX(y, p, x))
  }else{
    fit <- capture.output(refVARX(VARX(y, p, x), thres = t))
  }

  if (method == "AIC") {
    CriteriaText = fit[which(grepl("AIC: ", fit))][criteria_index]
  } else if (method == "BIC") {
    CriteriaText = fit[which(grepl("BIC: ", fit))][criteria_index]
  } else {
    return(Inf)
  }

  criteria = as.numeric(substr(CriteriaText, start = 5, stop = nchar(CriteriaText)))

  if (is.na(criteria)) {
    return(Inf)
  }
  return(criteria)
}

.calculate_varx_one_step_cov <- function(sigma_first, sigma_last, mat, s_next){
  sigma_first[1,1] = s_next ^ 2
  sigma_next = sigma_first + mat %*% sigma_last %*% t(mat)

  return(sigma_next)
}

varx_prediction_error <- function(varx_mod, GMSL_data){
  sigma_0 <- eye(dim(varx_mod$Sigma)[1]+1)
  sigma_0[2:(dim(varx_mod$Sigma)[1]+1),2:(dim(varx_mod$Sigma)[1]+1)] = varx_mod$Sigma
  sigma_0[1,1]=GMSL_data$GMSL_sigma[1]^2

  if(dim(varx_mod$beta)[2]==1){
    pi_mat <- rbind(c(0,rep(0,dim(varx_mod$Sigma)[1])),cbind(varx_mod$beta,varx_mod$Phi))

  }else{
    pi_mat <- rbind(c(0,rep(0,dim(varx_mod$Sigma)[1])),cbind(varx_mod$beta[,"GMSL"],varx_mod$Phi))

  }
  sigma <- sigma_0
  error <- c(sqrt(sigma_0[2,2]))
  for(i in 2:length(GMSL_data$GMSL_sigma)){
    sigma <- .calculate_varx_one_step_cov(sigma_0,sigma,pi_mat,GMSL_data$GMSL_sigma[i])
    error <- c(error,sqrt(sigma[2,2]))
  }

  return(error)


}

.forward_varx <- function(y, x, p, p.value, method, current_vars, candidate_vars, max_iter = 10) {
  candidate_vars <- setdiff(candidate_vars, current_vars)
  best_vars <- current_vars
  best_criteria <- calculate_varx_criteria(p, y[, current_vars, drop = F], x, p.value, method)

  improved <- TRUE
  iteration <- 0
  while (improved && iteration < max_iter) {
    improved <- FALSE
    iteration <- iteration + 1
    improved <- FALSE
    for(var in candidate_vars){
      temp_vars<-c(current_vars,var)
      temp_criteria <- calculate_varx_criteria(p,y[,temp_vars,drop=F],x,p.value, method)
      if(temp_criteria<best_criteria){
        best_criteria <- temp_criteria
        best_vars <- temp_vars
        improved <- TRUE
      }
      candidate_vars <- setdiff(candidate_vars,best_vars)
    }
  }

  if (iteration >= max_iter) {
    warning("Reached maximum iterations in .forward_varx")
  }

  if (iteration >= max_iter) {
    warning("Reached maximum iterations in .backward_varx")
  }

  if (iteration >= max_iter) {
    warning("Reached maximum iterations in .bidirectional_varx")
  }

  return(best_vars)
}

.backward_varx <- function(y, x, p, p.value, method, current_vars, max_iter = 10) {
  best_vars <- current_vars
  best_criteria <- calculate_varx_criteria(p, y[, current_vars, drop = F], x, p.value, method)

  improved <- TRUE
  iteration <- 0
  while (improved && iteration < max_iter) {
    improved <- FALSE
    iteration <- iteration + 1
    improved <- F
    for(var in current_vars[which(current_vars != "LSL")]){
      temp_vars <- setdiff(current_vars, var)
      if (length(temp_vars) > 0) {
        temp_criteria <- calculate_varx_criteria(p, y[, temp_vars, drop = FALSE], x, p.value, method)
        if (temp_criteria < best_criteria) {
          best_criteria <- temp_criteria
          best_vars <- temp_vars
          improved <- TRUE
        }
      }
    }
    current_vars<-best_vars
  }
  return(best_vars)
}

.bidirectional_varx <- function(y, x, p, p.value, method, current_vars, candidate_vars, max_iter = 10) {
  best_vars <- current_vars
  best_criteria <- calculate_varx_criteria(p, y[, current_vars, drop = F], x, p.value, method)

  improved <- TRUE
  iteration <- 0
  while (improved && iteration < max_iter) {
    improved <- FALSE
    iteration <- iteration + 1
    improved <- F
    temp_vars <- .forward_varx(y, x, p, p.value, method, current_vars, candidate_vars)
    temp_vars <- .backward_varx(y, x, p, p.value, method, temp_vars)
    temp_criteria <- calculate_varx_criteria(p,y[,temp_vars,drop=F],x,p.value,method)
    if(temp_criteria<best_criteria){
      best_criteria <- temp_criteria
      best_vars <- temp_vars
      improved <- TRUE
    }
    current_vars <- best_vars
  }
  return(best_vars)
}

stepwise_varx <- function(y, x, p, p.value, direction = "both", method = "AIC", max_iter = 10) {

  if(direction == "forward"){
    current_vars = "LSL"
    candidate_vars <- names(y)
    best_vars <- .forward_varx(y, x, p, p.value, method, current_vars, candidate_vars, max_iter)
  }
  else if(direction == "backward"){
    current_vars <- names(y)
    best_vars <- .backward_varx(y, x, p, p.value, method, current_vars, max_iter)
  }
  else if(direction == "both"){
    current_vars = "LSL"
    candidate_vars <- names(y)
    best_vars <- .bidirectional_varx(y, x, p, p.value, method, current_vars, candidate_vars, max_iter)
  }

  return(best_vars)
}

get_VARX <- function(y, x, min.p = 1, max.p = 1, p.value = 0.20, var_sel_method = "AIC", stepwise = FALSE, direction = "both", max_iter = 10) {
  best_aic <- Inf
  best_order <- 0

  if (stepwise) {
    criteria_index <- 1
  } else {
    criteria_index <- 2
  }

  selected_vars <- stepwise_varx(y, x, 1, p.value, direction, var_sel_method, max_iter)


  for (p in min.p:max.p) {
    current_aic <- tryCatch({
      if (stepwise) {
        selected_vars <- stepwise_varx(y, x, p, p.value, direction, var_sel_method, max_iter)
        if (length(selected_vars) == 0) {
          stop("No variables selected in stepwise selection.")
        }
        calculate_varx_criteria(p, y[, selected_vars, drop = FALSE], x, p.value, var_sel_method, criteria_index)
      } else {
        calculate_varx_criteria(p, y, x, p.value, var_sel_method, criteria_index)
      }
    }, error = function(e) { Inf })
    print(current_aic)
    if (current_aic < best_aic) {
      best_aic <- current_aic
      best_order <- p
    }
  }

  best_p <- best_order
  print(best_p)
  if (stepwise) {
    selected_vars <- stepwise_varx(y, x, best_p, p.value, direction, var_sel_method, max_iter)
    if (length(selected_vars) == 0) {
      stop("No variables selected in stepwise selection.")
    }
    best_varx_model <- VARX(y[, selected_vars, drop = FALSE], best_p, x)
  } else {
    mod <- VARX(y, best_p, x)
    df = length(x[,1]) - length(mod$se.coef)
    t = abs(qt(p.value/2, df))
    best_varx_model <- refVARX(mod, thres = t)
  }


  return(best_varx_model)
}

VARMA_stepAIC<- function(data, max.p = 2, max.q = 1){
  best_aic <- Inf
  best_order <- c(0, 0)

  for (p in 0:max.p) {
    for (q in 0:(max.q)) {
      current_aic <- tryCatch({
        calculate_varma_aic(p, q, data)
      },error = function(e){Inf})

      if (current_aic < best_aic) {
        best_aic <- current_aic
        best_order <- c(p, q)
      }
    }
  }

  best_p <- best_order[1]
  best_q <- best_order[2]

  best_varma_model <- refVARMA(VARMA(data, best_p, best_q))

  return(best_varma_model)
}

get_XVARMA <- function(data, max.p = 2, max.q = 1){
  variables = names(data)[!grepl("\\.l$",names(data))]
  varnames = variables[which(variables != "GMSL" & variables != "Year")]

  out = data.frame("GMSL"= data$GMSL)
  GMSL_coefs = rep(0,length(varnames))


  for(i in 1:length(varnames)){
    var = varnames[i]
    var_mod = lm(data[,var]~data$GMSL-1)
    best_var_mod = stepAIC(var_mod,direction="both")

    out = cbind(out, data.frame("NewVar"=best_var_mod$residuals))
    names(out)[i+1] <- var

    if(length(best_var_mod$coefficients)>0){
      GMSL_coefs[i] = best_var_mod$coefficients[1]
    }

  }

  resid_model = VARMA_stepAIC(out[,which(names(out)!="GMSL")], max.p = max.p, max.q = max.q)
  return(list(coef = GMSL_coefs, resmod = resid_model))
}


VARX_prediction <- function(location, measurement = "tg", start.year = -Inf,
                            stop.year = Inf, horizon = 2100, data_dir = DATA_DIR, restriction = "no",
                            selected_vars = c("LSL", "X10m_u_component_of_wind",
                                              "X10m_v_component_of_wind", "surface_pressure", "total_precipitation"
                                              , "X2m_temperature","X2m_dewpoint_temperature"
                                              ),
                            return_fitted_history_only = FALSE, model_only = FALSE, var_sel_method = "BIC", stepwise = TRUE, direction = "both", max_iter = 10, ...) {
  #' Generate VARX model predictions for sea level data.
  #'
  #' @param location Character. The location for which to generate predictions.
  #' @param measurement Character. The type of measurement (default is "tg").
  #' @param start.year Numeric. The start year for the data (default is -Inf).
  #' @param stop.year Numeric. The stop year for the data (default is Inf).
  #' @param horizon Numeric. The prediction horizon (default is 2100).
  #' @param data_dir Character. The directory containing the data (default is DATA_DIR).
  #' @param restriction Character. The restriction type for GMSL data (default is "no").
  #' @param selected_vars Character vector. The selected variables for the model (default includes several climate variables).
  #' @param return_fitted_history_only Logical. If TRUE, return only the fitted history (default is FALSE).
  #' @param model_only Logical. If TRUE, return only the fitted model (default is FALSE).
  #' @param var_sel_method Character. The variable selection method ("AIC" or "BIC", default is "AIC").
  #' @param stepwise Logical. If TRUE, use stepwise selection for variables (default is FALSE).
  #' @param direction Character. The direction for stepwise selection ("both", "forward", or "backward", default is "both").
  #' @param ... Additional arguments passed to the VARX model.
  #' @return A data frame with the VARX predictions. Column names are "Year", "LSL_mean", "LSL_var", "low_95", "up_95", "low_99", "up_99".
  #'
  #' @details This function preprocesses the sea level data, fits a VARX model with optional stepwise variable selection, and generates predictions for the specified horizon. The function can also return the fitted model or the fitted history only.
  #'
  #' @examples
  #' VARX_prediction("battery_park", start.year = 1950, stop.year = 2019, stepwise = TRUE)
  #'
  LSL_data <- preprocess_SEM_data(location, start.year = start.year,
                                  stop.year = stop.year, measurement.type = measurement, add_VLM = TRUE)
  recovery_data <- read.csv(file.path(data_dir, location,
                                      "era5data/recovery_data.csv"), header = TRUE)
  stop.year <- max(LSL_data$Year)

  if (grepl("ssp", restriction)) {
    GMSL_data <- read.csv('data/GMSL_prediction_SSP_perma.csv')[, c("Year",
                                                                    restriction)]
    GMSL_data <- cbind(GMSL_data, "sigma" = rep(NA,
                                                length(GMSL_data$Year)))
  } else if (restriction == "cop") {
    GMSL_data <- read.csv('data/GMSL_prediction_SSP_perma.csv')[, c("Year",
                                                                    "mean_cop", "sigma_cop")]
    names(GMSL_data) <- c("Year", "GMSL", "GMSL_sigma")
  } else {
    GMSL_data <- predict_GMSL(stop.year)[, c("Year", "mean", "sigma")]
    names(GMSL_data) <- c("Year", "GMSL", "GMSL_sigma")
  }

  vlm_rate <- get_ar6_vlm_rate(location)
  include_VLM <- abs(vlm_rate) > 0.002

  GMSL_data$GMSL <- (GMSL_data$GMSL - recovery_data$GMSL_Mean) /
    recovery_data$GMSL_SD
  GMSL_data$GMSL_sigma <- GMSL_data$GMSL_sigma / recovery_data$GMSL_SD

  if (include_VLM) {
    GMSL_data <- merge(GMSL_data, get_ar6_vlm(location, start.year =
                                                stop.year, end.year = horizon), by = "Year")
    GMSL_data$VLM <- (GMSL_data$VLM - recovery_data$VLM_Mean) /
      recovery_data$VLM_SD
  }

  variables <- names(LSL_data)[!grepl("\\.l$", names(LSL_data))]
  endog_vars <- variables[which(variables != "GMSL" & variables != "Year" &
                                  variables != "VLM")]
  if (!is.null(selected_vars)) {
    selected_vars <- intersect(endog_vars, selected_vars)
    if ("LSL" %in% selected_vars) {
      endog_vars <- selected_vars
    } else {
      endog_vars <- c("LSL", selected_vars)
    }
  }

  if (include_VLM) {
    exog_vars <- c("GMSL", "VLM")
  } else {
    exog_vars <- c("GMSL")
  }

  LSL_data = na.omit(LSL_data)

  mod_VARX <- get_VARX(LSL_data[endog_vars], LSL_data[exog_vars], var_sel_method = var_sel_method, stepwise = stepwise, direction = direction, max_iter = max_iter, ...)

  if(mod_VARX$aror == 0){
    stop("Model fitting error: aror = 0")
  }

  if(model_only){
    return(mod_VARX)
  }

  stop.year <- min(stop.year, max(LSL_data$Year))

  if (return_fitted_history_only) {
    fitted_values <- VARXpred(mod_VARX, newxt = LSL_data[which(LSL_data$Year<=stop.year), exog_vars], hstep = stop.year - min(LSL_data$Year))$pred %>% data.frame()
    fitted_values$LSL <- fitted_values$LSL * recovery_data$LSL_SD + recovery_data$LSL_Mean
    LSL_true <- LSL_data[which(LSL_data$Year<=stop.year), "LSL"]
    return(data.frame(Year = LSL_data[which(LSL_data$Year<=stop.year), "Year"][2:length(LSL_data[which(LSL_data$Year<=stop.year), "Year"])], fitted = fitted_values$LSL,
                      LSL_true = LSL_true[2:length(LSL_true)]))
  }

  LSL_pred <- data.frame("Year" = stop.year, "LSL_mean" =
                           LSL_data$LSL[which(LSL_data$Year == stop.year)], "LSL_var" = 0)

  predictions_VARX <- VARXpred(mod_VARX, newxt = GMSL_data[which(GMSL_data$Year > stop.year & GMSL_data$Year <= horizon), exog_vars], hstep = horizon - stop.year)
  errors_VARX <- varx_prediction_error(mod_VARX, GMSL_data[which(GMSL_data$Year > stop.year & GMSL_data$Year <= horizon),])
  if(is.null(dim(predictions_VARX$pred))){
    LSL_pred <- rbind(LSL_pred, data.frame("Year" = (stop.year + 1):horizon, "LSL_mean" = predictions_VARX$pred, "LSL_var" = errors_VARX^2))
  }
  else{
    LSL_pred <- rbind(LSL_pred, data.frame("Year" = (stop.year + 1):horizon, "LSL_mean" = predictions_VARX$pred[, "LSL"], "LSL_var" = errors_VARX^2))

  }

  LSL_pred$LSL_mean <- LSL_pred$LSL_mean * recovery_data$LSL_SD + recovery_data$LSL_Mean
  LSL_pred$LSL_var <- LSL_pred$LSL_var * recovery_data$LSL_SD^2
  LSL_pred$low_95 <- LSL_pred$LSL_mean - qnorm(0.975) * sqrt(LSL_pred$LSL_var)
  LSL_pred$up_95 <- LSL_pred$LSL_mean + qnorm(0.975) * sqrt(LSL_pred$LSL_var)
  LSL_pred$low_99 <- LSL_pred$LSL_mean - qnorm(0.995) * sqrt(LSL_pred$LSL_var)
  LSL_pred$up_99 <- LSL_pred$LSL_mean + qnorm(0.995) * sqrt(LSL_pred$LSL_var)

  return(LSL_pred)
}

XVARMA_prediction<-function(location,measurement = "tg", max.p = 2, max.q = 1, start.year = -Inf, stop.year = Inf, horizon = 2100,data_dir = DATA_DIR, restriction = "no",...){
  LSL_data = preprocess_SEM_data(location, start.year = start.year, stop.year = stop.year, measurement.type = measurement)
  recovery_data = read.csv(file.path(data_dir,location,"era5data/recovery_data.csv"),header = T)
  stop.year <- max(LSL_data$Year)

  if(grepl("ssp",restriction)){
    GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year",restriction)]
    GMSL_data = cbind(GMSL_data,"sigma"=rep(NA,length(GMSL_data$Year)))


  }
  else{
    GMSL_data = predict_GMSL(stop.year)[,c("Year","mean","sigma")]

  }



  GMSL_data[,2] = (GMSL_data[,2]-recovery_data$GMSL_Mean)/recovery_data$GMSL_SD
  GMSL_data[,3] = GMSL_data[,3]/recovery_data$GMSL_SD



  variables = names(LSL_data)[!grepl("\\.l$",names(LSL_data))]
  mod_vars = variables[which(variables != "GMSL" & variables != "Year")]


  mod_XVARMA = get_XVARMA(LSL_data, max.p = max.p, max.q = max.q)

  stop.year = min(stop.year, max(LSL_data$Year))


  predictions_VARX <- VARMApred(mod_XVARMA$resmod,h = horizon - stop.year)



  LSL_pred = data.frame("Year"=stop.year,"LSL_mean"=LSL_data$LSL[which(LSL_data$Year == stop.year)],"LSL_var"= 0)
  LSL_pred = rbind(LSL_pred,data.frame("Year"=(stop.year+1):horizon,"LSL_mean"=predictions_VARX$pred[,"LSL"],"LSL_var"= predictions_VARX$se.err[,1]^2))


  LSL_pred$LSL_var = LSL_pred$LSL_var + (c(0,mod_XVARMA$coef[1]*GMSL_data[which(GMSL_data$Year>stop.year),3]))^2
  LSL_pred$LSL_mean = LSL_pred$LSL_mean + (c(0,mod_XVARMA$coef[1]*GMSL_data[which(GMSL_data$Year>stop.year),2]))

  LSL_pred$LSL_mean = LSL_pred$LSL_mean * recovery_data$LSL_SD
  LSL_pred$LSL_var = LSL_pred$LSL_var * recovery_data$LSL_SD^2
  LSL_pred$low_95 = LSL_pred$LSL_mean - qnorm(0.975)*sqrt(LSL_pred$LSL_var)
  LSL_pred$up_95 = LSL_pred$LSL_mean + qnorm(0.975)*sqrt(LSL_pred$LSL_var)
  LSL_pred$low_99 = LSL_pred$LSL_mean - qnorm(0.995)*sqrt(LSL_pred$LSL_var)
  LSL_pred$up_99 = LSL_pred$LSL_mean + qnorm(0.995)*sqrt(LSL_pred$LSL_var)

  return(LSL_pred)

}

residual_analysis<- function(location,place = 1, measurement = "tg", max.p = 2, max.q = 1, start.year = -Inf, stop.year = Inf,data_dir = DATA_DIR, restriction = "no",...){
  data = preprocess_SEM_data(location, start.year = start.year, stop.year = stop.year, measurement.type = measurement)
  variables = names(data)[!grepl("\\.l$",names(data))]
  mod_vars = variables[which(variables != "GMSL" & variables != "Year" & variables != "VLM")]

  out = data.frame("SEM" = c(0,1),"VARX" = c(0,1),"XVARMA" = c(0,1))

  SEM_mod = auto_fit_SEM(data)
  out$SEM[1]=sort(get_adf_pvals(mod_vars,data,SEM_mod)$p.value,decreasing = T)[place]
  out$SEM[2]=sort(get_shapiro_pvals(mod_vars,data,SEM_mod)$p.value,decreasing = F)[place]


  VARX_mod = get_VARX(data[mod_vars],data$GMSL)
  VARX_res = VARX_mod$residuals
  VARX_adf = c()
  VARX_shapiro = c()
  for(i in 1:length(VARX_res[1,])){
    VARX_adf = c(VARX_adf,adf.test(VARX_res[,i])$p.value)
    VARX_shapiro = c(VARX_shapiro,shapiro.test(VARX_res[,i])$p.value)
  }
  out$VARX[1]=sort(VARX_adf,decreasing = T)[place]
  out$VARX[2]=sort(VARX_shapiro,decreasing = F)[place]


  XVARMA_mod = get_XVARMA(data)$resmod
  XVARMA_res = XVARMA_mod$residuals
  XVARMA_adf = c()
  XVARMA_shapiro = c()
  for(i in 1:length(XVARMA_res[1,])){
    XVARMA_adf = c(XVARMA_adf,adf.test(XVARMA_res[,i])$p.value)
    XVARMA_shapiro = c(XVARMA_shapiro,shapiro.test(XVARMA_res[,i])$p.value)
  }
  out$XVARMA[1]=sort(XVARMA_adf,decreasing = T)[place]
  out$XVARMA[2]=sort(XVARMA_shapiro,decreasing = F)[place]

  return(out)
}

test_VARX_variables <- function(station_list, measurement = "tg",start.year = -Inf, stop.year = Inf, horizon = 2100,data_dir = DATA_DIR, restriction = "no", selected_vars = c("X10m_u_component_of_wind","X10m_v_component_of_wind","X2m_temperature","surface_pressure","X2m_dewpoint_temperature","total_precipitation"),include_VLM=T,...){
  zero_list <- setNames(as.list(rep(0, length(selected_vars))), selected_vars)

  df <- as.data.frame(zero_list)

  for(location in station_list){
    station_df <- as.data.frame(zero_list)

    LSL_data = preprocess_SEM_data(location, start.year = start.year, stop.year = stop.year, measurement.type = measurement,add_VLM = include_VLM)
    recovery_data = read.csv(file.path(data_dir,location,"era5data/recovery_data.csv"),header = T)
    stop.year <- max(LSL_data$Year)

    if(grepl("ssp",restriction)){
      GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year",restriction)]
      GMSL_data = cbind(GMSL_data,"sigma"=rep(NA,length(GMSL_data$Year)))


    }
    else if(restriction == "cop"){
      GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year","mean_cop","sigma_cop")]
      names(GMSL_data) = c("Year","GMSL","GMSL_sigma")
    }else{
      GMSL_data = predict_GMSL(stop.year)[,c("Year","mean","sigma")]
      names(GMSL_data) = c("Year","GMSL","GMSL_sigma")
    }



    GMSL_data$GMSL = (GMSL_data$GMSL-recovery_data$GMSL_Mean)/recovery_data$GMSL_SD
    GMSL_data$GMSL_sigma = GMSL_data$GMSL_sigma/recovery_data$GMSL_SD

    if(include_VLM){
      GMSL_data = merge(GMSL_data, get_ar6_vlm(location,start.year = stop.year,end.year = horizon), by = "Year")
      GMSL_data$VLM = (GMSL_data$VLM-recovery_data$VLM_Mean)/recovery_data$VLM_SD
      exog_vars = c("GMSL","VLM")
    }
    else{
      exog_vars = c("GMSL")
    }



    for(var in selected_vars){
      mod_VARX = get_VARX(LSL_data[c("LSL",var)],LSL_data[exog_vars], max.p = 1)
      if(mod_VARX$Phi[1,2]!=0){
        station_df[var] = df[var]+1
      }

    }




    df = rbind(df,station_df)




  }

  return(df)
}


test_VARX_ensamble <- function(station_list, measurement = "tg",start.year = -Inf, stop.year = Inf, horizon = 2100,data_dir = DATA_DIR, restriction = "no", selected_vars = c("X10m_u_component_of_wind","X10m_v_component_of_wind","X2m_temperature","surface_pressure","X2m_dewpoint_temperature","total_precipitation"),include_VLM=T,...){
  zero_list <- setNames(as.list(rep(0, length(selected_vars))), selected_vars)

  df <- as.data.frame(zero_list)

  for(location in station_list){
    station_df <- as.data.frame(zero_list)

    LSL_data = preprocess_SEM_data(location, start.year = start.year, stop.year = stop.year, measurement.type = measurement,add_VLM = include_VLM)
    recovery_data = read.csv(file.path(data_dir,location,"era5data/recovery_data.csv"),header = T)
    stop.year <- max(LSL_data$Year)

    if(grepl("ssp",restriction)){
      GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year",restriction)]
      GMSL_data = cbind(GMSL_data,"sigma"=rep(NA,length(GMSL_data$Year)))


    }
    else if(restriction == "cop"){
      GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year","mean_cop","sigma_cop")]
      names(GMSL_data) = c("Year","GMSL","GMSL_sigma")
    }else{
      GMSL_data = predict_GMSL(stop.year)[,c("Year","mean","sigma")]
      names(GMSL_data) = c("Year","GMSL","GMSL_sigma")
    }



    GMSL_data$GMSL = (GMSL_data$GMSL-recovery_data$GMSL_Mean)/recovery_data$GMSL_SD
    GMSL_data$GMSL_sigma = GMSL_data$GMSL_sigma/recovery_data$GMSL_SD

    if(include_VLM){
      GMSL_data = merge(GMSL_data, get_ar6_vlm(location,start.year = stop.year,end.year = horizon), by = "Year")
      GMSL_data$VLM = (GMSL_data$VLM-recovery_data$VLM_Mean)/recovery_data$VLM_SD
      exog_vars = c("GMSL","VLM")
    }
    else{
      exog_vars = c("GMSL")
    }



    mod_VARX = get_VARX(LSL_data[c("LSL",selected_vars)],LSL_data[exog_vars], max.p = 1)


    sig_vars = mod_VARX$Phi[1,2:length(mod_VARX$Phi[1,])]!=0


    station_df[1,] = sig_vars

    df = rbind(df,station_df)





  }

  return(df)
}
