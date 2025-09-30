source("SEM_data_preprocessing.R")
library(lavaan)
library(regsem)



DATA_DIR = "data/plots"

fit_SEM <- function(data, lam = 0.1){
  #' Fit a Structural Equation Model (SEM) to the data
  #' @param data data.frame, the input data for fitting the SEM model.
  #' @param lam numeric, the lambda value for regularization, default is 0.1.
  #' @return object, the fitted SEM model.
  #' @export
  #' @examples
  #' fit_SEM(data = my_data, lam = 0.1)



  mod_formula <- '
    LSL ~ LSL.l + GMSL + X2m_dewpoint_temperature + X2m_temperature + surface_pressure + total_precipitation + 1
    total_precipitation ~ total_precipitation.l + X2m_dewpoint_temperature + X2m_temperature + surface_pressure + 1
    X2m_dewpoint_temperature ~ X2m_dewpoint_temperature.l + X2m_temperature +  surface_pressure + total_precipitation.l + 1
    surface_pressure ~ surface_pressure.l + X2m_temperature + 1
    X2m_temperature ~ X2m_temperature.l + 1
  '
  model <- sem(mod_formula, data=data)
  set.seed(123)
  out.reg <- multi_optim(model, max.try = 100, lambda = lam, type="alasso",
                         pars_pen = c(1:6,8:11,13:16,18:19,21:21))
  return(out.reg)
}

auto_fit_SEM <- function(LSL_data, starting.lam = 1,...){
  #' Automatically fit a SEM model by optimizing the lambda value
  #' @param LSL_data data.frame, the input data for fitting the SEM model.
  #' @param starting.lam numeric, the starting value for lambda in regularization, default is 1.
  #' @param ... additional arguments passed to other functions.
  #' @return object, the fitted SEM model with the optimal lambda.
  #' @export
  #' @examples
  #' auto_fit_SEM(LSL_data = my_data, starting.lam = 1)

  opt.lam = starting.lam
  opt.BIC = summary(fit_SEM(LSL_data,lam = opt.lam))$returnVals[5]
  new.BIC = opt.BIC
  new.lam = opt.lam
  max_iterations <- 10
  iterations <- 0

  while (new.BIC == opt.BIC && iterations < max_iterations) {
    iterations <- iterations + 1
    new.lam = new.lam/2
    new.BIC = summary(fit_SEM(LSL_data,lam = new.lam))$returnVals[5]
    if(new.BIC <= opt.BIC){
      opt.BIC = new.BIC
      opt.lam = new.lam
    }
  }

  if (iterations == max_iterations) {
    warning("Reached the maximum number of iterations.")
  }

  return(fit_SEM(LSL_data, lam = opt.lam))


}

get_mod_coef <- function(lhs_var,rhs_var,model){
  #' Get the coefficient for a given pair of variables from the SEM model
  #' @param lhs_var character, the left-hand side variable.
  #' @param rhs_var character, the right-hand side variable.
  #' @param model object, the fitted SEM model.
  #' @return numeric, the coefficient for the specified pair of variables.
  #' @export
  #' @examples
  #' get_mod_coef(lhs_var = "LSL", rhs_var = "GMSL", model = my_model)
  coef_name <- paste(rhs_var,"->",lhs_var)
  return(model$coefficients[[coef_name]])
}

SEM_A_matrix <- function(lhs_vars, model){
  #' Construct the A matrix for SEM prediction
  #' @param lhs_vars character vector, the left-hand side variables.
  #' @param model object, the fitted SEM model.
  #' @return matrix, the A matrix representing relationships among endogenous variables.
  #' @export
  #' @examples
  #' SEM_A_matrix(lhs_vars = c("LSL", "temperature"), model = fitted_model)
  A = matrix(nrow = length(lhs_vars), ncol = length(lhs_vars))

  for(row in 1:length(lhs_vars)){
    for(col in 1:length(lhs_vars)){
      if(row == col){
        A[row,col]=1.0
      }else{

        coef <- tryCatch(
          {
            get_mod_coef(lhs_vars[row],lhs_vars[col],model)
          },
          error = function(){return(0)}
        )
        if(is.null(coef)){
          coef <- 0
        }

        A[row, col] = -coef

      }

    }
  }

  return(A)
}

SEM_B_matrix <- function(lhs_vars, rhs_vars, model){
  #' Construct the B matrix for SEM prediction
  #' @param lhs_vars character vector, the left-hand side variables.
  #' @param rhs_vars character vector, the right-hand side variables.
  #' @param model object, the fitted SEM model.
  #' @return matrix, the B matrix representing relationships between endogenous and exogenous variables.
  #' @export
  #' @examples
  #' SEM_B_matrix(lhs_vars = c("LSL", "temperature"), rhs_vars = c("GMSL", "temperature.l"), model = fitted_model)
  B = matrix(nrow = length(lhs_vars), ncol = length(rhs_vars))

  for(row in 1:length(lhs_vars)){
    for(col in 1:length(rhs_vars)){
      coef <- tryCatch(
        {
          get_mod_coef(lhs_vars[row],rhs_vars[col],model)
        },
        error = function(){return(0)}
      )
      if(is.null(coef)){
        coef <- 0
      }

      B[row, col] = coef



    }
  }
  return(B)
}

SEM_prediction_matrix <- function(data, model){
  #' Create the prediction matrix for SEM
  #' @param data data.frame, the input data.
  #' @param model object, the fitted SEM model.
  #' @return matrix, the prediction matrix combining A and B matrices.
  #' @export
  #' @examples
  #' SEM_prediction_matrix(data = my_data, model = my_model)
  variables = names(data)[!grepl("\\.l$",names(data))]

  lhs_vars = variables[which(variables != "GMSL" & variables != "Year" & variables != "X10m_u_component_of_wind" & variables != "X10m_v_component_of_wind" & variables !="VLM")]
  rhs_vars = c("GMSL",paste0(lhs_vars, ".l"),"1")

  A = SEM_A_matrix(lhs_vars,model)

  B = SEM_B_matrix(lhs_vars,rhs_vars,model)

  return(inv(A)%*%B)

}

SEM_prediction <- function(location,measurement = "tg", starting.lam = 1, start.year = -Inf, stop.year = Inf, horizon = 2100,data_dir = DATA_DIR, restriction = "no",...){
  LSL_data = preprocess_SEM_data(location, start.year = start.year, stop.year = stop.year, measurement.type = measurement)
  recovery_data = read.csv(file.path(data_dir,location,"era5data/recovery_data.csv"),header = T)
  stop.year <- max(LSL_data$Year)

  LSL_data = LSL_data[,which(names(LSL_data)!= "X10m_u_component_of_wind" & names(LSL_data) != "X10m_v_component_of_wind" & names(LSL_data)!= "X10m_u_component_of_wind.l" & names(LSL_data) != "X10m_v_component_of_wind.l" & names(LSL_data)!="VLM")]

  if(grepl("ssp",restriction)){
    GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year",restriction)]
    GMSL_data = cbind(GMSL_data,"sigma"=rep(NA,length(GMSL_data$Year)))


  }
  else if(restriction == "cop"){
    GMSL_data = read.csv('data/GMSL_prediction_SSP_perma.csv')[,c("Year","mean_cop","sigma_cop")]

  }
  else{
    GMSL_data = predict_GMSL(stop.year)[,c("Year","mean","sigma")]

  }

  GMSL_data[,2] = (GMSL_data[,2]-recovery_data$GMSL_Mean)/recovery_data$GMSL_SD
  GMSL_data[,3] = GMSL_data[,3]/recovery_data$GMSL_SD


  mod_SEM <- auto_fit_SEM(LSL_data, starting.lam = starting.lam)





  variables = names(LSL_data)[!grepl("\\.l$",names(LSL_data))]
  lhs_vars = variables[which(variables != "GMSL" & variables != "Year" & variables != "X10m_u_component_of_wind" & variables != "X10m_v_component_of_wind" & variables != "VLM")]


  print(get_shapiro_pvals(lhs_vars,LSL_data,mod_SEM))

  pred_mat <- SEM_prediction_matrix(LSL_data, mod_SEM)


  pred_vcv <- rbind(c(1,rep(0,length(pred_mat[1,])-2)),pred_mat[,1:length(pred_mat[1,])-1])

  inv_A_mat <- eye(length(pred_vcv[1,]))
  inv_A_mat[2:(length(pred_vcv[1,])),2:(length(pred_vcv[1,]))] <- inv(SEM_A_matrix(lhs_vars, mod_SEM))

  init_vcv_mat <- get_vcv(LSL_data, mod_SEM)
  new_vcv_mat <- init_vcv_mat

  pred_seed <- unlist(c(GMSL_data[which(GMSL_data$Year == stop.year+1),2],LSL_data[which(LSL_data$Year == stop.year),lhs_vars],1))

  LSL_pred = data.frame("Year"=stop.year,"LSL_mean"=pred_seed[2],"LSL_var"= 0)


  while(stop.year < horizon){
    new_pred = pred_mat%*%pred_seed

    if(stop.year<2100){
      r_sigma <- GMSL_data[which(GMSL_data$Year==stop.year),3]/sqrt(new_vcv_mat[1,1])


      new_vcv_mat[,1] = r_sigma * new_vcv_mat[,1]
      new_vcv_mat[1,] = r_sigma * new_vcv_mat[1,]
    }


    new_vcv_mat <- pred_vcv%*%new_vcv_mat%*%t(pred_vcv) + inv_A_mat%*%init_vcv_mat%*%t(inv_A_mat)

    stop.year = stop.year+1
    LSL_pred = rbind(LSL_pred,data.frame("Year"=stop.year,"LSL_mean"=new_pred[1],"LSL_var"=new_vcv_mat[2,2]))


    pred_seed = c(GMSL_data[which(GMSL_data$Year == stop.year+1),2],new_pred,1)




  }


  LSL_pred$LSL_mean = LSL_pred$LSL_mean * recovery_data$LSL_SD + recovery_data$LSL_Mean
  LSL_pred$LSL_var = LSL_pred$LSL_var * recovery_data$LSL_SD^2
  LSL_pred$low_95 = LSL_pred$LSL_mean - qnorm(0.975)*sqrt(LSL_pred$LSL_var)
  LSL_pred$up_95 = LSL_pred$LSL_mean + qnorm(0.975)*sqrt(LSL_pred$LSL_var)
  LSL_pred$low_99 = LSL_pred$LSL_mean - qnorm(0.995)*sqrt(LSL_pred$LSL_var)
  LSL_pred$up_99 = LSL_pred$LSL_mean + qnorm(0.995)*sqrt(LSL_pred$LSL_var)

  return(LSL_pred)

}

get_var_res <- function(var_name, data, model){
  var_list = names(data)
  coefs = c()
  intercept = tryCatch(
    {
      get_mod_coef(var_name,1,model)
    },
    error = function(){return(0)}
  )

  if(is.null(coef)){
    coef <- 0
  }
  for(var in var_list){
    coef <- tryCatch(
      {
        get_mod_coef(var_name,var,model)
      },
      error = function(){return(0)}
    )
    if(is.null(coef)){
      coef <- 0
    }
    coefs = c(coefs,coef)

  }

  preds = apply(data,1,function(row){
    return(sum(coefs*row)+intercept)
  })

  return(data[,which(var_list==var_name)]-preds)

}

get_vcv <- function(data, model){
  variables = names(data)[!grepl("\\.l$",names(data))]
  lhs_vars = variables[which(variables != "GMSL" & variables != "Year")]

  return(cov(cbind(data$GMSL,sapply(lhs_vars,get_var_res,data,model))))
}


get_shapiro_pvals <- function(var_names, data, model){
  out = data.frame("Variable"=NULL,"p.value"=NULL)
  for(name in var_names){
    res = get_var_res(name,data,model)
    out = rbind(out,data.frame("Variable"=name,"p.value"=shapiro.test(res)$p.value))

  }
  return(out)
}

get_adf_pvals <- function(var_names, data, model){
  out = data.frame("Variable"=NULL,"p.value"=NULL)
  for(name in var_names){
    res = get_var_res(name,data,model)
    out = rbind(out,data.frame("Variable"=name,"p.value"=adf.test(res)$p.value))

  }
  return(out)
}

