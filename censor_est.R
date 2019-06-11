library(tidyverse)
library(magrittr)
library(rlist)
library(igraph)
library(VineCopula)
library(copula)
source("copula_cdfs.R")

# From Hofert et al., (2018)
# this wraps around functions to get timings
withTime <- function(expr, ...)
{
  st <- system.time(r <- expr, ...)
  list(value = r, sys.time = st)
}

censor <- function(x, method = "upper"){
  if(!(method %in% c("upper", "lower"))){
    stop("method has to be %in% c(\"upper\", \"lower\")")
  }
  if(method == "upper"){
    return(copula::pobs(x, ties.method = "max"))
  }else{
    return(copula::pobs(x, ties.method = "min"))
  }
}

check_rotations <- function(rot){
  # Rotation handling
  rotations <- c("0", "90", "180", "270")
  if(is.character(rot)){
    if(!(rot %in% rotations)){
      stop("Invalid rotation. Choose: 1 = \"0\", 2 = \"90\", 3 = \"180\", 4 = \"270\"" )
    }
    rot <- which(rotations == rot)
  }else if(is.numeric(rot)){
    if(!(rot %in% c(1,2,3,4))){
      stop("Invalid rotation. Choose: 1 = \"0\", 2 = \"90\", 3 = \"180\", 4 = \"270\"" )
    }
  }
  return(rot)
}

###################################### LOGLIK ###############################################

log_lik_censor <- function(u, v, cop, rot, check.pars = TRUE){
  cop_fun <- cop$cop_fun
  cop_dens <- cop$cop_dens
  cop_du <- cop$cop_du
  cop_dv <- cop$cop_dv
  theta <- cop$theta
  delta <- cop$delta
  loglik <- 0
  places <- c()
  # rot == 1 is 0 degrees, and no rotation
  if(rot == 2){
    # 90 degree
    u <- 1 - u
  }else if(rot == 3){
    # 180 degree
    u <- 1 - u
    v <- 1 - v
  }else if(rot == 4){
    # 270 degree
    v <- 1 - v
  }
  u_upper <- censor(u, method = "upper")
  u_lower <- censor(u, method = "lower")
  v_upper <- censor(v, method = "upper")
  v_lower <- censor(v, method = "lower")
  
  # Observation not tied in either margin
  places1 <- ((u_lower == u_upper) & (v_lower == v_upper))
  loglik <- sum(log(
    cop_dens(u = u_upper[places1],v = v_upper[places1], theta, delta, check.pars = check.pars)
                    ))
  
  # l smaller in both
  places2 <- ((u_lower < u_upper) & (v_lower < v_upper))
  loglik <- loglik + sum(log(
    cop_fun(u_upper[places2], v_upper[places2], theta, delta, check.pars = check.pars) - 
      cop_fun(u_upper[places2], v_lower[places2], theta, delta, check.pars = check.pars) - 
      cop_fun(u_lower[places2], v_upper[places2], theta, delta, check.pars = check.pars) +
      cop_fun(u_lower[places2], v_lower[places2], theta, delta, check.pars = check.pars)
  ))
  
  # u_lower smaller, v not tied
  places3 <- ((u_lower < u_upper) & (v_lower == v_upper))
  loglik <- loglik + sum(log(
    cop_dv(u_upper[places3], v_upper[places3], theta, delta, check.pars = check.pars) - 
      cop_dv(u_lower[places3], v_upper[places3], theta, delta, check.pars = check.pars)
  ))
  
  # v_lower smaller, u not tied
  places4 <- ((u_lower == u_upper) & (v_lower < v_upper))
  loglik <- loglik + sum(log(
    cop_du(u_upper[places4], v_upper[places4], theta, delta, check.pars = check.pars) - 
      cop_du(u_upper[places4], v_lower[places4], theta, delta, check.pars = check.pars)
  ))
  return(loglik)
}

########################### INTERNAL LOGLIK ###############################################

log_lik_censor.intern <- function(u_upper, u_lower, v_upper, v_lower, cop, check.pars){
  cop_fun <- cop$cop_fun
  cop_dens <- cop$cop_dens
  cop_du <- cop$cop_du
  cop_dv <- cop$cop_dv
  theta <- cop$theta
  delta <- cop$delta
  loglik <- 0
  
  # Observation not tied in either margin
  places1 <- ((u_lower == u_upper) & (v_lower == v_upper))
  loglik <- sum(log(
    cop_dens(u = u_upper[places1],v = v_upper[places1], theta, delta, check.pars = check.pars)
  ))
  
  # l smaller in both
  places2 <- ((u_lower < u_upper) & (v_lower < v_upper))
  loglik <- loglik + sum(log(
    cop_fun(u = u_upper[places2], v = v_upper[places2], theta, delta, check.pars = check.pars) - 
      cop_fun(u = u_upper[places2], v = v_lower[places2], theta, delta, check.pars = check.pars) - 
      cop_fun(u = u_lower[places2], v = v_upper[places2], theta, delta, check.pars = check.pars) +
      cop_fun(u = u_lower[places2], v = v_lower[places2], theta, delta, check.pars = check.pars)
  ))
  
  # u_lower smaller, v not tied
  places3 <- ((u_lower < u_upper) & (v_lower == v_upper))
  loglik <- loglik + sum(log(
    cop_dv(u = u_upper[places3], v = v_upper[places3], theta, delta, check.pars = check.pars) - 
      cop_dv(u = u_lower[places3],v = v_upper[places3], theta, delta, check.pars = check.pars)
  ))
  
  # v_lower smaller, u not tied
  places4 <- ((u_lower == u_upper) & (v_lower < v_upper))
  loglik <- loglik + sum(log(
    cop_du(u = u_upper[places4], v = v_upper[places4], theta, delta, check.pars = check.pars) - 
      cop_du(u = u_upper[places4], v = v_lower[places4], theta, delta, check.pars = check.pars)
  ))
  return(loglik)
}

######################################## ESTIMATION ############################################  

censor_est <- function(u, v, cop, rot = 1, method = "censor", check.pars = TRUE){
  # censor or regular method
  if(!is.list(cop)){
    cop <- get_cop(cop)
  }else if(7>sum(names(cop) == c("cop_name", "theta", "delta","cop_fun","cop_du","cop_dv","cop_dens"))){
    stop("copula mismatch or something")
  }
  # Rotation handling
  rotations <- c("0", "90", "180", "270")
  if(is.character(rot)){
    if(!(rot %in% rotations)){
      stop("Invalid rotation. Choose: 1 = \"0\", 2 = \"90\", 3 = \"180\", 4 = \"270\"" )
    }
    rot <- which(rotations == rot)
  }else if(is.numeric(rot)){
    if(!(rot %in% c(1,2,3,4))){
      stop("Invalid rotation. Choose: 1 = \"0\", 2 = \"90\", 3 = \"180\", 4 = \"270\"" )
    }
  }
  cop$rotation <- rot
  # rot == 1 is 0 degrees, and no rotation
  if(rot == 2){
    # 90 degree
    u <- 1 - u
  }else if(rot == 3){
    # 180 degree
    u <- 1 - u
    v <- 1 - v
  }else if(rot == 4){
    # 270 degree
    v <- 1 - v
  }
  
  if(!(method %in% c("censor", "regular"))){
    stop("Invalid method. Choose from:  \"censor\", \"regular\"" )
  }
  
  if(cop$n.param == 1){
    if(method == "censor"){
      u_upper <- censor(u, method = "upper")
      u_lower <- censor(u, method = "lower")
      v_upper <- censor(v, method = "upper")
      v_lower <- censor(v, method = "lower")
      
      ll <- function(theta){
        cop$theta <- theta
        return(log_lik_censor.intern(u_upper, u_lower, v_upper, v_lower, cop, check.pars = check.pars))
      }
    }else{
      ll <- function(theta){
        cop$theta <- theta
        return(sum(log(
          cop$cop_dens(u = u,v = v, theta = theta, delta, check.pars = check.pars)
        )))
      }
    }
    optimlist <- list()
    objectives <- c()
    for(i in 1:length(cop$optim.limits)){
      cop$theta <- cop$optim.limits[[i]]$start
      optimlist %<>% rlist::list.append(optimize(f = ll, maximum = TRUE,
                                                 interval = c(cop$optim.limits[[i]]$low,
                                                              cop$optim.limits[[i]]$up))) 
      objectives %<>% c(optimlist[[i]]$objective)
    }
    if(sum(is.na(objectives)) == length(objectives)){
      stop("All NA produced by optimize")
    }else{
      cop$theta <- optimlist[[which.max(objectives)]]$maximum
      cop$log.lik <- optimlist[[which.max(objectives)]]$objective
      cop$AIC <- 2*cop$n.param - 2*optimlist[[which.max(objectives)]]$objective
    }
    return(cop)
  }else{
    cop$theta <- cop$start.value[1]
    cop$delta <- cop$start.value[2]
    if(cop$cop_name %in% c("Tawn", "Tawn2")){
      tau <- cor(u,v,method = "kendall")
      cop$tau <- tau
      parlower <- c(1.001, max(tau - 0.1, 1e-04))
      parupper <- c(20, min(tau + 0.2, 0.99))
    }else{
      parlower <- cop$optim.lower
      parupper <- cop$optim.upper
    }
    
    if(method == "censor"){
      u_upper <- censor(u, method = "upper")
      u_lower <- censor(u, method = "lower")
      v_upper <- censor(v, method = "upper")
      v_lower <- censor(v, method = "lower")
      
      ll <- function(par){
        cop$theta <- par[1]
        cop$delta <- par[2]
        return(log_lik_censor.intern(u_upper, u_lower, v_upper, v_lower, cop, check.pars = check.pars))
      }
    }else{
      ll <- function(par){
        cop$theta <- par[1]
        cop$delta <- par[2]
        return(sum(log(
          cop$cop_dens(u = u,v = v, theta = par[1], delta = par[2], check.pars = check.pars)
        )))
      }
    }
    
    optimout <- optim(par = c(cop$theta, cop$delta), fn = ll,
                      method = "L-BFGS-B", 
                      lower = parlower, upper = parupper,
                      control = list(fnscale = -1, maxit = 500))
    cop$theta <- optimout$par[1]
    cop$delta <- optimout$par[2]
    cop$log.lik <- optimout$value
    cop$AIC <- 2*cop$n.param - 2*optimout$value
    return(cop)
  }
}

############################## INTERNAL ESTIMATOR ######################################

censor_est.intern <- function(u_upper, u_lower, v_upper, v_lower, cop, rot, check.pars = TRUE){
  # Internal functions unclude less checks, and recieve rotated data,
  # This one is solely to minimize calculations in parametric bootstrapping
  # NOT ADAPTED TO REGULAR ESTIMATION
  if(cop$n.param == 1){
    ll <- function(theta){
      cop$theta <- theta
      return(log_lik_censor.intern(u_upper = u_upper,u_lower = u_lower,
                                   v_upper = v_upper, v_lower = v_lower,cop = cop, check.pars = check.pars))
    }
    optimlist <- list()
    objectives <- c()
    for(i in 1:length(cop$optim.limits)){
      cop$theta <- cop$optim.limits[[i]]$start
      optimlist %<>% rlist::list.append(optimize(f = ll, maximum = TRUE,
                                          interval = c(cop$optim.limits[[i]]$low,
                                                       cop$optim.limits[[i]]$up))) 
      objectives %<>% c(optimlist[[i]]$objective)
    }
    if(sum(is.na(objectives)) == length(objectives)){
      #stop("All NA produced by optimize")
      warning("All NA produced by optimize")
    }else{
      cop$theta <- optimlist[[which.max(objectives)]]$maximum
      cop$log.lik <- optimlist[[which.max(objectives)]]$objective
      cop$AIC <- 2*cop$n.param - 2*optimlist[[which.max(objectives)]]$objective
    }
    return(cop)
  }else{
    cop$theta <- cop$start.value[1]
    cop$delta <- cop$start.value[2]
    parlower <- cop$optim.lower
    parupper <- cop$optim.upper
    
    ll <- function(par){
      cop$theta <- par[1]
      cop$delta <- par[2]
      return(log_lik_censor.intern(u_upper = u_upper,u_lower = u_lower,
                                   v_upper = v_upper, v_lower = v_lower,cop = cop, check.pars = check.pars))
    }
    optimout <- optim(par = c(cop$theta, cop$delta), fn = ll,#cop$start.value, fn = ll,
                      method = "L-BFGS-B", 
                      lower = parlower, upper = parupper,
                      control = list(fnscale = -1, maxit = 500))
    cop$theta <- optimout$par[1]
    cop$delta <- optimout$par[2]
    cop$log.lik <- optimout$value
    cop$AIC <- 2*cop$n.param - 2*optimout$value
    return(cop)
  }
}

censor_est.intern.tawn <- function(u_upper, u_lower, v_upper, v_lower, cop, rot, check.pars = TRUE){
  # Internal functions unclude less checks, and recieve rotated data,
  # This one is solely to minimize calculations in parametric bootstrapping
  # NOT ADAPTED TO REGULAR ESTIMATION
  cop$theta <- cop$start.value[1]
  cop$delta <- cop$start.value[2]
  tau <- cor(u,v,method = "kendall")
  parlower <- c(1.001, max(tau - 0.1, 1e-04))
  parupper <- c(20, min(tau + 0.2, 0.99))
  
  ll <- function(par){
    cop$theta <- par[1]
    cop$delta <- par[2]
    return(log_lik_censor.intern(u_upper = u_upper,u_lower = u_lower,
                                 v_upper = v_upper, v_lower = v_lower,cop = cop, check.pars = check.pars))
  }
  optimout <- optim(par = c(cop$theta, cop$delta), fn = ll,#cop$start.value, fn = ll,
                    method = "L-BFGS-B", 
                    lower = parlower, upper = parupper,
                    control = list(fnscale = -1, maxit = 500))
  cop$theta <- optimout$par[1]
  cop$delta <- optimout$par[2]
  cop$log.lik <- optimout$value
  cop$AIC <- 2*cop$n.param - 2*optimout$value
  return(cop)
}


####################################################################################################
############################################ COPULA SELECTION #######################################
####################################################################################################

censor_copula_select <- function(u,v, indeptest = FALSE, level = 0.05, method = "censor",
                                 include_tawn = TRUE, include_amh = TRUE, include_t = TRUE){
  tau <- cor(u, v, method = "kendall")
  N <- length(u)
  if(indeptest){
    f <- sqrt((9 * N * (N - 1))/(2 * (2 * N + 5))) * abs(tau)
    p.value = 2 * (1 - pnorm(f))
    if(p.value>level){
      return(get_cop("Independence"))
    }
  }
  
  copula_estimates <- list()
  copula_stats <- tibble(Name = character(), Theta = double(), Delta = double(), Rotation = integer(), AIC = double(), Time = double())
  copula_aics <- c()
  copula_names <- c()
  copula_timings <- c()
  ## Radially symmetric copulas
  copulas <- c("Frank", "normal", "t")#, "t2")
  if(!include_t){
    copulas <- copulas[-3]
  }
  for(cop_name in copulas){
    # (u, v, cop, rot = 1, method = "censor", check.pars = TRUE)
    fitted_copula <- withTime(try(censor_est(u,v,cop_name, rot = 1, method = method)))
    if(class(fitted_copula$value) == "try-error"){
      copula_estimates %<>% rlist::list.append(get_cop(cop_name))
      copula_aics %<>% c(Inf)
      copula_names %<>% c(paste("XX", cop_name, 1))
      copula_stats %<>% bind_rows(tibble(Name = cop_name, Theta = NA, Delta = NA, Rotation = 1,
                                         AIC = Inf, Time = fitted_copula$sys.time[3]))
    }else{
      copula_estimates %<>% rlist::list.append(fitted_copula$value)
      copula_aics %<>% c(fitted_copula$value$AIC)
      copula_names %<>% c(paste(fitted_copula$value$cop_name,1))
      copula_stats %<>% bind_rows(tibble(Name = cop_name, Theta = fitted_copula$value$theta, Delta = fitted_copula$value$delta,
                                         Rotation = 1, AIC = fitted_copula$value$AIC, Time = fitted_copula$sys.time[3]))
    }
  }
  
  ## 0 and 180 degree rotations
  if(tau>0){
    copulas <- c("AMH", "Clayton", "Gumbel", "Joe",
                 "BB1", "BB6", "BB7", "BB8", "Tawn", "Tawn2")
    if((tau>=1/3) | (!include_amh)){
      copulas <- copulas[-1]
    }
    if(!include_tawn){
      copulas <- copulas[!copulas %in% c("Tawn", "Tawn2")]
    }
  }else{
    copulas <- c("AMH", "Clayton")
    if((tau<=(5 - 8/log(2)) / 3) | (!include_amh)){
      copulas <- copulas[-1]
    }
  }
  for(cop_name in copulas){
    for(rot in c(1,3)){
      fitted_copula <- withTime(try(censor_est(u,v,cop_name, rot = rot, method = method)))
      if(class(fitted_copula$value) == "try-error"){
        copula_estimates %<>% rlist::list.append(get_cop(cop_name))
        copula_aics %<>% c(Inf)
        copula_names %<>% c(paste("XX", cop_name, rot))
        copula_stats %<>% bind_rows(tibble(Name = cop_name, Theta = NA, Delta = NA,
                                           Rotation = rot, AIC = Inf, Time = fitted_copula$sys.time[3]))
      }else{
        copula_estimates %<>% rlist::list.append(fitted_copula$value)
        copula_aics %<>% c(fitted_copula$value$AIC)
        copula_names %<>% c(paste(fitted_copula$value$cop_name,rot))
        copula_stats %<>% bind_rows(tibble(Name = cop_name, Theta = fitted_copula$value$theta, Delta = fitted_copula$value$delta,
                                           Rotation = rot, AIC = fitted_copula$value$AIC, Time = fitted_copula$sys.time[3]))
      }
    }
  }
  ## 90 and 270 degree rotations
  if(tau<0){
    copulas <- c("AMH", "Clayton", "Gumbel", "Joe",
                 "BB1", "BB6", "BB7", "BB8", "Tawn", "Tawn2")
    if((tau<=-1/3) | (!include_amh)){
      copulas <- copulas[-1]
    }
    if(!include_tawn){
      copulas <- copulas[!copulas %in% c("Tawn", "Tawn2")]
    }
  }else{
    copulas <- c("AMH", "Clayton")
    if((tau>=-(5 - 8/log(2)) / 3) | (!include_amh)){
      copulas <- copulas[-1]
    }
  }
  for(cop_name in copulas){
    for(rot in c(2,4)){
      fitted_copula <- withTime(try(censor_est(u,v,cop_name, rot = rot, method = method)))
      if(class(fitted_copula$value) == "try-error"){
        copula_estimates %<>% rlist::list.append(get_cop(cop_name))
        copula_aics %<>% c(Inf)
        copula_names %<>% c(paste("XX", cop_name,rot))
        copula_stats %<>% bind_rows(tibble(Name = cop_name, Theta = NA, Delta = NA,
                                           Rotation = rot, AIC = Inf, Time = fitted_copula$sys.time[3]))
      }else{
        copula_estimates %<>% rlist::list.append(fitted_copula$value)
        copula_aics %<>% c(fitted_copula$value$AIC)
        copula_names %<>% c(paste(fitted_copula$value$cop_name,rot))
        copula_stats %<>% bind_rows(tibble(Name = cop_name, Theta = fitted_copula$value$theta, Delta = fitted_copula$value$delta,
                                           Rotation = rot, AIC = fitted_copula$value$AIC, Time = fitted_copula$sys.time[3]))
      }
    }
  }
  cop <- copula_estimates[[which.min(copula_aics)]]
  cop$tau <- tau
  copula_stats %<>% arrange(AIC)
  return(list(Copula = cop,
         Copulas = copula_names,
         AICs = copula_aics, 
         Stats = copula_stats))
}

######################################################################################################################

censor_copula_select.intern <- function(u,v, indeptest = FALSE, level = 0.05, method = "censor",
                                        include_tawn = TRUE, include_amh = TRUE, include_t = TRUE){
  tau <- cor(u, v, method = "kendall")
  N <- length(u)
  if(indeptest){
    f <- sqrt((9 * N * (N - 1))/(2 * (2 * N + 5))) * abs(tau)
    p.value = 2 * (1 - pnorm(f))
    if(p.value>level){
      return(get_cop("Independence"))
    }
  }
  
  copula_estimates <- list()
  copula_aics <- c()
  ## Radially symmetric copulas
  copulas <- c("Frank", "normal", "t")
  if(!include_t){
    copulas <- copulas[-3]
  }
  for(cop_name in copulas){
    fitted_copula <- try(censor_est(u,v,cop_name, rot = 1, method = method))
    if(class(fitted_copula) == "try-error"){
      print(cop_name)
      copula_estimates %<>% rlist::list.append(get_cop(cop_name))
      copula_aics %<>% c(Inf)
    }else{
      copula_estimates %<>% rlist::list.append(fitted_copula)
      copula_aics %<>% c(fitted_copula$AIC)
    }
  }
  ## 0 and 180 degree rotations
  if(tau>0){
    copulas <- c("AMH", "Clayton", "Gumbel", "Joe",
                 "BB1", "BB6", "BB7", "BB8", "Tawn", "Tawn2")
    if((tau>=1/3)| (!include_amh)){
      copulas <- copulas[-1]
    }
    if(!include_tawn){
      copulas <- copulas[!copulas %in% c("Tawn", "Tawn2")]
    }
  }else{
    copulas <- c("AMH", "Clayton")
    if((tau<=(5 - 8/log(2)) / 3)| (!include_amh)){
      copulas <- copulas[-1]
    }
  }
  for(cop_name in copulas){
    for(rot in c(1,3)){
      fitted_copula <- try(censor_est(u,v,cop_name, rot = rot, method = method))
      if(class(fitted_copula) == "try-error"){
        copula_estimates %<>% rlist::list.append(get_cop(cop_name))
        copula_aics %<>% c(Inf)
      }else{
        copula_estimates %<>% rlist::list.append(fitted_copula)
        copula_aics %<>% c(fitted_copula$AIC)
      }
    }
  }
  ## 90 and 270 degree rotations
  if(tau<0){
    copulas <- c("AMH", "Clayton", "Gumbel", "Joe",
                 "BB1", "BB6", "BB7", "BB8", "Tawn", "Tawn2")
    if((tau<=-1/3)| (!include_amh)){
      copulas <- copulas[-1]
    }
    if(!include_tawn){
      copulas <- copulas[!copulas %in% c("Tawn", "Tawn2")]
    }
  }else{
    copulas <- c("AMH", "Clayton")
    if((tau>=-(5 - 8/log(2)) / 3) | (!include_amh)){
      copulas <- copulas[-1]
    }
  }
  for(cop_name in copulas){
    for(rot in c(2,4)){
      fitted_copula <- try(censor_est(u,v,cop_name, rot = rot, method = method))
      if(class(fitted_copula) == "try-error"){
        copula_estimates %<>% rlist::list.append(get_cop(cop_name))
        copula_aics %<>% c(Inf)
      }else{
        copula_estimates %<>% rlist::list.append(fitted_copula)
        copula_aics %<>% c(fitted_copula$AIC)
      }
    }
  }
  cop <- copula_estimates[[which.min(copula_aics)]]
  cop$tau <- tau
  return(cop)
}

###########################################################################################
################################## GOF ####################################################
###########################################################################################

qecdf <- function(x_ecdf,u){
  #breaks <- unique(c(-0.001,0,sort(x),1))
  breaks <- c(-0.001,x_ecdf)
  labels <- breaks[-1]
  return(labels[cut(u, breaks = breaks, labels = labels, include.lowest = F, right = T)])
}

censor_gof_test <- function(u, v, cop, rot = 1, N = 1000){
  if(cop$cop_name == "Independence"){
    return(list(
      p.value.CvM = NA,
      p.value.KS = NA
    ))
  }
  rot <- check_rotations(rot)
  if(rot != cop$rotation){
    warning("Mismatch between specified rotation and the stored copula rotation")
  }
  # rot == 1 is 0 degrees, and no rotation
  if(rot == 2){
    # 90 degree
    u <- 1 - u
  }else if(rot == 3){
    # 180 degree
    u <- 1 - u
    v <- 1 - v
  }else if(rot == 4){
    # 270 degree
    v <- 1 - v
  }
  u_upper <- censor(u, method = "upper")
  u_lower <- censor(u, method = "lower")
  v_upper <- censor(v, method = "upper")
  v_lower <- censor(v, method = "lower")
  
  # From here, rotation has been handeled, I think
  n <- length(u)
  u_upper0 <- sort(u_upper)
  v_upper0 <- sort(v_upper)
  u_lower0 <- sort(u_lower)
  v_lower0 <- sort(v_lower)
  #ecdf_u <- c(0,unique(rank(sort(u_upper), ties.method = "max"))/(n+1),1)
  #ecdf_v <- c(0,unique(rank(sort(v_upper), ties.method = "max"))/(n+1),1)
  
  if(cop$cop_name %in% c("AMH", "Clayton", "Gumbel", "Joe", "Frank")){
    cop_fun <- cop$cop_fun
    sample_copula <- archmCopula(cop$cop_name, param = cop$theta)
    sample_fun <- function(n){
      return(rCopula(n, sample_copula))
    }
    censor_est_fun <- function(ui_upper, ui_lower, vi_upper, vi_lower, cop){
      return(censor_est.intern(ui_upper, ui_lower, vi_upper, vi_lower, cop, check.pars = F))
    }
  }else if(cop$cop_name %in% c("normal", "t")){
    cop_fun <- cop$cop_fun
    sample_copula <- ellipCopula(cop$cop_name, param = c(cop$theta), df = cop$delta)
    sample_fun <- function(n){
      return(rCopula(n, sample_copula))
    }
    censor_est_fun <- function(ui_upper, ui_lower, vi_upper, vi_lower, cop){
      return(censor_est.intern(ui_upper, ui_lower, vi_upper, vi_lower, cop, check.pars = F))
    }
  }else if(cop$cop_name %in% c("BB1", "BB6", "BB7", "BB8")){
    cop_fun <- cop$cop_fun
    sample_copula <- BiCop(family = cop$fam, par = cop$theta, par2 = cop$delta)
    sample_fun <- function(n){
      return(BiCopSim(n, obj = sample_copula))
    }
    censor_est_fun <- function(ui_upper, ui_lower, vi_upper, vi_lower, cop){
      return(censor_est.intern(ui_upper, ui_lower, vi_upper, vi_lower, cop, check.pars = F))
    }
  }else if(cop$cop_name %in% c("Tawn", "Tawn2")){
    cop_fun <- cop$cop_fun
    sample_copula <- BiCop(family = BiCopName(cop$cop_name), par = cop$theta, par2 = cop$delta)
    sample_fun <- function(n){
      samples <- BiCopSim(n, obj = sample_copula)
    }
    censor_est_fun <- function(ui_upper, ui_lower, vi_upper, vi_lower, cop){
      return(censor_est.intern.tawn(ui_upper, ui_lower, vi_upper, vi_lower, cop, check.pars = F))
    }
  }
  cvm = sum((C.n(cbind(u_upper, v_upper), cbind(u_upper, v_upper), ties.method = "max") -
               cop_fun(u_upper, v_upper, cop$theta, cop$delta))^2) # squared ???
  ks = max(abs(C.n(cbind(u_upper, v_upper), cbind(u_upper, v_upper), ties.method = "max") -
                 cop_fun(u_upper, v_upper, cop$theta, cop$delta)))
  cvm_n <- vector("numeric", N)
  ks_n <- vector("numeric", N)
  thetas <- vector("numeric", N)
  deltas <- vector("numeric", N)
  
  for(i in 1:N){
    samples <- sample_fun(n)
    # Here the ties are generated
    # ties.method does not really matter here,
    # as long as it is either max or min
    ui_upper <- u_upper0[rank(samples[,1], ties.method = "max")]
    vi_upper <- v_upper0[rank(samples[,2], ties.method = "max")]
    ui_lower <- u_lower0[rank(samples[,1], ties.method = "min")]
    vi_lower <- u_lower0[rank(samples[,2], ties.method = "min")]
    # ui <- qecdf(ecdf_u, samples[,1])
    # vi <- qecdf(ecdf_v, samples[,2])
    # ui_upper <- censor(ui, method = "upper")
    # vi_upper <- censor(vi, method = "upper")
    # ui_lower <- censor(ui, method = "lower")
    # vi_lower <- censor(vi, method = "lower")
    
    cop_n <- censor_est_fun(ui_upper, ui_lower, vi_upper, vi_lower, cop)
    cvm_n[i] = sum((C.n(cbind(ui_upper, vi_upper), cbind(ui_upper, vi_upper), ties.method = "max") -
                      cop_fun(ui_upper, vi_upper, cop_n$theta, cop_n$delta))^2) # squared ???
    ks_n[i] = max(abs(C.n(cbind(ui_upper, vi_upper), cbind(ui_upper, vi_upper), ties.method = "max") -
                        cop_fun(ui_upper, vi_upper, cop_n$theta, cop_n$delta)))
    thetas[i] = cop_n$theta
    deltas[i] = cop_n$delta
  }
  return(list(
    CvM =cvm,
    p.value.CvM = sum(cvm_n >=cvm)/N,
    CvM.values = cvm_n,
    KS = ks,
    p.value.KS = sum(ks_n >=ks)/N,
    KS.values = ks_n,
    Thetas = thetas,
    Deltas = deltas
  ))
}

##########################################################################################################################
################################### RVINE STRUCTURE #################################################################
##########################################################################################################################

# Rotations and pobs handeled outside the transformations

transform <- function(u, v, cop){
  rot <- cop$rotation
  # rot == 1 is 0 degrees, and no rotation
  if(rot == 2){
    # 90 degree
    return(cop$cop_fun(u = 1, v = v, theta = cop$theta, delta = cop$delta) -
             cop$cop_fun(u = 1 - u, v = v, theta = cop$theta, delta = cop$delta))
  }else if(rot == 3){
    # 180 degree
    return(1 + cop$cop_fun(u = 1 - u, v = 1 - v, theta = cop$theta, delta = cop$delta) -
             cop$cop_fun(u = 1, v = 1 - v, theta = cop$theta, delta = cop$delta) -
             cop$cop_fun(u = 1 - u, v = 1, theta = cop$theta, delta = cop$delta))
  }else if(rot == 4){
    # 270 degree
    return(cop$cop_fun(u = u, v = 1, theta = cop$theta, delta = cop$delta) -
             cop$cop_fun(u = u, v = 1 - v, theta = cop$theta, delta = cop$delta))
  }
  return(cop$cop_fun(u = u, v = v, theta = cop$theta, delta = cop$delta))
}

transform_u <- function(u, v, cop){
  rot <- cop$rotation
  # rot == 1 is 0 degrees, and no rotation
  if(rot == 2){
    # 90 degree
    return(cop$cop_du(u = 1 - u, v = v, theta = cop$theta, delta = cop$delta))
  }else if(rot == 3){
    # 180 degree
    return(1 - cop$cop_du(u = 1 - u, v = 1 - v, theta = cop$theta, delta = cop$delta))
  }else if(rot == 4){
    # 270 degree
    return(1 - cop$cop_du(u = u, v = 1 - v, theta = cop$theta, delta = cop$delta))
  }
  return(cop$cop_du(u = u, v = v, theta = cop$theta, delta = cop$delta))
}

transform_v <- function(u, v, cop){
  rot <- cop$rotation
  # rot == 1 is 0 degrees, and no rotation
  if(rot == 2){
    # 90 degree
    return(1 - cop$cop_dv(u = 1 - u, v = v, theta = cop$theta, delta = cop$delta))
  }else if(rot == 3){
    # 180 degree
    return(1 - cop$cop_dv(u = 1 - u, v = 1 - v, theta = cop$theta, delta = cop$delta))
  }else if(rot == 4){
    # 270 degree
    return(cop$cop_dv(u = u, v = 1 - v, theta = cop$theta, delta = cop$delta))
  }
  return(cop$cop_dv(u = u, v = v, theta = cop$theta, delta = cop$delta))
}

cop_name2BiCop <- function(cop){
  family <- cop$fam
  rot <- cop$rotation
  theta <- cop$theta
  delta <- cop$delta
  if(is.na(delta)){
    delta <- 0
  }
  if(is.na(theta)){
    theta <- 0
  }
  # 90 degree
  if(rot == 2){
    theta <- -theta
    delta <- -delta
    family <- family + 20
  }else if(rot == 3){
    family <- family + 10
  }else if(rot == 4){
    theta <- -theta
    delta <- -delta
    family <- family + 30
  }
  return(list(Family = family,
              Par = theta,
              Par2 = delta))
}


rvine_translation <- function(tree, names){
  tree_0 <- tree
  d <- length(tree) + 1
  R_matrix <- matrix(0, nrow = d, ncol = d)
  cop_matrix <- matrix(0, nrow = d, ncol = d)
  theta_matrix <- matrix(0, nrow = d, ncol = d)
  delta_matrix <- matrix(0, nrow = d, ncol = d)
  aic_matrix <- matrix(0, nrow = d, ncol = d)
  
  # starting from top level of the tree
  for(i in (d-1):2){
    var <- tree[[i]]$t.cond[i,1]
    var %<>% c(tree[[i]]$t.cond[i,2])
    biCop <- cop_name2BiCop(tree[[i]]$Copulas[[1]])
    copulas <- c(biCop$Family)
    thetas <- c(biCop$Par)
    deltas <- c(biCop$Par2)
    aics <- c(tree[[i]]$Copulas[[1]]$AIC)
    for(j in (i-1):1){
      edge <- ceiling(which(tree[[j]]$t.cond[j,] == var[1])/2)
      ind <- which(tree[[j]]$t.cond[j,(edge*2 - 1):(2*edge)] != var[1]) - 1
      
      var %<>% c(tree[[j]]$t.cond[j,2*edge - 1 + ind])
      biCop <- cop_name2BiCop(tree[[j]]$Copulas[[edge]])
      copulas %<>% c(biCop$Family)
      thetas %<>% c(biCop$Par)
      deltas %<>% c(biCop$Par2)
      aics %<>% c(tree[[j]]$Copulas[[edge]]$AIC)
      
      # Removing used entries
      tree[[j]]$t.cond <- matrix(tree[[j]]$t.cond[,-c(edge*2 - 1, edge*2)], nrow=j)
      tree[[j]]$Copulas[[edge]] <- NULL
    }
    R_matrix[(d+1-length(var)):d,d - i] <- var
    cop_matrix[(d + 2 -length(var)):d,d - i] <- copulas
    theta_matrix[(d + 2 -length(var)):d,d - i] <- thetas
    delta_matrix[(d + 2 -length(var)):d,d - i] <- deltas
    aic_matrix[(d + 2 -length(var)):d,d - i] <- aics
  }
  ind <- !(tree[[2]]$t.cond[2,] %in% diag(R_matrix))
  var <- rev(tree[[2]]$t.cond[,ind])
  R_matrix[(d+1-length(var)):d,d - 1] <- var
  R_matrix[d,d] <- var[2]
  
  biCop <- cop_name2BiCop(tree[[1]]$Copulas[[1]])
  cop_matrix[d,d-1] <- biCop$Family
  theta_matrix[d,d-1] <- biCop$Par
  delta_matrix[d,d-1] <- biCop$Par2
  aic_matrix[d,d-1] <- tree[[1]]$Copulas[[1]]$AIC
  
  return(list(RVine = RVineMatrix(Matrix = R_matrix,
                     family = cop_matrix,
                     par = theta_matrix,
                     par2 = delta_matrix,
                     names = names),
              AICs = aic_matrix,
              Tree = tree_0))
}


censor_RVine_select <- function(data, indeptest = TRUE, level = 0.05, censor_level = 1,
                                include_tawn = TRUE, include_amh = FALSE, include_t = TRUE){
  # 0 censor_level corresponds to non-censored estimation
  d <- ncol(data)
  n <- nrow(data)
  tree <- list()
  copulas <- list()
  
  tau <- cor(data, method = "kendall")
  weights <- tau[upper.tri(tau)]
  sources <- matrix(rep(1:d, d), ncol = d)
  sources <- sources[upper.tri(sources)]
  destinations <- matrix(rep(1:d, d), ncol = d, byrow = TRUE)
  destinations <- destinations[upper.tri(destinations)]
  totals <- cbind(abs(weights), sources, destinations)
  totals <- totals[order(totals[,1], decreasing = TRUE),]
  graph <- graph_from_edgelist(cbind(as.character(totals[,2]), as.character(totals[,3])), directed = FALSE)
  E(graph)$weight <- -as.numeric(totals[,1])
  graph <- minimum.spanning.tree(graph)
  #plot(graph, edge.label=round(E(graph)$weight, 3))
  
  ################### transforming ########################
  # using matrix, not list, to remove data fast ...
  transformed <- matrix(0, nrow=n, ncol = 2*length(E(graph)))
  transformed_cond <- matrix(0, nrow=1, ncol = 2*(d-1))
  if(censor_level < 1){
    estimation_method <- "regular"
  }else{
    estimation_method <- "censor"
  }
  for(j in 1:length(E(graph))){
    ind1 <- as.numeric(tail_of(graph, j)$name)
    ind2 <- as.numeric(head_of(graph, j)$name)
    ui <- data[,ind1]#pobs(data[,ind1])
    vi <- data[,ind2]#pobs(data[,ind2])
    cop <- censor_copula_select.intern(ui, vi,
                                       indeptest = indeptest, level = level, method = estimation_method,
                                       include_tawn = include_tawn, include_amh = include_amh, include_t = include_t)
    #debug1 <- cop$cop_name
    transformed[,(j*2-1)] <- transform_u(u = ui, v = vi, cop = cop)
    transformed[,(j*2)] <- transform_v(u = ui, v = vi, cop = cop)
    transformed_cond[1,(j*2-1)] <- ind1
    transformed_cond[1,(j*2)] <- ind2
    copulas %<>% rlist::list.append(Copula = cop)
  }
  tree %<>% rlist::list.append(list(Copulas = copulas, Graph = graph,
                             t.data = transformed, t.cond = transformed_cond))
  
  ############## Transformation on higher levels
  for(i in 2:(d-1)){
    # Building structure for the next full tree
    tau <- cor(transformed, method = "kendall")
    sources <- c()
    destinations <- c()
    weights <- c()
    vertex_data <- matrix(nrow=n)
    vertex_cond <- matrix(nrow=i-1)
    vertex_var <- c()
    for(j in 1:(length(E(graph))-1)){
      for(k in (j+1):length(E(graph))){
        # Head is the second variable of "ends"
        first_ends <- ends(graph, j, names = FALSE)
        sec_ends <- ends(graph, k, names = FALSE)
        if(any(first_ends %in% sec_ends)){
          sources %<>% c(paste(ends(graph, j, names = TRUE), collapse = ","))
          destinations %<>% c(paste(ends(graph, k, names = TRUE), collapse = ","))
          
          if(i>=3){
            # First involved variables
            first_involved <- unique(c(transformed_cond[1:i-1,(j*2-1):(j*2)]))
            second_involved <- unique(c(transformed_cond[1:i-1,(k*2-1):(k*2)]))
            conditioned <- unique(first_involved[first_involved %in% second_involved])
            cond1 <- which((transformed_cond[i-1,(j*2-1):(j*2)] %in% conditioned)) +j*2-2
            cond2 <- which((transformed_cond[i-1,(k*2-1):(k*2)] %in% conditioned)) +k*2-2
            weights %<>% c(tau[cond1, cond2])
            
            vertex_data %<>% cbind(transformed[,cond1], transformed[,cond2])
            vertex_cond %<>% cbind(conditioned, conditioned)
            ind1 <- which(!(transformed_cond[i-1,(j*2-1):(j*2)] %in% conditioned)) +j*2-2
            ind2 <- which(!(transformed_cond[i-1,(k*2-1):(k*2)] %in% conditioned)) +k*2-2
            vertex_var %<>% c(transformed_cond[i-1,ind1] , transformed_cond[i-1,ind2])
          }else{
            ind1 <- which(first_ends %in% sec_ends) + j*2 - 2
            ind2 <- which(sec_ends %in% first_ends) + k*2 - 2
            cond1 <- which(!(first_ends %in% sec_ends)) + j*2 - 2 
            cond2 <- which(!(sec_ends %in% first_ends)) + k*2 - 2
            
            weights %<>% c(tau[ind1, ind2])
            vertex_data %<>% cbind(transformed[,ind1], transformed[,ind2])
            
            vertex_cond %<>% cbind(transformed_cond[,ind1], transformed_cond[,ind2])
            vertex_var %<>% c(transformed_cond[i-1,cond1], transformed_cond[i-1,cond2])
          }
        }
      }
    }
    vertex_data <- vertex_data[,-1]
    vertex_cond <- matrix(vertex_cond[,-1], nrow=i-1)
    graph <- graph_from_edgelist(cbind(sources, destinations), directed = FALSE)
    E(graph)$weight <- -abs(weights)
    new_graph <- minimum.spanning.tree(graph)
    #plot(graph, edge.label=round(E(graph)$weight, 3))
    #plot(new_graph, edge.label=round(E(new_graph)$weight, 3))
    
    ########################################################
    ####### tedious way of finding differing edges #########
    if(length(E(graph)) != length(E(new_graph))){
      graph_names <- vector("character", length = length(E(graph)))
      for(j in 1:length(E(graph))){
        graph_names[j] <- paste(ends(graph,j), collapse = " ")
      }
      new_graph_names <- vector("character", length = length(E(new_graph)))
      for(j in 1:length(E(new_graph))){
        new_graph_names[j] <- paste(ends(new_graph,j), collapse = " ")
      }
      data_removal <- which(!(graph_names %in% new_graph_names))
      vertex_data <- vertex_data[,-c(data_removal*2 - 1, data_removal*2)]
      vertex_cond <- vertex_cond[,-c(data_removal*2 - 1, data_removal*2)]
      vertex_var <- vertex_var[-c(data_removal*2 - 1, data_removal*2)]
    }
    #########################################################
    #########################################################
    
    if(i > censor_level){
      estimation_method <- "regular"
    }
    graph <- new_graph
    copulas <- list()
    transformed <- matrix(0, nrow = n, ncol = 2*length(E(graph)))
    for(j in 1:length(E(graph))){
      ind1 <- j*2 - 1
      ind2 <- j*2
      ui <- vertex_data[,ind1] #pobs(vertex_data[,ind1])
      vi <- vertex_data[,ind2]#pobs(vertex_data[,ind2])
      cop <- censor_copula_select.intern(ui, vi,
                                         indeptest = indeptest, level = level, method = estimation_method,
                                         include_tawn = include_tawn, 
                                         include_amh = include_amh,
                                         include_t = include_t)
      transformed[,ind1] <- transform_u(u = ui, v = vi, cop = cop)
      transformed[,ind2] <- transform_v(u = ui, v = vi, cop = cop)
      copulas %<>% rlist::list.append(Copula = cop)
    }
    # Kan flyttes høyere???
    transformed_cond <- rbind(vertex_cond, vertex_var)
    tree %<>% rlist::list.append(list(Copulas = copulas, Graph = graph,
                               t.data = transformed, t.cond = transformed_cond))
  }
  return(rvine_translation(tree = tree, names = colnames(data)))
}

