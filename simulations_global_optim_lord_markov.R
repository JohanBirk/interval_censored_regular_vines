library(tidyverse)
library(magrittr)
library(rlist)
library(igraph)
library(VineCopula)
library(copula)
library(doParallel)
library(foreach)

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

###################################### LOGLIK ###############################################

log_lik_censor <- function(u, v, cop, rot=1, check.pars = FALSE){
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

censor_est_sim <- function(u, v, cop, rot = 1, method = "censor", check.pars = FALSE){
  # censor or regular method
  if(!is.list(cop)){
    cop <- get_cop(cop)
    cop$rotation <- 1
    
  }else if(7>sum(names(cop) == c("cop_name", "theta", "delta","cop_fun","cop_du","cop_dv","cop_dens"))){
    stop("copula mismatch or something")
  }
  
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
      #cop$theta <- theta
      #cop_dens(u = u_upper[places1],v = v_upper[places1], theta, delta, check.pars = check.pars)
      return(sum(log(
        cop$cop_dens(u = u, v = v, theta = theta, delta = NULL, check.pars = check.pars)
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
    #cop$theta <- NA
    
    ## Using inverse kendall to find better parameters
    optimlist <- list()
    objectives <- c()
    delta <- 4
    inv_tau <- VineCopula::BiCopTau2Par(family = cop$fam, cor(u, v, method = "kendall"))
    for(i in 1:length(cop$optim.limits)){
      cop$theta <- cop$optim.limits[[i]]$start
      optimlist %<>% rlist::list.append(optimize(f = ll, maximum = TRUE,
                                                 interval = c(max(cop$optim.limits[[i]]$low, inv_tau - delta),
                                                              min(cop$optim.limits[[i]]$up, inv_tau + delta))))
      objectives %<>% c(optimlist[[i]]$objective)
    }
    if(sum(is.na(objectives)) == length(objectives)){
      #stop("All NA produced by optimize")
      cop$theta <- NA
    }else{
      cop$theta <- optimlist[[which.max(objectives)]]$maximum
      cop$log.lik <- optimlist[[which.max(objectives)]]$objective
      cop$AIC <- 2*cop$n.param - 2*optimlist[[which.max(objectives)]]$objective
    }
  }else{
    cop$theta <- optimlist[[which.max(objectives)]]$maximum
    cop$log.lik <- optimlist[[which.max(objectives)]]$objective
    cop$AIC <- 2*cop$n.param - 2*optimlist[[which.max(objectives)]]$objective
  }
  return(cop)
}

################################################################################

simulation_study_dvine <- function(data, est_cens = TRUE, optim_global = FALSE, delta = 0.2){
  N <- nrow(data)
  if(est_cens){
    estimation_method <- "censor"
    data <- pobs(data)
  }else{
    estimation_method <- "regular"
  }
  copulas_names <- list(
    list("frank", "gumbel", "clayton"),
    list("clayton", "joe"),
    list("normal")
  )
  ## initializing
  tree <- list()
  copulas <- list()
  cond_data <- matrix(ncol = 6, nrow=N)
  for(i in 1:3){
    data1 <- data[,i]
    data2 <- data[,i+1]
    cop <- censor_est_sim(data1, data2, cop = copulas_names[[1]][[i]], method = estimation_method)
    copulas %<>% list.append(cop)
    cond_data[,c(i*2-1,i*2)] <- cbind(data1,data2)
  }
  tree %<>% list.append(list(Copulas = copulas, CData = cond_data))
  
  
  for(i in 2:3){
    copulas <- list()
    cond_data <- matrix(ncol = (4-i)*2, nrow = N)
    for(j in 1:(4-i)){
      data1 <- transform_v(u = tree[[i-1]]$CData[,j*2-1], v = tree[[i-1]]$CData[,j*2], 
                           cop = tree[[i-1]]$Copulas[[j]])
      data2 <- transform_u(u = tree[[i-1]]$CData[,(j+1)*2-1], v = tree[[i-1]]$CData[,(j+1)*2], 
                           cop =tree[[i-1]]$Copulas[[j+1]])
      cop <- censor_est_sim(u = data1, v = data2, cop = copulas_names[[i]][[j]], method = estimation_method)
      copulas %<>% list.append(cop)
      cond_data[,c(j*2-1,j*2)] <- cbind(data1,data2)
    }
    tree %<>% list.append(list(Copulas = copulas, CData = cond_data))
  }
  
  
  if(!optim_global){
    ## returning the sequential result
    result <- c()
    for(i in 1:3){
      for(j in 1:(4-i)){
        result %<>% c(tree[[i]]$Copulas[[j]]$theta)
      }
    }
    result <- matrix(result, nrow = 1)
    colnames(result) <- c("Frank1", "Gumbel1", "Clayton1", "Clayton2", "Joe2", "Normal3")
    return(bind_cols(as_tibble(result), Fail_global = 0))
  }else{
    ## Joint optimization of full model
    fail_global <- 0
    #print(tree[[3]]$Copulas[[1]]$theta)
    ll_g <- function(par){
      logliks <- 0
      #print(par)
      cond_data <- matrix(ncol = 6, nrow=N)
      for(i in 1:3){
        data1 <- data[,i]
        data2 <- data[,i+1]
        
        tree[[1]]$Copulas[[i]]$theta <- par[i]
        if(est_cens){
          logliks <- logliks + log_lik_censor(data1, data2, tree[[1]]$Copulas[[i]], rot = 1)
        }else{
          logliks <- logliks + sum(log(
            tree[[1]]$Copulas[[i]]$cop_dens(u = data1,v = data2, theta = tree[[1]]$Copulas[[i]]$theta, check.pars = FALSE)
          ))
        }
        
        
        cond_data[,c(i*2-1,i*2)] <- cbind(data1,data2)
      }
      
      theta_count <- 4
      for(i in 2:3){
        cond_data_old <- cond_data
        cond_data <- matrix(ncol = (4-i)*2, nrow = N)
        for(j in 1:(4-i)){
          data1 <- transform_v(u = cond_data_old[,j*2-1], v = cond_data_old[,j*2], 
                               cop = tree[[i-1]]$Copulas[[j]])
          data2 <- transform_u(u = cond_data_old[,(j+1)*2-1], v = cond_data_old[,(j+1)*2], 
                               cop =tree[[i-1]]$Copulas[[j+1]])
          
          tree[[i]]$Copulas[[j]]$theta <- par[theta_count]
          theta_count <- theta_count + 1
          if(est_cens){
            logliks <- logliks + log_lik_censor(data1, data2, tree[[i]]$Copulas[[j]], rot = 1)
          }else{
            logliks <- logliks +sum(log(
              tree[[i]]$Copulas[[j]]$cop_dens(u = data1,v = data2, theta = tree[[i]]$Copulas[[j]]$theta, check.pars = FALSE)
            ))
          }
          
          cond_data[,c(j*2-1,j*2)] <- cbind(data1,data2)
        }
      }
      return(logliks)
    }
    
    if(tree[[1]]$Copulas[[1]]$theta > 0){
      tree[[1]]$Copulas[[1]]$optim.limits <- list(tree[[1]]$Copulas[[1]]$optim.limits[[2]])
    }else{
      tree[[1]]$Copulas[[1]]$optim.limits <- list(tree[[1]]$Copulas[[1]]$optim.limits[[1]])
    }
    
    ## Preparation of global optimization
    start_params <- c()
    upper_lims <- c()
    lower_lims <- c()
    for(i in 1:3){
      for(j in 1:(4-i)){
        cop <- tree[[i]]$Copulas[[j]]
        start_params %<>% c(cop$theta)
        lim_range <- cop$optim.limits[[1]]$up - cop$optim.limits[[1]]$low
        upper_lims %<>% c(min(cop$optim.limits[[1]]$up, cop$theta + delta*lim_range))
        lower_lims %<>% c(max(cop$optim.limits[[1]]$low, cop$theta - delta*lim_range))
      }
    }
    
    optimout <- try(optim(par = start_params, fn = ll_g,
                          method = "L-BFGS-B", 
                          lower = lower_lims, upper = upper_lims,
                          control = list(fnscale = -1, maxit = 500)))
    if("try-error" %in% class(optimout)){
      fail_global <- 1
      upper_lims <- c()
      lower_lims <- c()
      delta <- 0.025
      for(i in 1:3){
        for(j in 1:(4-i)){
          cop <- tree[[i]]$Copulas[[j]]
          lim_range <- cop$optim.limits[[1]]$up - cop$optim.limits[[1]]$low
          upper_lims %<>% c(min(cop$optim.limits[[1]]$up, cop$theta + delta*lim_range))
          lower_lims %<>% c(max(cop$optim.limits[[1]]$low, cop$theta - delta*lim_range))
        }
      }
      optimout <- optim(par = start_params, fn = ll_g,
                        method = "L-BFGS-B", 
                        lower = lower_lims, upper = upper_lims,
                        control = list(fnscale = -1, maxit = 100))
    }
    result <- matrix(optimout$par, nrow = 1)
    colnames(result) <- c("Frank1", "Gumbel1", "Clayton1", "Clayton2", "Joe2", "Normal3")
    return(bind_cols(as_tibble(result), Fail_global = fail_global))
  }
}

censor_est_sim_full_vine <- function(u_upper, u_lower, v_upper, v_lower, cop, check.pars = FALSE){
  # censor or regular method
  if(!is.list(cop)){
    cop <- get_cop(cop)
    cop$rotation <- 1
    
  }else if(7>sum(names(cop) == c("cop_name", "theta", "delta","cop_fun","cop_du","cop_dv","cop_dens"))){
    stop("copula mismatch or something")
  }
  
  ll <- function(theta){
    cop$theta <- theta
    return(log_lik_censor.intern(u_upper, u_lower, v_upper, v_lower, cop, check.pars = check.pars))
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
    ## Using inverse kendall to find better parameters
    
    optimlist <- list()
    objectives <- c()
    delta <- 1.5
    inv_tau <- VineCopula::BiCopTau2Par(family = cop$fam, cor(u_lower, v_lower, method = "kendall"))
    for(i in 1:length(cop$optim.limits)){
      cop$theta <- cop$optim.limits[[i]]$start
      optimlist %<>% rlist::list.append(optimize(f = ll, maximum = TRUE,
                                                 interval = c(max(cop$optim.limits[[i]]$low, inv_tau - delta),
                                                              min(cop$optim.limits[[i]]$up, inv_tau + delta))))
      objectives %<>% c(optimlist[[i]]$objective)
    }
    if(sum(is.na(objectives)) == length(objectives)){
      #stop("All NA produced by optimize")
      cop$theta <- NA
    }else{
      cop$theta <- optimlist[[which.max(objectives)]]$maximum
      cop$log.lik <- optimlist[[which.max(objectives)]]$objective
      cop$AIC <- 2*cop$n.param - 2*optimlist[[which.max(objectives)]]$objective
    }
  }else{
    cop$theta <- optimlist[[which.max(objectives)]]$maximum
    cop$log.lik <- optimlist[[which.max(objectives)]]$objective
    cop$AIC <- 2*cop$n.param - 2*optimlist[[which.max(objectives)]]$objective
  }
  return(cop)
}

simulation_study_dvine_censor_full <- function(data, optim_global = FALSE, delta = 0.2){
  N <- nrow(data)
  data_upper <- pobs(data, ties.method = "max")
  data_lower <- pobs(data, ties.method = "min")
  copulas_names <- list(
    list("frank", "gumbel", "clayton"),
    list("clayton", "joe"),
    list("normal")
  )
  ## initializing
  tree <- list()
  copulas <- list()
  cond_data_upper <- matrix(ncol = 6, nrow=N)
  cond_data_lower <- matrix(ncol = 6, nrow=N)
  for(i in 1:3){
    data1_upper <- data_upper[,i]
    data1_lower <- data_lower[,i]
    data2_upper <- data_upper[,i+1]
    data2_lower <- data_lower[,i+1]
    cop <- censor_est_sim_full_vine(u_upper = data1_upper, u_lower = data1_lower, v_upper = data2_upper, v_lower = data2_lower, cop = copulas_names[[1]][[i]])
    copulas %<>% list.append(cop)
    cond_data_upper[,c(i*2-1,i*2)] <- cbind(data1_upper,data2_upper)
    cond_data_lower[,c(i*2-1,i*2)] <- cbind(data1_lower,data2_lower)
  }
  tree %<>% list.append(list(Copulas = copulas, CData_upper = cond_data_upper, CData_lower = cond_data_lower))
  
  ## From here, we have to ensure that the upper limits are larger than the lower
  for(i in 2:3){
    copulas <- list()
    cond_data_upper <- matrix(ncol = (4-i)*2, nrow = N)
    cond_data_lower <- matrix(ncol = (4-i)*2, nrow = N)
    for(j in 1:(4-i)){
      data1_upper <- transform_v(u = tree[[i-1]]$CData_upper[,j*2-1], v = tree[[i-1]]$CData_upper[,j*2], 
                           cop = tree[[i-1]]$Copulas[[j]])
      data1_lower <- transform_v(u = tree[[i-1]]$CData_lower[,j*2-1], v = tree[[i-1]]$CData_lower[,j*2], 
                                      cop = tree[[i-1]]$Copulas[[j]])
      
      data2_upper <- transform_u(u = tree[[i-1]]$CData_upper[,(j+1)*2-1], v = tree[[i-1]]$CData_upper[,(j+1)*2], 
                           cop =tree[[i-1]]$Copulas[[j+1]])
      data2_lower <- transform_u(u = tree[[i-1]]$CData_lower[,(j+1)*2-1], v = tree[[i-1]]$CData_lower[,(j+1)*2], 
                                      cop =tree[[i-1]]$Copulas[[j+1]])
      
      cop <- censor_est_sim_full_vine(u_upper = pmax(data1_upper, data1_lower),
                                      u_lower = pmin(data1_upper, data1_lower),
                                      v_upper = pmax(data2_upper, data2_lower),
                                      v_lower = pmin(data2_upper, data2_lower), cop = copulas_names[[i]][[j]])
      #cop <- censor_est_sim(u = data1, v = data2, cop = copulas_names[[i]][[j]], method = estimation_method)
      copulas %<>% list.append(cop)
      #cond_data[,c(j*2-1,j*2)] <- cbind(data1,data2)
      cond_data_upper[,c(j*2-1,j*2)] <- cbind(data1_upper,data2_upper)
      cond_data_lower[,c(j*2-1,j*2)] <- cbind(data1_lower,data2_lower)
    }
    tree %<>% list.append(list(Copulas = copulas,  CData_upper = cond_data_upper, CData_lower = cond_data_lower))
  }
  
  
  if(!optim_global){
    ## returning the sequential result
    result <- c()
    for(i in 1:3){
      for(j in 1:(4-i)){
        result %<>% c(tree[[i]]$Copulas[[j]]$theta)
      }
    }
    result <- matrix(result, nrow = 1)
    colnames(result) <- c("Frank1", "Gumbel1", "Clayton1", "Clayton2", "Joe2", "Normal3")
    return(bind_cols(as_tibble(result), Fail_global = 0))
  }else{
    ## Joint optimization of full model
    fail_global <- 0
    #print(tree[[3]]$Copulas[[1]]$theta)
    ll_g <- function(par){
      logliks <- 0
      #print(par)
      cond_data_upper <- matrix(ncol = 6, nrow=N)
      cond_data_lower <- matrix(ncol = 6, nrow=N)
      for(i in 1:3){
        data1_upper <- data[,i]
        data1_lower <- data[,i]
        data2_upper <- data[,i+1]
        data2_lower <- data[,i+1]
        
        tree[[1]]$Copulas[[i]]$theta <- par[i]
        logliks <- logliks + 
          log_lik_censor.intern(u_upper = data1_upper, u_lower = data1_lower, v_upper = data2_upper, v_lower = data2_lower, tree[[1]]$Copulas[[i]], check.pars = FALSE)
          #log_lik_censor.intern(data1, data2, tree[[1]]$Copulas[[i]], rot = 1)
        
        
        cond_data_upper[,c(i*2-1,i*2)] <- cbind(data1_upper,data2_upper)
        cond_data_lower[,c(i*2-1,i*2)] <- cbind(data1_lower,data2_lower)
      }
      
      theta_count <- 4
      for(i in 2:3){
        cond_data_old_upper <- cond_data_upper
        cond_data_old_lower <- cond_data_lower
        cond_data_upper <- matrix(ncol = (4-i)*2, nrow = N)
        cond_data_lower <- matrix(ncol = (4-i)*2, nrow = N)
        for(j in 1:(4-i)){
          data1_upper <- transform_v(u = cond_data_old_upper[,j*2-1], v = cond_data_old_upper[,j*2], 
                               cop = tree[[i-1]]$Copulas[[j]])
          data1_lower <- transform_v(u = cond_data_old_lower[,j*2-1], v = cond_data_old_lower[,j*2], 
                                     cop = tree[[i-1]]$Copulas[[j]])
          data2_upper <- transform_u(u = cond_data_old_upper[,(j+1)*2-1], v = cond_data_old_upper[,(j+1)*2], 
                               cop =tree[[i-1]]$Copulas[[j+1]])
          data2_lower <- transform_u(u = cond_data_old_lower[,(j+1)*2-1], v = cond_data_old_lower[,(j+1)*2], 
                                     cop =tree[[i-1]]$Copulas[[j+1]])
          
          tree[[i]]$Copulas[[j]]$theta <- par[theta_count]
          theta_count <- theta_count + 1
          logliks <- logliks +
            log_lik_censor.intern(u_upper = data1_upper, u_lower = data1_lower, v_upper = data2_upper, v_lower = data2_lower, tree[[i]]$Copulas[[j]], check.pars = FALSE)
            #log_lik_censor(data1, data2, tree[[i]]$Copulas[[j]], rot = 1)
          
          cond_data_upper[,c(j*2-1,j*2)] <- cbind(data1_upper,data2_upper)
          cond_data_lower[,c(j*2-1,j*2)] <- cbind(data1_lower,data2_lower)
        }
      }
      return(logliks)
    }
    
    if(tree[[1]]$Copulas[[1]]$theta > 0){
      tree[[1]]$Copulas[[1]]$optim.limits <- list(tree[[1]]$Copulas[[1]]$optim.limits[[2]])
    }else{
      tree[[1]]$Copulas[[1]]$optim.limits <- list(tree[[1]]$Copulas[[1]]$optim.limits[[1]])
    }
    
    ## Preparation of global optimization
    start_params <- c()
    upper_lims <- c()
    lower_lims <- c()
    for(i in 1:3){
      for(j in 1:(4-i)){
        cop <- tree[[i]]$Copulas[[j]]
        start_params %<>% c(cop$theta)
        lim_range <- cop$optim.limits[[1]]$up - cop$optim.limits[[1]]$low
        upper_lims %<>% c(min(cop$optim.limits[[1]]$up, cop$theta + delta*lim_range))
        lower_lims %<>% c(max(cop$optim.limits[[1]]$low, cop$theta - delta*lim_range))
      }
    }
    
    optimout <- try(optim(par = start_params, fn = ll_g,
                          method = "L-BFGS-B", 
                          lower = lower_lims, upper = upper_lims,
                          control = list(fnscale = -1, maxit = 500)))
    if("try-error" %in% class(optimout)){
      fail_global <- 1
      upper_lims <- c()
      lower_lims <- c()
      delta <- 0.025
      for(i in 1:3){
        for(j in 1:(4-i)){
          cop <- tree[[i]]$Copulas[[j]]
          lim_range <- cop$optim.limits[[1]]$up - cop$optim.limits[[1]]$low
          upper_lims %<>% c(min(cop$optim.limits[[1]]$up, cop$theta + delta*lim_range))
          lower_lims %<>% c(max(cop$optim.limits[[1]]$low, cop$theta - delta*lim_range))
        }
      }
      optimout <- optim(par = start_params, fn = ll_g,
                        method = "L-BFGS-B", 
                        lower = lower_lims, upper = upper_lims,
                        control = list(fnscale = -1, maxit = 100))
    }
    result <- matrix(optimout$par, nrow = 1)
    colnames(result) <- c("Frank1", "Gumbel1", "Clayton1", "Clayton2", "Joe2", "Normal3")
    return(bind_cols(as_tibble(result), Fail_global = fail_global))
  }
}

###########################################################################################################
round_lower <- function(x, per_ties, place = 1, N){
  x[rank(x) <= per_ties*N] <- round(x[rank(x) <= per_ties*N],place)
  return(x)
}

###########################################################################################################
total_trials <- 1200

sample_sizes <- c(500,1000,5000)
errorhandling_ <- c("stop", "remove")[1]

experiments <- c("symm", "asymm", "symmG", "asymmG")

tree_names <- c("Frank1", "Gumbel1", "Clayton1", "Clayton2", "Joe2", "Normal3")
failed_pars <- matrix(rep(NA,7), nrow = 1)
colnames(failed_pars) <- c(tree_names, "Fail_global")
failed_pars %<>% as_tibble

for(experiment in experiments){
  ## Full ties experiment
  if(T){
    print("Full ties")
    print(experiment)
    if(experiment %in% c("symmG", "asymmG")){
      do_global <- TRUE
    }else{
      do_global <- FALSE
    }
    numCores <- parallel::detectCores()
    doParallel::registerDoParallel(min(numCores,25))
    b <- 15
    library(tictoc)
    tic()
    foreach (seed = 1:total_trials, .inorder = FALSE, .combine = rbind, .errorhandling = errorhandling_,
             .packages = c("tidyverse","magrittr", "VineCopula")) %dopar% {
               print(seed)
               source(file = 'copula_cdfs.R', local = TRUE)
               
               sampling_stats_out <- dplyr::tibble()
               for(tau in seq(from= 0.3, to = 0.9, by = 0.1)){ 
                 Matrix <- c(1, 4, 3, 2,
                             0, 2, 4, 3,
                             0, 0, 3, 4,
                             0, 0, 0, 4)
                 Matrix <- matrix(Matrix, 4, 4)
                 # define R-vine pair-copula family matrix
                 family <- c(0, 1, 3, 5, 
                             0, 0, 6, 4,
                             0, 0, 0, 3,
                             0, 0, 0, 0)
                 family <- matrix(family, 4, 4)
                 # define R-vine pair-copula parameter matrix
                 par <- c(0, BiCopTau2Par(family = family[2,1], tau = tau-0.2), BiCopTau2Par(family = family[3,1], tau = tau-0.1), BiCopTau2Par(family = family[4,1], tau = tau),
                          0, 0, BiCopTau2Par(family = family[3,2], tau = tau-0.1), BiCopTau2Par(family = family[4,2], tau = tau-0.1),
                          0, 0, 0, BiCopTau2Par(family = family[4,3], tau = tau),
                          0, 0, 0, 0)
                 par <- matrix(par, 4, 4)
                 # define second R-vine pair-copula parameter matrix
                 par2 <- c(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0)
                 
                 par2 <- matrix(par2, 4, 4)
                 ## define RVineMatrix object
                 RVM <- RVineMatrix(Matrix = Matrix, family = family,
                                    par = par, par2 = par2,
                                    names = c("V1", "V2", "V3", "V4"))
                 
                 for(N in sample_sizes){
                   set.seed(100 + seed)
                   data <- RVineSim(N, RVM)
                   # experiments <- c("symm", "asymm", "symmG", "asymmG")
                   if(experiment %in% c("symm", "symmG")){
                     censor_data <- (apply(data, 2, cut, breaks = 0:b/b, labels = FALSE) - 0.5) / b
                   }else{
                     data1 <- data[,1]
                     data2 <- (cut(data[,2], breaks = 0:b/b, labels = FALSE) - 0.5) / b
                     data3 <- (cut(data[,3], breaks = 0:b/b, labels = FALSE) - 0.5) / b
                     data4 <- (cut(data[,4], breaks = 0:(2*b)/(2*b), labels = FALSE) - 0.5) / (2*b)
                     censor_data <- cbind(data1,data2, data3, data4)
                   }
                   
                   parRef <- try(simulation_study_dvine(pobs(data), est_cens = FALSE, optim_global = do_global))
                   if("try-error" %in% class(parRef)){
                     parRef <- failed_pars #rep(NA,6)
                     #names(parRef) <- tree_names
                   }
                   sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                     Samples = N,
                     Tau = tau,
                     Method = "Reference",
                     parRef
                   ))
                   
                   ## Average
                   parAvg <- try(simulation_study_dvine(pobs(censor_data), est_cens = FALSE, optim_global = do_global))
                   if("try-error" %in% class(parAvg)){
                     parAvg <- failed_pars
                     #names(parAvg) <- tree_names
                   }
                   sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                     Samples = N,
                     Tau = tau,
                     Method = "Average",
                     parAvg
                   ))
                   
                   ## Random
                   parsRand <- tibble()
                   for(i in 1:100){
                     set.seed(seed*100 + i)
                     
                     pars <- try(simulation_study_dvine(pobs(censor_data, ties.method = "random"), est_cens = FALSE, optim_global = do_global))
                     if("try-error" %in% class(pars)){
                       pars <- failed_pars
                     }
                     parsRand %<>% bind_rows(pars)
                   }
                   parRand <- parsRand %>% map_df(~mean(.,na.rm = TRUE)) #%>% as.numeric
                   #names(parRand) <- tree_names
                   sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                     Samples = N,
                     Tau = tau,
                     Method = "Random",
                     parRand
                   ))
                   
                   ## Censor
                   # For censoring, pseudo-observations are handeled inside the function
                   parCens <- try(simulation_study_dvine(censor_data, est_cens = TRUE, optim_global = do_global))
                   if("try-error" %in% class(parCens)){
                     parCens <- failed_pars#rep(NA,6)
                     #names(parCens) <- tree_names
                   }
                   
                   sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                     Samples = N,
                     Tau = tau,
                     Method = "Censor",
                     parCens
                   ))
                   
                   ## Censor full
                   # For censoring, pseudo-observations are handeled inside the function
                   parCens_full <- try(simulation_study_dvine_censor_full(censor_data, optim_global = do_global))
                   if("try-error" %in% class(parCens_full)){
                     parCens_full <- failed_pars#rep(NA,6)
                     #names(parCens) <- tree_names
                   }
                   
                   sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                     Samples = N,
                     Tau = tau,
                     Method = "CensorFull",
                     parCens_full
                   ))
                 }
               }
               sampling_stats_out
             } -> sampling_stats_vine
    print(sampling_stats_vine)
    toc()
    # experiments <- c("symm", "asymm", "symmG", "asymmG")
    if(experiment == "symm"){
      sampling_stats_vine_symm <- sampling_stats_vine
      save(sampling_stats_vine_symm, file = "stored_data/sampling_stats_vine_symm_lord.RData")
    }else if(experiment == "asymm"){
      sampling_stats_vine_asymm <- sampling_stats_vine
      save(sampling_stats_vine_asymm, file = "stored_data/sampling_stats_vine_asymm_lord.RData")
    }else if(experiment == "symmG"){
      sampling_stats_vine_symm_global <- sampling_stats_vine
      save(sampling_stats_vine_symm_global, file = "stored_data/sampling_stats_vine_symm_global_lord.RData")
    }else if(experiment == "asymmG"){
      sampling_stats_vine_asymm_global <- sampling_stats_vine
      save(sampling_stats_vine_asymm_global, file = "stored_data/sampling_stats_vine_asymm_global_lord.RData")
    }
  }
  
  ## Lower tail experiment
  if(T){
    print("Tail ties")
    print(experiment)
    if(experiment %in% c("symmG", "asymmG")){
      do_global <- TRUE
    }else{
      do_global <- FALSE
    }
    numCores <- parallel::detectCores()
    doParallel::registerDoParallel(min(numCores,25))
    b <- 15
    library(tictoc)
    tic()
    foreach (seed = 1:total_trials, .inorder = FALSE, .combine = rbind, .errorhandling = errorhandling_,
             .packages = c("tidyverse","magrittr", "VineCopula")) %dopar% {
               print(seed)
               source(file = 'copula_cdfs.R', local = TRUE)
               
               sampling_stats_out <- dplyr::tibble()
               for(tau in c(0.25, 0.75)){ 
                 Matrix <- c(1, 4, 3, 2,
                             0, 2, 4, 3,
                             0, 0, 3, 4,
                             0, 0, 0, 4)
                 Matrix <- matrix(Matrix, 4, 4)
                 # define R-vine pair-copula family matrix
                 family <- c(0, 1, 3, 5, 
                             0, 0, 6, 4,
                             0, 0, 0, 3,
                             0, 0, 0, 0)
                 family <- matrix(family, 4, 4)
                 # define R-vine pair-copula parameter matrix
                 par <- c(0, BiCopTau2Par(family = family[2,1], tau = tau-0.2), BiCopTau2Par(family = family[3,1], tau = tau-0.1), BiCopTau2Par(family = family[4,1], tau = tau),
                          0, 0, BiCopTau2Par(family = family[3,2], tau = tau-0.1), BiCopTau2Par(family = family[4,2], tau = tau-0.1),
                          0, 0, 0, BiCopTau2Par(family = family[4,3], tau = tau),
                          0, 0, 0, 0)
                 par <- matrix(par, 4, 4)
                 # define second R-vine pair-copula parameter matrix
                 par2 <- c(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0)
                 
                 par2 <- matrix(par2, 4, 4)
                 ## define RVineMatrix object
                 RVM <- RVineMatrix(Matrix = Matrix, family = family,
                                    par = par, par2 = par2,
                                    names = c("V1", "V2", "V3", "V4"))
                 
                 for(N in sample_sizes){
                   for(per_ties in seq(from = 0.1, to = 0.5, by = 0.1)){
                     set.seed(100 + seed)
                     data <- RVineSim(N, RVM)
                     # experiments <- c("symm", "asymm", "symmG", "asymmG")
                     if(experiment %in% c("symm", "symmG")){
                       censor_data <- apply(data,2,round_lower,per_ties = per_ties, place = 1, N = N)
                     }else{
                       data1 <- data[,1]
                       data2 <- round_lower(data[,2],per_ties = per_ties, place = 1, N = N)
                       data3 <- round_lower(data[,3],per_ties = per_ties, place = 1, N = N)
                       data4 <- round_lower(data[,4],per_ties = per_ties, place = 2, N = N)
                       censor_data <- cbind(data1,data2, data3, data4)
                     }
                     
                     parRef <- try(simulation_study_dvine(pobs(data), est_cens = FALSE, optim_global = do_global))
                     if("try-error" %in% class(parRef)){
                       parRef <- failed_pars #rep(NA,6)
                       #names(parRef) <- tree_names
                     }
                     sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                       Samples = N,
                       Tau = tau,
                       Ties = per_ties,
                       Method = "Reference",
                       parRef
                     ))
                     
                     ## Average
                     parAvg <- try(simulation_study_dvine(pobs(censor_data), est_cens = FALSE, optim_global = do_global))
                     if("try-error" %in% class(parAvg)){
                       parAvg <- failed_pars
                       #names(parAvg) <- tree_names
                     }
                     sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                       Samples = N,
                       Tau = tau,
                       Ties = per_ties,
                       Method = "Average",
                       parAvg
                     ))
                     
                     ## Random
                     parsRand <- tibble()
                     for(i in 1:100){
                       set.seed(seed*100 + i)
                       
                       pars <- try(simulation_study_dvine(pobs(censor_data, ties.method = "random"), est_cens = FALSE, optim_global = do_global))
                       if("try-error" %in% class(pars)){
                         pars <- failed_pars
                       }
                       parsRand %<>% bind_rows(pars)
                     }
                     parRand <- parsRand %>% map_df(~mean(.,na.rm = TRUE)) #%>% as.numeric
                     #names(parRand) <- tree_names
                     sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                       Samples = N,
                       Tau = tau,
                       Ties = per_ties,
                       Method = "Random",
                       parRand
                     ))
                     
                     ## Censor
                     # For censoring, pseudo-observations are handeled inside the function
                     parCens <- try(simulation_study_dvine(censor_data, est_cens = TRUE, optim_global = do_global))
                     if("try-error" %in% class(parCens)){
                       parCens <- failed_pars#rep(NA,6)
                       #names(parCens) <- tree_names
                     }
                     
                     sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                       Samples = N,
                       Tau = tau,
                       Ties = per_ties,
                       Method = "Censor",
                       parCens
                     ))
                     
                     ## Censor full
                     # For censoring, pseudo-observations are handeled inside the function
                     parCens_full <- try(simulation_study_dvine_censor_full(censor_data, optim_global = do_global))
                     if("try-error" %in% class(parCens_full)){
                       parCens_full <- failed_pars#rep(NA,6)
                       #names(parCens) <- tree_names
                     }
                     
                     sampling_stats_out %<>% dplyr::bind_rows(bind_cols(
                       Samples = N,
                       Tau = tau,
                       Ties = per_ties,
                       Method = "CensorFull",
                       parCens_full
                     ))
                   }
                 }
               }
               sampling_stats_out
             } -> sampling_stats_tail_vine
    print(sampling_stats_tail_vine)
    toc()
    # experiments <- c("symm", "asymm", "symmG", "asymmG")
    if(experiment == "symm"){
      sampling_stats_tail_vine_symm <- sampling_stats_tail_vine
      save(sampling_stats_tail_vine_symm, file = "stored_data/sampling_stats_tail_vine_symm_lord.RData")
    }else if(experiment == "asymm"){
      sampling_stats_tail_vine_asymm <- sampling_stats_tail_vine
      save(sampling_stats_tail_vine_asymm, file = "stored_data/sampling_stats_tail_vine_asymm_lord.RData")
    }else if(experiment == "symmG"){
      sampling_stats_tail_vine_symm_global <- sampling_stats_tail_vine
      save(sampling_stats_tail_vine_symm_global, file = "stored_data/sampling_stats_tail_vine_symm_global_lord.RData")
    }else if(experiment == "asymmG"){
      sampling_stats_tail_vine_asymm_global <- sampling_stats_tail_vine
      save(sampling_stats_tail_vine_asymm_global, file = "stored_data/sampling_stats_tail_vine_asymm_global_lord.RData")
    }
  }
  
}

####

print(sampling_stats_vine)