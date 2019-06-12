library(tidyverse)
library(magrittr)
library(rlist)
library(igraph)
library(VineCopula)
library(copula)
library(doParallel)
library(foreach)
source(file = 'censor_est.R', local = TRUE)
source(file = 'copula_cdfs.R', local = TRUE)
source(file = "censor_gof_parallel.R")

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


censor_est.full <- function(u_upper, u_lower, v_upper, v_lower, cop_name, check.pars = FALSE){
  ## This function receives rotated data, and estimates the copula. Intended for use with fully censored R-vine.
  cop <- get_cop(cop_name)
  ## When conditional data are computed, the data computed from max rank estimated marginals are not always
  ## larger than when computed by min rank estimated marginals.
  u_up <- pmax(u_upper, u_upper)
  u_lo <- pmin(u_upper, u_lower)
  v_up <- pmax(v_upper, v_upper)
  v_lo <- pmin(v_upper, v_lower)
  
  if(cop$n.param == 1){
    ll <- function(theta){
      cop$theta <- theta
      return(log_lik_censor.intern(u_upper = u_up,u_lower = u_lo,
                                   v_upper = v_up, v_lower = v_lo, cop = cop, check.pars = check.pars))
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
      warning("All NA produced by optimize. Narrowing the parameter search.")
      ## Using inverse kendall to find better parameters
      
      optimlist <- list()
      objectives <- c()
      lim_adjust <- 3
      inv_tau <- VineCopula::BiCopTau2Par(family = cop$fam, cor(u_lower, v_lower, method = "kendall"))
      for(i in 1:length(cop$optim.limits)){
        cop$theta <- cop$optim.limits[[i]]$start
        optimlist %<>% rlist::list.append(optimize(f = ll, maximum = TRUE,
                                                   interval = c(max(cop$optim.limits[[i]]$low, inv_tau - lim_adjust),
                                                                min(cop$optim.limits[[i]]$up, inv_tau + lim_adjust))))
        objectives %<>% c(optimlist[[i]]$objective)
      }
      if(sum(is.na(objectives)) == length(objectives)){
        stop("All NA produced by optimize after narrowed search.")
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
  }else{
    cop$theta <- cop$start.value[1]
    cop$delta <- cop$start.value[2]
    parlower <- cop$optim.lower
    parupper <- cop$optim.upper
    
    ll <- function(par){
      cop$theta <- par[1]
      cop$delta <- par[2]
      return(log_lik_censor.intern(u_upper = u_up,u_lower = u_lo,
                                   v_upper = v_up, v_lower = v_lo,cop = cop, check.pars = check.pars))
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


censor_copula_select.full <- function(u_upper, u_lower, v_upper, v_lower, indeptest = FALSE, level = 0.05,
                                        include_tawn = TRUE, include_amh = TRUE, include_t = TRUE){
  tau <- cor(u_upper, v_upper, method = "kendall")
  N <- length(u_upper)
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
    fitted_copula <- try(censor_est.full(u_upper, u_lower, v_upper, v_lower, cop_name))
    if(class(fitted_copula) == "try-error"){
      print(cop_name)
      copula_estimates %<>% rlist::list.append(get_cop(cop_name))
      copula_aics %<>% c(Inf)
    }else{
      fitted_copula$rotation <- 1
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
      if(rot == 1){
        ## no rotation
        fitted_copula <- try(censor_est.full(u_upper = u_upper, u_lower = u_lower, v_upper = v_upper, v_lower = v_lower, cop_name))
      }else if(rot == 3){
        ## 180 degrees 1 - v, 1 - u
        fitted_copula <- try(censor_est.full(u_upper = 1 - u_lower, u_lower = 1 - u_upper, v_upper = 1 - v_lower, v_lower = 1 - v_upper, cop_name))
      }
      if(class(fitted_copula) == "try-error"){
        copula_estimates %<>% rlist::list.append(get_cop(cop_name))
        copula_aics %<>% c(Inf)
      }else{
        fitted_copula$rotation <- rot
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
      if(rot == 2){
        ## 90 degrees 1 - u
        fitted_copula <- try(censor_est.full(u_upper = 1 - u_lower, u_lower = 1 - u_upper, v_upper = v_upper, v_lower = v_lower, cop_name))
      }else if(rot == 4){
        ## 270 degrees 1 - v
        fitted_copula <- try(censor_est.full(u_upper = u_upper, u_lower = u_lower, v_upper = 1 - v_lower, v_lower = 1 - v_upper, cop_name))
      }
      if(class(fitted_copula) == "try-error"){
        copula_estimates %<>% rlist::list.append(get_cop(cop_name))
        copula_aics %<>% c(Inf)
      }else{
        fitted_copula$rotation <- rot
        copula_estimates %<>% rlist::list.append(fitted_copula)
        copula_aics %<>% c(fitted_copula$AIC)
      }
    }
  }
  cop <- copula_estimates[[which.min(copula_aics)]]
  cop$tau <- tau
  return(cop)
}

rvine_translation_pval <- function(tree, names){
  tree_0 <- tree
  d <- length(tree) + 1
  R_matrix <- matrix(0, nrow = d, ncol = d)
  cop_matrix <- matrix(0, nrow = d, ncol = d)
  theta_matrix <- matrix(0, nrow = d, ncol = d)
  delta_matrix <- matrix(0, nrow = d, ncol = d)
  aic_matrix <- matrix(0, nrow = d, ncol = d)
  cvm_matrix <- matrix(0, nrow = d, ncol = d)
  ks_matrix <- matrix(0, nrow = d, ncol = d)
  
  # starting from top level of the tree
  for(i in (d-1):2){
    var <- tree[[i]]$t.cond[i,1]
    var %<>% c(tree[[i]]$t.cond[i,2])
    biCop <- cop_name2BiCop(tree[[i]]$Copulas[[1]])
    copulas <- c(biCop$Family)
    thetas <- c(biCop$Par)
    deltas <- c(biCop$Par2)
    aics <- c(tree[[i]]$Copulas[[1]]$AIC)
    cvms <- c(tree[[i]]$Copulas[[1]]$p.value.CvM)
    kss <- c(tree[[i]]$Copulas[[1]]$p.value.KS)
    for(j in (i-1):1){
      edge <- ceiling(which(tree[[j]]$t.cond[j,] == var[1])/2)
      ind <- which(tree[[j]]$t.cond[j,(edge*2 - 1):(2*edge)] != var[1]) - 1
      
      var %<>% c(tree[[j]]$t.cond[j,2*edge - 1 + ind])
      biCop <- cop_name2BiCop(tree[[j]]$Copulas[[edge]])
      copulas %<>% c(biCop$Family)
      thetas %<>% c(biCop$Par)
      deltas %<>% c(biCop$Par2)
      aics %<>% c(tree[[j]]$Copulas[[edge]]$AIC)
      cvms %<>% c(tree[[j]]$Copulas[[edge]]$p.value.CvM)
      kss %<>% c(tree[[j]]$Copulas[[edge]]$p.value.KS)
      
      # Removing used entries
      tree[[j]]$t.cond <- matrix(tree[[j]]$t.cond[,-c(edge*2 - 1, edge*2)], nrow=j)
      tree[[j]]$Copulas[[edge]] <- NULL
    }
    R_matrix[(d+1-length(var)):d,d - i] <- var
    cop_matrix[(d + 2 -length(var)):d,d - i] <- copulas
    theta_matrix[(d + 2 -length(var)):d,d - i] <- thetas
    delta_matrix[(d + 2 -length(var)):d,d - i] <- deltas
    aic_matrix[(d + 2 -length(var)):d,d - i] <- aics
    cvm_matrix[(d + 2 -length(var)):d,d - i] <- cvms
    ks_matrix[(d + 2 -length(var)):d,d - i] <- kss
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
  cvm_matrix[d,d-1] <- tree[[1]]$Copulas[[1]]$p.value.CvM
  ks_matrix[d,d-1] <- tree[[1]]$Copulas[[1]]$p.value.KS
  
  return(list(RVine = RVineMatrix(Matrix = R_matrix,
                                  family = cop_matrix,
                                  par = theta_matrix,
                                  par2 = delta_matrix,
                                  names = names),
              AICs = aic_matrix,
              p.value.CvM = cvm_matrix,
              p.value.KS = ks_matrix,
              Tree = tree_0))
}

## These are modified since the upper limit becomes the lower under rotation.
transform_u.full <- function(u_upper, u_lower, v_upper, v_lower, cop, method){
  rot <- cop$rotation
  if(method == "upper"){
    # rot == 1 is 0 degrees, and no rotation
    if(rot == 2){
      # 90 degree
      return(cop$cop_du(u = 1 - u_lower, v = v_upper, theta = cop$theta, delta = cop$delta))
    }else if(rot == 3){
      # 180 degree
      return(1 - cop$cop_du(u = 1 - u_lower, v = 1 - v_lower, theta = cop$theta, delta = cop$delta))
    }else if(rot == 4){
      # 270 degree
      return(1 - cop$cop_du(u = u_upper, v = 1 - v_lower, theta = cop$theta, delta = cop$delta))
    }
    return(cop$cop_du(u = u_upper, v = v_upper, theta = cop$theta, delta = cop$delta))
  }else{
    # rot == 1 is 0 degrees, and no rotation
    if(rot == 2){
      # 90 degree
      return(cop$cop_du(u = 1 - u_upper, v = v_lower, theta = cop$theta, delta = cop$delta))
    }else if(rot == 3){
      # 180 degree
      return(1 - cop$cop_du(u = 1 - u_upper, v = 1 - v_upper, theta = cop$theta, delta = cop$delta))
    }else if(rot == 4){
      # 270 degree
      return(1 - cop$cop_du(u = u_lower, v = 1 - v_upper, theta = cop$theta, delta = cop$delta))
    }
    return(cop$cop_du(u = u_lower, v = v_lower, theta = cop$theta, delta = cop$delta))
  }
}

transform_v.full <- function(u_upper, u_lower, v_upper, v_lower, cop, method){
  rot <- cop$rotation
  if(method == "upper"){
    # rot == 1 is 0 degrees, and no rotation
    if(rot == 2){
      # 90 degree
      return(cop$cop_dv(u = 1 - u_lower, v = v_upper, theta = cop$theta, delta = cop$delta))
    }else if(rot == 3){
      # 180 degree
      return(1 - cop$cop_dv(u = 1 - u_lower, v = 1 - v_lower, theta = cop$theta, delta = cop$delta))
    }else if(rot == 4){
      # 270 degree
      return(1 - cop$cop_dv(u = u_upper, v = 1 - v_lower, theta = cop$theta, delta = cop$delta))
    }
    return(cop$cop_dv(u = u_upper, v = v_upper, theta = cop$theta, delta = cop$delta))
  }else{
    # rot == 1 is 0 degrees, and no rotation
    if(rot == 2){
      # 90 degree
      return(cop$cop_dv(u = 1 - u_upper, v = v_lower, theta = cop$theta, delta = cop$delta))
    }else if(rot == 3){
      # 180 degree
      return(1 - cop$cop_dv(u = 1 - u_upper, v = 1 - v_upper, theta = cop$theta, delta = cop$delta))
    }else if(rot == 4){
      # 270 degree
      return(1 - cop$cop_dv(u = u_lower, v = 1 - v_upper, theta = cop$theta, delta = cop$delta))
    }
    return(cop$cop_dv(u = u_lower, v = v_lower, theta = cop$theta, delta = cop$delta))
  }
}

censor_RVine_select_full <- function(data, indeptest = TRUE, level = 0.05,
                                              include_tawn = TRUE, include_amh = FALSE, include_t = TRUE, core_lim=20,
                                              calc_pVal = FALSE, N_bootstrap = NA){
  data_upper <- pobs(data, ties.method = "max")
  data_lower <- pobs(data, ties.method = "min")
  d <- ncol(data)
  n <- nrow(data)
  if(is.na(N_bootstrap)){
    N_bootstrap <- 10*n
  }
  tree <- list()
  copulas <- list()
  
  ## Constructing a matrix of sources, destinations and weights of each edge
  ## that can be converted to an igraph like graph
  tau <- cor(data_upper, method = "kendall")
  weights <- tau[upper.tri(tau)]
  sources <- matrix(rep(1:d, d), ncol = d)
  sources <- sources[upper.tri(sources)]
  destinations <- matrix(rep(1:d, d), ncol = d, byrow = TRUE)
  destinations <- destinations[upper.tri(destinations)]
  totals <- cbind(abs(weights), sources, destinations) ## Total of weights, sources and destinations
  totals <- totals[order(totals[,1], decreasing = TRUE),]
  graph <- graph_from_edgelist(cbind(as.character(totals[,2]), as.character(totals[,3])), directed = FALSE)
  E(graph)$weight <- -as.numeric(totals[,1])
  graph <- minimum.spanning.tree(graph)
  
  ################### transforming ########################
  transformed_upper <- matrix(0, nrow=n, ncol = 2*length(E(graph)))
  transformed_lower <- matrix(0, nrow=n, ncol = 2*length(E(graph)))
  transformed_cond <- matrix(0, nrow=1, ncol = 2*(d-1))
  for(j in 1:length(E(graph))){
    ind1 <- as.numeric(tail_of(graph, j)$name)
    ind2 <- as.numeric(head_of(graph, j)$name)
    ui_upper <- data_upper[,ind1]
    ui_lower <- data_lower[,ind1]
    vi_upper <- data_upper[,ind2]
    vi_lower <- data_lower[,ind2]
    cop <- censor_copula_select.full(ui_upper, ui_lower, vi_upper, vi_lower,
                                       indeptest = indeptest, level = level,
                                       include_tawn = include_tawn, include_amh = include_amh, include_t = include_t)
    if(calc_pVal){
      test <- censor_gof_test_parallel(ui_upper, vi_upper, cop = cop, N = N_bootstrap, core_lim = core_lim)
      cop$p.value.CvM <- test$p.value.CvM
      cop$p.value.KS <- test$p.value.KS
    }

    transformed_upper[,(j*2-1)] <- transform_u.full(u_upper = ui_upper, u_lower =  ui_lower, v_upper = vi_upper, v_lower = vi_upper,
                                                    cop = cop, method = "upper")
    transformed_lower[,(j*2-1)] <- transform_u.full(u_upper = ui_upper, u_lower =  ui_lower, v_upper = vi_upper, v_lower = vi_upper,
                                                    cop = cop, method = "lower")
    transformed_upper[,(j*2)] <- transform_v.full(u_upper = ui_upper, u_lower =  ui_lower, v_upper = vi_upper, v_lower = vi_upper,
                                                  cop = cop, method = "upper")
    transformed_lower[,(j*2)] <- transform_v.full(u_upper = ui_upper, u_lower =  ui_lower, v_upper = vi_upper, v_lower = vi_upper,
                                                  cop = cop, method = "lower")
    transformed_cond[1,(j*2-1)] <- ind1
    transformed_cond[1,(j*2)] <- ind2
    copulas %<>% rlist::list.append(Copula = cop)
  }
  tree %<>% rlist::list.append(list(Copulas = copulas, Graph = graph,
                                    t.data_upper = transformed_upper, t.data_lower = transformed_lower, t.cond = transformed_cond))
  
  ############## Transformation on higher levels
  for(i in 2:(d-1)){
    # Building structure for the next full tree
    tau <- cor(transformed_upper, method = "kendall")
    sources <- c()
    destinations <- c()
    weights <- c()
    vertex_data_upper <- matrix(nrow=n)
    vertex_data_lower <- matrix(nrow=n)
    vertex_cond <- matrix(nrow=i-1)
    vertex_var <- c()
    ## Looping over all edges in the previous tree to find all possible options for the next tree given this
    for(j in 1:(length(E(graph))-1)){
      for(k in (j+1):length(E(graph))){
        first_ends <- ends(graph, j, names = FALSE)
        sec_ends <- ends(graph, k, names = FALSE)
        if(any(first_ends %in% sec_ends)){
          sources %<>% c(paste(ends(graph, j, names = TRUE), collapse = ","))
          destinations %<>% c(paste(ends(graph, k, names = TRUE), collapse = ","))
          
          if(i>=3){
            ## This if/else is simply due to compatibility issues with the graph from the first tree
            ## For T_3 and so on, I have better control over the conditioned variables and the conditioning sets
            ## and which edges go where. This is a more tedious process for finding this.
            first_involved <- unique(c(transformed_cond[1:i-1,(j*2-1):(j*2)]))
            second_involved <- unique(c(transformed_cond[1:i-1,(k*2-1):(k*2)]))
            conditioned <- unique(first_involved[first_involved %in% second_involved])
            cond1 <- which((transformed_cond[i-1,(j*2-1):(j*2)] %in% conditioned)) +j*2-2
            cond2 <- which((transformed_cond[i-1,(k*2-1):(k*2)] %in% conditioned)) +k*2-2
            weights %<>% c(tau[cond1, cond2])
            
            vertex_data_upper %<>% cbind(transformed_upper[,cond1], transformed_upper[,cond2])
            vertex_data_lower %<>% cbind(transformed_lower[,cond1], transformed_lower[,cond2])
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
            vertex_data_upper %<>% cbind(transformed_upper[,ind1], transformed_upper[,ind2])
            vertex_data_lower %<>% cbind(transformed_lower[,ind1], transformed_lower[,ind2])
            
            
            vertex_cond %<>% cbind(transformed_cond[,ind1], transformed_cond[,ind2])
            vertex_var %<>% c(transformed_cond[i-1,cond1], transformed_cond[i-1,cond2])
          }
        }
      }
    }
    vertex_data_upper <- vertex_data_upper[,-1]
    vertex_data_lower <- vertex_data_lower[,-1]
    vertex_cond <- matrix(vertex_cond[,-1], nrow=i-1)
    graph <- graph_from_edgelist(cbind(sources, destinations), directed = FALSE)
    E(graph)$weight <- -abs(weights)
    new_graph <- minimum.spanning.tree(graph)
    
    ########################################################
    ## In this step, we find the index to which edges have been removed to construct the spanning tree
    ## I did not find a "natural" way of doing this with the igraph library, but this tedious solution works
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
      vertex_data_upper <- vertex_data_upper[,-c(data_removal*2 - 1, data_removal*2)]
      vertex_data_lower <- vertex_data_lower[,-c(data_removal*2 - 1, data_removal*2)]
      vertex_cond <- vertex_cond[,-c(data_removal*2 - 1, data_removal*2)]
      vertex_var <- vertex_var[-c(data_removal*2 - 1, data_removal*2)]
    }
    #########################################################
    #########################################################
    
    graph <- new_graph
    copulas <- list()
    transformed_upper <- matrix(0, nrow = n, ncol = 2*length(E(graph)))
    transformed_lower <- matrix(0, nrow = n, ncol = 2*length(E(graph)))
    for(j in 1:length(E(graph))){
      ind1 <- j*2 - 1
      ind2 <- j*2
      ui_upper <- vertex_data_upper[,ind1]
      ui_lower <- vertex_data_lower[,ind1]
      vi_upper <- vertex_data_upper[,ind2]
      vi_lower <- vertex_data_lower[,ind2]
      cop <- censor_copula_select.full(ui_upper, ui_lower, vi_upper, vi_lower,
                                         indeptest = indeptest, level = level, 
                                         include_tawn = include_tawn, 
                                         include_amh = include_amh,
                                         include_t = include_t)
      if(calc_pVal){
        test <- censor_gof_test_parallel(ui_upper, vi_upper, cop = cop, N = N_bootstrap, core_lim = core_lim)
        cop$p.value.CvM <- test$p.value.CvM
        cop$p.value.KS <- test$p.value.KS
      }
      transformed_upper[,ind1] <- transform_u.full(u_upper = ui_upper, u_lower = ui_lower, v_upper = vi_upper, v_lower = vi_lower, cop = cop, method = "upper")
      transformed_lower[,ind1] <- transform_u.full(u_upper = ui_upper, u_lower = ui_lower, v_upper = vi_upper, v_lower = vi_lower, cop = cop, method = "lower")
      transformed_upper[,ind2] <- transform_v.full(u_upper = ui_upper, u_lower = ui_lower, v_upper = vi_upper, v_lower = vi_lower, cop = cop, method = "upper")
      transformed_lower[,ind2] <- transform_v.full(u_upper = ui_upper, u_lower = ui_lower, v_upper = vi_upper, v_lower = vi_lower, cop = cop, method = "lower")
      copulas %<>% rlist::list.append(Copula = cop)
    }
    transformed_cond <- rbind(vertex_cond, vertex_var)
    tree %<>% rlist::list.append(list(Copulas = copulas, Graph = graph,
                                      t.data_upper = transformed_upper, t.data_lower = transformed_lower, t.cond = transformed_cond))
  }
  return(rvine_translation_pval(tree = tree, names = colnames(data)))
}