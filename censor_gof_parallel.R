library(tidyverse)
library(magrittr)
library(rlist)
library(igraph)
library(VineCopula)
library(copula)
library(doParallel)
library(foreach)


censor_gof_test_parallel <- function(u, v, cop, rot = 1, N = 1000, core_lim = 20){
  source("censor_est.R")
  if(cop$cop_name == "Independence"){
    return(list(
      p.value.CvM = NA,
      p.value.KS = NA
    ))
  }
  rot <- check_rotations(rot)
  if(rot != cop$rotation){
    if(!is.na(cop$theta)){
      rot <- cop$rotation
      warning("Mismatch between specified rotation and the stored copula rotation. Rotation from stored copula is used")
      print(rot)
    }else{
      warning("Mismatch between specified rotation and the stored copula rotation. rot = 1 used")
    }
  }
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(min(numCores, core_lim))
  
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
  # ecdf_u <- c(0,unique(rank(sort(u_upper), ties.method = "max"))/(n+1),1)
  # ecdf_v <- c(0,unique(rank(sort(v_upper), ties.method = "max"))/(n+1),1)
  
  if(cop$cop_name %in% c("AMH", "Clayton", "Gumbel", "Joe", "Frank")){
    cop_fun <- cop$cop_fun
    sample_copula <- copula::archmCopula(cop$cop_name, param = cop$theta)
    sample_fun <- function(n){
      return(copula::rCopula(n, sample_copula))
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
  #cvm_n <- vector("numeric", N)
  #ks_n <- vector("numeric", N)
  #thetas <- vector("numeric", N)
  #deltas <- vector("numeric", N)
  
  #bootstraps <- tibble(CvM = double(), KS = double(), Theta = double(), Delta = double())
  foreach(i = 1:N, .combine = rbind) %dopar% {
    source("censor_est.R", local = T)
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
    
    cop_n <- try(censor_est_fun(ui_upper, ui_lower, vi_upper, vi_lower, cop))
    if("try-error" %in% class(cop_n)){
      tibble(
        CvM = NA,
        KS = NA,
        Theta = NA,
        Delta = NA
      )
    }else{
      tibble(
        CvM = sum((C.n(cbind(ui_upper, vi_upper), cbind(ui_upper, vi_upper), ties.method = "max") -
                     cop_fun(ui_upper, vi_upper, cop_n$theta, cop_n$delta))^2),
        KS = max(abs(C.n(cbind(ui_upper, vi_upper), cbind(ui_upper, vi_upper), ties.method = "max") -
                       cop_fun(ui_upper, vi_upper, cop_n$theta, cop_n$delta))),
        Theta = cop_n$theta,
        Delta = cop_n$delta
      )
    }
    
    #cvm_n[i] = sum((C.n(cbind(ui_upper, vi_upper), cbind(ui_upper, vi_upper), ties.method = "max") -
    #                  cop_fun(ui_upper, vi_upper, cop_n$theta))^2) # squared ???
    #ks_n[i] = max(abs(C.n(cbind(ui_upper, vi_upper), cbind(ui_upper, vi_upper), ties.method = "max") -
    #                    cop_fun(ui_upper, vi_upper, cop_n$theta)))
    #thetas[i] = cop_n$theta
    #deltas[i] = cop_n$delta
  } -> bootstraps
  if(any(is.na(bootstraps$CvM))){
    print(sum(is.na(bootstraps$CvM)))
  }
  return(list(
    CvM = cvm,
    p.value.CvM = mean(bootstraps$CvM[!is.na(bootstraps$CvM)] >= cvm),
    KS = ks,
    p.value.KS = mean(bootstraps$KS[!is.na(bootstraps$KS)] >= ks),
    Results = bootstraps
  ))
}
