library(tidyverse)
library(magrittr)
library(VineCopula)
library(copula)
library(rlist)

############################ INDEPENDENCE ##########################################
cop_independence <- function(u, v, theta = NULL, delta = NULL, check.pars = FALSE){
  return(u*v)
}

cop_independence_du <- function(u, v, theta = NULL, delta = NULL, check.pars = FALSE){
  return(v)
}

cop_independence_dv <- function(u, v, theta = NULL, delta = NULL, check.pars = FALSE){
  return(u)
}

cop_independence_dens <- function(u, v, theta = NULL, delta = NULL, check.pars = FALSE){
  return(rep(0, length(u[,1])))
}

################### AMH ############################
cop_amh <- function(u,v,theta,delta=NULL, check.pars = TRUE){
  return(u*v/(1 - theta*(1-u)*(1-v)))
}

cop_amh_du <- function(u,v,theta,delta=NULL, check.pars = TRUE){
  return((theta*(v-1)*v + v)/(1 - theta*(u-1)*(v-1))^2)
}

cop_amh_dv <- function(u,v,theta,delta=NULL, check.pars = TRUE){
  return((theta*(u-1)*u + u)/(1 - theta*(v-1)*(u-1))^2)
}

cop_amh_dens <- function(u,v,theta,delta=NULL, check.pars = TRUE){
  return(-(theta*
             (theta*(u - 1)*(v - 1) + u*v + u + v - 2) + 1)/
           (theta*(u - 1)*(v - 1) - 1)^3)
}

################### CLAYTON ########################
cop_clayton <- function(u,v,theta,delta=NULL, check.pars = TRUE){
  return(pmax((u^(-theta) + v^(-theta) - 1), 0)^(-1/theta))
}

cop_clayton_du <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return((u^(-theta - 1))*
           (u^(-theta) + v^(-theta) - 1)^(-1/theta - 1)*
           (u^(-theta) + v^(-theta) - 1 >0))
}

cop_clayton_dv <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return((v^(-theta - 1))*
           (u^(-theta) + v^(-theta) - 1)^(-1/theta - 1)*
           (u^(-theta) + v^(-theta) - 1 >0))
}

cop_clayton_dens_2 <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(((theta + 1)*
            u^(theta - 1)*
            v^(theta - 1)*
            (u^(-theta) + v^(-theta) - 1)^(-1/theta))/
           (-u^theta + u^theta * v^theta - v^theta)^2*
           (u^(-theta) + v^(-theta) - 1 >0))
}

cop_clayton_dens <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return((theta + 1)*
           (u*v)^(-(theta+1))*
           (u^(-theta) + v^(-theta) - 1)^(-(2*theta + 1)/(theta))*
                                           (u^(-theta) + v^(-theta) - 1 >0))
}


################### Frank ##########################

cop_frank <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(-1/theta * log(1 + exp(-(-log((exp(-theta *u) - 1)/(exp(-theta) - 1)) + -log((exp(-theta * v) - 1)/(exp(-theta) - 1))))*(exp(-theta) - 1)))
}

cop_frank_du <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return((exp(theta)*(exp(theta*v) - 1))/(exp(theta)*
                                            (exp(theta*v) + exp(theta*u) - 1) - exp(theta*(v + u))))
}

cop_frank_dv <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return((exp(theta)*(exp(theta*u) - 1))/(exp(theta)*
                                            (exp(theta*u) + exp(theta*v) - 1) - exp(theta*(u + v))))
}

cop_frank_dens <- function(u,v, theta, delta=NULL, check.pars = TRUE){
  return(
    -1/theta * (exp(-(-log((exp(-theta * u) - 1)/(exp(-theta) - 1)) +
                        -log((exp(-theta * v) - 1)/(exp(-theta) - 1)))) * 
                  (exp(-theta * v) *
                     theta/(exp(-theta) - 1)/((exp(-theta * v) - 1)/(exp(-theta) - 1)))*
                  (exp(-theta * u) * theta/(exp(-theta) - 1)/
                     ((exp(-theta * u) - 1)/(exp(-theta) - 1))) *
                  (exp(-theta) - 1)/(1 + exp(-(-log((exp(-theta * u) - 1)/(exp(-theta) -1)) +
                                                 -log((exp(-theta * v) - 1)/(exp(-theta) - 1)))) *
                                       (exp(-theta) - 1)) - exp(-(-log((exp(-theta * u) - 1)/(exp(-theta) -1)) +
                                                                    -log((exp(-theta * v) - 1)/(exp(-theta) - 1)))) * 
                  (exp(-theta * u) *
                     theta/(exp(-theta) - 1)/((exp(-theta *u) - 1)/(exp(-theta) - 1))) *
                  (exp(-theta) - 1) * (exp(-(-log((exp(-theta *u) - 1)/(exp(-theta) - 1)) +
                                               -log((exp(-theta * v) - 1)/(exp(-theta) - 1)))) *
                                         (exp(-theta * v) * theta/(exp(-theta) - 1)/((exp(-theta * v) - 1)/(exp(-theta) - 1))) *
                                         (exp(-theta) - 1))/(1 + exp(-(-log((exp(-theta * u) - 1)/(exp(-theta) - 1)) + -log((exp(-theta * v) - 1)/(exp(-theta) - 1))))
                                                             * (exp(-theta) - 1))^2)
  )
}


################### Gumbel ##########################

cop_gumbel <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)))
}

cop_gumbel_du <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(((-log(v))^theta*exp(-((-log(v))^theta + (-log(u))^theta)^(1/theta))*
                   ((-log(v))^theta + (-log(u))^theta)^(1/theta - 1))/(u*log(u)) -
                  (exp(-((-log(v))^theta + (-log(u))^theta)^(1/theta))*
                     ((-log(v))^theta + (-log(u))^theta)^(1/theta))/(u*log(u)))
}

cop_gumbel_dv <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(((-log(u))^theta*exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta))*
            ((-log(u))^theta + (-log(v))^theta)^(1/theta - 1))/(v*log(v)) -
           (exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta))*
              ((-log(u))^theta + (-log(v))^theta)^(1/theta))/(v*log(v)))
}

cop_gumbel_dens <- function(u, v, theta, delta=NULL, check.pars = TRUE){
  exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)) * 
    (((-log(u))^theta + (-log(v))^theta)^((1/theta) - 1) * 
       ((1/theta) * ((-log(v))^(theta - 1) * (theta * (1/v))))) * 
    (((-log(u))^theta + (-log(v))^theta)^((1/theta) - 1) * 
       ((1/theta) * ((-log(u))^(theta - 1) * (theta * (1/u))))) - 
    exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta)) *
    (((-log(u))^theta + (-log(v))^theta)^(((1/theta) - 1) - 1) *
       (((1/theta) - 1) * ((-log(v))^(theta - 1) *
          (theta * (1/v)))) * ((1/theta) *
              ((-log(u))^(theta - 1) * (theta * (1/u)))))
}

################### Joe ##########################

cop_joe <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(1 - ((1-u)^theta + (1-v)^theta - (1-u)^theta*(1-v)^theta)^(1/theta))
}

cop_joe_du <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  #return(- ( (1-u)^theta + (1-v)^theta - (1-u)^theta*(1-v)^theta )^(1/theta - 1) *
  #            ( theta * (1-u)^theta*(1-v)^(theta - 1) - theta*(1-v*(theta - 1)) )/theta)
  return((1 - u)^(theta - 1)*
           (1 - (1 - v)^theta)*
           ((1 - v)^theta - (1 - u)^theta*
              ((1 - v)^theta - 1))^(1/theta - 1))
}

cop_joe_dv <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  #return(- ( (1-v)^theta + (1-u)^theta - (1-v)^theta*(1-u)^theta )^(1/theta - 1) *
  #         ( theta * (1-v)^theta*(1-u)^(theta - 1) - theta*(1-u*(theta - 1)) )/theta)
  return((1 - (1 - u)^theta)*
           (1 - v)^(theta - 1)*
           ((1 - v)^theta - (1 - u)^theta*
              ((1 - v)^theta - 1))^(1/theta - 1))
}

cop_joe_dens <- function(u, v, theta,delta=NULL, check.pars = TRUE){
  return(-(1 - u)^(theta - 1)*
           (1 - v)^(theta - 1)*
           ((1 - v)^theta - (1 - u)^theta*
              ((1 - v)^theta - 1))^(1/theta - 2)*
           (-theta + (1 - u)^theta *
              ((1 - v)^theta - 1) - (1 - v)^theta + 1))
}


################################################################################################
################################################################################################

################################### BB1 ###############################
# delta >= 1, theta > 0

cop_bb1 <- function(u,v,theta, delta, check.pars = TRUE){
  return(
    (1 + ( (u^(-theta) - 1)^delta + (v^(-theta) - 1)^delta)^(1/delta))^(-1/theta)
  )
}

cop_bb1_du <- function(u, v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 7, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb1_dv <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 7, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb1_dens <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 7, par = theta, par2 = delta, check.pars = check.pars))
}


################################### BB6 ###############################
# theta >= 1, delta >= 1

cop_bb6 <- function(u,v,theta, delta, check.pars = TRUE){
  return(
    1 - (1 - exp( - ((#log( ??????????????????
      (-log(1 - (1-u)^theta))^delta + (-log(1 - (1-v)^theta))^delta
    ))^(1/delta) ) )^(1/theta)
  )
}

cop_bb6_du <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 8, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb6_dv <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 8, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb6_dens <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 8, par = theta, par2 = delta, check.pars = check.pars))
}


################################### BB7 ###############################
# joe clayton theta >=1, delta > 0

cop_bb7 <- function(u,v,theta, delta, check.pars = TRUE){
  return(
    1 - (
      1 - (
        (1 - (1-u)^theta)^(-delta) + (1 - (1-v)^theta)^(-delta) - 1
      )^(-1/delta)
    )^(1/theta)
  )
}

cop_bb7_du <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 9, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb7_dv <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 9, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb7_dens <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 9, par = theta, par2 = delta, check.pars = check.pars))
}

################################### BB8 ###############################
# delta in (0,1]

cop_bb8 <- function(u,v,theta, delta, check.pars = TRUE){
  return(
    (1 - (1 - exp(-(
      -log((1 - (1-delta*u)^theta)/(1 - (1-delta)^theta)) -
        log((1 - (1-delta*v)^theta)/(1 - (1-delta)^theta))
    ))*(
      1-(1-delta)^theta
    ))^(1/theta))/delta
  )
}

cop_bb8_du <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 10, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb8_dv <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 10, par = theta, par2 = delta, check.pars = check.pars))
}

cop_bb8_dens <- function(u,v,theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 10, par = theta, par2 = delta, check.pars = check.pars))
}
################### Tawn Type 1 ##########################

cop_tawn1 <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopCDF(u1 = u,u2 = v, family = 104, par = theta, par2 = delta, check.pars = check.pars))
}

cop_tawn1_du <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 104, par = theta, par2 = delta, check.pars = check.pars))
}

cop_tawn1_dv <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 104, par = theta, par2 = delta, check.pars = check.pars))
}

cop_tawn1_dens <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 104, par = theta, par2 = delta, check.pars = check.pars))
}

################### Tawn Type 2 ##########################

cop_tawn2 <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopCDF(u1 = u,u2 = v, family = 204, par = theta, par2 = delta, check.pars = check.pars))
}

cop_tawn2_du <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 204, par = theta, par2 = delta, check.pars = check.pars))
}

cop_tawn2_dv <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 204, par = theta, par2 = delta, check.pars = check.pars))
}

cop_tawn2_dens <- function(u, v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 204, par = theta, par2 = delta, check.pars = check.pars))
}

##########################################################
##################### ELLIPTICAL #########################
##########################################################

############# NORMAL ################################################################################
cop_normal <- function(u,v, theta, delta = NULL, check.pars = TRUE){
  return(as.numeric(VineCopula::BiCopCDF(u1 = u,u2 = v, family = 1, par = theta, par2 = 0, check.pars = check.pars)))
}

cop_normal_du <- function(u,v, theta, delta = NULL, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 1, par = theta, par2 = 0, check.pars = check.pars))
}

cop_normal_dv <- function(u,v, theta, delta = NULL, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 1, par = theta, par2 = 0, check.pars = check.pars))
}

cop_normal_dens <- function(u,v, theta, delta = NULL, check.pars = TRUE){
  return((1/sqrt(1-(theta^2)))*exp(
  #return((1/sqrt(1-(min(abs(theta), 1)^2)))*exp(
    -((theta^2)*(qnorm(u)^2 + qnorm(v)^2) - 2*theta*qnorm(u)*qnorm(v))/(2*(1 - (theta^2)))
    ))
  #return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 1, par = theta, par2 = 0, check.pars = check.pars))
}


################# STUDENT-T ##########################################
cop_student <- function(u,v, theta, delta, check.pars = TRUE){
  return(as.numeric(VineCopula::BiCopCDF(u1 = u,u2 = v, family = 2, par = theta, par2 = delta, check.pars = check.pars)))
}

cop_student_du <- function(u,v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc1(u1 = u,u2 = v, family = 2, par = theta, par2 = delta, check.pars = check.pars))
}

cop_student_dv <- function(u,v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopHfunc2(u1 = u,u2 = v, family = 2, par = theta, par2 = delta, check.pars = check.pars))
}

cop_student_dens <- function(u,v, theta, delta, check.pars = TRUE){
  return(VineCopula::BiCopPDF(u1 = u,u2 = v, family = 2, par = theta, par2 = delta, check.pars = check.pars))
}

########################################################################################################

# Here is a simple copula construction based on a list. An object would be better
get_cop <- function(cop_name){
  cop_names <- c(amh = 1,
                 clayton = 2,
                 frank = 3,
                 gumbel = 4,
                 joe = 5,
                 bb1 = 6,
                 bb6 = 7,
                 bb7 = 8,
                 bb8 = 9,
                 tawn = 10,
                 tawn2 = 11,
                 normal = 12,
                 t = 13,
                 independence = 14,
                 t2 = 15)
  if(!(tolower(cop_name) %in% tolower(names(cop_names))) 
     & !(cop_name %in% cop_names)){
    stop(paste("copula not recognized! Choose between: ", paste(names(cop_names), collapse = ", ")))
  }else{
    if(is.character(cop_name)){
      cop_num <- cop_names[tolower(cop_name)]
    }else{
      cop_num <- cop_name
    }
  }
  copula_one_param <- list(
    cop_name = NULL,
    num = NULL,
    theta = NA,
    delta = NA,
    cop_fun = NULL,
    cop_du = NULL,
    cop_dv = NULL,
    cop_dens = NULL,
    optim.limits = list(list(low = NA, up = NA)),
    start.value = NA,
    n.param = 1,
    log.lik = NA,
    AIC = NA,
    rotation = NA,
    fam = NA,
    p.value.CvM = NA,
    p.value.KS = NA
  )
  copula_two_param <- list(
    cop_name = NULL,
    num = NULL,
    theta = NA,
    delta = NA,
    cop_fun = NULL,
    cop_du = NULL,
    cop_dv = NULL,
    cop_dens = NULL,
    start.value = c(NA,NA),
    optim.lower = c(NA,NA),
    optim.upper = c(NA,NA),
    n.param = 2,
    log.lik = NA,
    AIC = NA,
    rotation = NA,
    tau = NA,
    fam = NA,
    p.value.CvM = NA,
    p.value.KS = NA
  )
  if(cop_num == 1){
    copula_one_param$cop_name = "AMH"
    copula_one_param$num = 1
    copula_one_param$cop_fun = cop_amh
    copula_one_param$cop_du = cop_amh_du
    copula_one_param$cop_dv = cop_amh_dv
    copula_one_param$cop_dens = cop_amh_dens
    copula_one_param$optim.limits = list(list(low = -1, up = 1, start = 0.4))
    return(copula_one_param)
  }else if(cop_num == 2){
    copula_one_param$cop_name = "Clayton"
    copula_one_param$num = 2
    copula_one_param$fam = 3
    copula_one_param$cop_fun = cop_clayton
    copula_one_param$cop_du = cop_clayton_du
    copula_one_param$cop_dv = cop_clayton_dv
    copula_one_param$cop_dens = cop_clayton_dens
    # copula_one_param$optim.limits = list(list(low = -1, up = -0.001, start = -0.5),
    #                                      list(low = 0.001, up = 10, start = 1))
    # The VineCopula-package does not support negative values for clayton's
    # Also, limits affect the optimization of the clayton copula
    copula_one_param$optim.limits = list(list(low = 0.001, up = 22, start = 15))
    
    return(copula_one_param)
  }else if(cop_num == 3){
    copula_one_param$cop_name = "Frank"
    copula_one_param$num = 3
    copula_one_param$fam = 5
    copula_one_param$rotation = 1
    copula_one_param$cop_fun = cop_frank
    copula_one_param$cop_du = cop_frank_du
    copula_one_param$cop_dv = cop_frank_dv
    copula_one_param$cop_dens = cop_frank_dens
    copula_one_param$optim.limits = list(list(low = -45, up = -0.001, start = -5),
                                         list(low = 0.001, up = 45, start = 5))
    return(copula_one_param)
  }else if(cop_num == 4){
    copula_one_param$cop_name = "Gumbel"
    copula_one_param$num = 4
    copula_one_param$fam = 4
    copula_one_param$cop_fun = cop_gumbel
    copula_one_param$cop_du = cop_gumbel_du
    copula_one_param$cop_dv = cop_gumbel_dv
    copula_one_param$cop_dens = cop_gumbel_dens
    copula_one_param$optim.limits = list(list(low = 1.001, up = 25, start = 2))
    return(copula_one_param)
  }else if(cop_num == 5){
    copula_one_param$cop_name = "Joe"
    copula_one_param$num = 5
    copula_one_param$fam = 6
    copula_one_param$cop_fun = cop_joe
    copula_one_param$cop_du = cop_joe_du
    copula_one_param$cop_dv = cop_joe_dv
    copula_one_param$cop_dens = cop_joe_dens
    copula_one_param$optim.limits = list(list(low = 1.001, up = 25, start = 1.16))
    return(copula_one_param)
    ## The limits for these are the same as in VineCopula
  }else if(cop_num == 6){
    copula_two_param$cop_name = "BB1"
    copula_two_param$num = 6
    copula_two_param$fam = 7
    copula_two_param$cop_fun = cop_bb1
    copula_two_param$cop_du = cop_bb1_du
    copula_two_param$cop_dv = cop_bb1_dv
    copula_two_param$cop_dens = cop_bb1_dens
    copula_two_param$start.value = c(0.5, 1.5)
    copula_two_param$optim.lower = c(0.00001, 1)
    copula_two_param$optim.upper = c(7, 7)
    return(copula_two_param)
  }else if(cop_num == 7){
    copula_two_param$cop_name = "BB6"
    copula_two_param$num = 7
    copula_two_param$fam = 8
    copula_two_param$cop_fun = cop_bb6
    copula_two_param$cop_du = cop_bb6_du
    copula_two_param$cop_dv = cop_bb6_dv
    copula_two_param$cop_dens = cop_bb6_dens
    copula_two_param$start.value = c(1.5, 1.5)
    copula_two_param$optim.lower = c(1, 1)
    copula_two_param$optim.upper = c(6, 8)
    return(copula_two_param)
  }else if(cop_num == 8){
    copula_two_param$cop_name = "BB7"
    copula_two_param$num = 8
    copula_two_param$fam = 9
    copula_two_param$cop_fun = cop_bb7
    copula_two_param$cop_du = cop_bb7_du
    copula_two_param$cop_dv = cop_bb7_dv
    copula_two_param$cop_dens = cop_bb7_dens
    copula_two_param$start.value = c(1.5, 0.5)
    copula_two_param$optim.lower = c(1, 0.00001)
    copula_two_param$optim.upper = c(6, 75)
    return(copula_two_param)
  }else if(cop_num == 9){
    copula_two_param$cop_name = "BB8"
    copula_two_param$num = 9
    copula_two_param$fam = 10
    copula_two_param$cop_fun = cop_bb8
    copula_two_param$cop_du = cop_bb8_du
    copula_two_param$cop_dv = cop_bb8_dv
    copula_two_param$cop_dens = cop_bb8_dens
    copula_two_param$start.value = c(3, 0.2)
    copula_two_param$optim.lower = c(1.001, 0.0001)
    copula_two_param$optim.upper = c(7.999, 1-0.0001)
    return(copula_two_param)
  }else if(cop_num == 10){
    copula_two_param$cop_name = "Tawn"
    copula_two_param$num = 10
    copula_two_param$fam = 104
    copula_two_param$cop_fun = cop_tawn1
    copula_two_param$cop_du = cop_tawn1_du
    copula_two_param$cop_dv = cop_tawn1_dv
    copula_two_param$cop_dens = cop_tawn1_dens
    copula_two_param$start.value = c(1.2, 0.5)
    copula_two_param$optim.lower = c(1.0001, 0.00001)
    copula_two_param$optim.upper = c(20, 0.9999)
    return(copula_two_param)
  }else if(cop_num == 11){
    copula_two_param$cop_name = "Tawn2"
    copula_two_param$num = 11
    copula_two_param$fam = 204
    copula_two_param$cop_fun = cop_tawn2
    copula_two_param$cop_du = cop_tawn2_du
    copula_two_param$cop_dv = cop_tawn2_dv
    copula_two_param$cop_dens = cop_tawn2_dens
    copula_two_param$start.value = c(1.2, 0.5)
    copula_two_param$optim.lower = c(1.0001, 0.0001)
    copula_two_param$optim.upper = c(20, 0.9999)
    return(copula_two_param)
  }
  else if(cop_num == 12){
    copula_one_param$cop_name = "normal"
    copula_one_param$num = 12
    copula_one_param$fam = 1
    copula_one_param$rotation = 1
    copula_one_param$cop_fun = cop_normal
    copula_one_param$cop_du = cop_normal_du
    copula_one_param$cop_dv = cop_normal_dv
    copula_one_param$cop_dens = cop_normal_dens
    copula_one_param$optim.limits = list(list(low = -0.999, up = 0.999, start = 0.2))
    return(copula_one_param)
  }
  else if(cop_num == 13){
    copula_two_param$cop_name = "t"
    copula_two_param$num = 13
    copula_two_param$fam = 2
    copula_two_param$rotation = 1
    copula_two_param$cop_fun = cop_student
    copula_two_param$cop_du = cop_student_du
    copula_two_param$cop_dv = cop_student_dv
    copula_two_param$cop_dens = cop_student_dens
    copula_two_param$start.value = c(0.2, 4)
    copula_two_param$optim.lower = c(-0.999, 2.0001)
    copula_two_param$optim.upper = c(0.999, 30)
    return(copula_two_param)
  }
  else if(cop_num == 14){
    copula_one_param$cop_name = "Independence"
    copula_one_param$fam = 0
    copula_one_param$rotation = 1
    copula_one_param$cop_fun = cop_independence
    copula_one_param$cop_du = cop_independence_du
    copula_one_param$cop_dv = cop_independence_dv
    copula_one_param$cop_dens = cop_independence_dens
    copula_one_param$n.param = 0
    return(copula_one_param)
  }
  else if(cop_num == 15){
    copula_two_param$cop_name = "t"
    copula_two_param$num = 13
    copula_two_param$fam = 2
    copula_two_param$n.param = 1
    copula_two_param$rotation = 1
    copula_two_param$cop_fun = cop_student
    copula_two_param$cop_du = cop_student_du
    copula_two_param$cop_dv = cop_student_dv
    copula_two_param$cop_dens = cop_student_dens
    copula_two_param$start.value = c(0.2, 4)
    copula_two_param$optim.lower = c(-0.9999, 2.0001)
    copula_two_param$optim.upper = c(0.9999, 30)
    copula_two_param$optim.limits = list(list(low = -0.999, up = 0.999, start = 0.2))
    return(copula_two_param)
  }
}
