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

censor_RVine_select_pVal_parallel <- function(data, indeptest = TRUE, level = 0.05, censor_level = 1,
                                              include_tawn = TRUE, include_amh = FALSE, include_t = TRUE, core_lim=20,
                                              calc_pVal = FALSE, N_bootstrap = NA){
  d <- ncol(data)
  n <- nrow(data)
  if(is.na(N_bootstrap)){
    N_bootstrap <- 10*n
  }
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
    if(calc_pVal){
      test <- censor_gof_test_parallel(ui, vi, cop = cop, N = N_bootstrap, core_lim = core_lim)
      cop$p.value.CvM <- test$p.value.CvM
      cop$p.value.KS <- test$p.value.KS
    }
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
      if(calc_pVal){
        test <- censor_gof_test_parallel(ui, vi, cop = cop, N = N_bootstrap, core_lim = core_lim)
        cop$p.value.CvM <- test$p.value.CvM
        cop$p.value.KS <- test$p.value.KS
      }
      transformed[,ind1] <- transform_u(u = ui, v = vi, cop = cop)
      transformed[,ind2] <- transform_v(u = ui, v = vi, cop = cop)
      copulas %<>% rlist::list.append(Copula = cop)
    }
    # Kan flyttes høyere???
    transformed_cond <- rbind(vertex_cond, vertex_var)
    tree %<>% rlist::list.append(list(Copulas = copulas, Graph = graph,
                                      t.data = transformed, t.cond = transformed_cond))
  }
  return(rvine_translation_pval(tree = tree, names = colnames(data)))
}

