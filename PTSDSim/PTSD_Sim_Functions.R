install.packages('Metrics')

library(RCurl)
library(corpcor)
library(BDgraph)
library(BDgraph)
library(cvTools)
library(Metrics)
library(huge)
library(dplyr)
library(glasso)
library(qgraph)


#############################
# Cross-Validation Function #
#############################
l = function(Theta1,Sigma2,p){
  
  invTheta1 = solve(Theta1,diag(1,p))
  
  res = -log(det(invTheta1)) - sum(diag(Sigma2%*%Theta1)) 
  
  return(res)
  
}

crossvalidaatioGlassolle = function(Y,rho){
  
  n = nrow(Y)
  p = ncol(Y)
  
  I = diag(1,p)
  
  k = 5 # ns. ten-fold crossvalidation
  
  ind = sample(1:n, replace=F)
  
  
  #If the sample size is not divisible by five, the last section is 
  #larger than the one in the last section
  if(n %% k != 0){
    
    new = n - (n %% k) 
    
    ind1 = ind[1:(new/k)]
    ind2 = ind[((new/k)+1):(2*(new/k))]
    ind3 = ind[(2*(new/k)+1):(3*(new/k))]
    ind4 = ind[(3*(new/k)+1):(4*(new/k))]
    ind5 = ind[(4*(new/k)+1):(5*(new/k))]; ind5 = c(ind5,ind[(5*new/k + 1):n])
    
  }
  
  if(n %% k == 0){
    
    ind1 = ind[1:(n/k)]
    ind2 = ind[((n/k)+1):(2*(n/k))]
    ind3 = ind[(2*(n/k)+1):(3*(n/k))]
    ind4 = ind[(3*(n/k)+1):(4*(n/k))]
    ind5 = ind[(4*(n/k)+1):(5*(n/k))]
    
  }
  
  Y1 = Y[ind1, ]
  Y2 = Y[ind2, ]
  Y3 = Y[ind3, ]
  Y4 = Y[ind4, ]
  Y5 = Y[ind5, ]
  
  crossY = rbind(Y1,Y2,Y3,Y4,Y5)
  
  CV = rep(0,5)
  
  # Calculate with a glass the estimate for Theta and the Otoscovariance,
  # each one in turn drops out of turn Y1, then Y2, etc. at each Rho value::
  
  pitg = 1
  CVRho = 0
  
  while(pitg < length(rho) + 1){
    
    S         = var(crossY[-(1:nrow(Y1)),])
    hatTheta1 = glasso(S,rho[pitg])$wi
    S1        = var(Y1)
    
    CV[1] = l(hatTheta1,S1,p)
    
    S         = var(crossY[-((nrow(Y1)+1):(2*nrow(Y1))),])
    hatTheta2 = glasso(S,rho[pitg])$wi
    S2        = var(Y2)
    
    CV[2] = l(hatTheta2,S2,p)
    
    S         = var(crossY[-(2*(nrow(Y1) + 1):(3*nrow(Y1))),])
    hatTheta3 = glasso(S,rho[pitg])$wi
    S3        = var(Y3)
    
    CV[3] = l(hatTheta3,S3,p)
    
    S         = var(crossY[-((3*nrow(Y1) +1):(4*nrow(Y1))),])
    hatTheta4 = glasso(S,rho[pitg])$wi
    S4        = var(Y4)
    
    CV[4] = l(hatTheta4,S4,p)
    
    S         = var(crossY[-((4*nrow(Y1) + 1):(5*nrow(Y1))),])
    hatTheta5 = glasso(S,rho[pitg])$wi
    S5        = var(Y5)
    
    CV[5] = l(hatTheta5,S5,p)
    
    CVRho = c(CVRho,mean(CV))
    
    pitg = pitg + 1
    
  }
  
  # Eventually, Rho is chosen as the one with the highest CVRho
  
  CVRho = CVRho[-1]
  
  ind = which( CVRho == max(CVRho))
  
  ParasrhonArvo = rho[ind]
  
  return(ParasrhonArvo) }





########################
# Precision Matrix to ## 
# Partial Correlation ## 
######## Matrix ########
########################
prec2pc <- function(prec_mat){
  #argument : 
  #prec_mat = precision matrix 
  
  prec_mat <- as.matrix(prec_mat)
  
  #switch signs of off-diagonal 
  neg_prec_mat <- - prec_mat
  diag(neg_prec_mat) <- diag(prec_mat)
  
  #standardize values to get partial correlation matrix 
  pc_mat<- cov2cor(neg_prec_mat)
  return(pc_mat)
}



#########################
# True Matrix Generator # 
#########################

true_mat_gen <- function(n, p, sparsity){
  # arguments:
  # n = the sample size (rows)
  # p = the nodes (columns)
  # sparsity = connections
  
  # matrices for storage
  mat <- matrix(0, p, p)
  cor_mat <- matrix(NA, p, p)
  
  # number covariances
  n_off_diag <- (p * (p - 1)) / 2
  
  # number of zeroes in the graph
  n_zero = round(sparsity * n_off_diag, digits = 0)
  
  # number of non-zeroes in the graph
  n_nonzero <- n_off_diag - n_zero
  
  # vector of zeroes
  pcs_zero <- rep(0, n_zero)
  
  
  
  for(i in 1:1000){
    # sample partial correlation values 
    pcs_nonzero <- sample(dat, n_nonzero)
    mat[upper.tri(mat)] <- sample(c(pcs_nonzero, pcs_zero), size = n_off_diag, replace = F)
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    diag(mat) <- 1
    pcor_mat <- mat 
    
    # check positive definite of partial correlations
    # if true, convert to correlation matrix
    if(is.positive.definite(pcor_mat) == TRUE) {
      cor_mat <- pcor2cor(pcor_mat)
      
      # veryify 0 nan's in the cor_mat and no correlations > 1
      if(sum(is.nan(cor_mat))  == 0 & sum(abs(cor_mat[upper.tri(cor_mat)]) > 1) == 0 ) break 
    }
    
  }
  # adjacency matrix
  sig_mat <- ifelse(mat == 0, 0, 1)
  
  # returned objects
  list(cor_mat = cor_mat, pcor_mat = pcor_mat, sig_mat = sig_mat, n_attempts = i)
}



###########################
# Specificity/Sensitivity #
###########################
performance_measures <- function(true, est){
  # arguments: 
  # true = True/Population Matrix
  # est = Estimated matrix
  
  true <- as.matrix(true)
  est <- as.matrix(est)
  
  #Convert partial correlation matrix
  #to adjacency matrix 
  pc_sig <- ifelse(true == 0, 0, 1)
  est_sig <- ifelse(est == 0, 0, 1)
  
  #Flatten upper-tri of the matrix to vector
  pc_true <- pc_sig[upper.tri(pc_sig)]
  pc_est  <- est_sig[upper.tri(est_sig)]
  
  #Count number of true/false positive and negatives
  #in estimated matrix 
  
  TP <- sum(pc_true == 1   & pc_est == 1)
  TN <- sum(pc_true == 0   & pc_est == 0)
  FP <- sum(pc_true == 0   & pc_est == 1)
  FN <- sum(pc_true == 1   & pc_est == 0)
  
  # Calculate Specificity 
  Spec <- TN / (TN + FP)
  
  # Calculate Sensitivity 
  Sens <- TP / (TP + FN)
  
  #Calculate True Detection Rate 
  TDR <- TP/ (TP + FP)
  
  res <- rbind(TP,TN, FP, FN, Spec, Sens, TDR)
  return(res)
}



###########################
### Simulation Function ###
###########################

mapply_func <- function(n, p, sparsity, network){
  
  #Use population matrix to simulate data 
  d <- MASS::mvrnorm(n=n, mu = rep(0,p), Sigma = mat_list$cor_mat)
  
  
  ###############
  ## Fit Model ##
  ############### 
  
  h <- huge::huge(d, method = "glasso",scr = F, nlambda = 100, lambda.min.ratio = .01)
  
  ## select with ebic
  h_gl <- huge.select(h, criterion = "ebic")
  
  ## select with stars
  h_st  <- huge.select(h, criterion = "stars", stars.thresh = .1, 
                       stars.subsample.ratio = .8, rep.num = 20, verbose=F)
  ## select with ric
  h_ric <- huge.select(h, criterion = "ric", rep.num = 20, verbose=F)
  
  ## select with cv
  tune <- seq(0.001, 1, by = .01)
  rho_sel <- crossvalidaatioGlassolle(d, tune)
  cv <- glasso(cov(d), rho = rho_sel)
  
  #Compute Specificity + Sensitivity
  edges <- matrix(0, 7, 4)
  edges <- as.data.frame(edges, row.names=c("TP", "TN", "FP", "FN", 
                                            "Spec", "Sens", 
                                            "TDR"))
  colnames(edges) <- c("ebic", "stars", "ric", "cv")
  
  #ebic
  edges[, 1]<- performance_measures(mat_list$sig_mat, h_gl$opt.icov)
  
  #stars
  
  edges[, 2] <- performance_measures(mat_list$sig_mat, h_st$opt.icov)
  
  #ric
  
  edges[, 3] <- performance_measures(mat_list$sig_mat, h_ric$opt.icov)
  
  
  #cv
  
  edges[, 4] <- performance_measures(mat_list$sig_mat, cv$wi)
  edges <- as.matrix(edges)
  
  ##############
  #Save Lambdas#
  ##############
  
  lam_ebic <- h_gl$opt.lambda
  lam_stars <- h_st$opt.lambda
  lam_ric <- h_ric$opt.lambda
  lam_cv <- rho_sel
  lambdas <- cbind(lam_ebic, lam_stars, lam_ric, lam_cv)
  rownames(lambdas)<- "Lambda"
  colnames(lambdas)<- c("ebic", "stars", "ric", "cv")
  edges <- rbind(edges, lambdas)
  
  #####################
  # Network Densities #
  #####################
  
  #number of possible edges
  total <- (p * (p - 1)) / 2
  
  #estimated edges/ total edges
  #estimated edges = true positives + false positives
  ebic_den <- (edges[1,1] + edges[3, 1])/ total 
  stars_den <- (edges[1,2] + edges[3, 2])/ total
  ric_den <- (edges[1,3] + edges[3, 3])/ total
  cv_den <- (edges[1,4] + edges[3, 4])/ total
  
  density <- cbind(ebic_den, stars_den, ric_den, cv_den)
  rownames(density)<- "Density"
  colnames(density)<- c("ebic", "stars", "ric", "cv")
  edges <- rbind(edges, density)
  
  
  #######################
  # Global Correlations #
  #######################
  
  #take upper tri of estimated and true matrices
  #calculate correlation between each of the estimated matrices and the true matrix
  
  true_pc <- as.numeric(mat_list$pcor_mat[upper.tri(mat_list$pcor_mat)]) 
  #ebic
  ebic_pc <- as.numeric(prec2pc(h_gl$opt.icov)[upper.tri(h_gl$opt.icov)])
  ebic_g_cor <- cor(true_pc, ebic_pc)
  #stars
  stars_pc <- as.numeric(prec2pc(h_st$opt.icov)[upper.tri(h_st$opt.icov)])
  stars_g_cor <- cor(true_pc, stars_pc)
  #ric
  ric_pc <- as.numeric(prec2pc(h_ric$opt.icov)[upper.tri(h_ric$opt.icov)])
  ric_g_cor <- cor(true_pc, ric_pc)
  #cv
  cv_pc <- as.numeric(prec2pc(cv$wi)[upper.tri(cv$wi)])
  cv_g_cor <- cor(true_pc, cv_pc)
  
  global_cor <- cbind(ebic_g_cor, stars_g_cor, ric_g_cor, cv_g_cor)
  rownames(global_cor) <- "GlobalCor"
  colnames(global_cor)<- c("ebic", "stars", "ric", "cv")
  edges <- as.matrix(rbind(edges, global_cor))
  
  #######################
  ### Estimated Edge ####
  #### Correlations #####
  #######################
  
  # find correlation between true edge and corresponding estimated edge
  # subset edges that are nonzero in the true matrix or correspond to the 
  # that edge in the population matrix 
  
  #ebic
  ebic_est_cor <- cor(true_pc[true_pc != 0], ebic_pc[true_pc != 0])
  #stars
  stars_est_cor <- cor(true_pc[true_pc != 0], stars_pc[true_pc != 0])
  #ric
  ric_est_cor <- cor(true_pc[true_pc != 0], ric_pc[true_pc != 0])
  #cv
  cv_est_cor <- cor(true_pc[true_pc != 0], cv_pc[true_pc != 0])
  
  est_cor <- cbind(ebic_est_cor, stars_est_cor, ric_est_cor, cv_est_cor)
  rownames(est_cor) <- "EstCor"
  colnames(est_cor)<- c("ebic", "stars", "ric", "cv")
  edges <- as.matrix(rbind(edges, est_cor))
  
  ## reshape for storing nicely
  edges <- reshape2::melt(edges)
  ## Save PC matrices
  pc_df <- cbind(true_pc, ebic_pc, stars_pc, ric_pc, cv_pc)
  #save(ptsd_list, stars_list, ric_list, cv_list, ebic_list, 
  #file = paste("cov_out", p,  ,n, sparsity, ".Rda", sep = "."))
  return(list(edges=edges, pc_df=pc_df))
}




#############################
# Labeling Results Function #
#############################

sim_run <- function(n, p, sparsity, network){ 
  # mapply takes the function, mapply_fun, and the arguments
  # n and p. Simplify ensures things are returned nicely
  res_map <- mapply(mapply_func, n, p, sparsity, network, SIMPLIFY = F)
  pc_df <- res_map[[1]]$pc_df
  res_map <- list(res_map[[1]]$edges)
  # give condition columns to easily compute things later on
  for(i in 1:length(res_map)){
    
    res_map[[i]]$n <- n 
    res_map[[i]]$p <- p
    res_map[[i]]$sparsity<- sparsity
    res_map[[i]]$network <- network
  }
  return(list(res_map=res_map, pc_df=pc_df))
}

#############################
#     Mean + Variance       #
#       Calculation         #
#############################
mean_var_calc <- function(list){
  value_vec <- lapply(list, "[[", "value")
  value_df <- list.cbind(value_vec)
  mean <- apply(value_df, 1, mean, na.rm =TRUE)
  variance <- apply(value_df, 1, var, na.rm =TRUE)
  minus_value_df <- as.data.frame(list[1])[,-3]
  mean_res<- cbind(minus_value_df, mean, var)
}





