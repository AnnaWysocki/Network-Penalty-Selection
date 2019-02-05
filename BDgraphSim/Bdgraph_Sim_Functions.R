library(BDgraph)
library(cvTools)
library(huge)
library(glasso)
library(corpcor)


###############################
# Cross Validation  Functions #
###############################


l = function(Theta1,Sigma2,p){
  
  invTheta1 = solve(Theta1,diag(1,p))
  
  res = -log(det(invTheta1)) - sum(diag(Sigma2%*%Theta1)) 
  
  return(res)
  
}

crossvalidaatioGlassolle = function(Y,rho){
  
  n = nrow(Y)
  p = ncol(Y)
  
  I = diag(1,p)
  
  k = 5 # five-fold crossvalidation
  
  ind = sample(1:n, replace=F)
  
  
  ## If the sample size is not divisible by five, the last section is larger
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
  
  CVRho = CVRho[-1]
  
  ind = which( CVRho == max(CVRho))
  
  ParasrhonArvo = rho[ind]
  
  return(ParasrhonArvo) }


## Function to convert precision matrix to partial correlation matrix 
prec2pc <- function(prec_mat){
  prec_mat <- as.matrix(prec_mat)
  neg_prec_mat <- - prec_mat
  diag(neg_prec_mat) <- diag(prec_mat)
  pc_mat<- cov2cor(neg_prec_mat)
  return(pc_mat)
}

## Creates the population matrix used to create simulated data sets
true_mat_gen <- function(n, p, density, df){
  # Arguments:
  # n =  Sample Size 
  # p = Number of Variables 
  # density = Estimated Edges/ Number of Possible Edges
  # df = Degree of Freedom: Controls Partial Correlaton Range
  
  dat <- bdgraph.sim(p = p, n = n, prob = density, b = df)
  
  # Create sparse partial correlation matrix form precision matrix K
  
  pc_mat <- prec2pc(dat$K) 
  pcor_mat <- pc_mat * dat$G # G is the adjacency matrix of the true graph structure
  diag(pcor_mat) <- diag(pc_mat)
  
  # Create correlation matrix 
  cor_mat <- pcor2cor(pcor_mat)
  
  
  
  # Sigma is the generated covariance matrix
  # G is the genearted adjacency matrix 
  # Converting the covariance matrix into a correlation matrix and multiplying it by
  #  the adjacency matrix gives us the sparse correlation matrix 
  
  # Returned objects
  list(cor_mat = cor_mat, pcor_mat = pcor_mat, sig_mat = dat$G)
}

## Function to calculate edge detection performance metrics 
performance_measures <- function(true, est){
  true <- as.matrix(true)
  est <- as.matrix(est)
  pop_sig <- ifelse(true == 0, 0, 1)
  est_sig <- ifelse(est == 0, 0, 1)
  pc_true <- pop_sig[upper.tri(pop_sig)]
  pc_est  <- est_sig[upper.tri(est_sig)]
  
  # Calculate number of true + false positives/negatives 
  TP <- sum(pc_true == 1   & pc_est == 1)
  TN <- sum(pc_true == 0   & pc_est == 0)
  FP <- sum(pc_true == 0   & pc_est == 1)
  FN <- sum(pc_true == 1   & pc_est == 0)
  
  # Calculate Specificity 
  Spec <- TN / (TN + FP)
  
  # Calculate False Positive Rate
  FPR <- 1 - Spec
  
  # Calculate Sensitivity 
  Sens <- TP / (TP + FN)
  
  # Calculate True Detection Rate
  TDR <- TP / (TP + FP)
  
  res <- rbind(TP,TN, FP, FN, Spec, FPR, Sens, TDR)
  return(res)
}

## Function: Create population matrix, simulate dataset, estimate networks, compute performance
mapply_func <- function(n, p, density , df){
  
  # Creates a sparse population correlation, partial correlation, and adjacency matrix
  mat_list <- true_mat_gen(n, p, density, df)
  
  # Simulate a dataset from the population correlation matrix
  d <- scale(MASS::mvrnorm(n=n, mu = rep(0,p), Sigma = (mat_list$cor_mat)))
  
  # Fit Model 
  h <- huge::huge(d, method = "glasso",scr = F, nlambda = 100, lambda.min.ratio = .01)
  
  ## Select Lambda with EBIC
  h_gl <- huge.select(h, criterion = "ebic")
  
  ## Select with StARS
  h_st  <- huge.select(h, criterion = "stars", stars.thresh = .1, 
                       stars.subsample.ratio = .8, rep.num = 20, verbose=F)
  ## Select with RIC
  h_ric <- huge.select(h, criterion = "ric", rep.num = 20, verbose=F)
  
  ## Select with CV
  rho_sel <- crossvalidaatioGlassolle(d, tune)
  cv <- glasso(cov(d), rho = rho_sel)
  
  
  #Compute Specificity + Sensitivity #####
  edges <- matrix(0, 8, 4)
  
  rownames(edges) <- c("TP", "TN", "FP", "FN", 
                       "Spec", "FPR", "Sens", 
                       "TDR")
  colnames(edges) <- c("ebic", "stars", "ric", "cv")
  
  #ebic
  edges[, 1]<- performance_measures(mat_list$sig_mat, h_gl$opt.icov)
  
  #stars
  
  edges[, 2] <- performance_measures(mat_list$sig_mat, h_st$opt.icov)
  
  #ric
  
  edges[, 3] <- performance_measures(mat_list$sig_mat, h_ric$opt.icov)
  
  
  #cv
  
  edges[, 4] <- performance_measures(mat_list$sig_mat, cv$wi)
  
  
  # Save selected lambdas from each penalty selection method
  
  lam_ebic <- h_gl$opt.lambda
  lam_stars <- h_st$opt.lambda
  lam_ric <- h_ric$opt.lambda
  lam_cv <- rho_sel
  lambdas <- cbind(lam_ebic, lam_stars, lam_ric, lam_cv)
  rownames(lambdas)<- "Lambda"
  
  
  # Calculate estimated density for each network
  
  total <- (p * (p - 1)) / 2 # calculates number of possible edges
  
  # number of true positives + number of false positives = number of estimated edges
  # number estimated edges/ number possible edges = density 
  
  ebic_den <- edges["TP","ebic"] + edges["FP","ebic"]/ total 
  stars_den <- edges["TP","stars"] + edges["FP","stars"]/ total 
  ric_den <- edges["TP","ric"] + edges["FP","ric"]/ total 
  cv_den <- edges["TP","cv"] + edges["FP","cv"]/ total 
  
  dens <- cbind(ebic_den, stars_den, ric_den, cv_den)
  rownames(dens)<- "Density"
  
  
  # Calculate Global Correlation: correlation between edge weights in population 
  #   and corresponding estimated edge weight 
  
  true_pc <- as.numeric(mat_list$pcor_mat[upper.tri(mat_list$pcor_mat)]) 
  
  # EBIC
  # opt.icov is the precision matrix from the estimated networks 
  ebic_pc <- as.numeric(prec2pc(h_gl$opt.icov)[upper.tri(h_gl$opt.icov)])
  ebic_g_cor <- cor(true_pc, ebic_pc)
  
  # StARS
  stars_pc <- as.numeric(prec2pc(h_st$opt.icov)[upper.tri(h_st$opt.icov)])
  stars_g_cor <- cor(true_pc, stars_pc)
  
  # RIC
  ric_pc <- as.numeric(prec2pc(h_ric$opt.icov)[upper.tri(h_ric$opt.icov)])
  ric_g_cor <- cor(true_pc, ric_pc)
  
  # CV
  cv_pc <- as.numeric(prec2pc(cv$wi)[upper.tri(cv$wi)])
  cv_g_cor <- cor(true_pc, cv_pc)
  
  global_cor <- cbind(ebic_g_cor, stars_g_cor, ric_g_cor, cv_g_cor)
  
  
  # Estimated Edge Correlation: Correlation between true edge weights in the population 
  # and corresponding estimated edge weight 
  
  # EBIC
  ebic_est_cor <- cor(true_pc[true_pc != 0], ebic_pc[true_pc != 0])
  
  # StARS
  stars_est_cor <- cor(true_pc[true_pc != 0], stars_pc[true_pc != 0])
  
  # RIC
  ric_est_cor <- cor(true_pc[true_pc != 0], ric_pc[true_pc != 0])
  
  # CV
  cv_est_cor <- cor(true_pc[true_pc != 0], cv_pc[true_pc != 0])
  
  est_cor <- cbind(ebic_est_cor, stars_est_cor, ric_est_cor, cv_est_cor)
  
  
  
  # Combine all performance measures 
  edges <- rbind(edges, lambdas, dens, global_cor, est_cor)
  rownames(edges) <- c("TP", "TN", "FP", "FN", 
                       "Spec", "FPR", "Sens", 
                       "TDR", "Lambda", "Density", "Global_Cor", "Est_Cor")
  
  ## reshape for storing nicely
  edges <- reshape2::melt(edges)
  
  return(edges)
}


# Function to label results with condition parameters 
sim_run <- function(n, p, density, pc){ 
  # mapply takes the function, mapply_fun, and the arguments
  # n and p. Simplify ensures things are returned nicely
  res_map <- mapply(mapply_func, n, p, density, pc, SIMPLIFY = F)
  # give condition columns to easily compute things later on
  for(i in 1:length(res_map)){
    
    res_map[[i]]$n <- n 
    res_map[[i]]$density<- density 
    res_map[[i]]$pc<- pc_range
  }
  return(res_map)
}
