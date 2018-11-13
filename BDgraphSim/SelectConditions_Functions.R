# Necessary packages
library(BDgraph)

## Precision matrix to partial correlation matrix 

 # Inverting a precision matrix, reversing the signs
 # on the off-diagonals and standardizing the matrix
 # creates a partial correlation matrix 

prec2pc <- function(prec_mat){
  prec_mat <- as.matrix(prec_mat)
  neg_prec_mat <- - prec_mat
  diag(neg_prec_mat) <- diag(prec_mat)
  pc_mat<- cov2cor(neg_prec_mat)
  return(pc_mat)
}

## Calculates range that 90% of the partial correlations fall into 
calc_range <- function(p, dens, df){
  network <- bdgraph.sim(p=p, n = 100, prob= dens, b = df) 
  
  # K is the generated precision matrix, G is the adjacency matrix 
  # By converting K to a precision matrix and multiplying it by G
  # We get a sparse partial correlation matrix 
  pop_mat <- prec2pc(network$K) * network$G 
 
   # Remove 0s so we move forward with only the estimated edges 
  pop_mat <- pop_mat[pop_mat != 0] 
  
  interval <- quantile(pop_mat, probs = c(.05, .95)) 
  return(interval)}


## Calculates the *average* range that 90% of partial correlations fall into 
## for a specific degree of freedom at a specific level of density
av_interval <- function(p, dens, df, rep_num){
  
  # Calculates PC range for a single network and replicates this
  # process rep_num times with the same degree of freedom and 
  # density 
  int_list <- unlist(replicate(rep_num, calc_range(p, dens, df)))
  
  up_low <- apply(int_list, 1 , mean)
  return(up_low)
}
