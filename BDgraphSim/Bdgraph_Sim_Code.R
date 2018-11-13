#################################
###     BDGraph Simulation    ###
###  for Penalty Selection    ###
###          Project          ###
#################################

setwd("C:/Users/annawy/Box Sync/Research/Davis/Penalty Parameter Selection/R/BDgraph")

source("Bdgraph_Sim_Functions.R")

tune <- seq(0.001, 1, by = .01) # Create tuning interval for CV
results <- list() # Create list to save results
sim <- 1 # simulation iterations

# Load Condition Data Frame

load("Bdgraph_Conditions.RData") # loads in data frame named "conditions" 

# Carry out simulation
for(j in 1: (nrow(conditions)/10)){
  print(j)
  
  # Select condition parameters 
  n <- conditions$n[j]
  p <- conditions$p[j]
  density <- conditions$density[j]
  df <- conditions$df[j]
  pc_range <- conditions$pc_range[j]
  
  # Repeat sim_run function (simulate data, estimate networks, assess performance) sim times
  int_res <- replicate(sim, sim_run(n, p , density, pc))
  
  results<- c(results, int_res)
  
}

# List to Dataframe 
bdgraph_results <- do.call(rbind, results)

save(bdgraph_results,  file= 'BDgraph_Results.RData')


??cor2cov
