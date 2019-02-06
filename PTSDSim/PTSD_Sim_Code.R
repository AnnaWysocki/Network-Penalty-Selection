########################
##    Nested PTSD     ## 
##  Simulation Code   ##
########################



## Load packages and functions 

library(RCurl)

source('https://raw.githubusercontent.com/AnnaWysocki/Network-Penalty-Selection/master/PTSDSim/PTSD_Sim_Functions.R')


## Load Data
data_ptsd <- read.csv(text= getURL("https://raw.githubusercontent.com/AnnaWysocki/Network-Penalty-Selection/master/PTSDSim/PTSD_Data.csv"))



## Create partial correlation bank from PTSD dataset 
dat <- cor2pcor(cor(na.omit(data_ptsd)))[upper.tri(cor2pcor(cor(na.omit(data_ptsd))))]

## Remove partial correlations that are less than abs(.05)
dat <- subset(dat, abs(dat) > 0.05)

#number of simulations
sim <- 1000
nested_sim <- 1000

#Create conditions/condition matrix 
n <- c(100, 200, 250, 500, 1000, 2000, 3000)
p <- c(10, 20)
sparsity <- c (.50, .65, .80)
conditions <- expand.grid(n, p, sparsity)
iteration <- as.list(1:sim)


#Create list for storing results
results <- list()

#Simulation

for (j in 1:nrow(conditions)){
  
  n <- conditions$Var1[j]
  p <- conditions$Var2[j]
  sparsity <- conditions$Var3[j]
   
  
  for (i in 1: sim){
  
  # Create population matrix 
  mat_list <- true_mat_gen(n, p, sparsity)
  network <- i
  
  # Creates nested_sim number of datasets using the same population matrix
  # Fits network to each dataset and saves performance indices 
  int_res <- replicate(nested_sim, sim_run(n, p , sparsity, network))
  
  # Adds column with network number 
  int_res <- Map(cbind, int_res, iteration = iteration)
  results<- c(results, int_res)}
  
  
  }

#Create and save results dataframe 
ptsd_results <- do.call(rbind.data.frame, results)
save(ptsd_results,  file = 'NestedPTSD_Results.RData')

