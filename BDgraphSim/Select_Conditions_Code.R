#################################
###     Select Degrees of     ###
###  Freedom for  BDgraph     ###
###  Simulation Conditions    ###
#################################
library(RCurl)
source("https://raw.githubusercontent.com/AnnaWysocki/Network-Penalty-Selection/master/BDgraphSim/SelectConditions_Functions.R")

## Vector of possible degrees of freedom
df_vec <- seq(3, 50, by =1)

## Vector of levels of density included in simulation
dens_vec <- seq(.1, .9, by = .2)

## For loop computes the average partial correlation range 
##  for each degree of freedom at each level of density

## Data frame for storing results 
range_df <- data.frame(NULL)

for(i in 1:length(dens_vec)){
  print(i)
  dens <- dens_vec[i] # select a density
  for(j in 1:length(df_vec)){
    df <- df_vec[j] # select a degree of freedom
    
    # calculate the average range for that df and density across 1000 networks 
    r <- av_interval(p = 20, dens = dens, df= df, rep_num = 1000) 
    range_df <- rbind(range_df, cbind(r[1], r[2], df, dens)) 
  }
}


## Select df for abs(.25) range 

smallpc_conditions <- data.frame(NULL)

for(i in 1: length(dens_vec)){
  
  # Selects results for one level of sparsity
  dat <- range_df[range_df$dens == dens_vec[i],] 
  
  # selects degree of freedom that has the smallest average discrepancy
  # from .25 
  best_row <- which(abs(dat$V2-.25) == min(abs(dat$V2-.25)))
  smallpc_conditions<- rbind(smallpc_conditions, dat[best_row,])
}


## Select df for abs(.35) range

largepc_conditions <- data.frame(NULL)

for(i in 1: length(dens_vec)){
  dat <- range_df[range_df$dens == dens_vec[i],] 
  best_row <- which(abs(dat$V2-.35)==min(abs(dat$V2-.35)))
  largepc_conditions<- rbind(largepc_conditions, dat[best_row,])
}


########################################
# Create Condition Grid for Simulation #
########################################


## Sample Size Conditions
n <- c(100, 200, 250, 500, 1000, 1500)

## Number of Variables Conditions 
p <- 20


## Large Conditions
  ## Create dataframe where each row represents a simulation condition containing
  ## the necessary parameters
  large_conditions<- expand.grid(largepc_conditions$dens, n, p, pc_range= .35)
  
  # Add degrees of freedom to dataframe
  large_conditions <- cbind(largepc_conditions$df, large_conditions)
  
  
  colnames(large_conditions) <- c("df", "density", "n", "p", "pc_range")

## Small Conditions 
small_conditions<- expand.grid(smallpc_conditions$dens, n, p, pc_range= .25)
small_conditions <- cbind(smallpc_conditions$df, small_conditions)
names(small_conditions) <- names(large_conditions)

## Combine small and large conditions to one data frame 
conditions <- rbind(small_conditions, large_conditions)

save(conditions, file = "Bdgraph_Conditions")


