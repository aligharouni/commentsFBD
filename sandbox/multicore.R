## Goal: To speed up run time! 
## Run a loop  in parallel mode, or usingparallel::parLapply
## I am looking at this, but ask for better ref? https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745

library(parallel)
numCores <- detectCores() # get the number of cores available


# simple_example1 ---------------------------------------------------------

# Generate data
data <- 1:1e3
data_list <- list("1" = data,
                  "2" = data,
                  "3" = data,
                  "4" = data)

# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores())
# Run parallel computation
time_parallel <- system.time(
  mean_data <- parallel::parLapply(cl,
                                   data_list,
                                   mean)
)
# Close cluster
parallel::stopCluster(cl)


# simple _example2 --------------------------------------------------------
## a pilot example to select trajs from a whole traj matrix and compute a pairwise dist between the subset and the whole set
n_traj0 <- 4 ## 3 trajectories, 2 data points
n_point0 <- 2
traj <- matrix(NA,n_point0,n_traj0) ## trajectories are in columns
traj[,1] <- c(0,0)
traj[,2] <- c(1,1)
traj[,3] <- c(0,1)
traj[,4] <- c(0,3)


## only select col 1:2 and col 3
data_list <- list("sub1" = subset.matrix(traj,select = c(1:n_traj0)) ,
                  "sub2" = subset.matrix(traj,select = c(2:n_traj0)),
                  "sub3" = subset.matrix(traj,select = c(3:n_traj0))
                  )




## single core
system.time(out_list <- sapply(data_list,dtraj))

## 8 cores
system.time({
  # Detect the number of available cores and create cluster
  cl <- parallel::makeCluster(detectCores())
  clusterExport(cl, 'Frechet') ##to export the variables you need to the clusters
  # Run parallel computation
  out_list_p <- parallel::parSapply(cl,data_list,dtraj)
})

# Close cluster
parallel::stopCluster(cl)

## Ali; Hummm that is odd, the multicore computation is slower than the sapply! what am I missing?

# automating_data_list ----------------------------------------------------
## automate making the data_list as a list with # of trajs elements
n_traj <- ncol(traj) # number of trajectories stored in columns of traj matrix
test_data <- sapply(1:(n_traj-1),FUN=function(i) subset.matrix(traj,select = c(i:n_traj)))
system.time(out_list_test <- sapply(test_data,dtraj))


