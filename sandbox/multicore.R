## Goal: To speed up run time! 
## Run a loop  in parallel mode, or usingparallel::parLapply
## I am looking at this, but ask for better ref? https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745

library(parallel)
numCores <- detectCores() # get the number of cores available

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
  parallel::parLapply(cl,
                      data_list,
                      mean)
)

mean_data <- parallel::parLapply(cl,
                    data_list,
                    mean)

# Close cluster
parallel::stopCluster(cl)



