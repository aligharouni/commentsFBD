

source('functions.R')
library(tidyverse)
library(readr)

###################################
## Juul's work replicating panel b;
###################################
## upload Fig.2 ensemble  
ensemble_J <- (read_csv("./data/juul1.csv", col_names = FALSE)) ## the columns are the trajs
tvec <- seq(0,nrow(ensemble_J)-1) ## time domain

###################################
## Pointwise L2 Norm on Juul's data
###################################
## The distance between 2 curves is the area between them

n_traj <- ncol(ensemble_J)
t <- 1:nrow(ensemble_J)

## Memory allocation for distance matrix
dmatl2 <- matrix(NA, n_traj, n_traj) ## Euclidean L2 norm
dmatl2 <- as.matrix(dist(t(ensemble_J), method = "euclidean", upper = TRUE, p = 2))
stopifnot(isSymmetric(dmatl2))

## l2 mean corresponding to arithmetic mean: minimize sum of squares of distances
md2l2 <- rowSums(dmatl2^2)  
m2l2 <- which.min(md2l2)
## corresponding to median: minimize sum of distances 
md1l2 <- rowSums(dmatl2)
m1l2 <- which.min(md1l2)
## Correlated mean and median?
cor.test(md1l2,md2l2,method="spearman")$estimate

plot(sort(md1))

central_curves_l2 <- which(md1l2<quantile(md1l2,0.5))
## Find the envelope, i.e., the POINTWISE min/max of the central curves 
fminl2 <- apply(ensemble_J[,central_curves_l2],1,min)
fmaxl2 <- apply(ensemble_J[,central_curves_l2],1,max)







