
source('functions.R')
library(tidyverse)
library(readr)
library(latex2exp)

g_labels <- c("trajectories", "fixed-time 25th-75th percentiles")
g_cols <- c("gray","black")
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

par(las=1, bty="l")
matplot(x=tvec,y=ensemble_J,type = "l", col = g_cols[1], lty=1,
        ylim=c(0,820),
        main = TeX("Curve boxplot with $L_2$ norm"),
        xlab = "Day",
        ylab = "Newly hospitalized"
)

polygon(c(tvec,rev(tvec)),c(fminl2,rev(fmaxl2)), col=adjustcolor("blue",alpha.f=0.2),border=NA)
legend("topleft", legend = g_labels,
       col= g_cols, lty=1, bg="white", cex=0.6, bty='n') 
mpt <- t(apply(ensemble_J,1,quantile, c(0.25,0.75)))
matlines(mpt,col="black",lty=1,lwd=2)

###################################
## Functional Boxplot (FDA) on Juul's data
###################################
## I think this is similar to Juul's et al. method but instead the sample size is 2 at each draw 
## Data depth is a well-known and useful nonparametric tool for analyzing functional data. It provides a novel way of ranking a sample of curves from the center outwards and defining robust statistics, such as the median or trimmed means.
library("fda")
library("roahd")

fD <- fData(tvec,as.matrix(t(ensemble_J))) ## note that the trajs are on rows in fD object
# dev.new()
ff <- fbplot(fD,method='MBD', plot=TRUE,
             outliercol = 2, 
             fullout = FALSE,
             outline=FALSE,
             main = "functional boxplot, fbplot()",
             xlab = "Day",
             ylab = "Newly hospitalized" 
             ) 
## if plot=false, it outputs the calculations and scores


###################################
## mahalanobison probes on Juul's data
###################################
class(ensemble_J)

n_probes <- 4
# probes.fun <- function(x) {
#   return(c(min = min(x), mean = median(x), max = max(x), curvemean=rowMeans(x)))
# }

probes <- matrix(NA,n_traj,n_probes) 
probes[,1] <- apply(t(ensemble_J), 1, FUN = min)
probes[,2] <- apply(t(ensemble_J), 1, FUN = max)   
probes[,3] <- apply(t(ensemble_J), 1, FUN = median)   
probes[,4] <- rowMeans(t(ensemble_J))


Sigma <- cov(probes)
mahalmat <- matrix(NA, n_traj, n_traj)
diag(mahalmat) <- 0

for (i in 1:(n_traj-1)) {
  for (j in (i+1):n_traj) {
    ## calculating pairwise M distance of both trajectories from their mean
    ## (identical, take the first one)
    mahalmat[i,j] <- mahalanobis(probes[c(i,j),], #taj i and j
                                  center=colMeans(probes[c(i,j),]),
                                  cov=Sigma)[1]
  }
}

mahalmat <- make_symm(mahalmat)
md_mah <- rank(rowMeans(mahalmat))
central_curves_mah <- which(md_mah<quantile(md_mah,0.5))

## central envelope, ie POINTWISE min & max of these curves
fmin_mah <- apply(ensemble_J[,central_curves_mah],1,min)
fmax_mah <- apply(as.matrix(ensemble_J[,central_curves_mah]),1,max)

par(las=1, bty="l")
matplot(x=tvec,y=ensemble_J,type = "l", col = g_cols[1], lty=1,
        ylim=c(0,820),
        main = "Mahalanobis norm on the probes",
        xlab = "Day",
        ylab = "Newly hospitalized"
)
polygon(c(tvec,rev(tvec)),c(fmin_mah,rev(fmax_mah)), col=adjustcolor("blue",alpha.f=0.2),
        border=NA)
legend("topleft", legend = g_labels,
       col= g_cols, lty=1, bg="white", cex=0.6, bty='n') 

mpt <- t(apply(ensemble_J,1,quantile, c(0.25,0.75)))
matlines(mpt,col="black",lty=1,lwd=2)



