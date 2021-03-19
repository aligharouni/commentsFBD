## setwd("/home/ag/projects/AliMac_scripts/scripts/")
## To replicate centrality scoring work by Juul et al. 2021, NATURE PHYSICS VOL17
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
## Juul et al., Fig. 2 panel b; 
## the centrality of curves was ranked in different ways: all-or-nothing ranking for the full predicted time interval, Ncurves = 50 and Nsamples = 100
###################################
set.seed(1234)
EnRkTemp_J <- EnsRank_all(ensemble_J,numSample=100,sizeSample=50) ##  as Juul's et al did.
##saveRDS(EnRkTemp_J, file = "Juul_replicate.rds")

rnkmat <- EnRkTemp_J[["ensembRank"]] ## take the resulted rank mat

md1 <- rank(colSums(rnkmat)) ## , ties="random")
plot(sort(md1)) ## Note the discontinuity!?
central_curves <- which(md1>quantile(md1,0.5))

## central envelope, ie POINTWISE min & max of these curves
fmin <- apply(ensemble_J[,central_curves],1,min)
fmax <- apply(as.matrix(ensemble_J[,central_curves]),1,max)

g_labels <- c("trajectories", "fixed-time 25th-75th percentiles")
g_cols <- c("gray","black")

par(las=1, bty="l")
matplot(x=tvec,y=ensemble_J,type = "l", col = g_cols[1], lty=1,
        ylim=c(0,820),
        main = "Curve boxplot (Juul's et al. Fig. 2(b))",
        xlab = "Day",
        ylab = "Newly hospitalized"
        )
## lines(tvec,rowMeans(ensemble_J), col = g_cols[2], lwd = 2) ## fixed-time mean
## lines(tvec,apply(ensemble_J,1,median), col = g_cols[3], lwd = 2) ## fixed-time median
polygon(c(tvec,rev(tvec)),c(fmin,rev(fmax)), col=adjustcolor("blue",alpha.f=0.2),
        border=NA)
legend("topleft", legend = g_labels,
       col= g_cols, lty=1, bg="white", cex=0.6, bty='n') 

mpt <- t(apply(ensemble_J,1,quantile, c(0.25,0.75)))
matlines(mpt,col="black",lty=1,lwd=2)



# Note: A test run can be found at sandbox/pre_Juul_replicate.R



