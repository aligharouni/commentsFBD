# setwd("/home/ag/projects/AliMac_scripts/scripts/")
source('functions.R')
###################################
## Juul's work
###################################

M <- (read_csv("./data/juul1.csv", col_names = FALSE)) ## the columns are the trajs
ensemble_J <- M
tvec <- seq(0,nrow(ensemble_J)-1)

set.seed(1234)
EnRkTemp_J <- EnsRank_all(ensemble_J,numSample=100,sizeSample=50) ## sizeSample=50 (Juul's N_curves), numSample=100 (N_samples)

##saveRDS(EnRkTemp_J, file = "Juul_replicate.rds")

rnkmat <- EnRkTemp_J[["ensembRank"]] ## rank mat

barplot(colSums(rnkmat)) 
md1 <- rank(colSums(rnkmat))

central_curves <- which(md1<=quantile(md1,0.5)) ## Question: md1<?

## central envelope, ie POINTWISE min/max of these curves
fmin <- apply(ensemble_J[,central_curves],1,min)
fmax <- apply(as.matrix(ensemble_J[,central_curves]),1,max)


g_labels <- c("trajectories", "fixed-time mean",
              "fixed-time median", "Juul central")
g_cols <- c("gray","blue","purple","red")

matplot(x=tvec,y=ensemble_J,type = "l", col = g_cols[1], lty=1,ylim=c(0,800))
lines(tvec,rowMeans(ensemble_J), col = g_cols[2], lwd = 2)
lines(tvec,apply(ensemble_J,1,median), col = g_cols[3], lwd = 2)
polygon(c(tvec,rev(tvec)),c(fmin,rev(fmax)), col=adjustcolor("black",alpha.f=0.2),
        border=NA)
legend("topright", legend = g_labels,
       col= g_cols, lty=1, bg="white") 




