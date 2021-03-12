
# setwd("/home/ag/projects/AliMac_scripts/scripts/")
source('functions.R')

###################################
## Stage 1; Given an envelope, rank the ensemble (matrix of curves) starting from weight=0
###################################
## simple example of 3 trajectories 
tvec <- seq(0,10,length=11)
y0 <- replicate(length(tvec),0)
y1 <- replicate(length(tvec),1)
y2 <- replicate(length(tvec),2)

ensemble <- cbind(y0,y1,y2) ## ensemble: each column stores a traj
subens <- ensemble[,1:2] ##subensemble
matplot(tvec,ensemble,type="l")

## defining the envelop of an ensemble as pointwise min max
pmin_subens <- apply(subens, 1, FUN=min)
pmax_subens <- apply(subens, 1, FUN=max)
envelope <- cbind(pmin_subens,pmax_subens) ## temporary envelope

## test
IsInEnv(as.matrix(y2),envelope) ## a single trajectory
test <- IsInEnv(ensemble,envelope) ## a ensemble
names(test)[which(test=="TRUE")]

## output the rank of the input an ensemble, given 1 envelope
RnkEns_1sample(ensemble,envelope,weight=0)
##test
EnRkTemp <- EnsRank_all(ensemble,numSample=20,sizeSample=2)

rnkmat <- EnRkTemp[["ensembRank"]] ## rank mat

barplot(colSums(rnkmat)) 
md1 <- rank(colSums(rnkmat))

central_curves <- which(md1<quantile(md1,0.5))

## central envelope, ie POINTWISE min/max of these curves
fmin <- apply(as.matrix(ensemble[,central_curves]),1,min)
fmax <- apply(as.matrix(ensemble[,central_curves]),1,max)


g_labels <- c("trajectories", "fixed-time mean",
              "fixed-time median", "Juul central")
g_cols <- c("gray","blue","purple","red")

matplot(x=tvec,y=ensemble,type = "l", col = g_cols[1], lty=1,ylim=c(0,3))
lines(tvec,rowMeans(ensemble), col = g_cols[2], lwd = 2)
lines(tvec,apply(ensemble,1,median), col = g_cols[3], lwd = 2)
polygon(c(tvec,rev(tvec)),c(fmin,rev(fmax)), col=adjustcolor("black",alpha.f=0.2),
        border=NA)
legend("topright", legend = g_labels,
       col= g_cols, lty=1, bg="white") 

