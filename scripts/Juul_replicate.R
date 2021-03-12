
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
RnkEns(ensemble,envelope,weight=0)
###################################
## Stage 2: Form the subensemble by sampling from the whole ensemble
###################################
## take the ensemble, sample from it and call it subensemble, create the envelope of the subensemble (i.e., pointwise min-max), rank the curves in the ensemble (rank=1 if a curve is entirely in the envelope and 0 otherwise), for each sample store the ranks in a row of matrix "EnsRank" 
nTrajs <- ncol(ensemble) ## number of trajectories
numSample <- 4 ## number of sampling from the ensemble
sizeSample <- 2 ## sample size (must be >=2), i.e., how many curve to be sampled in each draw.

##form the matrix of samples
smat <- t(replicate(numSample,sample(1:ncol(ensemble), sizeSample, replace=FALSE)) )

## Allocate list and Matrices
Envlist <- vector(mode="list",length = numSample)
Envlist <- lapply(Envlist, FUN =function(i) matrix(NA,nrow=nrow(ensemble),ncol =2 ))

## Calculate the Envelope for each subensemble in smat:
Envlist <- EnvOut(ensemble, smat)

## initialize the rank matrix, EnsRank(i,j): rank of j'th curve of the ensemble wrt i'the envelope (subensemble or sample)
EnsRank <- matrix(NA,nrow = numSample,ncol = ncol(ensemble)) 
temp <- lapply(Envlist,FUN= function(envelope_in)RnkEns(ensemble,envelope_in,weight=0))
EnsRank <- matrix(unlist(temp),nrow = numSample,byrow = TRUE)
rm("temp")

matplot(1:ncol(ensemble),t(EnsRank))



###################################






