

## Make a matrix, symetric
## input: matrix x
## symmetrize (is there something built-in? Matrix::forceSymmetric
##  converts to a fancy Matrix object ...
make_symm <- function(x){
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  return(x)
}

## Definition of log Euclidean norm in R following the 
## paper: arsigny2006log
logedu <- function(x,y){
  if(x<=0 | y<=0) return(NA)
  else
    return(abs(log10(x)-log10(y)))
} 

## pairwise distance
## to do pass the distance function such as Frechet, DTW, etc

## dtraj: Calculate all of the pairwise distances between first traj (col1) and all others 
dtraj <- function(submatrix){
  ## compare the first column with the other columns
  out <- apply(as.matrix(submatrix[,-1]),2,FUN=function(x) SimilarityMeasures::Frechet(subset.matrix(submatrix,select = 1),matrix(x)) )
  return(out)
}

##########################################
## FUNCTIONs FOR REPLICATING JUUL'S WORK
##########################################
## output a list of envelopes (pointwise min-max) of a subensemble 
EnvOut <- function(ensemble,subEnsMat){
  ##ensemble: the whole ensemble, set of curves, with the curves on columns 
  ##subEnsMat: subensemble matrix with the (i,j) element as i'th sample & j'th curve of the ensemble
  envelopes <- vector(mode="list",length = nrow(subEnsMat)) ## initializing
  for(i in 1:(nrow(subEnsMat))){
    subens <- ensemble[,subEnsMat[i,]]
    pmin_subens <- apply(subens, 1, FUN=min)
    pmax_subens <- apply(subens, 1, FUN=max)
    envelopes[[i]] <- cbind(pmin_subens,pmax_subens) ## temporary envelope
  }
  return(envelopes)
}

## A  function to determine whether the curves of an ensemble fall entirely in an envelope: 
IsInEnv <- function(ensemble,envelope){
  ## isinenv() returns  a logi vector; TRUE if the input curve is entirely inside the envelope.
  ## the envelope is an n by 2 matrix, determining the lower/upper bounds of an envelope. 
  o <- apply(ensemble,2,FUN = function(traj){traj >= envelope[,1] && traj <= envelope[,2]})
  return(o)
}

## rank the curves in an ensemble based on the whether it is in an subensemble
RnkEns <- function(ensemble,envelope,weight=0){
 ##weight <- 0 ## initial rank assigned to each curve in the ensemble
  inc <- 1 # increase the rank of a curve by inc (increment) if it is outside of the evelope
  nTrajs <- ncol(ensemble) ## number of trajectories in the ensemble 
  rank <- replicate(nTrajs,weight) ## initialize the ranks of the ensemble
  fall <- IsInEnv(ensemble,envelope)
  rank <- ifelse(fall,weight,weight+inc)
  return(rank)
}

lapply(Envlist,FUN= function(envelope_in)RnkEns(ensemble,envelope_in,weight=0))


##########################################




