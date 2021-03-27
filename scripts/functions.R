

## Make a matrix, symetric
## input: matrix x
## symmetrize (is there something built-in? Matrix::forceSymmetric
##  converts to a fancy Matrix object ...
make_symm <- function(x){
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  return(x)
}

## Definition of log Euclidean norm in R following this paper: arsigny2006log
## For the purpose of disproving our conjecture about the central curve.
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

get_envelope <- function(ensemble, indices) {
  data.frame(tvec=seq(nrow(ensemble))-1,
             lwr=apply(ensemble[,indices],1,min),
             upr=apply(ensemble[,indices],1,max))
}

##########################################
## FUNCTIONs FOR REPLICATING JUUL'S WORK
##########################################
## output a list of envelopes (pointwise min-max) of a subensemble 
EnvOut <- function(ensemble,subEnsMat){
  ## ensemble: the whole ensemble, set of curves, with the curves on columns 
  ## We sample n times from the ensemble and the info of all samples are stored in subEnsMat as follows.  
  ## subEnsMat: subensemble matrix with the (i,j) element as i'th sample & j'th curve of the ensemble
  envelopes <- vector(mode="list",length = nrow(subEnsMat)) ## initializing
  for(i in 1:(nrow(subEnsMat))){
    subens <- ensemble[,subEnsMat[i,]]
    pmin_subens <- apply(subens, 1, FUN=min)
    pmax_subens <- apply(subens, 1, FUN=max)
    envelopes[[i]] <- cbind(pmin_subens,pmax_subens) ## temporary envelope
  }
  return(envelopes)
}

## A  function to determine whether the curves of an ensemble fall entirely in an given envelope: 
IsInEnv <- function(ensemble,envelope){
  ## isinenv() returns  a logi vector; TRUE if the input curve is entirely inside the envelope.
  ## the envelope is a n by 2 matrix, determining the lower bound (envelope[,1]) and the upper bound (envelope[,2]) of an envelope. 
    ## o <- apply(ensemble,2,FUN = function(traj){traj >= envelope[,1] && traj <= envelope[,2]})
  o <- apply(ensemble,2,FUN = function(traj){ all(traj >= envelope[,1] & traj <= envelope[,2]) } )
  return(o)
}

## rank the curves in an ensemble based on the whether it is in an subensemble
RnkEns_1sample <- function(ensemble,envelope) {
  ## ensemble: matrix of whole ensemble with curves on the column
  ## envelope: list of envelopes corresponding to each subensemble
  ## could calculate a weighted value instead?
  ## would need to add weights back in
  ## Juul et al.; they compute a 'centrality score', i.e.
  ##  curves that fall within the envelope get 1, outside -> 0
  fall <- IsInEnv(ensemble,envelope)
  return(ifelse(fall,1,0))
  }

## take the ensemble, sample from it and call it subensemble, create the envelope of the subensemble (i.e., pointwise min-max), rank the curves in the ensemble (rank=1 if a curve is entirely in the envelope and 0 otherwise), for each sample store the ranks in a row of matrix "EnsRank" 
EnsRank_all <- function(ensemble,numSample,sizeSample){
  ## Inputs: 
  ## ensemble: matrix of whole ensemble with curves on the column.
  ## numSample: number of sampling from the ensemble, integer. 
  ## sizeSample: sample size (must be >=2), i.e., how many curve to be sampled in each draw, integer.
  ## Output: 
  ## matrix EnsRank, with EnsRank(i,j): rank of j'th curve of the ensemble wrt i'the envelope (subensemble or sample)
  nTrajs <- ncol(ensemble)
  ## form the matrix of samples
  smat <- t(replicate(numSample,sample(1:ncol(ensemble), sizeSample, replace=FALSE)) )
  ## Allocate a list of envelopes for each subensemble
  Envlist <- vector(mode="list",length = numSample)
  Envlist <- lapply(Envlist, FUN =function(i) matrix(NA,nrow=nrow(ensemble),ncol =2 ))
  
  ## Calculate the Envelope for each subensemble in smat:
  Envlist <- EnvOut(ensemble, smat)
  
  ## initialize the rank matrix, EnsRank(i,j): rank of j'th curve of the ensemble wrt i'the envelope (subensemble or sample)
  EnsRank <- matrix(NA,nrow = numSample,ncol = ncol(ensemble)) 
  ## Rank the ensemble across all subensembles
  temp <- lapply(Envlist,FUN=function(envelope_in) RnkEns_1sample(ensemble,envelope_in))
  ## output a list of sample matrix (smat) and ranks
  EnsRank <- matrix(unlist(temp),nrow = numSample,byrow = TRUE)
  return(list(sampMat=smat,ensembRank=EnsRank))
}

##########################################
## Probes FUNCTIONs
##########################################

## epidemic duration (e.g. time between 10% and 90% cumulative cases)
epi_du <- function(curve){
  q <- quantile(cumsum(curve),c(0.1,0.9))
  return(q[["90%"]]-q[["10%"]])
}

## estimating r (epidemic initial growth rate) will fit a exponential to a curve between t0 and t1 time interval
growth_rate <- function(curve,t){
  t0 <-t[1]
  t1 <-t[2]
  dat <- data.frame(y=curve[t0:t1], time=seq(0,t1-t0))
  m <- lm(log(y)~time,data=dat)
  r <- coef(m)[2]
  return(r)
}

probes <- function(x) {
  ## x is a vector
  c(min = min(x), 
    peak = max(x), 
    peak_time=which.max(x), ## return the index of the first max it hits. Q: what to do with multiple peaks?
    mean = mean(x), 
    epi_duration = epi_du(x),
    r_init = growth_rate(x,c(1,50)) ## growth rate in initial 50 days of the pandemic
  )
}

