

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




