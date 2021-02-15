

## Make a matrix, symetric
## input: matrix x
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
