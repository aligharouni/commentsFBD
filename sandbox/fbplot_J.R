## How to pick J>2 in fbplot?

library(fda)
## Creating toy example with trajectories in columns  
tvec <- seq(1,10,length=10)
y1 <- ifelse(tvec<5,0,1)
y2 <- ifelse(tvec<6,0,1)
y3 <- ifelse(tvec<7,0,1)

ensemble_toy <- cbind(y1,y2,y3)
# data <- ensemble_toy
matplot(tvec,ensemble_toy,type="l")

## Digging into fda::fbplot from url: https://rdrr.io/cran/fda/src/R/fbplot.R

## what I don't understand is why exp(log)?
combinat=function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))} ##this is exp(log(C(n,p)))
}

# #BD2
fBD2=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=apply(rmat,1,min)-1
  up=n-apply(rmat,1,max)
  (up*down+n-1)/combinat(n,2)

}
## modified fBD2, combinant(n,2) to combinant(n,J) with J to be determined eg, J=50 as in Juul's work 
fBD2_modified=function(data,J){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=apply(rmat,1,min)-1
  up=n-apply(rmat,1,max)
  (up*down+n-1)/combinat(n,J)
  
}

fBD2(ensemble_toy)
fBD2_modified(ensemble_toy,J=3)
