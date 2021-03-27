
###################################
## Stage 1; 
###################################
y <- 1:100

## cumulative sums
cumsum(y)

## peak
max(y)

## peak time
which(y==max(y))
which.max(y)

## epidemic duration (e.g. time between 10% and 90% cumulative cases)
y_temp <- c(0, replicate(10,1)) 
q_temp <- quantile(cumsum(y_temp),c(0.1,0.9))
epi_duration <- q_temp[["90%"]]-q_temp[["10%"]]

## epidemic duration (e.g. time between 10% and 90% cumulative cases)
epi_du.test <- function(curve){
  q <- quantile(cumsum(curve),c(0.1,0.9))
  return(q[["90%"]]-q[["10%"]])
  }

epi_du.test(y_temp)

t <- 0:100 
curve <- exp(t)
## estimating r (epidemic initial growth rate) will fit a exponential to a curve between t0 and t1 time interval
growth_rate.test <- function(curve,t){
  t0 <-t[1]
  t1 <-t[2]
  dat <- data.frame(y=curve[t0:t1], time=seq(0,t1-t0))
  m <- lm(log(y)~time,data=dat)
  r <- coef(m)[2]
  return(r)
}

growth_rate.test(curve,c(1,50))

## example 2, implementation of probes on an ensemble (matrix with curves on columns)
tvec <- seq(0,10,length=11)
y0 <- replicate(length(tvec),0)
y1 <- replicate(length(tvec),1)
y2 <- replicate(length(tvec),2)
ens_temp <- cbind(y0,y1,y2) ## ensemble: each column stores a trajectory

apply(ens_temp,2,FUN = epi_du.test)

apply(matrix(curve),2,FUN = function(curve_in)growth_rate.test(curve_in,c(1,5)))
###################################
## apply multiple functions to an ensemble
###################################

probes.test <- function(x) {
  ## x is a vector
  c(min = min(x), 
    peak = max(x), 
    peak_time=which.max(x), ## return the index of the first max it hits. Q: what to do with multiple peaks?
    mean = mean(x), 
    epi_duration = epi_du(x)
  )
}

probes_temp <- apply(ens_temp,2, probes.test)
t(probes_temp)

