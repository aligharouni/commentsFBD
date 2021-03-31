
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
which(cumsum(y_temp)==q_temp[["10%"]])
epi_duration <- which(cumsum(y_temp)==q_temp[["90%"]])-which(cumsum(y_temp)==q_temp[["10%"]])

## epidemic duration (e.g. time between 10% and 90% cumulative cases)
epi_du.test <- function(curve){
  total_sz <- sum(curve)
  cdf <- cumsum(curve)
  t1 <- which(cdf >= 0.1*total_sz)[1]
  t2 <- which(cdf >= 0.9*total_sz)[1]
  ## the quantile approach didn't work properly
  # cdf <- cumsum(curve)
  # q <- quantile(cdf,prob=c(0.1,0.9))
  # t2 <- which(cdf>=q[["90%"]])[1] 
  # t1 <- which(cdf==q[["10%"]])
  delta <- t2-t1
  return(delta)  
}

epi_du.test(y_temp)
epi_du.test(replicate(20,0))

t <- 0:100 
# curve <- t
curve <- exp(t)
# curve <- c(replicate(10,1),exp(1*t))
plot(curve)
epi_du.test(curve)


## estimating r (epidemic initial growth rate) will fit a exponential to a curve between t0 and t1 time interval
# growth_rate.test <- function(curve,t){
#   t0 <-t[1]
#   t1 <-t[2]
#   dat <- data.frame(y=curve[t0:t1], time=seq(0,t1-t0))
#   m <- lm(log(y)~time,data=dat)
#   r <- coef(m)[2]
#   return(r)
# }

## Estimating the initial growth rate when the prevalence is 10%;
prevalence <- 0.1
growth_rate.test <- function(curve,prevalence){
  total_sz <- sum(curve)
  cdf <- cumsum(curve)
  t1 <- which(cdf > 0)[1] ## start the fit from the first day of nonzero prevalence 
  t2 <- which(cdf >= prevalence*total_sz)[1]
  ## Quantile approach which didn't work properly
  # q <- quantile(cumsum(curve),prob=prevalence,type = 9)
  # t <- which(cumsum(curve)>=unname(q))[1] 
  dat <- data.frame(y=curve[t1:t2], time=seq(0,t2-t1))
  m <- lm(log(y)~time,data=dat) ## 
  r <- coef(m)[2]
  return(r)
}

growth_rate.test(curve,prevalence=0.1)

## example 2, implementation of probes on an ensemble (matrix with curves on columns)
tvec <- seq(0,100,length=101)
# y0 <- replicate(length(tvec),0)
y1 <- tvec
y2 <- tvec^2
y3 <- exp(0.1*tvec)
ens_temp <- cbind(y1,y2,y3) ## ensemble: each column stores a trajectory
matplot(ens_temp)

apply(ens_temp,2,FUN = epi_du.test)

apply(matrix(curve),2,FUN = function(curve_in)growth_rate.test(curve_in,prevalence=0.1))
apply(ens_temp,2,FUN = function(curve_in)growth_rate.test(curve_in,prevalence=0.1))

##Juul's example
curve <- as.vector(unlist(ensemble_J[,1]))
plot(curve)
plot(cumsum(curve))
growth_rate.test(curve,prevalence=0.1)
epi_du.test(curve)
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

