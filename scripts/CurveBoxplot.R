
source('functions.R')
library(tidyverse)
library(latex2exp)
library(fda)
library(roahd)

envelope_list <- list()


###################################
## Juul's work replicating panel b;
###################################
## upload Fig.2 ensemble  
ensemble_J <- read_csv("./data/juul1.csv", col_names = FALSE) ## the columns are the trajs
tvec <- seq(nrow(ensemble_J))-1 ## time domain

get_envelope <- function(ensemble, indices) {
    data.frame(tvec=seq(nrow(ensemble))-1,
               lwr=apply(ensemble[,indices],1,min),
               upr=apply(ensemble[,indices],1,max))
}


###################################
## Juul et al., Fig. 2 panel b; 
## the centrality of curves was ranked in different ways: all-or-nothing ranking for the full predicted time interval, Ncurves = 50 and Nsamples = 100
###################################
set.seed(1234)
EnRkTemp_J <- EnsRank_all(ensemble_J,numSample=100,sizeSample=50) ##  as Juul's et al did.
##saveRDS(EnRkTemp_J, file = "Juul_replicate.rds")

rnkmat <- EnRkTemp_J[["ensembRank"]] ## take the resulted rank mat

md1 <- rank(colSums(rnkmat)) ## , ties="random")
plot(sort(md1)) ## Note the discontinuity!?
central_curves <- which(md1>quantile(md1,0.5))

envelope_list <- c(envelope_list,
                   list(juul=get_envelope(ensemble_J, central_curves)))

###################################
## Pointwise L2 Norm on Juul's data
###################################
## The distance between 2 curves is the area between them

n_traj <- ncol(ensemble_J)
t <- 1:nrow(ensemble_J)

## Memory allocation for distance matrix
dmatl2 <- matrix(NA, n_traj, n_traj) ## Euclidean L2 norm
dmatl2 <- as.matrix(dist(t(ensemble_J), method = "euclidean", upper = TRUE, p = 2))
stopifnot(isSymmetric(dmatl2))

## l2 mean corresponding to arithmetic mean: minimize sum of squares of distances
md2l2 <- rowSums(dmatl2^2)  
m2l2 <- which.min(md2l2)
## corresponding to median: minimize sum of distances 
md1l2 <- rowSums(dmatl2)
m1l2 <- which.min(md1l2)
## Correlated mean and median?
cor.test(md1l2,md2l2,method="spearman")$estimate

central_curves_l2 <- which(md1l2<quantile(md1l2,0.5))
## Find the envelope, i.e., the POINTWISE min/max of the central curves


envelope_list <- c(envelope_list,
                   list(L2norm=get_envelope(ensemble_J, central_curves_l2)))

###################################
## Functional Boxplot (FDA) on Juul's data
###################################

## I think this is similar to Juul's et al. method but instead the sample size is 2 at each draw 
## Data depth is a well-known and useful nonparametric tool for analyzing functional data. It provides a novel way of ranking a sample of curves from the center outwards and defining robust statistics, such as the median or trimmed means.

fD <- fData(tvec,as.matrix(t(ensemble_J))) ## note that the trajs are on rows in fD object
# dev.new()
ff <- fbplot(fD,method='MBD', plot=FALSE,
             outliercol = 2, 
             fullout = FALSE,
             outline=FALSE,
             main = "functional boxplot, fbplot()",
             xlab = "Day",
             ylab = "Newly hospitalized" 
             ) 
## if plot=false, it outputs the calculations and scores

envelope_list <- c(envelope_list,
                   list(fda=get_envelope(ensemble_J, ff$Depth>median(ff$Depth))))

###################################
## Mahalanobis probes on Juul's data
###################################

n_probes <- 4
# probes.fun <- function(x) {
#   return(c(min = min(x), mean = median(x), max = max(x), curvemean=rowMeans(x)))
# }

probes <- matrix(NA,n_traj,n_probes) 
probes[,1] <- apply(t(ensemble_J), 1, FUN = min)
probes[,2] <- apply(t(ensemble_J), 1, FUN = max)   
probes[,3] <- apply(t(ensemble_J), 1, FUN = median)   
probes[,4] <- rowMeans(t(ensemble_J))


Sigma <- cov(probes)
mahalmat <- matrix(NA, n_traj, n_traj)
diag(mahalmat) <- 0

for (i in 1:(n_traj-1)) {
  for (j in (i+1):n_traj) {
    ## calculating pairwise M distance of both trajectories from their mean
    ## (identical, take the first one)
    mahalmat[i,j] <- mahalanobis(probes[c(i,j),], #taj i and j
                                  center=colMeans(probes[c(i,j),]),
                                  cov=Sigma)[1]
  }
}

mahalmat <- make_symm(mahalmat)
md_mah <- rank(rowMeans(mahalmat))
central_curves_mah <- which(md_mah<quantile(md_mah,0.5))

envelope_list <- c(envelope_list,
                   list(mahal=get_envelope(ensemble_J, central_curves_mah)))


envdat <- dplyr::bind_rows(envelope_list, .id="method")

long_ensemble <- ensemble_J %>% as.matrix() %>% reshape2::melt() %>%
    rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))

theme_set(theme_bw())
ggplot(envdat, aes(tvec)) + geom_ribbon(aes(ymin=lwr,ymax=upr),colour=NA, alpha=0.4, fill="blue") +
    facet_wrap(~method) + geom_line(data=long_ensemble, aes(y=value,group=grp), alpha=0.05)
