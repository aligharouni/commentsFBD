
# This is meant to compare and to verify whether 
# 1- Juul's algorithm with J=2 and fda match?
# 2- FBD + J=50 in the comparison
# 3- do roahd and fda give similar curves? YES
# rm(list=ls())
source('functions.R')
library(tidyverse)
library(latex2exp)
library(fda)
library(roahd)
library(directlabels)
library(colorspace)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)

## Placeholder to save the method specific envelopes including Juul's work, L2, FBP, Mahalonobis on probes.
envelope_list <- list()

###################################
## Juul's work replicating panel b;
###################################
## upload Fig.2 ensemble
ensemble_J <- read_csv("./data/juul1.csv", col_names = FALSE) ## the columns are the trajs
tvec <- seq(nrow(ensemble_J))-1 ## time domain
n_traj <- ncol(ensemble_J)
###################################
## Juul et al., Fig. 2 panel b;
## the centrality of curves was ranked in different ways: all-or-nothing ranking for the full predicted time interval, Ncurves = 50 and Nsamples = 100
###################################
set.seed(1234)
nSample = c(100,1000)
sSample =c(2,50)
## NON-PARALLELIZED APPROACH
# for (n in nSample){
#   for (j in sSample){
#     print(c(n,j))
#     ##  This is AG algorithm to simulate Juul's et al algorithm.
#     EnRkTemp_J <- EnsRank_all(ensemble_J, numSample=n, sizeSample=j) 
#     ##saveRDS(EnRkTemp_J, file = "Juul_replicate.rds")
#     rnkmat <- EnRkTemp_J[["ensembRank"]] ## take the resulted rank mat
#     md1 <- rank(colSums(rnkmat)) ## , ties="random")
#     central_curves <- which(md1>quantile(md1,0.5))
#     envelope_list[[paste("AG_J",j,"_",n, sep = "")]] <- 
#       cbind(get_envelope(ensemble_J, central_curves), 
#             nsample=rep.int(n,length(tvec))
#       )
#   } }

## setup parallelized compute double R loop
## see https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()

envelope_list_temp <- list()
envelope_list_temp <- 
  foreach(j=sSample, .inorder = FALSE) %:% 
  foreach(n=nSample,.inorder = FALSE) %dopar% {
    ##  This is AG algorithm to simulate Juul's et al algorithm.
    EnRkTemp_J <- EnsRank_all(ensemble_J, numSample=n, sizeSample=j)
    rnkmat <- EnRkTemp_J[["ensembRank"]] ## take the resulted rank mat
    md1 <- rank(colSums(rnkmat)) ## , ties="random")
    central_curves <- which(md1>quantile(md1,0.5))
    ## return a list
    c(data = get_envelope(ensemble_J, central_curves), 
      nsample = n,
      J = j,
      name = paste("AG_J",j,"_",n, sep = "") ## name the list here
    )
  }

## Note sure if there is a better way to do this: 
## Fill envelope_list from envelope_list_temp 
for (j1 in 1:length(sSample)){
  for (n1 in 1:length(nSample)){
    envelope_list[[envelope_list_temp[[j1]][[n1]][["name"]] ]] <- 
      envelope_list_temp[[j1]][[n1]][-6] ## remove the name
   }}
# rm(envelope_list_temp)

parallel::stopCluster(cl = my.cluster)                   

## Delete?
###################################
## ckeck: 1- Juul's algorithm with J=2 and fda match?
# set.seed(1234)
# ## BMB: parallelize??
# nsamp <- c(10000,100000)
# i <- 1
# EnRkTemp_J2 <- EnsRank_all2(ensemble_J, numSample=nsamp[i], sizeSample=2)
# ## note numSample in fbplot is choose(500,2)
# rnkmat2 <- EnRkTemp_J2[["ensembRank"]]
# md12 <- rank(colSums(rnkmat2))
# central_curves2 <- which(md12>quantile(md12,0.5))
# envelope_list <- c(envelope_list,
#                      list(AG_J2_1=cbind(get_envelope(ensemble_J, central_curves2),
#                                      nsample=rep.int(nsamp[i],length(tvec)))
#                           )) 

###################################
## Functional Boxplot (FDA) on Juul's data
###################################
## This is similar to Juul's et al. method but instead the sample size is 2 at each draw
## Data depth is a well-known and useful nonparametric tool for analyzing functional data. It provides a novel way of ranking a sample of curves from the center outwards and defining robust statistics, such as the median or trimmed means.
fD <- fData(tvec,as.matrix(t(ensemble_J))) ## note that the trajs are on rows in fD object
# dev.new()
ff <- roahd::fbplot(fD,method='MBD', plot=FALSE, 
             outliercol = 2,
             fullout = FALSE,
             outline=FALSE,
             main = "functional boxplot, fbplot()",
             xlab = "Day",
             ylab = "Newly hospitalized"
)

envelope_list[["fda_roahd"]] <- c( data = get_envelope(ensemble_J, ff$Depth>median(ff$Depth)),
                                   nsample = choose(500,2),
                                   J = 2 )

###################################
## From Juul's Python curveboxplot() function
###################################
J_vals <- c(2, 50)
ss <- c(100,10000) ## sample size
# bnds_juul <- data_frame(read_csv( paste("./data/juul_boundary_ss",ss[1],".csv", sep = ""), col_names = TRUE))

populate_envlist <- function(J, sampsize){
  # read the csv file created in Juul_work.Rmd in python
  bnds_juul <- data_frame(read_csv(paste("./data/juul_boundary_ss",sampsize,".csv", sep = ""), col_names = TRUE))
  name <- paste("Juul_J",J,"_",sampsize, sep = "")
  ## note adding 0 to the end of the lwr and upr, since mismatch with length of tvec
  b_lwr <- paste(J,"b_lwr50",sep = "_")
  b_upr <- paste(J,"b_upr50",sep = "_")
  out <- c( data = data.frame(tvec=tvec, 
            lwr=c(as.numeric(bnds_juul %>% pull(b_lwr)),0), 
            upr=c(as.numeric(bnds_juul %>% pull(b_upr)),0)),
            nsample = sampsize,
            J = J)
  return((out))
} 

# Maybe FIXME: smarter way? lapply or something else??
for(sampsize in ss){
  for(J in J_vals){
    print(c(sampsize,J))
    name <- paste("Juul_J",J,"_",sampsize, sep="")
    print(name)
    envelope_list[[name]] <- populate_envlist(J,sampsize)
  }
}


###################################
## Settings/plotting
###################################
envdat <- dplyr::bind_rows(envelope_list, .id="method")
saveRDS(envdat, file = "envdat.rds")

