
# This is meant to parallelize sampling compare and to verify whether 
# 1- Juul's algorithm with J=2 and fda match?

source('functions.R')
library(tidyverse)
library(latex2exp)
library(fda)
library(roahd)
library(directlabels)
library(colorspace)

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
## ckeck: 1- Juul's algorithm with J=2 and fda match?
set.seed(1234)
## BMB: parallelize??

library(doParallel) # setting up parallel
library(foreach)
cl <- makeCluster(detectCores())
registerDoParallel(cl)




nsamp <- c(10000,100000)
i <- 1
EnRkTemp_J2 <- EnsRank_all2(ensemble_J, numSample=nsamp[i], sizeSample=2)
## note numSample in fbplot is choose(500,2)
rnkmat2 <- EnRkTemp_J2[["ensembRank"]]
md12 <- rank(colSums(rnkmat2))
central_curves2 <- which(md12>quantile(md12,0.5))
envelope_list <- c(envelope_list,
                     list(AG_J2_1=cbind(get_envelope(ensemble_J, central_curves2),
                                     nsample=rep.int(nsamp[i],length(tvec)))
                          )) 

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

envelope_list <- c(envelope_list,
                   list(fda_roahd=get_envelope(ensemble_J, ff$Depth>median(ff$Depth))))

###################################
## From Juul's Python curveboxplot() function
###################################
# J_vals <- c(2, 5, 10, 20, 30, 40, 50)
bnds_juul <- data_frame(read_csv("./data/juul_boundary.csv", col_names = TRUE)) 

## note adding 0 to the end of the lwr and upr, since mismatch with length of tvec

envelope_list <- c(envelope_list,
                   list(juul_J2=data.frame(tvec=tvec, 
                                           lwr=c(as.numeric(bnds_juul %>% pull("2_b_lwr50")),0), 
                                           upr=c(as.numeric(bnds_juul %>% pull("2_b_upr50")),0) ))
                   )

###################################
## Settings/plotting
###################################
envdat <- dplyr::bind_rows(envelope_list, .id="method")
saveRDS(envdat, file = "envdat.rds")

