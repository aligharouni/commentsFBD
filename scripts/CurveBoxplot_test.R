
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
EnRkTemp_J <- EnsRank_all(ensemble_J,numSample=100,sizeSample=50) ##  as Juul's et al did.
##saveRDS(EnRkTemp_J, file = "Juul_replicate.rds")

rnkmat <- EnRkTemp_J[["ensembRank"]] ## take the resulted rank mat

md1 <- rank(colSums(rnkmat)) ## , ties="random")
central_curves <- which(md1>quantile(md1,0.5))

envelope_list <- c(envelope_list,
                   list(juul=get_envelope(ensemble_J, central_curves)))
###################################
## ckeck: 1- Juul's algorithm with J=2 and fda match?
set.seed(110)
EnRkTemp_J2 <- EnsRank_all2(ensemble_J,numSample=124750,sizeSample=2) 
## note numSample in fbplot is choose(500,2)
rnkmat2 <- EnRkTemp_J2[["ensembRank"]]
md12 <- rank(colSums(rnkmat2))
central_curves2 <- which(md12>quantile(md12,0.5))
envelope_list <- c(envelope_list,
                   list(juul_J2=get_envelope(ensemble_J, central_curves2)))

###################################
## Functional Boxplot (FDA) on Juul's data
###################################
## This is similar to Juul's et al. method but instead the sample size is 2 at each draw
## Data depth is a well-known and useful nonparametric tool for analyzing functional data. It provides a novel way of ranking a sample of curves from the center outwards and defining robust statistics, such as the median or trimmed means.
fD <- fData(tvec,as.matrix(t(ensemble_J))) ## note that the trajs are on rows in fD object
# dev.new()
ff <- roahd::fbplot(fD,method='MBD', plot=FALSE, #method='MBD'
             outliercol = 2,
             fullout = FALSE,
             outline=FALSE,
             main = "functional boxplot, fbplot()",
             xlab = "Day",
             ylab = "Newly hospitalized"
)
# ff <- roahd::fbplot(fD,method='BD2', plot=FALSE, #method='MBD'
#                     main = "functional boxplot, fbplot()",
#                     xlab = "Day",
#                     ylab = "Newly hospitalized"
# )

envelope_list <- c(envelope_list,
                   list(fda_roahd=get_envelope(ensemble_J, ff$Depth>median(ff$Depth))))
###################################
## Settings/plotting
###################################
envdat <- dplyr::bind_rows(envelope_list, .id="method")

long_ensemble <- ensemble_J %>% as.matrix() %>% reshape2::melt() %>%
  rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))

theme_set(theme_bw())

labvec <- c(juul = "AG\n($J = 50$)",
            fda_roahd = "FBD\n($J = 2$)",
            juul_J2 = "AG\n($J = 2$)"
)

envdat <- mutate(envdat, across(method, ~ labvec[.]))

envdat_temp <- envdat 
# envdat_temp <- envdat %>% filter(method==c("FBD\n($J = 2$)","Juul\n($J = 2$)"))
# envdat_temp <- envdat %>% filter(method==c("FBD\n($J = 50$)","Juul\n($J = 2$)"))
cent_plot2 <- (ggplot(envdat_temp, aes(tvec))
               + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=method,colour=method),alpha=0.4,lwd=1)
)
print(cent_plot2)

###################################
# 3- do roahd and fda give similar curves?
###################################
envelope_list <- list()

## fda
ff_fda <-fda::fbplot(as.matrix(ensemble_J),method='MBD',ylim=c(0,1000))

envelope_list <- c(envelope_list,
                   list(fda_fda=get_envelope(ensemble_J, ff_fda$depth>median(ff_fda$depth))))
## roahd
fD <- fData(tvec,as.matrix(t(ensemble_J))) 
ff_roahd <-roahd::fbplot(fD,method='MBD', plot=FALSE, 
                         outliercol = 2,
                         fullout = FALSE,
                         outline=FALSE
)
## The following gives a uniroot() error! not sure at this point.
# ff_roahd <-roahd::fbplot(fD,method='MBD', plot=FALSE, 
#                          adjust = list( N_trials = 10,
#                                         trial_size = 2 * fD$N,
#                                         VERBOSE = TRUE ),
#                          # outliercol = 2,
#                          fullout = FALSE,
#                          outline=FALSE
# )

envelope_list <- c(envelope_list,
                   list(fda_roahd=get_envelope(ensemble_J, ff_roahd$Depth>median(ff_roahd$Depth))))

envdat <- dplyr::bind_rows(envelope_list, .id="method")
theme_set(theme_bw())
labvec <- c(fda_fda = "FBD\n(fda)",
            fda_roahd = "FBD\n(roahd)")

envdat <- mutate(envdat, across(method, ~ labvec[.]))

envdat_temp <- envdat 
cent_plot2 <- (ggplot(envdat_temp, aes(tvec))
               + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=method,colour=method),alpha=0.4,lwd=1)
)
print(cent_plot2)
## Conclusion: the 2 packages fda and roah gives the same results.

#################################################
## Extra Exploration: more test (#1) on a simple and known example
##i.e., explore the differences between fda() and Juul's algorithm (fir J=2) 
#################################################
envelope_list <- list()
## Creating toy example with trajectories in columns  
tvec <- seq(1,10,length=10)
y1 <- ifelse(tvec<5,0,1)
y2 <- ifelse(tvec<6,0,1)
y3 <- ifelse(tvec<7,0,1)
y4 <- ifelse(tvec<8,0,1)


ens_toy <- cbind(y1,y2,y3,y4)
matplot(tvec,ens_toy,type="l")

## fda:
fD <- fData(tvec,as.matrix(t(ens_toy))) 

ff_roahd <-roahd::fbplot(fD,method='MBD', plot=FALSE, 
                         outliercol = 2,
                         fullout = FALSE,
                         outline=FALSE
)

envelope_list <- c(envelope_list,
                   list(fda_roahd=get_envelope(ens_toy, ff_roahd$Depth>median(ff_roahd$Depth))))

# envelope_list$fda_roahd

## Juul's algorithm with J=2 on the ens_toy:
set.seed(1234)
# EnRkTemp_J2 <- EnsRank_all(ens_toy,numSample=6,sizeSample=2)
EnRkTemp_J2 <- EnsRank_all2(ens_toy,numSample=6,sizeSample=2)
rnkmat2 <- EnRkTemp_J2[["ensembRank"]]
md12 <- rank(colSums(rnkmat2))
central_curves2 <- which(md12>quantile(md12,0.5))
envelope_list <- c(envelope_list,
                   list(juul_J2=get_envelope(ens_toy, central_curves2)))

envdat <- dplyr::bind_rows(envelope_list, .id="method")
long_ensemble <- ens_toy %>% as.matrix() %>% reshape2::melt() %>%
  rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))

ggplot(envdat, aes(tvec)) + 
  geom_ribbon(aes(ymin=lwr,ymax=upr), colour=NA,alpha=0.4, fill="blue") +
  facet_wrap(~method) + 
  geom_line(data=long_ensemble, aes(y=value,group=grp), alpha=0.05)
 



