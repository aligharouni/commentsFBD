# This is meant to compare and to verify whether 
# 1- Juul's algorithm with J=2 and fda match?
# 2- FBD + J=50 in the comparison
# 3- do roahd and fda give similar curves?
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
plot(sort(md1)) ## Note the discontinuity!?
## Discontinuity is due to ties; if ties were broken with 'random' this
## would go away
central_curves <- which(md1>quantile(md1,0.5))

envelope_list <- c(envelope_list,
                   list(juul=get_envelope(ensemble_J, central_curves)))
###################################
## ckeck: Does juul algorithm converges to fda for ncurves=2?
#set.seed(1234)
EnRkTemp_J2 <- EnsRank_all(ensemble_J,numSample=1000,sizeSample=2)
rnkmat2 <- EnRkTemp_J[["ensembRank"]]
md12 <- rank(colSums(rnkmat2))
central_curves2 <- which(md12>quantile(md1,0.5))
envelope_list <- c(envelope_list,
                   list(juul_J2=get_envelope(ensemble_J, central_curves2)))

###################################
## Functional Boxplot (FDA) on Juul's data
###################################

## This is similar to Juul's et al. method but instead the sample size is 2 at each draw
## Data depth is a well-known and useful nonparametric tool for analyzing functional data. It provides a novel way of ranking a sample of curves from the center outwards and defining robust statistics, such as the median or trimmed means.

fD <- fData(tvec,as.matrix(t(ensemble_J))) ## note that the trajs are on rows in fD object
# dev.new()
ff <- fbplot(fD,method='BD', plot=FALSE, #method='MBD'
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
envdat <- dplyr::bind_rows(envelope_list, .id="method")

long_ensemble <- ensemble_J %>% as.matrix() %>% reshape2::melt() %>%
  rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))

theme_set(theme_bw())
## ggplot(envdat, aes(tvec)) + geom_ribbon(aes(ymin=lwr,ymax=upr),colour=NA, alpha=0.4, fill="blue") +
##    facet_wrap(~method) + geom_line(data=long_ensemble, aes(y=value,group=grp), alpha=0.05)

labvec <- c(juul = "FBD\n($J = 50$)",
            fda = "FBD\n($J = 2$)",
            mahal = "Mahalanobis\n(6 features)",
            juul_J2 = "Juul\n($J = 2$)"
)

envdat_temp <- envdat %>% filter(method==c("FBD\n($J = 2$)","Juul\n($J = 2$)"))
cent_plot2 <- (ggplot(envdat_temp, aes(tvec))
               + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=method,colour=method),alpha=0.4,lwd=1)
)
print(cent_plot2)



