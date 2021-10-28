## Dr. Ben Bolker and Ali Gharouni
## Spring 2021, UNiversity of McMaster
## Comparison of Juul's work with alternative centrality scoring of an epidemic ensemble.
## We used L2 norm, FBP and Mahalanobis on probes.
# setwd("projects/commentsFBD/scripts/")
source('functions.R')
library(tidyverse)
library(latex2exp)
library(fda)
library(roahd)
library(directlabels)
library(colorspace)
library(ggrepel)
library(Cairo)
library(ggrastr)
library(tikzDevice)

do_mahal <- FALSE

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
## we want *RANKS* > 0.1 quantile to get 90% central region
central_curves <- which(md1 > quantile(md1,0.1))

envelope_list <- c(envelope_list,
                   list(juul=get_envelope(ensemble_J, central_curves)))
###################################
## Pointwise L2 Norm on Juul's data
###################################
## The distance between 2 curves is the area between them

## Memory allocation for distance matrix
dmatl2 <- matrix(NA, n_traj, n_traj) ## Euclidean L2 norm
dmatl2 <- as.matrix(dist(t(ensemble_J), method = "euclidean", upper = TRUE, p = 2))
stopifnot(isSymmetric(dmatl2))

## l2 mean corresponding to arithmetic mean: minimize sum of squares of distances
md2l2 <- rowSums(dmatl2^2)
m2l2 <- which.min(md2l2)
## corresponding to median: minimize sum of distances
md1l2 <- rowSums(dmatl2) ## these are depths
m1l2 <- which.min(md1l2)
## Correlated mean and median?
cor.test(md1l2,md2l2,method="spearman")$estimate

central_curves_l2 <- which(md1l2<quantile(md1l2,0.9))
## Find the envelope, i.e., the POINTWISE min/max of the central curves

envelope_list <- c(envelope_list,
                   list(L2norm=get_envelope(ensemble_J, central_curves_l2)))

###################################
## Functional Boxplot (FDA) on Juul's data
###################################

## This is similar to Juul's et al. method but instead the sample size is 2 at each draw
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
## get 90% boundaries
fda_new <- get_envelope(ensemble_J, ff$Depth>quantile(ff$Depth, 0.1))
write.csv(fda_new, file = "data/fda_new.csv", row.names=FALSE)

envelope_list <- c(envelope_list,
                   list(fda=fda_new))
## get_envelope(ensemble_J, ff$Depth>median(ff$Depth))))

library(fda)
fbplot_env <- list()
tvec <- 0:149
for (m in c("BD2", "MBD", "Both")) {
    ## fda::fbplot can't handle tibbles
    ff <- fda::fbplot(as.data.frame(ensemble_J), x = tvec, plot=FALSE, method=m)
    fbplot_env[[m]] <- get_envelope(ensemble_J, ff$depth>quantile(ff$depth, 0.1))
}
fbplot_env <- dplyr::bind_rows(fbplot_env, .id="method")
write.csv(fbplot_env, file = "data/fbplot.csv", row.names=FALSE)

###################################
## Mahalanobis probes on Juul's data
###################################

probes_mat <- t(apply(ensemble_J,2, probes))
pdf("probes_hist.pdf", width=8,height=5)
op <- par(mfrow=c(2,3))
for (i in 1:ncol(probes_mat)) {
    hist(probes_mat[,i],main=colnames(probes_mat)[i])
}
par(op)
dev.off()

Sigma <- cov(probes_mat)
mahalmat <- matrix(NA, n_traj, n_traj)
diag(mahalmat) <- 0

for (i in 1:(n_traj-1)) {
  for (j in (i+1):n_traj) {
    ## calculating pairwise M distance of both trajectories from their mean
    ## (identical, take the first one)
    mahalmat[i,j] <- mahalanobis(probes_mat[c(i,j),], #taj i and j
                                  center=colMeans(probes_mat[c(i,j),]),
                                  cov=Sigma)[1]
  }
}

mahalmat <- make_symm(mahalmat)
md_mah <- rank(rowMeans(mahalmat))
central_curves_mah <- which(md_mah < quantile(md_mah, 0.9)) 

envelope_list <- c(envelope_list,
                   list(mahal=get_envelope(ensemble_J, central_curves_mah)))

envdat <- dplyr::bind_rows(envelope_list, .id="method")

long_ensemble <- ensemble_J %>% as.matrix() %>% reshape2::melt() %>%
    rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))



theme_set(theme_bw())
## ggplot(envdat, aes(tvec)) + geom_ribbon(aes(ymin=lwr,ymax=upr),colour=NA, alpha=0.4, fill="blue") +
##    facet_wrap(~method) + geom_line(data=long_ensemble, aes(y=value,group=grp), alpha=0.05)

labvec <- c(juul = "FBD\n($J = 50$)",
            fda = "FBD\n($J = 2$)",
            mahal = "Mahalanobis\n(5 features)",
            L2norm="$\\ell_2$"
            )

envdat <- mutate(envdat, across(method, ~ labvec[.]))
## dat for plotting in the main body of the manuscript;  without L_2 norm
envdat_sub <- envdat %>% filter(method != "$\\ell_2$") 

## MAIN PLOT #######################################################
## manual label positioning
maxdat <- (envdat_sub
    %>% group_by(method)
    %>% filter(upr == max(upr))
    %>% ungroup()
    ## juul (J=50) / roahd(fda) (J=2)/ mahal
    %>% mutate(xoff = c(+40, +10, -30),
               yoff = c(50, 70, 50))
)

cent_plot <- (ggplot(envdat_sub, aes(tvec))
  + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=method,colour=method),alpha=0.4,lwd=1)
  + rasterise(geom_line(data=long_ensemble, aes(y=value,group=grp), alpha=0.1),
              dpi = 300)
  ## + geom_dl(aes(label=method,y=upr,colour=method),
  ## method=list(cex=2,"top.points"))
  ## + geom_label_repel(data = maxdat, aes(label=method, y = upr, colour=method),
  ## force=50, force_pull = 0.1)
  + geom_segment(data = maxdat, aes(x = tvec + xoff, xend = tvec,
                                    y = upr + yoff, yend = upr,
                                    colour = method),
                 width = 2)
  ## labels after segments so segments are masked appropriately
  + geom_label(data = maxdat, aes(label = method,
                                  x = tvec + xoff,
                                  y = upr + yoff,
                                  colour = method),
               fill = "white",
               label.padding = unit(0.3, "lines"))
  + scale_colour_discrete_qualitative()
  + scale_fill_discrete_qualitative()
  + labs(x="time", y="number newly hospitalized")
  + theme(legend.position = "none")
)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}","\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))
 ggsave(cent_plot,
       device = tikz,
       filename = "cent_plot.tex",
       standAlone = TRUE,
       width = 6, height = 6, units = "in")

 # ggsave(cent_plot,
 #        device = "png",
 #        filename = "cent_plot.png",
 #        width = 6, height = 6, units = "in")

 ## APPENDIX PLOT #######################################################
 
# maxdat2 <- (envdat
#              %>% group_by(method)
#              %>% filter(upr == max(upr))
#              %>% ungroup()
#              ## juul (J=50) / roahd(fda) (J=2)/ mahal / L_2
#              %>% mutate(xoff = c(+40, +10, -30, -40),
#                         yoff = c(50, 70, 50, 90))
#  )
 
 maxdat2 <- (envdat
             %>% filter(method== "$\\ell_2$")
             %>% group_by(method)
             %>% filter(upr == max(upr))
             %>% ungroup()
             ## juul (J=50) / roahd(fda) (J=2)/ mahal / L_2
             %>% mutate(xoff = c(30),
                        yoff = c(90))
 )
 
 cent_plot2 <- ((cent_plot %+% envdat)
                + geom_segment(data = maxdat2, aes(x = tvec + xoff, xend = tvec,
                                                   y = upr + yoff, yend = upr,
                                                   colour = method),
                               width = 2)
                ## labels after segments so segments are masked appropriately
                + geom_label(data = maxdat2, aes(label = method,
                                                 x = tvec + xoff,
                                                 y = upr + yoff,
                                                 colour = method),
                             fill = "white",
                             label.padding = unit(0.3, "lines"))
 )
 
 options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}","\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))
 ggsave(cent_plot2,
        device = tikz,
        filename = "cent_plot2.tex",
        standAlone = TRUE,
        width = 6, height = 6, units = "in")
 
 
 
 
 