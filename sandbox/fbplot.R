## fda::fbplot and roahd::fbplot gives the same depth? Yes

## url: https://rdrr.io/cran/fda/man/fbplot.html
source('../scripts/functions.R')
library(fda)
library(roahd)

## Placeholder to save the method specific envelopes including Juul's work, L2, FBP, Mahalonobis on probes.
envelope_list <- list()


## 5.  model 2 with outliers
##
#magnitude
cov.fun=function(d,k,c,mu){
  k*exp(-c*d^mu)
}
k=6
n=50
p=30
t=seq(0,1,len=p)
d=dist(t,upper=TRUE,diag=TRUE)
d.matrix=as.matrix(d)
#covariance function in time
t.cov=cov.fun(d.matrix,1,1,1)
# Cholesky Decomposition
L=chol(t.cov)
#randomly introduce outliers
set.seed(1234)
C=rbinom(n,1,0.1)
s=2*rbinom(n,1,0.5)-1
cs.m=matrix(C*s,p,n,byrow=TRUE)
mu=4*t
e<-matrix(rnorm(n*p),p,n)
ensemble_tmp<-mu+t(L)%*%e

## functional boxplot
ff_fda <-fda::fbplot(ensemble_tmp,method='MBD',ylim=c(-2,6))

envelope_list <- c(envelope_list,
                   list(fda_fda=get_envelope(ensemble_tmp, ff_fda$depth>median(ff_fda$depth))))


## functional boxplot
tvec<-seq(1,p)
ydata_fdat <- fData(tvec,as.matrix(t(ensemble_tmp))) ## note that the trajs are on rows in fD object
ff_roahd <-roahd::fbplot(ydata_fdat,method='MBD')

envelope_list <- c(envelope_list,
                   list(fda_roahd=get_envelope(ensemble_tmp, ff_roahd$Depth>median(ff_roahd$Depth))))


envdat <- dplyr::bind_rows(envelope_list, .id="method")

long_ensemble <- ensemble_tmp %>% as.matrix() %>% reshape2::melt() %>%
  rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))

theme_set(theme_bw())

ggplot(envdat, aes(tvec)) + geom_ribbon(aes(ymin=lwr,ymax=upr),colour=NA, alpha=0.4, fill="blue") +
   facet_wrap(~method) + geom_line(data=long_ensemble, aes(y=value,group=grp), alpha=0.05)

labvec <- c(fda_fda = "FBD\n(fda)",
            fda_roahd = "FBD\n(roahd)")

envdat <- mutate(envdat, across(method, ~ labvec[.]))

envdat_temp <- envdat 
cent_plot2 <- (ggplot(envdat_temp, aes(tvec))
               + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=method,colour=method),alpha=0.4,lwd=1)
)
print(cent_plot2)


