


# A single box plot gives a graphical representation of the range of variation in a numerical variable, based on five numbers:
#   
#   The minimum and maximum values
# The median (or "middle") value
# Two intermediate values called the lower and upper quartiles
# In addition, the standard box plot computes a nominal data range from three of these numbers and flags points falling outside this range as outliers, representing them as distinct points.
# 
# The bag plot extends this representation to two numerical variables, showing their relationship, both within two-dimensional "bags" corresponding to the "box" in the standard boxplot, and indicating outlying points outside these limits.

library(aplpack)

## examples from RDocumentation https://www.rdocumentation.org/packages/aplpack/versions/1.3.3/topics/bagplot
# example: 100 random points and one outlier
dat<-cbind(rnorm(100)+100,rnorm(100)+300)
dat<-rbind(dat,c(105,295))
bagplot(dat,factor=2.5,create.plot=TRUE,approx.limit=300,
        show.outlier=TRUE,show.looppoints=TRUE,
        show.bagpoints=TRUE,dkmethod=2,
        show.whiskers=TRUE,show.loophull=TRUE,
        show.baghull=TRUE,verbose=FALSE)
# example of Rousseeuw et al., see R-package rpart
cardata <- structure(as.integer( c(2560,2345,1845,2260,2440,
                                   2285, 2275, 2350, 2295, 1900, 2390, 2075, 2330, 3320, 2885,
                                   3310, 2695, 2170, 2710, 2775, 2840, 2485, 2670, 2640, 2655,
                                   3065, 2750, 2920, 2780, 2745, 3110, 2920, 2645, 2575, 2935,
                                   2920, 2985, 3265, 2880, 2975, 3450, 3145, 3190, 3610, 2885,
                                   3480, 3200, 2765, 3220, 3480, 3325, 3855, 3850, 3195, 3735,
                                   3665, 3735, 3415, 3185, 3690, 97, 114, 81, 91, 113, 97, 97,
                                   98, 109, 73, 97, 89, 109, 305, 153, 302, 133, 97, 125, 146,
                                   107, 109, 121, 151, 133, 181, 141, 132, 133, 122, 181, 146,
                                   151, 116, 135, 122, 141, 163, 151, 153, 202, 180, 182, 232,
                                   143, 180, 180, 151, 189, 180, 231, 305, 302, 151, 202, 182,
                                   181, 143, 146, 146)), .Dim = as.integer(c(60, 2)), 
                     .Dimnames = list(NULL, c("Weight", "Disp.")))
bagplot(cardata,factor=3,show.baghull=TRUE,
        show.loophull=TRUE,precision=1,dkmethod=2)
title("car data Chambers/Hastie 1992")
# points of y=x*x
bagplot(x=1:30,y=(1:30)^2,verbose=FALSE,dkmethod=2)
# one dimensional subspace
bagplot(x=1:100,y=1:100)

