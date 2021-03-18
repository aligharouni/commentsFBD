
## consider 2 vectors, or 2 functions with the same domain
traj1 <- c(0,0,1,0,0)
traj2 <- c(0,0,0,1,0) ##shifted traj1
# traj2 <- c(0,0,0,0,0)
tvec <- 1:length(traj1)
m <-cbind(traj1,traj2)
matplot(x=tvec,y=m,type = 'l')

dist(t(m), method = "euclidean", diag = TRUE, upper = TRUE, p = 2) ## by default gives L2 
dmat <- as.matrix(dist(t(m), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))



