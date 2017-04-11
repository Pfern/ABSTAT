
# Challenge 06 Gibbs Sampler 
# --------------------------


# 1 Multinomial Distribution
# --------------------------
 ?rmultinom

# c - rmultinom(n,size,prob)
theta<-c(0.2,0.3,0.5)
data<-t(rmultinom(1,1000,theta))
data

# d - dmultinom(x,size,prob)
dmultinom(c(220,350,430),1000,theta)


# 2 Dirichlet Distribution
# --------------------------
# install package gtools
library(gtools)
?rdirichlet

# c - rdirichlet(n,alpha)
a.0<-c(1,1,1)
theta.0<-rdirichlet(1,a.0)
theta.0

# d - ddirichlet(x,alpha)
ddirichlet(c(0.15,0.25,0.6),a.0) # Note: prob density function > 0


# 3 Simulate vector theta
# ------------------------

# rdirichlet(x,alpha_updated)

# a<-a.0 # Initialize in order to check the function
rdirichlet(1,a+data) # Note: prob density function > 0


# 4 Gibbs Sampler
# ----------------

# a/b - Gibbs cicle

# using rdirichlet

a<-a.0
theta.all.iter<-matrix(0,5000,3)
for(i in 1:5000)
{
  theta.all.iter[i,]<-rdirichlet(1,a) 
  a<-a+data
  # cat(i,"\n")
}
tail(theta.all.iter)


# 5 Trace 7 Parameters estimation
# --------------------------------

# a - Trace

par(mfrow=c(3,1))

plot(1:5000,theta.all.iter[,1],main=expression(paste("Trace for ",theta[1])),ylab=expression(theta),xlab="iteration")
plot(1:5000,theta.all.iter[,2],main=expression(paste("Trace for ",theta[2])),ylab=expression(theta),xlab="iteration")
plot(1:5000,theta.all.iter[,3],main=expression(paste("Trace for ",theta[3])),ylab=expression(theta),xlab="iteration")

# b - Burn-in
# Consider a burn-in period of 1000 (dependes on the runs)

# c - Estimates
theta.final<-apply(theta.all.iter[1001:5000,],2,mean) # 2: 'by column'
round(theta.final,2)
