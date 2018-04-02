
# Example 06 Gibbs Sampling
# -------------------------


# Simulate Data y=(y_1,...,y_1000)
# --------------------------------

# Using simulates data we can check wether the algorithm 
# is estimating the parametrs correctly.

# 1000 observarions from Y~N(3,9)

n<-1000
mu<-3
sigma<-3  # sigma^2=9
y<-rnorm(n,mu,sigma);head(y)
summary(y)
var(y)
hist(y,col=rainbow(10))


# Store Results
# -------------
ns<-5000              # nr. of simulations
Mu<-Sigma2<-rep(0,ns) # store parameters' updates


# Gibbs  Sampler
# --------------

m<-mean(y); m
# initialize tau: tau^(0)
tau<-rgamma(1,n/2,sum((y-m)^2)/2); tau 
sigma2<-1/tau; sigma2

for (i in 1:ns) {
  mu<- Mu[i]<-rnorm(1,m,1/(n*tau))     # mu^(1), mu^(2), ...
  tau<-rgamma(1,n/2,sum((y-mu)^2)/2)   # tau^(1), tau^(2), ...
  Sigma2[i]<-1/tau
}


# Posterior summaries 
# -------------------


plot(1:ns,Mu[1:ns],type="l",col="lightblue",xlab="iteration",ylab="mu")
mean(Mu[1:ns])	
sd(Mu[1:ns])

plot(1:ns,Sigma2[1:ns],type="l",col="lightgreen",xlab="iteration",ylab="sigma^2")
mean(Sigma2[1:ns])
sd(Sigma2[1:ns])

# Discard first 1500 as burn-in
# The estimates do not differ from the previous ones

plot(1:ns,Mu[1:ns],type="l",col="lightblue")
abline(v=1500)
mean(Mu[1501:ns])	
sd(Mu[1501:ns])

plot(1:ns,Sigma2[1:ns],type="l",col="lightgreen")
abline(v=1500)
mean(Sigma2[1501:ns])
sd(Sigma2[1501:ns])



