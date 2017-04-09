
# Example 06 Gibbs Sampling
# -------------------------


# Simulate Data
# -------------

# 1000 observarions from Y~N(3,9)

n<-1000
mu<-3
sigma<-3
sigma2<-sigma^2
y<-rnorm(n,mu,sigma);head(y)
summary(y)
hist(y,col=rainbow(10))


# Store Results
# -------------
ns<-5000              # nr. of simulations
Mu<-Sigma2<-rep(0,ns) # store parameters' updates


# Gibbs  Sampler
# --------------

tau<-1     # initialize tau to calculate mu^(0)
m<-mean(y) # in order to calculate mu^(t)

for (i in 1:ns) {
  mu<- Mu[i]<-rnorm(1,m,1/(n*tau))
  tau<-rgamma(1,n/2,sum((y-mu)^2)/2)
  Sigma2[i]<-1/tau
}


# Posterior summaries 
# -------------------


plot(1:ns,Mu[1:ns],type="l",col="lightblue",xlab="iteration",ylab="mu")
plot(1:ns,Sigma2[1:ns],type="l",col="lightgreen",xlab="iteration",ylab="sigma^2")

mean(Mu[1:ns])	
sd(Mu[1:ns])
mean(Sigma2[1:ns])
sd(Sigma2[1:ns])


plot(1:ns,Mu[1:ns],type="l",col="lightblue")
abline(v=1500)
plot(1:ns,Sigma2[1:ns],type="l",col="lightgreen")
abline(v=1500)

mean(Mu[1501:ns])	# Discard first 1500 as burn-in
sd(Mu[1501:ns])
mean(Sigma2[1501:ns])
sd(Sigma2[1501:ns])



