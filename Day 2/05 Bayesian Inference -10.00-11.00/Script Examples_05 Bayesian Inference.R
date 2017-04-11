
# Examples 05 Bayesian Inference
# ------------------------------


# Occurence of Nucleotide A in a sequence
# ---------------------------------------

# X: nr. nucleotides A in a certain position of the sequence
# X ~ Bernoulli(theta)


# Prior distribution: theta ~ Beta(a,b) 
# initial hyperparameters: a=1 and b=1
a<-1; b<-1

# Sample 1  (t=sum_xi)
n<-100; x<-43

# Posterior distribution: theta|x ~ Beta(a*,b*)
# Updated hyperparametrs: a*=t+a and b*=n-t+b 
a.star<-x+a; a.star   
b.star<-n-x+b; b.star

s<-seq(0,1,by=0.01)
y<-dbeta(s,a,b)
y.star<-dbeta(s,a.star,b.star)
plot(s,y,type="n",xlim=c(0,1),ylim=c(0,2.8),xlab=expression(theta),ylab="density function: beta(5,7)")
lines(s,y)
lines(s,y.star,col="red")
text(0.8,2,"prior")
text(0.8,1.8,"posterior",col="red")

# theta_mean = a*/(a*+b*)  (posterior mean)
a.star/(a.star+b.star)

# P(theta<=theta_median|x)=0.5 --> theta_median=quantil_0.5  (posterior median)
qbeta(0.5,a.star,b.star)

# theta_mode = (a*-1)/(a*+b*-2)  (posterior mode)
(a.star-1)/(a.star+b.star-2)

# credible interval for theta with alpha=0.05
# lower bound = quantil_0.025
# upper bound = quantil_0.975
qbeta(c(0.025,0.975),a.star,b.star)

