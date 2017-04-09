#=================================#
#=====   MC and Bootstrap     ====#
#=================================#



#===========================================================
# Exercise 1  MC simulation
#===========================================================


set.seed(3) #  is used to set the random number seed, 
#  if we want the results reproducible we use 
#  set.seed before generate the number

S <- 1000   #  number of MC simulations
n <- 15     #  sample dimension
mu <- 1     # true mean value
sigma <- sqrt(5/3) # standard deviation

out<-t(replicate(S,rnorm(n,mu,sigma))) # 1000 MC replicates
str(out)
class(out)
trimmean <- function(Y){mean(Y,0.2)} # function to calculate 20% trimmean

outsampmean <- apply(out,1,mean)   
outtrimmean <- apply(out,1,trimmean)
outmedian <- apply(out,1,median)
summary.sim <- data.frame(mean=outsampmean,trim=outtrimmean,
                          median=outmedian)

meanMC<-function(Y){   # MC mean estimate 
  sum(Y)/length(Y)
}

MCbias<-function(Y,T){ # MC bias  T=true value of a parameter
  meanMC(Y)-T
}



#labels:
#mean=1,  20% trimmean=2,  median=3

MC_mean.1<-meanMC(outsampmean)
MC_mean.2<-meanMC(outtrimmean)
MC_mean.3<-meanMC(outmedian)

MCsd.1<-sd(outsampmean)
MCsd.2<-sd(outtrimmean)
MCsd.3<-sd(outmedian)


MCbias.1<-MCbias(outsampmean,1)
MCbias.2<-MCbias(outtrimmean,1)
MCbias.3<-MCbias(outmedian,1)

MC_MSE.1<-(MCsd.1)^2-(MCbias.1)^2
MC_MSE.2<-(MCsd.2)^2-(MCbias.2)^2
MC_MSE.3<-(MCsd.3)^2-(MCbias.3)^2

results<-matrix(c(1,1,1,1000,1000,1000,MC_mean.1,MC_mean.2,MC_mean.3,
                  MCsd.1,MCsd.2,MCsd.3,MCbias.1,MCbias.2,MCbias.3,MC_MSE.1,MC_MSE.2,MC_MSE.3),
                ncol=3,byrow=T,dimnames=list(c("True value","# replicates","MC mean","MC sd","MC bias","MC MSE"),
                                             c("Sample mean","Trimmed mean","Median")))

print(round(results,3))

RE12<-MCsd.1^2/MCsd.2^2
RE13<-MCsd.1^2/MCsd.3^2
RE23<-MCsd.2^2/MCsd.2^3

print(RE12)
print(RE13)
print(RE23)


#===========================================================
# Exercise 2 Parametric bootstrap CI
#===========================================================


#95% CI for exponential parameter lambda
# Given 300 data points with mean 2. 
# Assume the data is exp(lambda) 

# We are given the number of data points and mean 
n<-300  
xbar<-2  
# The MLE for lambda is 1/xbar  
lambda_hat<-1/xbar 
# Generate the bootstrap samples 
# Each column is one bootstrap sample (of 300 resampled values) 
nboot = 1000 
# Here's the key difference with the empirical 
# bootstrap. We draw the bootstrap sample from Exponential(hat lambda) 

x<-rexp(n*nboot,lambda_hat) 
bootstrap_sample<-matrix(x, nrow=n, ncol=nboot)  

# Compute the bootstrap lambda star 
lambda_star = 1./colMeans(bootstrap_sample) 

# Find the .025 and .975 quantile for delta star 
d = quantile(lambda_star, c(.025,.975))  

cat(d) 

#===========================================================
# Exercise 3 Parametric bootstrap CI for mean value 
#===========================================================

x<-c(38.43, 38.43, 38.39, 38.83, 38.45, 38.35, 38.43, 38.31, 38.32, 38.48, 38.50)

xbar<-mean(x)
sd<-0.57
nx<-11
B<-1000

boot_sample<-replicate(B,rnorm(nx,xbar,sd))
dim(boot_sample)
str(boot_sample)
boot_sample[,1]  # 1st bootsample
mean(boot_sample[,1]) # 1st mean 
boot_means<-apply(boot_sample,2,mean)
length(boot_means)

ci_mu<-quantile(boot_means,c(0.025,0.975))
print(ci_mu)

hist(boot_means, col="blue", nclass=30)

# Classical parametric CI

t.test(x)


#==================================================
# Exercise 4 - Non-parametric Bootstrap simulation
#=================================================

set.seed(1234)
x<-c(3.12,0.00,1.57,19.67, 0.22, 2.20)

x.star<-matrix(NA,1000,6)
i<-1
for (i in 1:1000){
  x.star[i,]<-sample(x,length(x),replace=T)
}

class(x.star)
dim(x.star)

mean.1000<-apply(x.star,1,mean)

meanB<-mean(mean.1000)  #bootstrap estimate of the mean

biasB<-meanB-mean(x)   #bias

hist(mean.1000)
qqnorm(mean.1000)

num<-NULL
i<-1
for (i in 1:1000){
  num[i]<-(mean.1000[i]-meanB)^2
}

se.B<-sqrt(sum(num)/999)



#==================================================================
#Exercise 5 non-parametric Bootstrap CI for correlation coefficien
#==================================================================

source("https://bioconductor.org/biocLite.R")
biocLite("multtest")

library(multtest);
data(golub)
x <- golub[2289,];
y <- golub[2430,]
cor(x,y)

xy<-matrix(c(x,y),ncol=2) #matrix containig the data

B <- 1000 #number of bootstrap samples
cor.star <- 0 #matrix containing correlation coefficients
for (i in 1:B){
  z<-sample(1:nrow(xy),replace=TRUE)
  cor.star[i] <- cor(xy[z,1],xy[z,2])
}

mean(cor.star) #bootstrap estimate of the correlation coefficient

plot(density(cor.star)) #bootstrap sampling distribution 

quantile(cor.star,c(0.025,0.975))#bootstrap 95% confidence interval

#parametric approach

install.packages("psychometric") # install package with function
library(psychometric) # load package with function

# The following command calculates lower and upper
# 95% confidence intervals (level)
# sample size (n) is 

CIr(r=cor(x,y), n = 38, level = .95)

#===========================================================================
#Exercise 6 - non-parametric bootstrap  test of hypothesis for mean value
#===========================================================================

library(multtest)
data(golub)
?golub
#golub.cl is a numeric vector indicating the tumor class,
#27 acute lymphoblastic leukemia (ALL) cases (code 0) 
#and 11 acute myeloid leukemia (AML) cases (code 1). 

gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

x <- golub[2058,golub.cl==0] # getting ALL expression levels


#non-parametric approach:bootstrap
n<-length(x)
mu0<-0

t.obs<-(mean(x)-mu0)*sqrt(n)/sd(x)
t.obs

z<-x-mean(x)+mu0
z
hist(z)


m <- 10000

t.star <- 0
for(j in 1:m)
{
  z.star <- sample(z,n,replace=T)
  t.star[j] <- (mean(z.star)-mu0)*sqrt(n)/sd(z.star)
}

p<-length(t.star[t.star>t.obs])/m  # one-sided
p


#parametric approach: t-test

hist(x,probability=T)
lines(density(x))
shapiro.test(x)  #normality test
t.test(x,alternative = "greater",mu = 0)



#================================================================
#Exercise 7 - Non parametric bootstrap -  Equality of two means
#================================================================

golub.cl
x1<-golub[1042,golub.cl==0]
x2<-golub[1042,golub.cl==1]


plot(density(x1),xlim=c(-2,3))  #empirical densities
lines(density(x2),col=2)
legend(-1,0.6,legend=c("ALL","AML"),col=1:2,lty=1)

n1<-length(x1)
n2<-length(x2)
xb1<-mean(x1)
xb2<-mean(x2)
vb1<-var(x1)
vb2<-var(x2)
t.obs<-(xb1-xb2)/sqrt(vb1/n1+vb2/n2)
t.obs
xb<-mean(c(x1,x2)) #combined mean of the two samples
z1<-x1-xb1+xb
z2<-x2-xb2+xb
mean(z1)
mean(z2)
xb
t.star<-0
B<-1000
for(i in 1:B){
  z1.star<-sample(z1,n1,replace=T)
  z2.star<-sample(z2,n2,replace=T)
  zb1<-mean(z1.star)
  zb2<-mean(z2.star)
  vz1<-var(z1.star)
  vz2<-var(z2.star)
  t.star[i]<-(zb1-zb2)/sqrt(vz1/n1+vz2/n2)
}
pvalue<-2*(sum(abs(t.star)>abs(t.obs))/B)
pvalue

#parametric approach
shapiro.test(x1)
shapiro.test(x2)
var.test(x1,x2)
test<-t.test(golub[1042,] ~ golub.cl, var.equal=T,alternative="two.sided")
test

#===============================================
#Exercise 8:  Permutation test 
#===============================================

X<-c(x1,x2) #combined sample
n<-n1+n2
B=1000
t.perm<-0
for(i in 1:B){
  indixes<-sample(1:n,n,replace=F)
  x1.star<-X[indixes[1:n1]]
  x2.star<-X[indixes[(n1+1):n]]
  m1.star<-mean(x1.star)
  m2.star<-mean(x2.star)
  v1.star<-var(x1.star)
  v2.star<-var(x2.star)
  t.perm[i]<-(m1.star-m2.star)/sqrt(v1.star/n1+v2.star/n2)
}
pvalor.perm<-2*(sum(abs(t.perm)>abs(t.obs))/B)
pvalor.perm



