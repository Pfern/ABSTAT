---
title: "Advanced Biostatistics 2018 - Practical Exercises"
#author: "carina"
#date: '2018-04-01'
output:
  html_document:
    keep_md: true
---

x <- rmarkdown::render("file.RMD", run_pandoc = FALSE, clean = FALSE)
knit_meta <- attr(x, "knit_meta") 
rmarkdown::render(input = 'file.knit.md', knit_meta = knit_meta )

<style type="text/css"> body, td { font-size: 18px; } code.r{ font-size: 18px; } pre { font-size: 16px } </style> 




## Probability Review


The follow exercises will cover the concepts of:

* Random variables
* Parameters
* Discrete Distributions
* Continuous Distributions



The exercises will be solved using R. If you have any questions press red light of your monitor.


###*Introduction*

Every random variable has an associated probability distribution function. This
function is called a probability mass function in the case of a discrete random
variable or probability density function in the case of a continuous random
variable. The distribution function is used to model how the outcome
probabilities are associated with the values of the random variable.



For practical computations R has built-in-functions for the binomial,
normal, t-Student, F, etc.,  where _d_ stands for density, _p_ for (cumulative) probability distribution, _q_ for quantiles, and _r_ for drawing
random samples.


In <https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html> you can find the R functions for several distributions.

> **Bernoulli Distribution $$X \sim Ber(1,\theta)$$**

We call _Bernoulli trial_ a single trial with
two possible outcomes: _success_ and _failure_, where $$\theta$$ is the probability of success.

* Parameter $$\theta\in[0,1]$$.


* Support $$\mathcal{X}=\{0,1\}$$.

* Probability mass function

$$P[X=x|\theta]=p(x|\theta)=\left\{
  \begin{array}{ll}
    \theta^x(1-\theta)^{1-x}, & \hbox{for $x=0,1$;} \\
    0, & \hbox{otherwise.}
  \end{array}
\right.
$$

* Population Moments
$$\mu=E[X]=\theta; \quad \sigma^2=var[X]=\theta(1-\theta).$$



**Exercise 1**   Let's say it is known that 2.6% of the population has a certain genetic disease. Then when we randomly
select one person and observe their disease status, we define X to be 1 if they have it and 0 if they don't. 

**1.1** Using R, what is the probablity of the selected person has no disease.



```r
dbinom(0,1,0.026)
```

**1.2** Plot the probability mass function (pmf) of this distribution


```r
  prob<-c(dbinom(0,1,0.026),dbinom(1,1,0.026))
  barplot(prob,ylab="P(X=k)",names.arg=c(0,1), width=1,xlim=c(0,4),ylim=c(0,1), main="Probability   mass function Ber(0.026)")
```




> **Binomial Distribution  $$X\sim Bin(n,\theta)$$**

The random variable which counts the number of successes in $$n$$ independent and identical
Bernoulli trials is called Binomial random variable.

The binomial distribution describes the behavior of a count variable $X$ if the following conditions apply:


1 The number of events $$n$$ is fixed.\\
2 Each observation is independent.\\
3 Each observation represents one of two outcomes (_success_ or _failure_).\\
4: The probability of _success_ $$\theta$$ is the same for each outcome.\\


* Parameter $$\theta\in[0,1]$$

* Support $$\mathcal{X}=\{0,1,...,n\}$$


* Probability mass function 
$$P[X=x|\theta]=p(x|\theta)=\left\{
  \begin{array}{ll}
   {n\choose x}\theta^x(1-\theta)^{n-x}, & \hbox{for $x\in S$;} \\
    0, & \hbox{otherwise.}
  \end{array}
\right.
$$


* Population Moments $$\mu=E[X]=n\theta,\quad \sigma^2=var[X]=n\theta(1-\theta).$$


**Exercise 2.** Suppose that for certain microRNA of size 20 the probability of a purine is binomially distributed with probability 0.7.

**2.1** What is the probability of having 14 purines?

**2.2** What is the probability of less than or equal to 14 purines?

**2.3** What is the probability of strictly more than 10 purines?

**2.4** How many purines do you expect? In other words: What is the mean of the distribution?

**2.5** What is the standard deviation of the distribution?



```r
#a)
dbinom(14,20,0.7)

#b)
pbinom(14,20,0.7)

#c)
1-pbinom(10,20,0.7)

#d)
mean_10<-20*0.7

#e)
sqrt(20*0.7*0.3)
```

**Exercise 3.** What probability does this R code represent: pbinom(15,20,0.7)-pbinom(10,20,0.7)+dbinom(10,20,0.7)




```r
# P(10<=X<=20)
```


> **Normal Distribution $$X \sim N(\mu,\sigma)$$**



* Parameters $$\mu\in \Re, \sigma>0$$, where $$\mu$$ is a location parameter and $$\sigma$$ is a scale parameter.

* Support $$x \in \Re$$

* Density distribution $$f_X(x)=\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{1}{2\sigma^2}(x-\mu)^2},\quad -\infty<x<+\infty.$$


* Population Moments
$$E[X]=\mu,\quad var[X]=\sigma^2.$$


**Exercise 4.** The distribution of the expression values of the ALL patients on the Zyxin gene are distributed according to N(1.6; 0.42).

** 4.1** Compute the probability that the expression values are smaller than 1.2?

** 4.2** What is the probability that the expression values are between 1.2 and 2.0?

** 4.3** What is the 15th quantile of that distribution. What does it mean?



```r
#a)
pnorm(1.2,1.6,0.42)

#b)
pnorm(2,1.6,0.42)-pnorm(1.2,1.6,0.42)

#c)
qnorm(0.15,1.6,0.42)
```


All normal distributions can be converted into the standard normal curve by subtracting the mean and dividing by the standard deviation: 

$$ Z=\frac{X-\mu}{\sigma}$$


  
    

then $Z$ has a Normal(0, 1) distribution. Such $Z$ is called a __standard normal random variable__.  The scale of $Z$ has no units and it is called the standardized scale. 


**Exercise 5.** Given a sample 0.12, 0.24, 0.01, 0.16, 0.18, 0.55,0.89, 1.00, 1.45 and 2.5 corresponding to intensity levels of one gene from 5 DNA chips.

Analyze the distribution considering normality using density function
in R and construct a normal QQ-plot to check for normality. Apply a
log2 transformation and repeat the analysis. Use the qqnorm() and
qqline() functions to get the normal QQ plot. Compare your results.


```r
set<-c(0.12, 0.24, 0.01, 0.16, 0.18, 0.55,0.89, 1.00, 1.45,
2.5)

qqnorm(set)
qqline(set)
```



## **Maximum Likelihood Estimation**


>Likelihood function is the function of the parameters given the data. The principle of maximum likelihood estimation is that somewhere on the range of parameter values, the likelihood function hits a maximum value. The parameter values for which the likelihood function obtains its maximum are called the maximum likelihood estimates for the parameters.


**Exercise 6**  A common use of maximum likelihood estimation in genetics is to estimate parameters for the [Hardy Weinberg Equilibrium](https://www.nature.com/scitable/definition/hardy-weinberg-equation-299) model for the distribution of genotypes in a population at equilibrium. In a given population gene pool, for a two allele gene locus, alleles A and a are at frequencies p and q, respectively. According to Hardy Weinberg equilibrium the genetic frequencies of genotypes are modeled by the simple equation: $p^2+2pq+q^2=1$.


  AA  |  Aa  |  aa
------|------|------
  p^2 |  2pq |  q^2
  
  The likelihood for this follows a [multinomial probability distribution](https://en.wikipedia.org/wiki/Multinomial_distribution) that we can model, where parameter theta,$\theta=p$ (and therefore $1-\theta=q$, theta between 0 and 1).


**6.1** Obtain the log likelihood function.

$$L(\theta) \propto \theta^{2n_{AA}}(2\theta(\theta))^{n_{Aa}}(1-\theta)^{2n_{aa}} $$

Letting n1=nAA and n2=nAa and n3=naa and logarithms (does not matter which base):

$$ log L(\theta)=2n_1log(\theta)+n_2log(2\theta(1-\theta))+2n_3log(1-\theta)$$

Eliminating the multiplicative term:

$$Log(\theta)=2n_1log(\theta)+n_2log(2)+n_2log(\theta)+n_2log(1-\theta)+2n_3log(1-\theta)$$


To analytically solve for the maximum likelihood take the first derivative of this (which we designate with a lowercase letter "l" prime):

$$  I'(\theta)=2n_1/\theta+n_2/\theta-n_2/(1-\theta) -2n_3/(1-\theta)$$
To find the value of theta that maximizes the above we could solve it
analytically (which is not hard in this case, but with other models is often analytically impossible) by setting $l'(\theta) =0$ and solving for theta, or we can solve it using R and finding the maximum value of the likelihood function over a grid of theta values, which is very simple.

**6.2** Suppose we collect data from a population sample of 500 (assume the
population obeys Hardy Weinberg equilibrium) and obtain the data in Table:

Genotype  |  Data  | Variable
----------|--------|----------
AA        |   134  | n_1
Aa        |   266  | n_2
aa        |   100  | n_3


Program R to solve for the maximum likelihood estimate of theta given
this data.



```r
#Create a grid of theta values
theta<-1:1000/1000

#Create a data vector to store (log)likelihood values
lik<-vector(length=1000)

#Enter data
n1<-134
n2<-266
n3<-100

#Given data, evaluate log likelihood
for(i in 1:1000){
 lik[i]<-2*n1*log(theta[i])+n2*log(2)+n2*log(theta[i])+
  n2*log(1-theta[i])+2*n3*log(1-theta[i])}


# Use which function to determine max value of lik

which(lik==max(lik))

lik[534]

#MLE value for theta (corresponding vector index to lik[534])
theta[534]
```

**6.3** Plot the results


```r
plot(theta,lik,xlab="theta",ylab="log likelihood",
main="MLE estimation for theta")
abline(v=theta[534],lty=2)
legend(x=0.54,y=-2000,legend="MLE theta=0.534")
```
## Simulation


The follow exercises will cover the concepts of:

* Simulation
* Parameters 

### Introduction

One of the greatest powers of using a computer lies in the ability to simulate things. What if you don't have any data, but you know the probability model that you want to work with? Through the power of simulation, you can use the computer to generate sample values for you. It's like doing an experiment, only in the virtual world instead of the real world. Simulation is an essential technique in computational statistics.


* Some distributions  has a corresponding function in R that starts with "r" for random. For example, if you whant to simulate a sample with normal distribution with mean value 2 and standard deviation 1 with size 20, you can use the above R code :



```r
x<-rnorm(20,mean=2,sd=1)
```




### What is simulation

* (Pseudo) random numbers generated from a computer, since
the method is completely deterministic. Thus, this numbers
have a similar behaviour from the random numbers.

* Random number generator is an algorithm that can generate $x_{i+1}$ from $x_i$.

* Require a start called "seed", i.e., a number that initiates the
deterministic/iterative process.

* The "seed" associated to a generator method (algorithm),
always produce the same sequence.


```r
set.seed(1234)
```

* This is an important characteristic, because that gives to the
user the possibility to reproduce exactly the same results.


* Basically the uniform distribution is simulated in this way.



> **Uniform Distribution $$X \sim U(a,b), a<b$$**

The Uniform distribution assigns equal probabilities to all possible ranges of equal length within which the r.v. can fall.
This distributions is widely used in bioinformatics, Bayesian analysis, quantitative genetics and so on.


One of the most important applications of the uniform distribution is in the generation of random numbers. That is, almost all random number generators generate random numbers on the (0,1) interval.

* Parameters $$ a<b \in (-\infty,\infty)$$

* Support  $$x \in (a,b)$$


* Density function $$f_X(x)=\frac{1}{b-a},\quad x\in [a,b]$$

* Population Moments $$E[X]=\frac{b+a}{2},\quad var[X]=\frac{(b-a)^2}{12}$$


**Exercise 6.** Generate random numbers from an Uniform
Distribution  with different size samples and construct
histograms and analyze them.



```r
# Get a vector of 4 numbers
runif(4)  #by default is generated an uniform in (0,1)


# Get a vector of 3 numbers from 0 to 100
runif(3, min=0, max=100)


# Get 3 integers from 0 to 100
# Use max=101 because it will never actually equal 101
floor(runif(3, min=0, max=101))  

# This will do the same thing
sample(1:100, 3, replace=T)

# To generate integers WITHOUT replacement:
sample(1:100, 3, replace=F)


#comparison of several samples generated from U(0,1)

x<-runif(100)
y<-runif(200)
z<-runif(300)
k<-runif(400)
w<-runif(500)
j<-runif(600)

par(mfrow=c(3,2))

hist(x,main="n=100",prob=T,xlab="")
abline(h=1,col="red")
hist(y,main="n=200",prob=T,xlab="")
abline(h=1,col="red")
hist(z,main="n=300",prob=T,xlab="")
abline(h=1,col="red")
hist(k,main="n=400",prob=T,xlab="")
abline(h=1,col="red")
hist(w,main="n=500",prob=T,xlab="")
abline(h=1,col="red")
hist(y,main="n=600",prob=T,xlab="")
abline(h=1,col="red")
```

**Exercise 7** Generate one sample of size 200 from each of the discrete Binomial and Poisson distributions. Compare the distributions with suitable plots. Change the parameters in each of the models to see how they infuence the distributions.



```r
bin<-rbinom(200,20,0.5) # B(20,0.5)
pois<-rpois(200,2)      # P(2)
bin2<-rbinom(200,200,0.9) # B(20,0.9)
pois2<-rpois(200,30)      # P(30)


par(mfrow=c(2,2))
hist(bin, main="B(20,0.5)",col="blue",prob=T)
lines(density(bin),col="red")

hist(pois, main="P(2)",col="blue",prob=T)
lines(density(pois),col="red")

hist(bin2, main="B(200,0.9)",col="blue",prob=T)
lines(density(bin2),col="red")

hist(pois2, main="P(30)",col="blue",prob=T)
lines(density(pois2),col="red")
```


