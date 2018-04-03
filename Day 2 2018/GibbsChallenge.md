## Gibbs Sampling

At the end of this challenge, you should be able to:

(1) understand the Gibbs Sampler mechanism,

(2) how to run it for the Multinomial-Dirichlet distributions with known **x**,

(3) implement the main functions in ![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7BR%7D%24).
<br/>

**Step 1.** Simulate the data - Multinomial Distribution

* Learn about [Multinomial ditribution](https://en.wikipedia.org/wiki/Multinomial_distribution). For the three-dimensional case we have:

![eq1](http://latex.codecogs.com/gif.latex?%7B%5Cbf%20X%7D%3D%28X_1%2CX_2%2CX_3%29%5Cfrown%20Multinomial%28n%3B%5Ctheta_1%2C%5Ctheta_2%2C%5Ctheta_3%29)

where ![](http://latex.codecogs.com/gif.latex?X_i%5Cin%5C%7B0%2C1%2C2%2C%5Cldots%5C%7D) and ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta_1&plus;%5Ctheta_2&plus;%5Ctheta_3%3D1%5Cquad%20%28%5Ctheta_i%3E0%29%24)

![](http://latex.codecogs.com/gif.latex?%24%24P%5B%28X_1%2CX_2%2CX_3%29%3D%28x_1%2Cx_2%2Cx_3%29%5C%20%7C%5C%20%5Ctheta_1%2C%5Ctheta_2%2C%5Ctheta_3%5D%3D%5Cdfrac%7Bn%21%7D%7Bx_1%21x_2%21x_3%21%7D%5C%20%5Ctheta_1%5E%7B%5C%3Bx_1%7D%5Ctheta_2%5E%7B%5C%3Bx_2%7D%5Ctheta_3%5E%7B%5C%3Bx_3%7D%24%24)

* Search for the ![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7BR%7D%24) function which generates multinomially distributed random number vectors and computes multinomial probabilities.

<details><summary>Click Here to see the answer</summary><p>

```r
?rmultinom
```
</p></details>
<br/>

* Simulate 1 random vector, ![](http://latex.codecogs.com/gif.latex?%24%7B%5Cbf%20x%7D%3D%28x_%7B1%7D%2Cx_%7B2%7D%2Cx_%7B3%7D%29%24), following a Multinomial distribution with parameters n=1000 and ![](http://latex.codecogs.com/gif.latex?%7B%5Cboldsymbol%5Ctheta%7D%3D%28%5Ctheta_1%2C%5Ctheta_2%2C%5Ctheta_3%29%3D%280.2%2C0.3%2C0.5%29). Store the simulated data in an object named ![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7Bdata%7D%24).

<details><summary>Click Here to see the answer</summary><p>

```r
theta<-c(0.2,0.3,0.5)
data<-rmultinom(1,1000,theta)
data
```
</p></details>
<br/>

* Calculate the probability of observing the vector (220,350,430), that is, ![](http://latex.codecogs.com/gif.latex?P%5B%28X_1%2CX_2%2CX_3%29%3D%28220%2C350%2C430%29%5C%20%7C%5C%20%5Cboldsymbol%7B%5Ctheta%7D%5D).

<details><summary>Click Here to see the answer</summary><p>

```r
dmultinom(c(220,350,430),1000,theta)
```
</p></details>
<br/>
<br/>

**Step 2.** Dirichlet Distribution - Prior / Posterior Distribution

* Learn about [Dirichlet ditribution](https://en.wikipedia.org/wiki/Dirichlet_distribution). For the three-dimensional case we have:

![](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%5Ctheta%3D%28%5Ctheta_1%2C%5Ctheta_2%2C%5Ctheta_3%29%5Cfrown%20Dirichlet%28a_1%2Ca_2%2Ca_3%29)

where ![](http://latex.codecogs.com/gif.latex?a_1%2Ca_1%2Ca_3%3E0) and ![](http://latex.codecogs.com/gif.latex?a%3Da_1&plus;a_2&plus;a_3),

![](http://latex.codecogs.com/gif.latex?p_%7B%5Cboldsymbol%5Ctheta%7D%28%5Cboldsymbol%5Ctheta%29%3D%5Cdfrac%7B%5CGamma%28a%29%7D%7B%5CGamma%28a_1%29%5CGamma%28a_2%29%5CGamma%28a_3%29%7D%5C%20%5Ctheta_1%5E%7B%5C%3B%7B%5Ccolor%7Bblue%7Da_1%7D-1%7D%5Ctheta_2%5E%7B%5C%3B%7B%5Ccolor%7Bblue%7Da_2%7D-1%7D%5Ctheta_3%5E%7B%5C%3B%7B%5Ccolor%7Bblue%7Da_3%7D-1%7D)

* Search for the ![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7BR%7D%24) function which generates Dirichlet distributed random number vectors and computes Dirichlet probabilities.

<details><summary>Click Here to see the answer</summary><p>

```r
?rdirichlet
library(gtools)
```
</p></details>
<br/>


* Simulate a random vector ![](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%5Ctheta%5E%7B%280%29%7D%3D%28%5Ctheta_1%5E%7B%280%29%7D%2C%5Ctheta_2%5E%7B%280%29%7D%2C%5Ctheta_3%5E%7B%280%29%7D%29), from a Dirichlet distribution with hyperparameters ![](http://latex.codecogs.com/gif.latex?%7B%5Cbf%20a%7D%5E%7B%280%29%7D%3D%28a_1%5E%7B%280%29%7D%2Ca_2%5E%7B%280%29%7D%2Ca_3%5E%7B%280%29%7D%29%3D%281%2C1%2C1%29). Store the simulated vector in ![](http://latex.codecogs.com/gif.latex?%5Ctexttt%7Btheta.0%7D) and vector ![](http://latex.codecogs.com/gif.latex?%7B%5Cbf%20a%7D%5E%7B%280%29%7D) in ![](http://latex.codecogs.com/gif.latex?%5Ctexttt%7Ba.0%7D).

<details><summary>Click Here to see the answer</summary><p>

```r
a.0<-c(1,1,1)
theta.0<-rdirichlet(1,a.0)
theta.0
```
</p></details>
<br/>

* Calculate the probability density for the vector (0.15,0.25,0.6), that is, ![](http://latex.codecogs.com/gif.latex?p_%7B%5Cboldsymbol%5Ctheta%7D%5B%280.15%2C0.25%2C0.6%29%5D).

<details><summary>Click Here to see the answer</summary><p>

```r
ddirichlet(c(0.15,0.25,0.6),a.0) # Note: prob density function > 0
```
</p></details>
<br/>

* The Dirichlet distribution is the conjugate prior of the Multinomial distribution, by achiving the following result:   ![](http://latex.codecogs.com/gif.latex?p_%7B%5Cboldsymbol%5Ctheta%7C%5Cbf%20x%7D%28%5Cboldsymbol%5Ctheta%29%5Cpropto%20%5Ctheta_1%5E%7B%5C%3B%7B%5Ccolor%7Bblue%7Da_1&plus;x_1%7D-1%7D%5Ctheta_2%5E%7B%5C%3B%7B%5Ccolor%7Bblue%7Da_2&plus;x_2%7D-1%7D%5Ctheta_3%5E%7B%5C%3B%7B%5Ccolor%7Bblue%7Da_3&plus;x_3%7D-1%7D).

Simulate ![](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%5Ctheta%5E%7B%281%29%7D) knowing **x** and ![](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%5Ctheta%5E%7B%280%29%7D), directly from de posterior distribution (Dirichlet).

<details><summary>Click Here to see the answer</summary><p>

```r
data<-t(data) # need a row vector
updated.theta<-rdirichlet(1,a.0+data) # Note: prob density function > 0
updated.theta
```
</p></details>
<br/>
<br/>

**Step 3.** The Gibbs Sampler.

* Develop a function in ![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7BR%7D%24) for the Gibbs Sampler. Consider 5000 iterations.

* Store the updated parameters for each iteration in a  matrix of order ![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7Bnr.iter%7D%5Ctimes%203%24).

<details><summary>Click Here to see the answer</summary><p>

```r
nr.iter<-5000
theta.all.iter<-matrix(0,5000,3)

a<-a.0
for(i in 1:nr.iter)
{
  theta.all.iter[i,]<-rdirichlet(1,a) 
  a<-a+data
  # cat(i,"\n")
}
tail(theta.all.iter)  # last 6 rows
```
</p></details>
<br/>
<br/>

**Step 4.** Trace / Parameters Estimation

* Represent graphically the trace for each parameter along the 5000 iterations.

<details><summary>Click Here to see the answer</summary><p>

```r
par(mfrow=c(1,3))

plot(1:nr.iter,theta.all.iter[,1],ylim=c(0,0.5),main=expression(paste("Trace for ",theta[1])),ylab=expression(theta),xlab="iteration")
plot(1:nr.iter,theta.all.iter[,2],ylim=c(0,0.5),main=expression(paste("Trace for ",theta[2])),ylab=expression(theta),xlab="iteration")
plot(1:nr.iter,theta.all.iter[,3],ylim=c(0,0.5),main=expression(paste("Trace for ",theta[3])),ylab=expression(theta),xlab="iteration")
```
</p></details>
<br/>

* Evaluate the need of setting a period of burn-in.

* Find estimates of the parameters according to the decisions made in (b).

<details><summary>Click Here to see the answer</summary><p>

```r
theta.final<-apply(theta.all.iter[1001:5000,],2,mean) # 2: 'by column'
round(theta.final,2)
```
</p></details>
<br/>
<br/>
