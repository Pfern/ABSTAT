
## Practical Exercise - EM Algorithm

The main locus for the blood type of mice is called Ag-B (B). Several alleles are associated to this locus but for some crossovers Mendel's laws do not seem to hold. A mating $AaBb\times AaBb\equiv F_1\times F_1$, originated a $F_2$ progeny, yielding


  | Genotype  | Frequency | Probability                     |
  |-----------|-----------|---------------------------------|
  | AABB      |    11     | ![](http://latex.codecogs.com/gif.latex?%24%281-%5Ctheta%29%5E2/4%24)               |
  | AABb      |    14     | ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%281-%5Ctheta%29/2%24)           |
  | AAbb      |     1     | ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%5E2/4%24)                    |
  | AaBB      |    10     | ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%281-%5Ctheta%29/2%24)            |
  | AaBb      |    27     | ![](http://latex.codecogs.com/gif.latex?%24%28%5Ctheta%5E2/2%29&plus;%5B%281-%5Ctheta%29%5E2%5D/2%24) |
  | Aabb      |    12     | ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%281-%5Ctheta%29/2%24)             |
  | aaBB      |     3     | ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%5E2/4%24)                    |
  | aaBb      |    13     | ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%281-%5Ctheta%29/2%24)             |
  | aabb      |    11     | ![](http://latex.codecogs.com/gif.latex?%24%281-%5Ctheta%29%5E2/4%24)                 |

Estimate the recombination fraction, $\theta$, from these data by the EM algorithm.

**Step 1**

* Read the data and state ao many recombinant gametes are there for each genotype.

<details><summary>Click Here to see the answer</summary><p>

```r
nAABB<-11  # 0 recombinant gametes
nAABb<-14  # 1 recombinant gamete
nAAbb<-1   # 2 recombinant gametes
nAaBB<-10  # 1 recombinant gamete
nAaBb<-27  # 0 or 2 recombinant gametes
nAabb<-12  # 1 recombinant gamete
naaBB<-3   # 2 recombinant gametes
naaBb<-13  # 1 recombinant gamete
naabb<-11  # 0 recombinant gametes
```
</p></details>
<br/>

* Calculate $n_1$, the number of individuals from 1 recombinant gametes ($\texttt{n1}$).

<details><summary>Click Here to see the answer</summary><p>

```r
n1 <- nAABb + nAaBB + nAabb + naaBb
n1
```
</p></details>
<br/>

* Calculate $n_2$, the number of individuals from 2 recombinant gametes ($\texttt{n2}$).

Note that $n_{AaBb}=n_2^*+n_0^*$.

<details><summary>Click Here to see the answer</summary><p>

```r
n2.star <- NULL
n2 <- nAAbb + naaBB + n2.star
```
</p></details>
<br/>

* Calculate $n$, the total number of individuals ($\texttt{n}$).

<details><summary>Click Here to see the answer</summary><p>

```r
n <- n1 + nAAbb + naaBB + nAABB + nAaBb + naabb
n
```
</p></details>
<br/>
<br/>

**Step 2**

* Initialize $\theta\in]0,0.5[$ ($\texttt{r}$).

<details><summary>Click Here to see the answer</summary><p>

```r
r <- 0.3
```
</p></details>
<br/>
<br/>

**Step 3 - E (Expectation)**

* Create function \texttt{expected} in order to calculate the expected value for $N_2^*$:

$N_2^*$: random variable representing the number of individuals from 2 recombinant gametes, among $n_{AaBb}$ individuals.

$N_2^*\frown Binomial(n_{AaBb},p)\quad$ with $\quad p=\dfrac{\theta^2}{\theta^2+(1-\theta)^2}$

then, $n_2^*=E(N_2^*)=n_{AaBb}\times p$


```r
expected <- function(r)
{
 n2.star <- nAaBb*r^2/(r^2+(1-r)^2)
 n2.star
}
```
</p></details>
<br/>
<br/>

**Step 4 - M (Maximization)**

* Create function \texttt{update.theta} in order to update $\theta$ according to:

$\theta=\dfrac{n_1+2(n_{AAbb}+n_{aaBB}+n_2^*)}{2n}$

meaning that the proportion of recombinant gametes is calculated as the total number of recombinant gametes (0, 1 or 2 for each individual) over the total number of gametes for $n$ individuals.

<details><summary>Click Here to see the answer</summary><p>

```r
update.theta <- function(n2.star)
{
r <- (n1+2*(nAAbb+naaBB+n2.star))/(2*n)
r
}
```
</p></details>
<br/>
<br/>

**Step 5 - Iterative procedure**

* Compute the cycle.

<details><summary>Click Here to see the answer</summary><p>

```r
i<-0
er<-1
error<-10^(-5)
while(er>=error)
{
 # Step E
 n2.star<-expected(r)
 # Step M
 r.updated<-update.theta(n2.star)
 # Stop criteria
 er<-abs(r-r.updated)
 i<-i+1
 r<-r.updated
 cat(i,r,"$\backslash$n") 
} 
```
</p></details>
<br/>
<br/>

**Step 6 - Print the results**

<details><summary>Click Here to see the answer</summary><p>

```r
cat("\nThe final solution, after",i,"iterations, is r* =",r,"\n")
```
</p></details>
<br/>
