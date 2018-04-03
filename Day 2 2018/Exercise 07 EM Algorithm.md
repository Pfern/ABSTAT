
## Practical Exercise - EM Algorithm

The main locus for the blood type of mice is called Ag-B (B). Several alleles are associated to this locus but for some crossovers Mendel's laws do not seem to hold. A mating AaBb x AaBb![](http://latex.codecogs.com/gif.latex?%5Cequiv%20F_1%5Ctimes%20F_1), originated a ![](http://latex.codecogs.com/gif.latex?%24F_2%24) progeny, yielding


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

Estimate the recombination fraction,![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%24), from these data by the EM algorithm.

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

* Calculate ![](http://latex.codecogs.com/gif.latex?%24n_1%24), the number of individuals from 1 recombinant gametes (![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7Bn1%7D%24)).

<details><summary>Click Here to see the answer</summary><p>

```r
n1 <- nAABb + nAaBB + nAabb + naaBb
n1
```
</p></details>
<br/>

* Calculate ![](http://latex.codecogs.com/gif.latex?%24n_2%24), the number of individuals from 2 recombinant gametes (![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7Bn2%7D%24)).

Note that ![](http://latex.codecogs.com/gif.latex?%24n_%7BAaBb%7D%3Dn_2%5E*&plus;n_0%5E*%24).

<details><summary>Click Here to see the answer</summary><p>

```r
n2.star <- NULL
n2 <- nAAbb + naaBB + n2.star
```
</p></details>
<br/>

* Calculate n, the total number of individuals (![](http://latex.codecogs.com/gif.latex?%24%5Ctexttt%7Bn%7D%24)).

<details><summary>Click Here to see the answer</summary><p>

```r
n <- n1 + nAAbb + naaBB + nAABB + nAaBb + naabb
n
```
</p></details>
<br/>
<br/>

**Step 2**

* Initialize ![](http://latex.codecogs.com/gif.latex?%5Ctheta%5Cin%5D0%2C0.5%5B%5Cquad%20%28%5Ctexttt%7Br%7D%29).

<details><summary>Click Here to see the answer</summary><p>

```r
r <- 0.3
```
</p></details>
<br/>
<br/>

**Step 3 - E (Expectation)**

* Create function ![](http://latex.codecogs.com/gif.latex?%5Ctexttt%7Bexpected%7D) in order to calculate the expected value for ![](http://latex.codecogs.com/gif.latex?%24N_2%5E*%24),

![](http://latex.codecogs.com/gif.latex?%5Ctexttt%7Bexpected%7D): random variable representing the number of individuals from 2 recombinant gametes, among ![](http://latex.codecogs.com/gif.latex?%24n_%7BAaBb%7D%24) individuals.

![](http://latex.codecogs.com/gif.latex?%24N_2%5E*%5Cfrown%20Binomial%28n_%7BAaBb%7D%2Cp%29%24) with ![](http://latex.codecogs.com/gif.latex?p%3D%5Cfrac%7B%5Ctheta%5E2%7D%7B%5Ctheta%5E2&plus;%281-%5Ctheta%29%5E2%7D)

then, ![](http://latex.codecogs.com/gif.latex?%24n_2%5E*%3DE%28N_2%5E*%29%3Dn_%7BAaBb%7D%5Ctimes%20p%24).

<details><summary>Click Here to see the answer</summary><p>

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

* Create function \texttt{update.theta} in order to update ![](http://latex.codecogs.com/gif.latex?%24%5Ctheta%24) according to:

![](http://latex.codecogs.com/gif.latex?%5Ctheta%3D%5Cfrac%7Bn_1&plus;2%28n_%7BAAbb%7D&plus;n_%7BaaBB%7D&plus;n_2%5E*%29%7D%7B2n%7D)

meaning that the proportion of recombinant gametes is calculated as the total number of recombinant gametes (0, 1 or 2 for each individual) over the total number of gametes for n individuals.

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
