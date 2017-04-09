
# Exercise 07 Expectation Maximization
# -------------------------------------

# Step 1
# Read the data

nAABB<-11  # 0 recombinant gametes
nAABb<-14  # 1 recombinant gamete
nAAbb<-1   # 2 recombinant gametes
nAaBB<-10  # 1 recombinant gamete
nAaBb<-27  # 0 or 2 recombinant gametes
nAabb<-12  # 1 recombinant gamete
naaBB<-3   # 2 recombinant gametes
naaBb<-13  # 1 recombinant gamete
naabb<-11  # 0 recombinant gametes

n1 <- nAABb + nAaBB + nAabb + naaBb
n1

n2.star<-NULL
n2 <- nAAbb + naaBB + n2.star
# nAaBb = n2.star + n0.star

n <- n1 + nAAbb + naaBB + nAABB + nAaBb + naabb
n

# STEP 2
# Initialize theta (r)

r <- 0.3

# STEP 3 (E)
# Expected value for N2.star

expected <- function(r){
n2.star <- nAaBb*r^2/(r^2+(1-r)^2)
n2.star}

# STEP 4 (M)
# Update r

update.theta <- function(n2.star){
r <- (n1+2*(nAAbb+naaBB+n2.star))/(2*n)
r}

# STEPS 5
# EM Algoritm: iterative process

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
 cat(i,r,"\n") 
}

# STEP 6
# Print results

cat("\nThe final solution, after",i,"iterations, is r^ =",r,"\n")

