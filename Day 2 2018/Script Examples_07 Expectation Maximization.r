
# Examples 07 EM Algorithm
# ------------------------


# MLE
# ---

# This script estimates the allelic freq from genotype freq
# using Maximum Likelihood Estimation

nAA<-12
nAB<-35
nA0<-21
nBB<-15
nB0<-30
n00<-9
n<-nAA+nAB+nA0+nBB+nB0+n00 

p<-(nAA+nA0/2+nAB/2)/n
q<-(nBB+nB0/2+nAB/2)/n
r<-(n00+nA0/2+nB0/2)/n
p;q;r


# EM Algorithm
# ------------

# This script estimates the allelic freq from phenotype freq
# using EM Algorithm

nA<-33  # nAA and nA0 unknown: nA=nAA+nA0
nB<-45  # nBB and nB0 unknown: nB=nBB+nB0
nAB<-35 # nAB known
n0<-9   # n00=n0 known

n<-nA+nB+nAB+n0

# STEP 1

# Parameters' vector par=(p,q). Note that r=1-p-q

# HYP 1 - Uniform Distribution
p<-1/3
q<-1/3
r<-1/3

# HYP 2 - Bernstein
p<-1-sqrt((n0+nB)/n)
q<-1-sqrt((n0+nA)/n)
r<-sqrt(n0/n)
p;q;r

# STEP 2 (E)

# Calculate approximate values for nAA, nA0, nBB, nB0
expected<-function(p,q,r){
nAA<-nA*p/(p+2*r)
nA0<-nA*2*r/(p+2*r)
nBB<-nB*q/(q+2*r)
nB0<-nB*2*r/(q+2*r)
c(nAA,nA0,nBB,nB0)
}
expected(p,q,r)

# STEP 3 (M) 

# Update p, q and r
param<-function(nAA,nA0,nBB,nB0){
p<-(2*nAA+nA0+nAB)/(2*n)
q<-(2*nBB+nB0+nAB)/(2*n)
r<-(2*n0+nA0+nB0)/(2*n)
c(p,q,r)
}
param(nAA,nA0,nBB,nB0)

# STEP 4 e 5

# EM ALGORITHM - iterative procedure

i<-0
er<-1
erro<-10^(-5)

while(sum(er>=erro)>0)
{
 # Step E
 e<-expected(p,q,r)
 # Step M
 par<-param(e[1],e[2],e[3],e[4])
 # STOP criteria
 er<-abs(c(p,q,r)-c(par[1],par[2],par[3]))
 i<-i+1
 p<-par[1]; q<-par[2]; r<-par[3]
 cat(i,p,q,r,"\n") # print all the values obtained in each iterate.
}

cat("\nObtained solution after ",i," iterations:\n
p* =",p,
"\nq* =",q,
"\nr* =",r,"\n")


# Example of Application - nudge
# ------------------------------

# Download the datafiles from the web

# http://bioinf.wehi.edu.au/limma/
# go to the end of the webpage
# download the Zebrafish data - zip.file
# working directory - swirl
# Chip: 4x4 blocks, 22 rows and 24 col, each

# Change directory to folder 'Swirl'
# install and load package limma in order to read the data
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)

# Don't forget to download the data, as described previously
swirl <- readTargets("SwirlSample.txt")
swirl
RG <- read.maimages(swirl$FileName, source="spot")
RG
RG$genes <- readGAL("fish.gal")
# Get some information about the files
head(RG$genes)
RG$printer <- getLayout(RG$genes)
RG$printer
show(RG)
summary(RG$R)

# Create the data matrix
lr <- matrix(nrow=8448,ncol=4)
lr[,1] <- log2(RG$R[,1])-log2(RG$G[,1])
lr[,2] <- log2(RG$G[,2])-log2(RG$R[,2])
lr[,3] <- log2(RG$R[,3])-log2(RG$G[,3])
lr[,4] <- log2(RG$G[,4])-log2(RG$R[,4])
li <- log2(RG$R)+log2(RG$G)

# Install and load package nudge
biocLite("nudge")
library(nudge)

# Run packade nudge
result <- nudge1(logratio = lr,logintensity = li, dye.swap = T)

names(result)
result$mu
result$sigma
result$mixprob
result$a
result$b
result$iter

# Top 20 genes with the highest probability of being DE
s <- sort(result$pdiff, decreasing = T, index.return = T)
rownames(lr) <- RG$genes$Name
cbind(rownames(lr)[s$ix[1:20]], round(s$x[1:20], 2))

# number of genes with a probability of DE, greater than 0.5
thresh <- 0.5
sum(result$pdiff >= thresh)
