

# 09 Multiple Testing
# -------------------

# Example 1
# ---------

raw.p.values<-c(0.042,0.001,0.031,0.014,0.007)
raw.p.values

# Bonferroni's Method

corr.p.values.Bonf<-p.adjust(p=raw.p.values,"bonferroni")
corr.p.values.Bonf

global.alpha<-0.05
which(corr.p.values.Bonf<global.alpha)

# Holm's Method

corr.p.values.Holm<-p.adjust(p=raw.p.values,"holm")
corr.p.values.Holm

which(corr.p.values.Holm<global.alpha)

# Benjamini and Hochbergs's Method

corr.p.values.BH<-p.adjust(p=raw.p.values,"BH")
corr.p.values.BH

which(corr.p.values.BH<global.alpha)


# Application
# -----------

# Read the data

library(limma)
# change dir - folder 'Swirl' - ctr+shift+H
# or
# setwd("C:/.../swirl")
swirl<-readTargets("SwirlSample.txt")
RG <- read.maimages(swirl$FileName, source="spot")
RG$genes <- readGAL("fish.gal")
RG$printer <- getLayout(RG$genes)

MA<-backgroundCorrect(RG,method="sub")
MA <- normalizeWithinArrays(MA,method="printtiploess")
MA <- normalizeBetweenArrays(MA,method="scale")

---
# RANKPROD

# Install and load package RankProd
source("http://bioconductor.org/biocLite.R")
biocLite("RankProd")
library(RankProd)

# Normalize data as in limma

MA$M[,2]<- (-1)*MA$M[,2]
MA$M[,4]<-(-1)*MA$M[,4]
data<-MA$M

# Calculate the number of replicates
k <- dim(data)[2]
k

# Associate to each replicate the number 1, meaning
# that the data is paires 
v<-c(rep(1,k))
v

# Run RankProd through function RP()
RP.out<- RP(data,v,num.perm=1000,logged=TRUE)
length(RP.out$pfp[RP.out$pfp<0.001])

# Visualize the genes with FDR (pfp) lower than 0.001
table<-topGene(RP.out,cutoff=0.001,logged=TRUE,logbase=2,method="pfp")
table

# Graphical representation of the FDR (pfp)
# plotRP(RP.out,cutoff=0.001)

# Visualize top40 genes (20 up and 20 down)
# table<-topGene(RP.out,num.gene=20,logged=TRUE,logbase=2,method="pfp")
# table


# ---------------------
# VOLCANO PLOT (not in the slides)

table.all<-topGene(RP.out,cutoff=1,logged=TRUE,logbase=2,method="pfp")
dim(table.all$Table1)
dim(table.all$Table2)
# column 3: FC:(class1/class2) 
# column 4: pfp
# column 5: P.value
# There are still some genes with FDR>1, which will not be represented,
# but they are for sure not DE

#The volcano plot can then be generated using the standard plot command:

# The par command sets "nice" graphical defaults
# Sum 0.0001 to pfp because of the zeros
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
plot(log2(table.all$Table1[,3]),-log10(table.all$Table1[,4]+0.00001),
     xlim=c(-2,2), ylim=c(0,5), #Set limits
     xlab="log2 fold change", ylab="-log10 FDR")
points(log2(table.all$Table2[,3]),-log10(table.all$Table2[,4]+0.00001)) 

# separate genes with FDR smaller than 0.001=10^(-3), ie,
# -log10(FDR)<-log10(10^(-3))=3
abline(3,0)

# mark those genes 
fdr1<-table.all$Table1[table.all$Table1[,4]<0.001,4]+0.00001 #pfp
fc1<-table.all$Table1[table.all$Table1[,4]<0.001,3] #FC
points(log2(fc1),-log10(fdr1),col="blue")
length(fdr1)

fdr2<-table.all$Table2[table.all$Table2[,4]<0.001,4]+0.00001 #pfp
fc2<-table.all$Table2[table.all$Table2[,4]<0.001,3] #FC
points(log2(fc2),-log10(fdr2),col="blue")
length(fdr2)
