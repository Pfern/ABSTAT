################################################################################
#########             Principal Component Analysis                    ##########
################################################################################



#=========================================================
#       Toy Example 
#=========================================================

X<-matrix(c(2.5, 2.4,   #data: 2 vectors with dimension 10, variables are in the same scale
            0.5, 0.7,
            2.2, 2.9,
            1.9, 2.2,
            3.1, 3.0,
            2.3, 2.7,
            2, 1.6,
            1, 1.1,
            1.5, 1.6,
            1.1, 0.9),ncol=2,byrow=T)

print(X)

#plot the data

plot(X)

#analyze the correlation between the two variables

cor(X[,1],X[,2])

#center the data 

Xred<-matrix(0,nrow=nrow(X),ncol=ncol(X))

#subtract the mean (if the variables are on different units you should divide
# by the sd also

Xred[,1]<-X[,1]-mean(X[,1])
Xred[,2]<-X[,2]-mean(X[,2])
Xred

cor(Xred[,1],Xred[,2])


#get two plots side by side with the raw data and the centered data

par(mfrow=c(1,2))
plot(X,main="raw data")
text(X[,1]+0.1,X[,2]+0.1,seq(1:nrow(X)) ,col="blue", cex=1,font=3)

plot(Xred,main="centered data")
abline(h=0)
abline(v=0)
text(Xred[,1]+0.1,Xred[,2]+0.1,seq(1:nrow(X)) ,col="blue", cex=1,font=3)
dev.off()


#calculate the covariance (if the variables are on different scale you should use
# correlation instead)

S<-cov(Xred)   #function cor() for correlation matrix
S


#compute the eigenvalues and eigenvectors based on the covariance matrix

ev<-eigen(S)
names(ev)

ev$values   #eigenvalues

sum(ev$values) #sum of the eigenvalues equals the total of the variance

# in the centered data set
var(Xred[,1])+var(Xred[,2])

#or in raw data
var(X[,1])+var(X[,2])

sum(diag(S))

#scree plot - bar plot of the eigenvalues

barplot(ev$values)

# The eigenvectors are the principal components.

ev$vector   #eigenvectors=loadings

# the PC are ortogonal (are not correlated), means the inner produc is equal zero 

t(ev$vector[,1])%*%ev$vector[,2]


#the eigenvectors have norm 1

sum(ev$vector[,1]^2)
sum(ev$vector[,2]^2)


# Let's see to what extent each variable contributes

PC<-ev$vector
rownames(PC)<-c("X1","X2")
colnames(PC)<-c("PC1","PC2")

# Let's plot them on the centered data plot

plot(Xred[,1],Xred[,2])

pc1.slope = ev$vectors[1,1]/ev$vectors[2,1]  #slope of the eigenvectors
pc2.slope = ev$vectors[1,2]/ev$vectors[2,2]

abline(0,pc1.slope,col="red")   #plot the eigenvectors on the existing data (centered or standardized)
abline(0,pc2.slope,col="green")  # to see the relationship and the eigenvetors

text(1,0.5,"(0.678,0.735)",cex=0.5,col="red")
text(-1,1,"(-0.735,0.678)",cex=0.5,col="green")


# See how much variation each eigenvector accounts for

pc1.var = 100*round(ev$values[1]/sum(ev$values),digits=2)
pc2.var = 100*round(ev$values[2]/sum(ev$values),digits=2)
pc1.var
pc2.var

#To get the data expressed in terms of principal components
#Multiply the centered data by the eigenvectors (principal components)

scores<-t(ev$vector)%*%t(Xred) 
t(scores)


# BiPlot (see how variables correlate)

xlab=paste("PC1 - ",pc1.var," % of variation")
ylab=paste("PC2 - ",pc2.var," % of variation")

sd = sqrt(ev$values)
u<-t(scores)[,1]/sd[1]
v<-t(scores)[,2]/sd[2]
plot(u,v,main="My First BiPlot",xlab=xlab,ylab=ylab,type="p")
abline(0,0,col="red")
abline(0,90,col="green")


# First plot the variables as vectors
loadings = ev$vectors
arrows(0,0,loadings[,1]*sd[1],loadings[,2]*sd[2],length=0.1, lwd=2,angle=20, col="red")
text(loadings[,1]*sd[1],loadings[,2]*sd[2],c("var1","var2"), col="red", cex=0.9)

cor(X[,1],X[,2])

# Second plot the scores as points (cases)
points(t(scores)[,1]/sd[1],t(scores)[,2]/sd[2],pch=16,col="blue")
text(t(scores)[,1]/sd[1]+0.1,t(scores)[,2]/sd[2]+0.1,seq(1:nrow(t(scores))) ,col="black", cex=1,font=3)
X

X[5,]
X[1,]
X[3,]
X[2,]
X[10,]

#PCA analysis using princomp function

z<-princomp(X)
z
names(z)


z$sdev    #he standard deviations of the principal components. square root of eigenvalues
z$loadings  #the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors)
z$center  #the means that were subtracted.
z$scale   #the scale applied to each variable
z$scores # the scores of the supplied data on the principal components.


plot(z)   #screeplot

biplot(z)  


#PCA analysis using prcomp function

?prcomp

g<-prcomp(X)
names(g)

eigenvalues<-g$sdev^2  #eigenvalues 

sum(eigenvalues) #variance of the data

g$rotation
z$loading


#primcomp()  vs. prcomp()


par(mfrow=c(1,2))
plot(g)
biplot(g)


