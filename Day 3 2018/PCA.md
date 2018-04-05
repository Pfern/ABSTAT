


##  Principal Component Analysis


Consider the iris dataset (included with R) which gives the petal
width, petal length, sepal width, sepal length and species for 150
irises. To view more information about the dataset, enter
help(iris).


**1** Wich are the variables that you can use for PCA?

<details><summary>Click Here to see the answer</summary><p>

Just the quantitative ones: petal width, petal length, sepal width, sepal length.

</p></details>
<br/>
<br/>

**2.** Plot the data and analyze the correlation between the variables.

<details><summary>Click Here to see the answer</summary><p>
 
 ```{r}
  data<-iris[,1:4]  # object with the variables that we are going to do PCA
  
  plot(data)
  cor(data) 
  
 ```
</p></details>
<br/>
<br/>

**3.** Are you going to use the covariance or the correlation matrix to do your PCA analysis?

<details><summary>Click Here to see the answer</summary><p>
  
The covariance matrix, because the data are on the same units.

</p></details>
<br/>
<br/>

**4.** Does your data need to be centered or standardized?

<details><summary>Click Here to see the answer</summary><p>

Just centered, because the variables are on the same units.

```{r}
Data.red<-matrix(0,nrow=nrow(data),ncol=ncol(data))
dim(Data.red)


i<-1
for(i in 1:4){

Data.red[,i]<-data[,i]-mean(data[,i])

}
```
</p></details>
<br/>
<br/>

**5.** Compute the variance or correlation matrix acoordingly to the previous answer.

<details><summary>Click Here to see the answer</summary><p>

```{r}
V<-cov(data) #original data
print(V)

V2<-cov(Data.red) #centered data
print(V2)

# observe that the two matrices are the same
```

</p></details>
<br/>
<br/>

**6.** Compute the eigenvalues and eigenvectors based on the covariance/correlation matrix.

<details><summary>Click Here to see the answer</summary><p>
  
 ```{r}
ev<-eigen(V2)
names(ev)

ev$values   #eigenvalues

sum(ev$values) #sum of the eigenvalues equals the total of the variance
```
  
</p></details>
<br/>
<br/>

**7.** Obtain the scree plot and decide how many PC do you want to keep.

<details><summary>Click Here to see the answer</summary><p>
  
 ```{r}
totvar<-sum(ev$values) #sum of the eigenvalues equals the total of the variance

j<-1
var.ev<-0
for (j in 1:4){
  
var.ev[j]  <-(ev$values[j]/sum(ev$values))*100 
  
}

print(var.ev)

barplot(ev$values)
```

Two PC

</p></details>
<br/>
<br/>

**8.** Plot the Biplot for variables and cases on the same plot.

<details><summary>Click Here to see the answer</summary><p>
  
 ```{r}
pc1.var<-round(var.ev[1])
pc2.var<-round(var.ev[2])

xlab=paste("PC1 - ",pc1.var," % of variation")
ylab=paste("PC2 - ",pc2.var," % of variation")

sd = sqrt(ev$values)
u<-t(scores)[,1]/sd[1]
v<-t(scores)[,2]/sd[2]
plot(u,v,main="BiPlot",xlab=xlab,ylab=ylab,type="p")
abline(0,0,col="red")
abline(0,90,col="green")

# First plot the variables as vectors
loadings = ev$vectors
arrows(0,0,loadings[,1]*sd[1],loadings[,2]*sd[2],length=0.1, lwd=2,angle=20, col="red")
text(loadings[,1]*sd[1]+0.2,loadings[,2]*sd[2]+0.2,c("Sepal.Length", "Sepal.Width",  "Petal.Length", "Petal.Width" ), col="red", cex=0.9)

# Second plot the scores as points (cases)
points(t(scores)[,1]/sd[1],t(scores)[,2]/sd[2],pch=16,col="blue")
text(t(scores)[,1]/sd[1]+0.1,t(scores)[,2]/sd[2]+0.1,seq(1:nrow(t(scores))) ,col="black", cex=1,font=3)



```
  
</p></details>
<br/>
<br/>




