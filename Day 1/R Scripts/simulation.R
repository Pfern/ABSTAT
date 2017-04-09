
                        #=================================#
                        #=====   Simluation           ====#
                        #=================================#


#===========================================================
# Exercise 1 Random numbers generation
#===========================================================


set.seed(1234)  # random number seed for sampling repeatability



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


#------------------------------------------
# Exercise 2 Multinomial simulation
#----------------------------------------

pa.star<-0.21
pb.star<-0.08
p0.star<-0.71



pA<-pa.star*(2-pa.star-2*pb.star)
pB<-pb.star*(2-pb.star-2*pa.star)
pAB<-2*pa.star*pb.star
p0<-(1-pa.star-pb.star)^2

pA+pB+pAB+p0



Result<-matrix(,5,4)
colnames(Result)<-c("A","B","AB","0")
rownames(Result)<-c("run1","run2","run3","run4","run5")
j<-1

for (j in 1:5) {
  
u<-runif(2128)  #simulation of uniform random numbers
A<-0
AB<-0
B<-0
n0<-0
i<-1
for(i in 1:2128){   #loop 
  
  if (u[i]>0 & u[i]<=pA)  #subinterval 1
    
    A<-A+1  
  
  if (u[i]>pA  & u[i]<=pA+pB) #subinterval 2
    B<-B+1
  
  if (u[i]>pA+pB  & u[i]<=pA+pB+pAB)  #subinteral 3
    AB<-AB+1
  
  if (u[i]>pA+pB+pAB & u[i]<=1)                 
    n0<-n0+1
  
}

Result[j,]<-c(A,B,AB,n0)
j<-j+1

}

print(Result)
print("A=725, B=258, AB=72, 0=1073")


#------------------------------------------
# Exercise 3
#----------------------------------------



bin<-rbinom(200,20,0.5) # B(20,0.5)
pois<-rpois(200,2)      # P(2)
bin2<-rbinom(200,200,0.9) # B(20,0.5)
pois2<-rpois(200,30)      # P(2)


par(mfrow=c(2,2))
hist(bin, main="B(20,0.5)",col="blue",prob=T)
lines(density(bin),col="red")

hist(pois, main="P(2)",col="blue",prob=T)
lines(density(pois),col="red")

hist(bin2, main="B(200,0.9)",col="blue",prob=T)
lines(density(bin2),col="red")

hist(pois2, main="P(30)",col="blue",prob=T)
lines(density(pois2),col="red")



#------------------------------------------
# Exercise 4
#----------------------------------------

mean<-120*(1/3)
mean
var<-mean*(1-(1/3))
var

bino_storm<-rbinom(100,120,1/3)

hist(bino_storm,main="Distribution of the number of stormy days in a season",col="cyan",prob=T)

mean_100<-100*mean

#------------------------------------------
# Exercise 5
#----------------------------------------

a.11<-pnorm(1.2,1.6,0.42)
b.11<-pnorm(2,1.6,0.42)-pnorm(2,1.6,0.42)
c.11<-pnorm(2.4,1.6,0.42)-pnorm(0.8,1.6,0.42)
d.11<-rnorm(1000,1.6,0.42)
sample_1000<-mean(d.11)






