                    ########################################################
                    ######        Project PCA and Bootstrap        #########
                    ########################################################
                    
library(ade4)
library(xtable)
# Preparing data
                    
data(meaudret)
names(meaudret)
mode(meaudret)
str(meaudret)

names(meaudret$env)
names(meaudret$design)
table(meaudret$design$season)
table(meaudret$design$site)
names(meaudret$spe)
mode(meaudret$spe.names)

meaudret$design
river<-meaudret$env
names(river)





#  descriptive analysis

#descriptive analysis
library(psych)
describe.by(river,season)
summary(river)


#create factors containing each individual's season and site and generate a colour vector to distinguish them.

site <- meaudret$design$site 
sitecol <- ifelse(site == c("S1","S2","S3","S4","S5"),c("blue","red","green","orange","cyan"))
season<-meaudret$design$season
seasoncol<-c("green","green","green","green","green","cyan","cyan","cyan","cyan","cyan","brown","brown","brown","brown","brown",
             "black","black","black","black","black")

#analysis of correlation
plot(river,col=seasoncol,pch=19)
cor(river)



# ------------------ Exercise 2 
#  a)Bootstrap: non-parametric approach for CI for the median

#normality of variables  (check for normality)

shapiro.test(river$Ammo)
shapiro.test(river$Nitr)
shapiro.test(river$Phos)

library(boot)

#we need a function of the parameter that we whant to estimate
#second argument "indices" is the indices of the observations for bootstrap sample



my.median = function(x, indices){
  return( median( x[indices] ) )
}

##Ammo

data.boot.ammo<-boot(river$Ammo,my.median,1000)  
mode(data.boot.ammo)
str(data.boot.ammo)
ci.boots.ammo<-boot.ci(data.boot.ammo,0.90,type="perc")
print(ci.boots.ammo)

##Nitr

data.boot.nitro<-boot(river$Nitr,my.median,1000)  
mode(data.boot.nitro)
str(data.boot.nitro)
ci.boots.nitro<-boot.ci(data.boot.nitro,0.90,type="perc")
print(ci.boots.nitro)

#Phos

data.boot.phos<-boot(river$Phos,my.median,1000)  
mode(data.boot.phos)
str(data.boot.nitro)
ci.boots.phos<-boot.ci(data.boot.phos,0.90,type="perc")
print(ci.boots.phos)


# b)  Hypothesis Test      H0: mu=3.5   vs H1: mu>3.5

#---------------- parametric approach
shapiro.test(river$Oxyd)
hist(river$Oxyd,probability=T)
lines(density(river$Oxyd))
t.test(x,alternative = "greater",mu = 3.5)
#------------------


#---bootstrap

n<-length(river$Oxyd)
mu0<-3.5

t.obs<-(mean(river$Oxyd)-mu0)*sqrt(n)/sd(river$Oxyd)
t.obs

z<-river$Oxyd-mean(river$Oxyd)+mu0
z
hist(z,prob=T)


m <- 10000

t.star <- 0
for(j in 1:m)
{
  z.star <- sample(z,n,replace=T)
  t.star[j] <- (mean(z.star)-mu0)*sqrt(n)/sd(z.star)
}

p<-length(t.star[t.star>t.obs])/m  # one-sided
p


# c)
my.data<-data.frame(meaudret$env$Oxyd,meaudret$design$season)
names(my.data)
names(my.data)<-c("Oxyd","season")

x1<-my.data[my.data$season=="winter",]
x2<-my.data[my.data$season=="autumn",]
x1<-x1[,1]
x2<-x2[,1]

my.data[my.data$site,]


plot(density(x1[,1]))  #empirical densities
lines(density(x2[,1]),col=2)
legend("topright",c("Series1","Series2","Series3"))
legend("topright",col=c("black","red"),c("winter","autumn"),lwd=1)

n1<-length(x1)
n2<-length(x2)
xb1<-mean(x1)
xb2<-mean(x2)
vb1<-var(x1)
vb2<-var(x2)
t.obs2<-(xb1-xb2)/sqrt(vb1/n1+vb2/n2)
t.obs2
xb<-mean(c(x1,x2)) #combined mean of the two samples
z1<-x1-xb1+xb
z2<-x2-xb2+xb
mean(z1)
mean(z2)
xb
t.star<-0
B<-1000
for(i in 1:B){
  z1.star<-sample(z1,n1,replace=T)
  z2.star<-sample(z2,n2,replace=T)
  zb1<-mean(z1.star)
  zb2<-mean(z2.star)
  vz1<-var(z1.star)
  vz2<-var(z2.star)
  t.star[i]<-(zb1-zb2)/sqrt(vz1/n1+vz2/n2)
}
pvalue<-(sum(abs(t.star)>abs(t.obs))/B)
pvalue

#parametric approach
shapiro.test(x1)
shapiro.test(x2)
var.test(x1,x2)
test<-t.test(x1,x2, var.equal=F,alternative="two.sided")
?t.test






#PCA analysis witg ade4


pca1 <- dudi.pca(river, scann = FALSE, nf = 3) 
names(pca1) 

#selection of the number of pc
#screeplot

pca1$eig
barplot(pca1$eig)

#How much variance is explained by each PC

(kip<-100*pca1$eig/sum(pca1$eig))
cumsum(kip)

#Scatter plot of individuals

# first plane (axes 1 and 2):
s.label(pca1$li, xax = 1, yax = 2)

par(mfrow=c(1,3))
s.label(pca1$li, xax = 1, yax = 2)
s.label(pca1$li, xax = 1, yax = 3)
s.label(pca1$li, xax = 2, yax = 3)

dev.off()

#add supplementary variable 

s.class(dfxy = pca1$li, fac =season , col =seasoncol , xax=2,yax=3)

#Scatter plot of variables

s.corcircle(pca1$co, xax = 1, yax = 2)

par(mfrow = c(3, 3))
for (i in 1:9) {
  plot(x = pca1$li[, 1], y = river[, i], pch = 19, col = gcol,
       xlab = "First axis of the PCA", las = 1, ylab = colnames(river)[i])
}

dev.off()

plot(river[,1],river[,2])



#Simultaneous plot of individuals and variables

scatter(pca1)

scatter(pca1, clab.row = 0, posieig = "none")

s.class(pca1$li, sex, col = gcol, add.plot = TRUE, cstar = 0, clabel = 0,
       cellipse = 0)





