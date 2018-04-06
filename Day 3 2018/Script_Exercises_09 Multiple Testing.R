
# Exercises 09 Multiple Testing
# -----------------------------


# Ex. 1
# ------

# Run these commands only once
install.packages("BiocInstaller",repos="https://bioconductor.org/packages/3.3/bioc")
source("http://bioconductor.org/biocLite.R")
biocLite("multtest")
#########################


library(multtest)
?mt.teststat

data(golub)
?golub

head(golub)
golub.cl # class of each of the 38 columns (indiv.)
ncol(golub) # 38 columns (individuals)
nrow(golub) # 3051 rows (genes)


# Ex. 2
# ------

? mt.teststat
welch_t<-mt.teststat(golub,golub.cl,test="t")  
head(welch_t)
# values of the welch t-statitic (by default).
# Others: ="t.equalvar", "wilcoxon", "f", "pairt", "blockf" (two-way anova)

qqnorm(welch_t)
qqline(welch_t)

wilk<-mt.teststat(golub,golub.cl,test="wilcoxon")  #values of the mann-whitney statistic 

qqnorm(wilk)    #normal QQ plot
qqline(wilk)



# Ex. 3
# ------

raw_p_t<-2*(1-pnorm(abs(welch_t)))   #raw p-values of the welcht t statistic
hist(raw_p_t,ylim = c(0,1500))
plot(sort(raw_p_t))

length(raw_p_t)*0.01 # 30.51 expected rejected hypothesis in 3051
length(raw_p_t[raw_p_t<0.01])  # 815 efectively rejected!


# Ex. 4
# ------

?mt.rawp2adjp 

procs = c("Bonferroni", "Holm", "BH")   
res = mt.rawp2adjp(raw_p_t, procs)
names(res)
adjp = res$adjp[order(res$index),]   #getting the adjust p-values
round(adjp,3)

hist(adjp[,2]) # adjp from Bonferroni
hist(adjp[,3]) # adjp from Holm
hist(adjp[,4],ylim = c(0,1500)) # adjp from Bonferroni


# Ex. 5
# ------

?mt.reject

mt.reject(adjp, c(0.01,0.05))$r 
mt.reject(adjp, seq(0,1, 0.05))$r  #table with several alphas and the corresponding number of HT rejected


# Ex. 6
# ------

#source("http://bioconductor.org/biocLite.R")
biocLite("fdrtool")
library(fdrtool)

fdr<-fdrtool(raw_p_t,statistic="pvalue")
names(fdr)
head(fdr$pval)
head(fdr$qval)     # vector with local q-values for each case
head(fdr$lfdr)     # vector with local fdr values for each case

# or

biocLite("qvalue")
library(qvalue)
fdr.q <- qvalue(raw_p_t)
names(fdr.q)
head(fdr.q$pvalues)
head(fdr.q$qvalues)
head(fdr.q$lfdr)

# or

fdr.p <- p.adjust(raw_p_t,"BH")
head(fdr.p)


