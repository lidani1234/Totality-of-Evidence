#Loading data
library("sas7bdat")
all7 = read.sas7bdat(file="all7new.sas7bdat")
all20 = read.sas7bdat(file="all20new.sas7bdat")

#Study 7 treatment arm and outcome differences by patient
arm7 = all7$arm
WD7 = all7$ch6mwd
walk7 = all7$chwalk
climb7 = all7$chclimb
descend7 = all7$chdescend

#Study 20 treatment arm and outcome differences by patient
arm20 = all20$ARM
WD20 = all20$ch6mwd
walk20 = all20$chWALK
climb20 = all20$chclimb
descend20 = all20$chdescend

#Observed Z/t statistic (signs set so positive indicates improvement with Ataluren)
z1 = as.numeric(t.test(WD7 ~ arm7)$statistic)
z2 = -as.numeric(t.test(walk7 ~ arm7)$statistic)
z3 = -as.numeric(t.test(climb7 ~ arm7)$statistic)
z4 = -as.numeric(t.test(descend7 ~ arm7)$statistic)

z5 = as.numeric(t.test(WD20 ~ arm20)$statistic)
z6 = -as.numeric(t.test(walk20 ~ arm20)$statistic)
z7 = -as.numeric(t.test(climb20 ~ arm20)$statistic)
z8 = -as.numeric(t.test(descend20 ~ arm20)$statistic)

#Constants and storage for running permutation test
iter = 1000000
Z = numeric(length = iter)
Z.ind = matrix(nrow=iter, ncol=8)
colnames(Z.ind) = c("study07.6MWD", "study07.walk", "study07.climb", 
                    "study07.descend", "study20.6MWD", "study20.walk", 
                    "study20.climb", "study20.descend")
n7 = length(arm7); ni7 = n7/2
n20 = length(arm20); ni20 = n20/2
arm7.perm = c(rep("Ataluren", ni7), rep("Placebo", ni7) )
arm20.perm = c(rep("Ataluren", ni20), rep("Placebo", ni20) )

#Observed test statistics - included as one observation in permutation test
Z[1] = z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8
Z.ind[1,] = c(z1, z2, z3, z4, z5, z6, z7, z8)

#Permutation test storing total Z score and all 8 individual z scores
set.seed(91571)
for(i in 2:iter){
  #Permutation indecies
  perm7 = sample(1:n7, ni7, replace=FALSE)
  perm20 = sample(1:n20, ni20, replace=FALSE)
  
  #Permuted group assignments
  WD7.perm = c(WD7[perm7], WD7[-perm7])
  walk7.perm = c(walk7[perm7], walk7[-perm7])
  climb7.perm = c(climb7[perm7], climb7[-perm7])
  descend7.perm = c(descend7[perm7], descend7[-perm7])
  
  WD20.perm = c(WD20[perm20], WD20[-perm20])
  walk20.perm = c(walk20[perm20], walk20[-perm20])
  climb20.perm = c(climb20[perm20], climb20[-perm20])
  descend20.perm = c(descend20[perm20], descend20[-perm20])
  
  #Z/t statistics
  Z.ind[i,1] = as.numeric(t.test(WD7.perm ~ arm7.perm)$statistic)
  Z.ind[i,2] = -as.numeric(t.test(walk7.perm ~ arm7.perm)$statistic)
  Z.ind[i,3] = -as.numeric(t.test(climb7.perm ~ arm7.perm)$statistic)
  Z.ind[i,4] = -as.numeric(t.test(descend7.perm ~ arm7.perm)$statistic)
  
  Z.ind[i,5] = as.numeric(t.test(WD20.perm ~ arm20.perm)$statistic)
  Z.ind[i,6] = -as.numeric(t.test(walk20.perm ~ arm20.perm)$statistic)
  Z.ind[i,7] = -as.numeric(t.test(climb20.perm ~ arm20.perm)$statistic)
  Z.ind[i,8] = -as.numeric(t.test(descend20.perm ~ arm20.perm)$statistic)
  
  Z[i] = sum(Z.ind[i,])
}

#Summary of permutation Z-score sums
summary(Z)

#One-sided p-value
length(which(Z >= Z[1])) / iter

#Two-sided p-value
length(which(abs(Z) >= abs(Z[1]))) / iter

#Histogram of permutation mean Z-scores
h = hist(Z/8, breaks = seq(from=-2.8, to=2.8, by=0.2))
h$density = h$counts/sum(h$counts)
plot(h, freq=FALSE, main="",
     xlab="Simulated Average Z Scores", ylab="Simulated Proportion")
abline(v = Z[1]/8, col="blue", lty=2)