library("readxl")
library("RSQLite")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(caret)
library(plotrix)
library(glmnet)
library(meta)
library(reshape2)
library(psych)
library(tableone)
library(wCorr)
library(rms)
library(boot)
library(leaps)
library(car)
library(reticulate)
library(rootSolve)
library(pROC)
library(wCorr)
library(MALDIquant)
library(tidyverse)      
library(lubridate)     
library(fpp2)          
library(zoo)            
library(latex2exp)
library(forestplot)
library(ggplot2)
library(caret)
library("nnet")
library("DBI")

# Goal - Process the CNQ1 Variants to get single unique variants with total counts, adopting from Japan processing script
# STATUS: Wildly inaccurate counts for final cohort - setting all to 0 pre-processing improved, but incomplete fix

# Load data and briefly process
setwd("G:/My Drive/KronckeLab/RyR2")
lit <- read_xlsx("RYR2_20241129.xlsx")

# Need to define total carriers to include gnomAD
names(lit)[names(lit)=="unaffected (mutation positive)"]<-"unaff"
lit$unaff <- as.integer(lit$unaff)
lit$unaff[is.na(lit$unaff)] <- 0
lit$CPVT <- as.integer(lit$CPVT)
lit$CPVT[is.na(lit$CPVT)] <- 0
lit$total_carriers <- lit$unaff + lit$CPVT 

lit_clean <- lit[,c("natAA", "resnum", "mutAA", "total_carriers", "unaff", "CPVT")]
lit_clean <- lit_clean[!is.na(lit_clean$resnum),]
lit_clean$var = paste(lit_clean$natAA, lit_clean$resnum, lit_clean$mutAA, sep = "")

res <- c(sum(lit_clean$total_carriers),sum(lit_clean$CPVT))
vars<-unique(lit_clean$var)
tmp<-data.frame(vars,"NA", "NA", "NA", 0, 0, 0)
tmp<-setNames(tmp, c("var", "natAA", "resnum", "mutAA", "total_carriers", "unaff", "CPVT"))
for (variant in vars){
  all_them <- lit_clean[lit_clean$var==variant,]
  cpvt_a <- sum(all_them$CPVT)
  unaff_a <- sum(all_them$unaff)
  total_carriers_a <- cpvt_a + unaff_a
  natAA_a <- all_them$natAA[1]
  resnum_a <- all_them$resnum[1]
  mutAA_a <- all_them$mutAA[1]
  tmp[tmp$var==variant, c("natAA", "resnum", "mutAA", "total_carriers", "unaff", "CPVT")]<-
    c(natAA_a, resnum_a, mutAA_a, as.integer(total_carriers_a), as.integer(unaff_a), as.integer(cpvt_a))
}

# Include distances between residue centroids from structure of KCNQ1 PDB-ID: 6UZZ
ryr2dist <- read.csv(file = "Human_RyR2_tetramer_AF2_preliminary.dists.csv", header = FALSE)
source('func_dist_seq.R')

d<-tmp
d$total_carriers<-as.integer(d$total_carriers)
d$unaff<-as.integer(d$unaff)
d$CPVT<-as.integer(d$CPVT)
d$resnum<-as.integer(d$resnum)
d<-d[d$total_carriers>0,]

# set initial weighting and penetrance
d$weight = 1-1/(0.01+d$total_carriers)
d$penetrance_cpvt <- d$CPVT/d$total_carriers
d[d$total_carriers < 1,"weight"] <- 0.000 # This is changed to "< 2" when evaluating ROC-AUC of n=1 variants from the literature

# Mean squared error
mse <- function(sm) {
  mean((sm$residuals)^2*(sm$weights))
}

# Weighted mean to determine LQT1 penetrance empirical prior
newdata = data.frame(wt=1)
model <- lm(penetrance_cpvt ~ 1, data=d, weights = d$weight)
summary(model)
p<-predict(model, newdata)
dev<- mse(model)#p*(1-p)

# Derive alpha and beta from weighted mean and MSE (estimated variance)
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# Estimated shape parameters for LQT1 empirical prior
alpha0 = estBetaParams(p,dev)$alpha
beta0 = estBetaParams(p,dev)$beta
print(paste("alpha0 = ", alpha0, "  beta0 = ", beta0))

# Bayesian LQT1 penetrance estimates from empirical priors
# and observed affected/unaffected counts:
d$cpvt_penetranceBayesian_initial <- (alpha0 + d[,"CPVT"])/((alpha0 + beta0 + d[,"total_carriers"]))
d$cpvt_penetranceBayesian<-d$cpvt_penetranceBayesian_initial

# Plot literature observed LQT1 penetrance versus residue number
fit <- loess(d[,"cpvt_penetranceBayesian_initial"]~as.numeric(d[,"resnum"]), span = 0.1)
plot(d$resnum, d$cpvt_penetranceBayesian_initial, xlab ="Residue", ylab = "CPVT Penetrance Estimate")

xrange <- seq(min(fit$x), max(fit$x), length.out = 100)
ps <- predict(fit, xrange, se=T)
lines(xrange, ps$fit*1, lwd=5)
lines(xrange, (ps$fit+1.96*ps$se.fit)*1, lty=2, lwd=4)
lines(xrange, (ps$fit-1.96*ps$se.fit)*1, lty=2, lwd=4)

d_resnums<-data.frame(seq(1,4965,1))
names(d_resnums)[names(d_resnums)=="seq.1..4965..1."]<-"resnum"
d_full<-merge(d_resnums,d,all = TRUE)
d_full[, "cpvt_dist"]<-NA
d_full[, "cpvt_dist_weight"]<-NA
for(resi in seq(1,4967,1)){
  l<-funcdist(resi, "var", d, ryr2dist, "penetrance_cpvt", "sigmoid", 7)
  d_full[d_full$resnum == resi, "cpvt_dist"]<-l[1]
  d_full[d_full$resnum == resi, "cpvt_dist_weight"]<-l[2]
}

# Plot literature observed LQT1 penetrance versus residue number
fit <- loess(d_full[,"cpvt_dist"]~as.numeric(d_full[,"resnum"]), span = 0.1)
plot(d_full$resnum, d_full$cpvt_dist, xlab ="Residue", ylab = "CPVT Penetrance Density")
xrange <- seq(min(fit$x), max(fit$x), length.out = 100)
ps <- predict(fit, xrange, se=T)
lines(xrange, ps$fit*1, lwd=5)
lines(xrange, (ps$fit+1.96*ps$se.fit)*1, lty=2, lwd=4)
lines(xrange, (ps$fit-1.96*ps$se.fit)*1, lty=2, lwd=4)

out<-d_full[,c("resnum","cpvt_dist")]
out<-unique(out)
out$cpvt_dist<-round(100*out$cpvt_dist)
write.table(out,file = "test.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.csv(d_full, file = "d_full_new.csv", row.names = FALSE)

#definition of "solab"
solab <- function(mean, variance){
  alpha <- (mean^2 * (1 - mean) - variance * mean) / variance
  beta <- alpha * (1 / mean - 1)
  return(c(alpha, beta))
}

#Calculate an expectation algorithm to predict the unknown pathogenicity score from known data.
#Use an EM algorithm to estimate LQT2 PPV for each variant (or LQT2 diagnosis probability)
f <- read.csv("G:/My Drive/KronckeLab/RyR2/d_full_Model1_AM_REVEL_ClinVar.csv")
covariates <- c("cpvt_dist_Model1", 'Alpha.missense.value',"REVEL")
complete_rows_in_f <- complete.cases(f[, covariates])

tmp<-f[complete_rows_in_f,]

# Initialize variables
delta <- 10
count <- 0

# Set initial values for p_mean_w using cpvt_penetranceBayesian_initial
tmp$p_mean_w <- tmp$cpvt_penetranceBayesian_initial  # Initialize p_mean_w with Bayesian initial penetrance estimates

# EM Loop
while(delta > 0.01 & count < 20) {
  count <- count + 1

  # Fit the model using 'f' with complete data for current covariates
  regression_formula <- as.formula(paste("p_mean_w ~", paste(covariates, collapse = " + ")))
  model <- glm(regression_formula, data = tmp, weights = tmp$weight)

  # Make predictions and update only those rows in 'tmp'
  prediction_output <- predict(model, newdata = tmp, se.fit = TRUE, type = "response")
  predictions <- pmin(pmax(prediction_output$fit, 0.0005), 1)
  variance_f <- (prediction_output$se.fit)^2

  # Assuming predictions and variance_f are vectors
  alpha_beta <- solab(predictions, variance_f)
  # Length of each segment (assuming alpha_beta is a single column)
  segment_length <- length(alpha_beta) / 2

  # Extract alpha values (first half of alpha_beta)
  alpha_f <- alpha_beta[1:segment_length]
  beta_f <- alpha_beta[(segment_length + 1):(2 * segment_length)]
  new_mean <- (alpha_f + tmp$CPVT) / (alpha_f + beta_f + tmp$total_carriers)

  delta <- 100 * sum(abs(new_mean - tmp$p_mean_w)) / length(tmp$var)
  tmp$p_mean_w <- new_mean
  tmp[tmp$p_mean_w<0,"p_mean_w"]<-0.005
  mean_alpha <- mean(alpha_f)
  mean_beta <- mean(beta_f)
  print(paste("Mean alpha:", mean_alpha, "Mean beta:", mean_beta))
  print(paste("Delta:", delta, "Count:", count))
}

# when tuning parameter is 11, predictions are equivalent to "nu" variant heterozygote phenotypes
nu <- 10
prior_mean <- tmp$p_mean_w
tmp$prior_mean <- tmp$p_mean_w
variance <- prior_mean*(1-prior_mean)
variance <- variance / (1 + nu)
ind_a <- seq(1, length(variance),1)
ind_b <- seq(length(variance)+1, length(variance)*2,1)
alpha <- solab(prior_mean,variance)[ind_a]
beta <- solab(prior_mean,variance)[ind_b]

new_mean <- (alpha + tmp$CPVT)/(alpha + beta + tmp$total_carriers)
tmp$p_mean_w <- new_mean
tmp$alpha <- alpha
tmp$beta <- beta

sub.tmp<-tmp
write.csv(sub.tmp, file = "sub.tmp_Model1_AM_REVEL_ClinVar.csv", row.names = FALSE)

##reload the data, so that you can restart from here
setwd("G:/My Drive/KronckeLab/RyR2")
sub.tmp <- read_csv("sub.tmp_Model1_AM_REVEL_ClinVar.csv")

#check datas
View(sub.tmp)
plot(sub.tmp$resnum,sub.tmp$p_mean_w)
plot(sub.tmp$Alpha.missense.value,sub.tmp$p_mean_w)
plot(sub.tmp$cpvt_dist_Model1,sub.tmp$p_mean_w)
plot(sub.tmp$REVEL,sub.tmp$p_mean_w)
plot(sub.tmp$Annotation_score,sub.tmp$p_mean_w)

# Calculate correlations of covariates with empirical posteriors
#We next evaluate our model using Spearman correlations, Pearson correlations, and Brier scores (see below). 
## Spearman Lit and Cohort merged 
#Data for text and Supplemental Table 1 
d <- read.csv("G:/My Drive/KronckeLab/RyR2/sub.tmp_Model1_AM_REVEL_ClinVar.csv")
calcPval=function(xName,yName,weightName,nPerms,new.mat2){
  # Pulls out variables
  x=new.mat2[,xName] 
  y=new.mat2[,yName] 
  w=new.mat2[,weightName]
  x2=x[!is.na(x)]
  y2=y[!is.na(x)]
  w2=w[!is.na(x)]
  
  # Calculate the real correlation
  realCorr=weightedCorr(x2,y2,method='spearman',weights=w2)
  # Do permutations, calculate fake correlations
  permutedCorrList=c()
  for(permNum in 1:nPerms){
    permutedX=sample(x2,length(x2),replace=FALSE)
    wCorrSim=weightedCorr(permutedX,y2,method='spearman',weights=w2)
    permutedCorrList=c(permutedCorrList,wCorrSim)
  }
  permutedCorrList2=abs(permutedCorrList)
  realCorr2=abs(realCorr)
  
  # Calculate pvalue
  summ=sum(realCorr2<permutedCorrList2)
  pValue=summ/nPerms
  return(list(realCorr,pValue,length(x2)))
}

calcAllPvals=function(yList,xList,nPerms,weightName,new.mat2){
  i=0
  resultTable=data.frame()
  for(yName in yList){
    for(xName in xList){
      i=i+1
      result=calcPval(xName,yName,weightName,nPerms,new.mat2)
      resultTable[i,'x']=xName
      resultTable[i,'y']=yName
      resultTable[i,'nPerms']=nPerms
      resultTable[i,'weightedCorr']=result[[1]]
      resultTable[i,'pValue']=result[[2]]
      resultTable[i,'n']=result[[3]]
      print(resultTable[i,'pValue'])
    }
  }
  print(resultTable)
  return(resultTable)
}

yList=c('cpvt_penetranceBayesian')
xList=c("cpvt_dist_Model1","Alpha.missense.value","REVEL","Annotation_score", 'prior_mean', 'p_mean_w')
tmp<-d[!is.na(d$cpvt_penetranceBayesian),] 

resultTable<-calcAllPvals(yList, xList, 1000, 'weight', tmp)

## Spearman Forest Plots for All Variant Heterozygotes
d2 <- d[!is.na(d$total_carriers),]
yList=c('cpvt_penetranceBayesian')
xList=c("cpvt_dist_Model1","Alpha.missense.value","REVEL","Annotation_score", 'prior_mean', 'p_mean_w')
resultTable<-calcAllPvals(yList, xList, 1000, 'weight', d2)

rm(d2)
i=0
tmp<-data.frame()
for (x in xList){
  i=i+2
  tmp[i-1,"Feature"]<-x
  t<-d[!is.na(d[,x]) & d$total_carriers>0,]
  t<-t[!is.na(t[,"var"]),]
  foo <- boot(t, function(data,indices)
    weightedCorr(t[indices,x],t$cpvt_penetranceBayesian[indices], method="spearman", weights = t$weight[indices]), R=1000)
  tmp[i-1,"Spearman"]<-foo$t0
  tmp[i-1,"Spearman_low"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[1][[1]] 
  tmp[i-1,"Spearman_high"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[2][[1]]
  tmp[i-1,"n"]<-length(t[,x])
}

# Data for Figure 2 
forestplot(tmp$Feature,tmp$Spearman,tmp$Spearman_low,tmp$Spearman_high)

# Part 4: ROC and AUC plots
## ROCs of observed cohort LQT2 diagnosis probability with 0.5 as cutoff all variants. 
setwd("G:/My Drive/KronckeLab/RyR2")
sub.tmp <- read_csv("sub.tmp_Model1_AM_REVEL_ClinVar.csv")

#sub.tmp <- sub.tmp[sub.tmp$total_carriers==1,]
#sub.tmp$cpvt_patho <- sub.tmp$CPVT
sub.tmp$cpvt_patho <- 0
sub.tmp$cpvt_patho[sub.tmp$p_mean_w > 0.5] <- 1
sub.tmp$cpvt_patho[sub.tmp$p_mean_w <= 0.5] <- 0

fglm <- sub.tmp[!is.na(sub.tmp$Alpha.missense.value), ]
#fglm <- sub.tmp[!is.na(sub.tmp$Alpha.missense.value) & sub.tmp$total_carriers == 1, ]
colrs <- c("green","blue","pink","purple","orange")
modfunc1 <- glm(cpvt_patho~cpvt_penetranceBayesian, data = fglm, family = 'binomial')
modfunc2 <- glm(cpvt_patho~cpvt_dist_Model1, data = fglm, family = 'binomial')
modfunc3 <- glm(cpvt_patho~Alpha.missense.value, data = fglm, family = 'binomial')
modfunc4 <- glm(cpvt_patho~REVEL, data = fglm, family = 'binomial')
modfunc5 <- glm(cpvt_patho~Annotation_score, data = fglm, family = 'binomial')
modfunc6 <- glm(cpvt_patho~prior_mean, data = fglm, family = 'binomial')

funs <- list(modfunc2, modfunc3,modfunc4,modfunc5,modfunc6)

i=0
par(pty="s")
tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(modfunc1)))], predict(modfunc1), ci=T)
plot.roc(tmp, col = "red")

for(li in funs){
  i=i+1
  tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(li)))], predict(li), ci=T)
  plot.roc(tmp,add = T, col = colrs[i])
}

## AUCs from ROCs observed cohort LQT2 diagnosis probability at multiple cutoffs
#sub.tmp <- read.csv("G:/My Drive/KronckeLab/RyR2/sub.tmp.csv")
# Wrapper for convenience of making GLMs and outputting AUCs
glm.mod=function(fglm,independent){
  in_string <- paste(independent, collapse=" + ")
  regression_formula <- as.formula(paste("cpvt_patho", in_string, sep=" ~ "))
  mod <- glm(regression_formula, data = fglm, family = 'binomial')
  tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(mod)))], predict(mod), ci=T)
  return(tmp$auc)
}

#sub.tmp$cpvt_patho <- sub.tmp$CPVT
#fglm <- sub.tmp[!is.na(sub.tmp$Alpha.missense.value) & sub.tmp$total_carriers == 1, ]
fglm <- unique(fglm)

colrs <- c("green", "blue","pink","purple","orange")
modfunc1 <- glm(cpvt_patho~cpvt_penetranceBayesian, data = fglm, family = 'binomial')
modfunc2 <- glm(cpvt_patho~cpvt_dist_Model1, data = fglm, family = 'binomial')
modfunc3 <- glm(cpvt_patho~Alpha.missense.value, data = fglm, family = 'binomial')
modfunc4 <- glm(cpvt_patho~REVEL, data = fglm, family = 'binomial')
modfunc5 <- glm(cpvt_patho~Annotation_score, data = fglm, family = 'binomial')
modfunc6 <- glm(cpvt_patho~prior_mean, data = fglm, family = 'binomial')

funs <- list(modfunc2, modfunc3,modfunc4,modfunc5,modfunc6)

i=0
par(pty="s")
tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(modfunc1)))], predict(modfunc1), ci=T)
plot.roc(tmp, col = "red")
for(li in funs){
  i=i+1
  tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(li)))], predict(li), ci=T)
  plot.roc(tmp,add = T, col = colrs[i])
  print(tmp$auc)
}

cutoffs <- seq(0.1,0.8,0.05)

cpvt_penetranceBayesian <- 0
cpvt_dist_Model1 <- 0
Alpha.missense.value <- 0
REVEL <-0
Annotation_score <-0
prior <-0
i=0
for (co in cutoffs){
  fglm$cpvt_patho[fglm$CPVT>=co] <- 1
  fglm$cpvt_patho[fglm$CPVT<co] <- 0
  print(paste(length(fglm$cpvt_patho), " ", sum(fglm$cpvt_patho)))
  if (!sum(fglm$cpvt_patho)<6){
    i=i+1
    
    cpvt_penetranceBayesian[i]<-glm.mod(fglm,"cpvt_penetranceBayesian")
    cpvt_dist_Model1[i]<-glm.mod(fglm,"cpvt_dist_Model1")
    Alpha.missense.value[i]<-glm.mod(fglm,"Alpha.missense.value")
    REVEL[i]<-glm.mod(fglm,"REVEL")
    Annotation_score[i]<-glm.mod(fglm,"Annotation_score")
    prior[i]<-glm.mod(fglm,"prior_mean")
  }
}

par(cex=1, bty='l', lwd=2)
plot(cutoffs[1:i],cpvt_penetranceBayesian,col="red",type = "l",ylim = c(0.5,1), ylab = "AUC")#, 

lines(cutoffs[1:i],cpvt_dist_Model1,col=colrs[1]) # green
lines(cutoffs[1:i],Alpha.missense.value,col=colrs[2]) # blue
lines(cutoffs[1:i],Alpha.missense.value,col=colrs[3]) # pink
lines(cutoffs[1:i],Annotation_score,col=colrs[4]) # purple
lines(cutoffs[1:i],prior,col=colrs[5]) # orange


#N>0
#fglm<-sub.tmp[!is.na(sub.tmp$Alpha.missense.value) & sub.tmp$total_carriers >0, ]
fglm<-fglm[!is.na(fglm$var),]
fglm <- unique(fglm)
fglm$cpvt_patho[fglm$cpvt_penetranceBayesian > 0.5] <- 1
fglm$cpvt_patho[fglm$cpvt_penetranceBayesian <= 0.5] <- 0

modfunc1 <- glm(cpvt_patho~cpvt_penetranceBayesian, data = fglm, family = 'binomial')
modfunc2 <- glm(cpvt_patho~cpvt_dist_Model1, data = fglm, family = 'binomial')
modfunc3 <- glm(cpvt_patho~Alpha.missense.value, data = fglm, family = 'binomial')
modfunc4 <- glm(cpvt_patho~REVEL, data = fglm, family = 'binomial')
modfunc5 <- glm(cpvt_patho~Annotation_score, data = fglm, family = 'binomial')
modfunc6 <- glm(cpvt_patho~prior_mean, data = fglm, family = 'binomial')

funs <- list(modfunc2, modfunc3,modfunc4,modfunc5,modfunc6)

i=0
par(pty="s")
tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(modfunc1)))], predict(modfunc1), ci=T)
plot.roc(tmp, col = "red")
for(li in funs){
  i=i+1
  tmp<-roc(fglm$cpvt_patho[row.names(fglm) %in% as.numeric(names(predict(li)))], predict(li), ci=T)
  plot.roc(tmp,add = T, col = colrs[i])
  print(tmp$auc)
}

cpvt_penetranceBayesian <- 0
cpvt_dist_Model1 <- 0 
Alpha.missense.value <- 0
Annotation_score<-0
prior <-0
i=0
for (co in cutoffs){
  fglm$cpvt_patho[fglm$cpvt_penetranceBayesian>=co] <- 1
  fglm$cpvt_patho[fglm$cpvt_penetranceBayesian<co] <- 0
  print(paste(length(fglm$cpvt_patho), " ", sum(fglm$cpvt_patho)))
  if (!sum(fglm$cpvt_patho)<6){
    i=i+1
    
    cpvt_penetranceBayesian[i]<-glm.mod(fglm,"cpvt_penetranceBayesian")
    cpvt_dist_Model1[i]<-glm.mod(fglm,"cpvt_dist_Model1")
    Alpha.missense.value[i]<-glm.mod(fglm,"Alpha.missense.value")
    REVEL[i]<-glm.mod(fglm,"REVEL")
    Annotation_score[i]<-glm.mod(fglm,"Annotation_score")
    prior[i]<-glm.mod(fglm,"prior_mean")
  }
}

par(cex=1, bty='l', lwd=2)
plot(cutoffs[1:i],cpvt_penetranceBayesian,col="red",type = "l",ylim = c(0.5,1), ylab = "AUC")#, 

lines(cutoffs[1:i],cpvt_dist_Model1,col=colrs[1]) # green
lines(cutoffs[1:i],Alpha.missense.value,col=colrs[2]) # blue
lines(cutoffs[1:i],REVEL,col=colrs[3]) # pink
lines(cutoffs[1:i],Annotation_score,col=colrs[4]) # purple
lines(cutoffs[1:i],prior,col=colrs[5]) # orange

# 5-fold Cross Validation for Spearman
d <- read.csv("sub.tmp_Model1_AM_REVEL_ClinVar.csv")
source("G:/My Drive/KronckeLab/RyR2/func_dist_seq.R")
ryr2dist <- read.csv(file = "G:/My Drive/KronckeLab/RyR2/Human_RyR2_tetramer_AF2_preliminary.dists.csv", header = FALSE)

# overlap of variants - remove non-uniques
d <- distinct(d, var, .keep_all = TRUE)
d <- d[d$total_carriers >0 , ]

# index rows randomly - need to randomize rows first in the final 
d$ID <- seq.int(nrow(d))
d <- d[sample(nrow(d)),]

# set prior and posterior to 0 when using saved or loaded data
d$prior_mean <- NULL
d$p_mean_w <- NULL
d$alpha <- NULL
d$beta <- NULL
d$cpvt_dist_Model1<- NULL
d$cpvt_dist_weight_Model1<- NULL
d$cpvt_dist<- NULL
d$cpvt_dist_weight<- NULL

# Turn off warnings to speed computation
options(warn = -1)

# Functions for 10-fold CV
regression <- function(dv, pivs, nivs, data) {
  # run a linear model with text arguments for dv and ivs - positive input variable, negative input variable 
  piv_string <- paste(pivs, collapse=" + ")
  niv_string <- paste(nivs, collapse=" - ")
  if(niv_string!="") iv_string <- paste(piv_string, " - ", niv_string, sep = "")
  if(niv_string=="") iv_string <- paste(piv_string)
  #print(iv_string)
  regression_formula <- as.formula(paste(dv, iv_string, sep=" ~ "))
  #print(regression_formula)
  glm(regression_formula, data, family = quasibinomial(link = "logit"), weights = data[,"weight"])
}

# solve for alpha and beta in Beta distribution
solab <- function(mean, variance){
  alpha <- (mean^2 * (1-mean) - variance * mean)/variance
  beta <- alpha * (1 / mean - 1)
  return(c(alpha,beta))
}

# make 5folds 
d2 <- d
k = 5
kfold = k
d3 <- createFolds(d2$ID, k=kfold)

# Make a vector to host Spearman correlations for ith fold 
spearman_cor <- vector("numeric", 5)
j = 0

# Working loop to set ith fold to weight 0! 
for(i in d3){
  j = j + 1
  
  # Part 1 - split data
  test <- d2[i,]
  test$weight <- 0
  train <- d2[-i,]
  full <- rbind(test, train) # 10 dataframes where ith variants have weight = 0 
  
  # part 1.5 - make cpvt_dist on each train data
  full$cpvt_dist <- NA
  full$cpvt_dist_weight <- NA
  for (resi in unique(train$resnum)) {
    l <- funcdist(resi, "var", train, ryr2dist, "penetrance_cpvt", "sigmoid", 7)
    full[full$resnum == resi, "cpvt_dist"] <- l[1]
    full[full$resnum == resi, "cpvt_dist_weight"] <- l[2]
  }
  
  # Part 2 - run EM algorithm
  full$p_mean_w <- full$cpvt_penetranceBayesian_initial
  covariates <- c("cpvt_dist", "Alpha.missense.value")
  delta <- 10
  count <- 0
  tmp <- full
  options(warn = -1)
  
  while(delta > 1 & count < 25){ # delta = 5 is roughly a change of 0.5%
    print(paste(delta, count))
    count <- count + 1
    alpha_f <- NULL
    beta_f <- NULL
    
    for(i in 1:nrow(tmp)){
      newdata = data.frame(var=tmp[i,"var"])
      newdata[covariates] <- tmp[i,covariates]
      model <- regression("p_mean_w", covariates, 
                          colnames(newdata)[colSums(is.na(newdata))>0], tmp)
      mean_f <- predict(model, newdata, type = "response")
      variance_f <- (predict(model, newdata,se.fit = T, type = "response")$se.fit)^2
      alpha <- solab(mean_f,variance_f)[1]
      beta <- solab(mean_f,variance_f)[2]
      tmp[i,"prior_mean"] <- mean_f
      if(alpha<0.01 | beta<0.01){
        alpha_f[i]=alpha0
        beta_f[i]=beta0
      }else{
        alpha_f[i]=alpha
        beta_f[i]=beta
      }
    }
    new_mean <- (alpha_f + tmp$CPVT)/(alpha_f + beta_f + tmp$total_carriers)
    full$prior_mean<- tmp$prior_mean
    delta <- sum(abs(new_mean-tmp$p_mean_w))
    tmp$p_mean_w <- new_mean
    full$p_mean_w <- tmp$p_mean_w
    print(delta)
  }
  
  full2 <- full[full$weight == 0, ]
  
  spearman_cor[j] <- weightedCorr(full2$cpvt_penetranceBayesian, full2$p_mean_w, method = 'spearman', weights = (1 - 1 / (0.01 + full2$total_carriers)))
}

print(spearman_cor)