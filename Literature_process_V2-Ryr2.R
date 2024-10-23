---
title: "Bayes KCNQ1 Model Evaluation"
author: "Matthew O'Neill"
date: "7/6/2022"
output: 
  html_document:
    toc: true
    toc_depth: 4
    toc_float: false
      # collapsed: true
    smooth_scroll: true
    code_folding: hide
    highlight: zenburn #textmate
    theme: flatly
    # number_sections: true
editor_options: 
  chunk_output_type: console
---
# Matthew O'Neill - 5/13/2021

```{r}
library("readxl")
# Goal - Process the KCNQ1 Variants to get single unique variants with total counts, adopting from Japan processing script

# STATUS: Wildly inaccurate counts for final cohort - setting all to 0 pre-processing improved, but incomplete fix 

# Load data and briefly process 
setwd("C:\\Users/KRONCKE/Box Sync/Kroncke_Lab/RyR2/VariantDatabase/")
lit <- read_xlsx("RYR2_20220711.xlsx") # have resnum data for 434 of the total positions, seems like good coverage, not including BRAVO

# Need to define total carriers to include gnomAD

names(lit)[names(lit)=="asymptomatic (unaffected, mutation positive)"]<-"unaff"
lit$unaff <- as.integer(lit$unaff)
lit$unaff[is.na(lit$unaff)] <- 0
lit$CPVT <- as.integer(lit$CPVT)
lit$CPVT[is.na(lit$CPVT)] <- 0

lit$total_carriers <- lit$unaff + lit$CPVT # LOSE IT HERE! Ambiguous phenotypes too?! BRETT

lit_clean <- lit[,c("natAA", "resnum", "mutAA", "total_carriers", "unaff", "CPVT")] # removed unaff
lit_clean <- lit_clean[!is.na(lit_clean$resnum),]
lit_clean$var = paste(lit_clean$natAA, lit_clean$resnum, lit_clean$mutAA, sep = "")

# QC Check - total carriers much larger than LQT1 status here! 

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

source('func_dist_seq.R')

# Include distances between residue centroids from structure of KCNQ1 PDB-ID: 6UZZ
ryr2dist <- read.csv(file = "../Human_RyR2_tetramer_AF2_preliminary.dists.csv", header = FALSE)

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
for(resi in seq(1,4965,1)){
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
```

### Calculate EM priors and posteriors for all variants
Use an EM algorithm to estimate LQT2 PPV for each variant (or LQT2 diagnosis probability)

```{r, include=FALSE}

covariates <- c("traff_score", "lqt2_dist", "am_pathogenicity", "pore") # 
complete_rows_in_f <- complete.cases(f[, covariates])

tmp<-f[complete_rows_in_f,]

# Initialize variables
delta <- 10
count <- 0

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
  new_mean <- (alpha_f + tmp$lqt2) / (alpha_f + beta_f + tmp$total_carriers)
  
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

new_mean <- (alpha + tmp$lqt2)/(alpha + beta + tmp$total_carriers)
tmp$p_mean_w <- new_mean
tmp$alpha <- alpha
tmp$beta <- beta

sub.tmp<-tmp # reassign tmp so I can still use it in the next chunk

```

# Pattern Mixture Model, adding Posteriors for variants without trafficking

```{r}

covariates <- c("lqt2_dist", "am_pathogenicity", "pore") 
new_complete_rows_in_f <- complete.cases(f[, covariates]) 

tmp<-f[new_complete_rows_in_f,]

delta <- 10
count <- 0

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
  new_mean <- (alpha_f + tmp$lqt2) / (alpha_f + beta_f + tmp$total_carriers)
  
  delta <- 100 * sum(abs(new_mean - tmp$p_mean_w)) / length(tmp$var)
  tmp$p_mean_w <- new_mean 
  tmp[tmp$p_mean_w<0,"p_mean_w"]<-0.005
  mean_alpha <- mean(alpha_f)
  mean_beta <- mean(beta_f)
  print(paste("Mean alpha:", mean_alpha, "Mean beta:", mean_beta))
  print(paste("Delta:", delta, "Count:", count))
}

# when tuning parameter is 11, predictions are equivalent to 10 variant heterozygote phenotypes 
prior_mean <- tmp$p_mean_w
tmp$prior_mean <- tmp$p_mean_w
variance <- prior_mean*(1-prior_mean)
variance <- variance / (1 + nu)
ind_a <- seq(1, length(variance),1)
ind_b <- seq(length(variance)+1, length(variance)*2,1)
alpha <- solab(prior_mean,variance)[ind_a]
beta <- solab(prior_mean,variance)[ind_b]

new_mean <- (alpha + tmp$lqt2)/(alpha + beta + tmp$total_carriers)
tmp$p_mean_w <- new_mean
tmp$alpha <- alpha
tmp$beta <- beta

sub.tmp1 <- subset(tmp, !var %in% sub.tmp$var)

```

# Pattern Mixture Model, adding Posteriors for variants without trafficking or AlphaMissense score

```{r}

covariates <- c("lqt2_dist", "pore") 
new_complete_rows_in_f <- complete.cases(f[, covariates]) 

tmp<-f[new_complete_rows_in_f,]

delta <- 10
count <- 0

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
  new_mean <- (alpha_f + tmp$lqt2) / (alpha_f + beta_f + tmp$total_carriers)
  
  delta <- 100 * sum(abs(new_mean - tmp$p_mean_w)) / length(tmp$var)
  tmp$p_mean_w <- new_mean 
  tmp[tmp$p_mean_w<0,"p_mean_w"]<-0.005
  mean_alpha <- mean(alpha_f)
  mean_beta <- mean(beta_f)
  print(paste("Mean alpha:", mean_alpha, "Mean beta:", mean_beta))
  print(paste("Delta:", delta, "Count:", count))
}

# when tuning parameter is 11, predictions are equivalent to 10 variant heterozygote phenotypes 
prior_mean <- tmp$p_mean_w
tmp$prior_mean <- tmp$p_mean_w
variance <- prior_mean*(1-prior_mean)
variance <- variance / (1 + nu)
ind_a <- seq(1, length(variance),1)
ind_b <- seq(length(variance)+1, length(variance)*2,1)
alpha <- solab(prior_mean,variance)[ind_a]
beta <- solab(prior_mean,variance)[ind_b]

new_mean <- (alpha + tmp$lqt2)/(alpha + beta + tmp$total_carriers)
tmp$p_mean_w <- new_mean
tmp$alpha <- alpha
tmp$beta <- beta

sub.tmp2 <- subset(tmp, !var %in% c(sub.tmp$var, sub.tmp1$var))

```

#combine subsets

```{r}
tmp <- rbind(sub.tmp,sub.tmp1,sub.tmp2)

#adjust penetrance estimates to same scale with trafficking and peak tail
tmp$p_mean_w <- tmp$p_mean_w*4

```

# Calculate correlations of covariates with empirical posteriors
We next evaluate our model using Spearman correlations, Pearson correlations, and Brier scores (see below). 

## Spearman Lit and Cohort merged 

```{r}

# Data for text and Supplemental Table 1 
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

yList=c('lqt1_penetranceBayesian')
xList=c("lqt1_dist", "cardiacboost", 
        'hm_peak','hm_Vhalfact','hm_tauact','hm_taudeact',
        'ht_peak','ht_Vhalfact','ht_tauact','ht_taudeact',
        'pph2_prob', 'provean_score', 'revel_score', "blast_pssm", "prior_mean", "p_mean_w"
        ) 
tmp<-d[!is.na(d$penetrance_lqt1) & !is.na(d$revel_score),] 

resultTable<-calcAllPvals(yList, xList, 1000, 'weight', tmp)


```

## Pearson Lit and Cohort Merged 

```{r}

# Data for text and Supplemental Table 1 

yList=c('lqt1_penetranceBayesian')
xList=c("lqt1_dist", "cardiacboost", 
        'hm_peak','hm_Vhalfact','hm_tauact','hm_taudeact',
        'ht_peak','ht_Vhalfact','ht_tauact','ht_taudeact',
        'pph2_prob', 'provean_score', 'revel_score', "blast_pssm", "prior_mean", "p_mean_w"
        ) 
tmp<-d[!is.na(d$penetrance_lqt1) & !is.na(d$revel_score),] 
resultTable<-calcAllPvals(yList, xList, 1000, 'weight', tmp)

```


## Spearman Forest Plots for All Variant Heterozygotes

```{r}

d2 <- d[!is.na(d$total_carriers) & d$mut_type == "missense",]
tmp<-d2[!is.na(d2$provean_score) & !is.na(d2$revel_score),] 
yList=c('lqt1_penetranceBayesian')
xList=c( 
        'hm_peak','hm_Vhalfact','hm_tauact','hm_taudeact',
        'ht_peak','ht_Vhalfact','ht_tauact','ht_taudeact',
        "cardiacboost", 'pph2_prob', 'revel_score', 'provean_score', 
        "blast_pssm", "lqt1_dist", 'prior_mean_w', 'p_mean_w')
# Heterozygous data is empty in this set! Find out what happened
resultTable<-calcAllPvals(yList, xList, 1000, 'weight', tmp)

rm(tmp)
rm(t)
i=0
tmp<-data.frame()
for (x in xList){
  i=i+2
  tmp[i-1,"Feature"]<-x
  t<-d[!is.na(d[,x]) & d$total_carriers>0,]
  t<-t[!is.na(t[,"var"]),]
  foo <- boot(t, function(data,indices)
  weightedCorr(t[indices,x],t$penetrance_lqt1[indices], method="spearman", weights = t$weight[indices]), R=1000)
  tmp[i-1,"Spearman"]<-foo$t0
  tmp[i-1,"Spearman_low"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[1][[1]] 
  tmp[i-1,"Spearman_high"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[2][[1]]
  tmp[i-1,"n"]<-length(t[,x])
  
}

# Data for Figure 2 
forestplot(tmp$Feature,tmp$Spearman,tmp$Spearman_low,tmp$Spearman_high)


```


# 10-fold Cross Validation for Spearman, Pearson, and Brier scores 
We perform 10-fold cross validation to assess the optimism of our estimates across the dataset. We do this using 3 different metrics. 

## Import Data 

```{r}

# overlap of variants - remove non-uniques
d <- distinct(d, var, .keep_all = TRUE)
d <- d[d$total_carriers >0 , ]

# index rows randomly - need to randomize rows first in the final 
d$ID <- seq.int(nrow(d))
d <- d[sample(nrow(d)),]

# set prior and posterior to 0 when using saved or loaded data
d$prior_mean_w <- NULL
d$p_mean_w <- NULL
d$alpha <- NULL
d$beta <- NULL

# Only test on data for which we have carriers
d <- d[d$total_carriers > 0, ]

# Turn off warnings to speed computation
options(warn = -1)

```


## Functions for 10-fold CV

```{r}

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

# calculate Brier score
brier_calc <- function(covariate, outcome, count){
  score <- sum((covariate - outcome)^2)/count
  print(score)
}

```
