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

# Load data and briefly process
setwd("G:/My Drive/KronckeLab/RyR2")
d <- read.csv("C:/Users/kohei/Downloads/Updated_REVEL_AlphaMissense_RYR2_unaff_renamed.csv")

# Include distances between residue centroids from structure of KCNQ1 PDB-ID: 6UZZ
ryr2dist <- read.csv(file = "Human_RyR2_tetramer_AF2_preliminary.dists.csv", header = FALSE)
source('func_dist_seq.R')

d$total_carriers<-as.integer(d$total_carriers)
d$unaff<-as.integer(d$unaff)
d$CPVT<-as.integer(d$CPVT)
colnames(d)[colnames(d) == "Residue_Number"] <- "resnum"
d$resnum<-as.integer(d$resnum)

# set initial weighting and penetrance
d$total_carriers[is.na(d$total_carriers)] <- 0
d$CPVT[is.na(d$CPVT)] <- 0
d$unaff[is.na(d$unaff)] <- 0

d$weight = 1-1/(0.01+d$total_carriers)
d$penetrance_cpvt <- d$CPVT/d$total_carriers
d[d$total_carriers < 1,"weight"] <- 0.000

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
d$cpvt_penetranceBayesian_initial <- (alpha0)/(alpha0 + beta0)
d$cpvt_penetranceBayesian <- (alpha0 + d[,"CPVT"])/((alpha0 + beta0 + d[,"total_carriers"]))

# Plot literature observed LQT1 penetrance versus residue number
fit <- loess(d[,"cpvt_penetranceBayesian_initial"]~as.numeric(d[,"resnum"]), span = 0.1)
plot(d$resnum, d$cpvt_penetranceBayesian_initial, xlab ="Residue", ylab = "CPVT Penetrance Estimate")

xrange <- seq(min(fit$x), max(fit$x), length.out = 100)
ps <- predict(fit, xrange, se=T)
lines(xrange, ps$fit*1, lwd=5)
lines(xrange, (ps$fit+1.96*ps$se.fit)*1, lty=2, lwd=4)
lines(xrange, (ps$fit-1.96*ps$se.fit)*1, lty=2, lwd=4)

d_resnums<-data.frame(seq(1,4967,1))
names(d_resnums)[names(d_resnums)=="seq.1..4967..1."]<-"resnum"
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

write.csv(d_full, file = "d_full_AM_REVEL_ClinVar_againstALL.csv", row.names = FALSE)

#definition of "solab"
solab <- function(mean, variance){
  alpha <- (mean^2 * (1 - mean) - variance * mean) / variance
  beta <- alpha * (1 / mean - 1)
  return(c(alpha, beta))
}

#Calculate an expectation algorithm to predict the unknown pathogenicity score from known data.
#Use an EM algorithm to estimate LQT2 PPV for each variant (or LQT2 diagnosis probability)
f <- read.csv("G:/My Drive/KronckeLab/RyR2/d_full_AM_REVEL_ClinVar_againstALL.csv")
covariates <- c("cpvt_dist", 'Alpha.missense.value',"REVEL")
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
write.csv(sub.tmp, file = "sub.tmp_AM_REVEL_againstALL.csv", row.names = FALSE)

##reload the data, so that you can restart from here
setwd("G:/My Drive/KronckeLab/RyR2")
sub.tmp <- read_csv("sub.tmp_AM_REVEL_againstALL.csv")

#check datas
View(sub.tmp)
plot(sub.tmp$resnum,sub.tmp$p_mean_w)
plot(sub.tmp$Alpha.missense.value,sub.tmp$p_mean_w)
plot(sub.tmp$cpvt_dist,sub.tmp$p_mean_w)
plot(sub.tmp$REVEL,sub.tmp$p_mean_w)