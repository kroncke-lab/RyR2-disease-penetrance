# Matthew O'Neill - 5/13/2021

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
