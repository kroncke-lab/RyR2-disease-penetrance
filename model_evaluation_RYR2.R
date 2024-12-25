
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

#precision recall curve, comes from Precision recall curve.R
# 必要なライブラリの読み込み
library(PRROC)
library(ROCR)
library(ggplot2)
library(readr)

# 作業ディレクトリの設定
setwd("G:/My Drive/KronckeLab/RyR2")
sub.tmp <- read_csv("sub.tmp_Model1_AM_REVEL_ClinVar.csv")

# データのフィルタリング
sub.tmp <- sub.tmp[sub.tmp$total_carriers == 1, ]
sub.tmp$cpvt_patho <- sub.tmp$CPVT

fglm <- sub.tmp[!is.na(sub.tmp$Alpha.missense.value), ]

# モデルの作成
modfunc1 <- glm(cpvt_patho ~ cpvt_dist_Model1, data = fglm, family = 'binomial')
modfunc2 <- glm(cpvt_patho ~ cpvt_penetranceBayesian, data = fglm, family = 'binomial')
modfunc3 <- glm(cpvt_patho~Alpha.missense.value, data = fglm, family = 'binomial')
modfunc4 <- glm(cpvt_patho~REVEL, data = fglm, family = 'binomial')
modfunc5 <- glm(cpvt_patho~Annotation_score, data = fglm, family = 'binomial')
modfunc6 <- glm(cpvt_patho~prior_mean, data = fglm, family = 'binomial')

# モデルリストの作成
funs <- list(
  list(name = "Model 2 (cpvt_penetranceBayesian)", model = modfunc2, color = "red"),
  list(name = "Model 3 (Alpha.missense.value)", model = modfunc3, color = "blue"),
  list(name = "Model 4 (REVEL)", model = modfunc4, color = "pink"),
  list(name = "Model 5 (Annotation_score)", model = modfunc5, color = "purple"),
  list(name = "Model 6 (prior_mean)", model = modfunc6, color = "orange")
)

# 各モデルのAUCを保存するリストを初期化
auc_values <- list()

# PRCの描画とAUCの計算
plot(0, 0, type = "n", xlab = "Recall", ylab = "Precision", 
     xlim = c(0, 1), ylim = c(0, 1), main = "Precision-Recall Curves for Multiple Models")

# Model 1 (modfunc1) の計算
predictions <- predict(modfunc1, type = "response")
labels <- fglm$cpvt_patho
pred <- prediction(predictions, labels)
perf <- performance(pred, "prec", "rec")

recall <- perf@x.values[[1]]  # Recall
precision <- perf@y.values[[1]]  # Precision

# NAを除外
valid_idx <- !is.na(recall) & !is.na(precision)
recall <- recall[valid_idx]
precision <- precision[valid_idx]

# 台形公式でAUCを計算
prc_auc <- sum(diff(recall) * (head(precision, -1) + tail(precision, -1)) / 2)
auc_values[["Model 1 (cpvt_dist_Model1)"]] <- prc_auc

# PRCのプロット
lines(recall, precision, col = "green", lwd = 2)

# 他のモデルを計算
for (li in funs) {
  predictions <- predict(li$model, type = "response")
  pred <- prediction(predictions, labels)
  perf <- performance(pred, "prec", "rec")
  
  recall <- perf@x.values[[1]]
  precision <- perf@y.values[[1]]
  
  # NAを除外
  valid_idx <- !is.na(recall) & !is.na(precision)
  recall <- recall[valid_idx]
  precision <- precision[valid_idx]
  
  # AUCの計算
  prc_auc <- sum(diff(recall) * (head(precision, -1) + tail(precision, -1)) / 2)
  auc_values[[li$name]] <- prc_auc
  
  # PRCをプロット
  lines(recall, precision, col = li$color, lwd = 2)
}

# AUCを水平線で表示
# AUCを水平線で表示
colors <- c("green", "red", "blue","pink", "purple", "orange")
plot(1:length(auc_values), unlist(auc_values), type = "n", xlab = "Models", ylab = "AUC",
     xlim = c(0.5, length(auc_values) + 0.5), 
     ylim = c(0, max(unlist(auc_values)) + 0.1), 
     main = "AUC for Each Model")

# 各モデルのAUCを水平線としてプロット
i <- 1
for (model_name in names(auc_values)) {
  abline(h = auc_values[[model_name]], col = colors[i], lwd = 2)  # 水平線を描画
  text(i, auc_values[[model_name]] + 0.02,  # 値の上にラベルを表示
       labels = paste0(round(auc_values[[model_name]], 3)), col = colors[i], cex = 0.8)
  points(i, auc_values[[model_name]], col = colors[i], pch = 16)  # ポイントを追加
  i <- i + 1
}

#Brier score
#definition of brier_calc
brier_calc <- function(p_i, y_i, N) {
  score <- sum((p_i - y_i)^2) / N
  return(score)
}

d <- read.csv("G:/My Drive/KronckeLab/RyR2/sub.tmp_Model1_AM_REVEL_ClinVar.csv")
d <- d[d$total_carriers == 1, ]

#y_i = actual outcome (disease status)
d$y_i <- 0
d$y_i <- d$CPVT
#d$y_i[d$p_mean_w > 0.5] <- 1
#d$y_i[d$p_mean_w <= 0.5] <- 0

#GLM logistic regression
fglm <- d[!is.na(d$Alpha.missense.value), ]
modfunc1 <- glm(y_i ~ cpvt_penetranceBayesian, data = fglm, family = 'binomial')
modfunc2 <- glm(y_i ~ cpvt_dist_Model1, data = fglm, family = 'binomial')
modfunc3 <- glm(y_i ~ Alpha.missense.value, data = fglm, family = 'binomial')
modfunc4 <- glm(y_i ~ REVEL, data = fglm, family = 'binomial')
modfunc5 <- glm(y_i ~ Annotation_score, data = fglm, family = 'binomial')
modfunc6 <- glm(y_i ~ prior_mean, data = fglm, family = 'binomial')

# make list
funs <- list(modfunc1, modfunc2, modfunc3, modfunc4, modfunc5, modfunc6)
names(funs) <- c("cpvt_penetranceBayesian", "cpvt_dist_Model1", "Alpha.missense.value", 
                 "REVEL", "Annotation_score", "prior_mean")

#definition of N
N <- sum(fglm$total_carriers, na.rm = TRUE)

#calcurate
i <- 0
for (name in names(funs)) {
  i <- i + 1
  model <- funs[[name]]
  
  #p_i=predicted probability for individual i
  p_i <- predict(model, type = "response")
  
  brier_score <- brier_calc(p_i, fglm$y_i, N)
  print(paste("Model:", name, "Brier Score:", brier_score))
  points(i, brier_score, col = colrs[i], pch = 19)
}
