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
