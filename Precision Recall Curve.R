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