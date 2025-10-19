library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
data <- read.csv("CHG-CMD_imputed.csv", na.strings = c("NA", ""))

# 数据预处理
# 将目标变量转换为二分类数值变量（1 = "Yes", 0 = "No"）
data$disease_status <- ifelse(data$cardiometabolic_diseases == "Yes", 1, 0)

# 去除包含NA值的行
data_clean <- data[complete.cases(data[, c("BMI", "CHG", "TyG", "disease_status")]), ]

# 初始化结果存储
results <- data.frame()
roc_list <- list()
auc_ci_results <- data.frame()

# 定义要评估的指标
predictors <- c("CHG", "TyG")

# 对每个指标进行ROC分析
for (predictor in predictors) {
  # 拟合ROC曲线
  roc_obj <- roc(response = data_clean$disease_status, 
                 predictor = data_clean[[predictor]],
                 levels = c(0, 1), direction = "<")
  
  # 存储ROC对象
  roc_list[[predictor]] <- roc_obj
  
  # 计算AUC的95%置信区间
  auc_ci <- ci.auc(roc_obj)
  
  # 存储AUC和CI结果
  auc_ci_results <- rbind(auc_ci_results, data.frame(
    Predictor = predictor,
    AUC = round(auc_ci[2], 3),
    CI_lower = round(auc_ci[1], 3),
    CI_upper = round(auc_ci[3], 3),
    AUC_CI = sprintf("%.3f (%.3f-%.3f)", auc_ci[2], auc_ci[1], auc_ci[3])
  ))
  
  # 计算最佳阈值（Youden指数）
  coords <- coords(roc_obj, "best", best.method = "youden", transpose = TRUE)
  
  # 计算混淆矩阵
  predictions <- ifelse(data_clean[[predictor]] > coords["threshold"], 1, 0)
  cm <- table(Predicted = predictions, Actual = data_clean$disease_status)
  
  # 计算性能指标
  if (ncol(cm) == 2 && nrow(cm) == 2) {
    TP <- cm[2, 2]
    TN <- cm[1, 1]
    FP <- cm[2, 1]
    FN <- cm[1, 2]
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    accuracy <- (TP + TN) / sum(cm)
    youden_index <- sensitivity + specificity - 1
  } else {
    sensitivity <- specificity <- accuracy <- youden_index <- NA
  }
  
  # 添加到结果表
  results <- rbind(results, data.frame(
    Predictor = predictor,
    AUC_CI = sprintf("%.3f (%.3f-%.3f)", auc_ci[2], auc_ci[1], auc_ci[3]),
    Cutoff = round(coords["threshold"], 3),
    Sensitivity = round(sensitivity * 100, 1),
    Specificity = round(specificity * 100, 1),
    Accuracy = round(accuracy * 100, 1),
    Youden_Index = round(youden_index * 100, 1)
  ))
}

# 进行DeLong检验比较两个AUC
if (length(roc_list) == 2) {
  delong_test <- roc.test(roc_list[["CHG"]], roc_list[["TyG"]], method = "delong")
  
  # 输出DeLong检验结果
  cat("\n=== DeLong Test for AUC Comparison ===\n")
  cat("CHG AUC:", auc_ci_results$AUC_CI[auc_ci_results$Predictor == "CHG"], "\n")
  cat("TyG AUC:", auc_ci_results$AUC_CI[auc_ci_results$Predictor == "TyG"], "\n")
  cat("DeLong test statistic:", round(delong_test$statistic, 3), "\n")
  cat("P-value:", format.pval(delong_test$p.value, digits = 3), "\n")
  
  # 判断结果
  if (delong_test$p.value > 0.05) {
    cat("Conclusion: No statistically significant difference between CHG and TyG AUCs (P > 0.05)\n")
  } else {
    cat("Conclusion: Statistically significant difference between CHG and TyG AUCs (P < 0.05)\n")
  }
  cat("========================================\n\n")
}

# 创建包含AUC值和置信区间的图例标签
auc_labels <- sapply(predictors, function(p) {
  auc_val <- round(auc(roc_list[[p]]), 3)
  ci_lower <- round(ci.auc(roc_list[[p]])[1], 3)
  ci_upper <- round(ci.auc(roc_list[[p]])[3], 3)
  paste0(p, " (AUC = ", auc_val, ", 95% CI: ", ci_lower, "-", ci_upper, ")")
})

# 设置图例名称
names(roc_list) <- auc_labels

# 绘制ROC曲线
p <- ggroc(roc_list, size = 1.0) +
  geom_abline(intercept = 1, slope = 1, 
              linetype = "dashed", color = "black", size = 0.8) +
  labs(title = "",
       x = "1-Specificity",
       y = "Sensitivity") +
  scale_color_discrete(name = "Model") +
  theme_bw(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = c(0.7, 0.3),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    legend.background = element_rect(fill = "white"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_fixed(ratio = 1)

# 显示图形
print(p)

# 保存图形
ggsave("ROC_Curve_with_AUC_CI.png", plot = p,
       width = 14, height = 10, dpi = 300, bg = "white")

# 打印详细的AUC和CI结果
cat("\n=== Detailed AUC and Confidence Intervals ===\n")
print(auc_ci_results)

# 打印性能指标结果表
cat("\n=== Performance Metrics ===\n")
print(results)

# 保存结果到CSV文件
write.csv(auc_ci_results, "AUC_CI_results.csv", row.names = FALSE)
write.csv(results, "ROC_performance_metrics.csv", row.names = FALSE)

# 保存DeLong检验结果
if (exists("delong_test")) {
  delong_results <- data.frame(
    Test = "DeLong test for AUC comparison",
    Statistic = round(delong_test$statistic, 3),
    P_value = delong_test$p.value,
    Conclusion = ifelse(delong_test$p.value > 0.05, 
                        "No significant difference between AUCs", 
                        "Significant difference between AUCs")
  )
  write.csv(delong_results, "DeLong_test_results.csv", row.names = FALSE)
}
