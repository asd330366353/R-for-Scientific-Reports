rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/Scientific Reports")

# 加载必要的包
library(pROC)
library(ggplot2)
library(boot)
library(ResourceSelection)

# 读取数据
data <- read.csv("CHG-CMD_imputed.csv", stringsAsFactors = TRUE)

# 创建二分类结局变量
data$cmd_outcome <- ifelse(data$cardiometabolic_diseases == "Yes", 1, 0)
table(data$cmd_outcome)

# 设置随机种子
set.seed(123)

# 构建三个逻辑回归模型
base_model <- glm(cmd_outcome ~ Age + Gender + Smoking + Drinking + Education + 
                    Residence + Married_status + Sleep + BMI, 
                  data = data, family = binomial)

chg_model <- glm(cmd_outcome ~ Age + Gender + Smoking + Drinking + Education + 
                   Residence + Married_status + Sleep + BMI + CHG, 
                 data = data, family = binomial)

tyg_model <- glm(cmd_outcome ~ Age + Gender + Smoking + Drinking + Education + 
                   Residence + Married_status + Sleep + BMI + TyG, 
                 data = data, family = binomial)

# 获取预测概率
data$pred_base <- predict(base_model, type = "response")
data$pred_chg <- predict(chg_model, type = "response")  
data$pred_tyg <- predict(tyg_model, type = "response")

# 1. 计算Brier分数
brier_base <- mean((data$pred_base - data$cmd_outcome)^2)
brier_chg <- mean((data$pred_chg - data$cmd_outcome)^2)
brier_tyg <- mean((data$pred_tyg - data$cmd_outcome)^2)

cat("Brier Scores:\n")
cat("Base model:", round(brier_base, 4), "\n")
cat("Base + CHG:", round(brier_chg, 4), "\n") 
cat("Base + TyG:", round(brier_tyg, 4), "\n\n")

# 2. Hosmer-Lemeshow检验
hl_base <- hoslem.test(data$cmd_outcome, data$pred_base, g = 10)
hl_chg <- hoslem.test(data$cmd_outcome, data$pred_chg, g = 10)
hl_tyg <- hoslem.test(data$cmd_outcome, data$pred_tyg, g = 10)

cat("Hosmer-Lemeshow Test P-values:\n")
cat("Base model:", round(hl_base$p.value, 4), "\n")
cat("Base + CHG:", round(hl_chg$p.value, 4), "\n")
cat("Base + TyG:", round(hl_tyg$p.value, 4), "\n\n")

# 3. 修正后的校准图函数 - 添加hl_pvalue参数
create_calibration_plot <- function(observed, predicted, model_name, color, hl_pvalue) {
  cal_data <- data.frame(obs = observed, pred = predicted)
  loess_fit <- loess(obs ~ pred, data = cal_data)
  
  cal_curve <- data.frame(
    predicted_prob = seq(0, 1, length = 100),
    ideal = seq(0, 1, length = 100)
  )
  cal_curve$observed_prob <- predict(loess_fit, newdata = cal_curve$predicted_prob)
  
  ggplot(cal_curve, aes(x = predicted_prob)) +
    geom_line(aes(y = observed_prob, color = model_name), size = 1.2) +
    geom_line(aes(y = ideal), linetype = "dashed", color = "gray50", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    # 在右上角添加HL检验P值标签
    annotate("text", x = 0.1, y = 0.9, 
             label = paste("H-L test: p =", 
                           ifelse(hl_pvalue < 0.001, "< 0.001", 
                                  format(round(hl_pvalue, 3), nsmall = 3))),
             size = 4.5, color = "black", hjust = 0, fontface = "bold") +
    labs(x = "Predicted Probability", 
         y = "Observed Proportion",
         title = paste("Calibration Plot -", model_name),
         color = "Model") +
    scale_color_manual(values = color) +
    theme_minimal() +
    xlim(0, 1) + ylim(0, 1) +
    theme(legend.position = "bottom")
}

# 生成校准图 - 现在传入HL检验的p值
p1 <- create_calibration_plot(data$cmd_outcome, data$pred_base, "Base Model", "blue", hl_base$p.value)
p2 <- create_calibration_plot(data$cmd_outcome, data$pred_chg, "Base + CHG", "green", hl_chg$p.value)
p3 <- create_calibration_plot(data$cmd_outcome, data$pred_tyg, "Base + TyG", "orange", hl_tyg$p.value)

print(p1)
print(p2)
print(p3)

# 4. 替代confusionMatrix的自定义函数
calculate_confusion_matrix <- function(observed, predicted, cutoff) {
  pred_class <- ifelse(predicted > cutoff, 1, 0)
  
  # 计算混淆矩阵
  TP <- sum(pred_class == 1 & observed == 1)
  TN <- sum(pred_class == 0 & observed == 0)
  FP <- sum(pred_class == 1 & observed == 0)
  FN <- sum(pred_class == 0 & observed == 1)
  
  # 计算性能指标
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  ppv <- TP / (TP + FP)
  npv <- TN / (TN + FN)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  return(list(
    table = matrix(c(TN, FP, FN, TP), nrow = 2, 
                   dimnames = list(Predicted = c(0, 1), Actual = c(0, 1))),
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    accuracy = accuracy
  ))
}

# 修改calculate_metrics函数
calculate_metrics <- function(observed, predicted, cutoff) {
  cm <- calculate_confusion_matrix(observed, predicted, cutoff)
  
  # 计算95% CI
  sens_ci <- binom.test(sum(predicted > cutoff & observed == 1), sum(observed == 1))$conf.int
  spec_ci <- binom.test(sum(predicted <= cutoff & observed == 0), sum(observed == 0))$conf.int
  ppv_ci <- binom.test(sum(predicted > cutoff & observed == 1), sum(predicted > cutoff))$conf.int
  npv_ci <- binom.test(sum(predicted <= cutoff & observed == 0), sum(predicted <= cutoff))$conf.int
  
  return(list(
    sensitivity = c(estimate = cm$sensitivity, lower = sens_ci[1], upper = sens_ci[2]),
    specificity = c(estimate = cm$specificity, lower = spec_ci[1], upper = spec_ci[2]),
    ppv = c(estimate = cm$ppv, lower = ppv_ci[1], upper = ppv_ci[2]),
    npv = c(estimate = cm$npv, lower = npv_ci[1], upper = npv_ci[2])
  ))
}

# 选择临床切点
roc_base <- roc(data$cmd_outcome, data$pred_base)
roc_chg <- roc(data$cmd_outcome, data$pred_chg)
roc_tyg <- roc(data$cmd_outcome, data$pred_tyg)

cutoff_base <- coords(roc_base, "best", ret = "threshold", transpose = TRUE)[1]
cutoff_chg <- coords(roc_chg, "best", ret = "threshold", transpose = TRUE)[1]  
cutoff_tyg <- coords(roc_tyg, "best", ret = "threshold", transpose = TRUE)[1]

cat("Optimal cutoffs (Youden index):\n")
cat("Base model:", round(cutoff_base, 3), "\n")
cat("Base + CHG:", round(cutoff_chg, 3), "\n")
cat("Base + TyG:", round(cutoff_tyg, 3), "\n\n")

# 计算各模型在最佳切点的性能指标
metrics_base <- calculate_metrics(data$cmd_outcome, data$pred_base, cutoff_base)
metrics_chg <- calculate_metrics(data$cmd_outcome, data$pred_chg, cutoff_chg)
metrics_tyg <- calculate_metrics(data$cmd_outcome, data$pred_tyg, cutoff_tyg)

# 打印结果函数
print_metrics <- function(metrics, model_name) {
  cat(model_name, "Model Metrics at Optimal Cutoff:\n")
  cat("Sensitivity:", round(metrics$sensitivity["estimate"], 3), 
      "95% CI:", round(metrics$sensitivity["lower"], 3), "-", 
      round(metrics$sensitivity["upper"], 3), "\n")
  cat("Specificity:", round(metrics$specificity["estimate"], 3),
      "95% CI:", round(metrics$specificity["lower"], 3), "-", 
      round(metrics$specificity["upper"], 3), "\n")
  cat("PPV:", round(metrics$ppv["estimate"], 3),
      "95% CI:", round(metrics$ppv["lower"], 3), "-", 
      round(metrics$ppv["upper"], 3), "\n")
  cat("NPV:", round(metrics$npv["estimate"], 3),
      "95% CI:", round(metrics$npv["lower"], 3), "-", 
      round(metrics$npv["upper"], 3), "\n\n")
}

print_metrics(metrics_base, "Base")
print_metrics(metrics_chg, "Base + CHG") 
print_metrics(metrics_tyg, "Base + TyG")

# 5. 汇总结果
results_summary <- data.frame(
  Model = c("Base", "Base + CHG", "Base + TyG"),
  Brier_Score = c(brier_base, brier_chg, brier_tyg),
  HL_P_value = c(hl_base$p.value, hl_chg$p.value, hl_tyg$p.value),
  Optimal_Cutoff = c(cutoff_base, cutoff_chg, cutoff_tyg),
  Sensitivity = c(metrics_base$sensitivity["estimate"], 
                  metrics_chg$sensitivity["estimate"],
                  metrics_tyg$sensitivity["estimate"]),
  Sensitivity_CI = paste0(round(c(metrics_base$sensitivity["lower"], 
                                  metrics_chg$sensitivity["lower"],
                                  metrics_tyg$sensitivity["lower"]), 3), "-",
                          round(c(metrics_base$sensitivity["upper"],
                                  metrics_chg$sensitivity["upper"], 
                                  metrics_tyg$sensitivity["upper"]), 3)),
  Specificity = c(metrics_base$specificity["estimate"],
                  metrics_chg$specificity["estimate"],
                  metrics_tyg$specificity["estimate"]),
  Specificity_CI = paste0(round(c(metrics_base$specificity["lower"],
                                  metrics_chg$specificity["lower"],
                                  metrics_tyg$specificity["lower"]), 3), "-",
                          round(c(metrics_base$specificity["upper"],
                                  metrics_chg$specificity["upper"],
                                  metrics_tyg$specificity["upper"]), 3))
)

print(results_summary)

# 保存结果
write.csv(results_summary, "model_performance_metrics1.csv", row.names = FALSE)

# 保存校准图
ggsave("calibration_base.png", p1, width = 6, height = 5, dpi = 300)
ggsave("calibration_chg.png", p2, width = 6, height = 5, dpi = 300) 
ggsave("calibration_tyg.png", p3, width = 6, height = 5, dpi = 300)

cat("分析完成！结果已保存到当前工作目录。\n")
