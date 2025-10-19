# 加载必要的包
library(MatchIt)
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.csv("CHG-CMD_imputed.csv", stringsAsFactors = TRUE)

# 将因变量转换为二元变量（Yes=1, No=0）
data$cardiometabolic_diseases <- ifelse(data$cardiometabolic_diseases == "Yes", 1, 0)

# 查看CHG变量的分布并创建二分类变量
summary(data$CHG)
median_chg <- median(data$CHG, na.rm = TRUE)
data$CHG_binary <- ifelse(data$CHG > median_chg, 1, 0)

# 将分类变量转换为因子
data$Gender <- as.factor(data$Gender)
data$Education <- as.factor(data$Education)
data$Smoking <- as.factor(data$Smoking)
data$Drinking <- as.factor(data$Drinking)
data$Residence <- as.factor(data$Residence)
data$Married_status <- as.factor(data$Married_status)
data$CHG_binary <- as.factor(data$CHG_binary)

# 构建倾向得分模型并进行匹配
ps_model <- glm(CHG_binary ~ Age + Gender + Education + Smoking + Drinking + 
                  Residence + Married_status + Sleep + BMI,
                family = binomial(), data = data)

matched_data <- matchit(CHG_binary ~ Age + Gender + Education + Smoking + Drinking + 
                          Residence + Married_status + Sleep + BMI,
                        data = data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = 0.2)

# 提取匹配后的数据
matched_df <- match.data(matched_data)

# 将因子转换为数值
matched_df$CHG_binary_num <- as.numeric(as.character(matched_df$CHG_binary))
matched_df$cardiometabolic_diseases_num <- as.numeric(as.character(matched_df$cardiometabolic_diseases))

# 进行 McNemar 检验
mcnemar_table <- table(matched_df$CHG_binary_num, matched_df$cardiometabolic_diseases_num)
print("McNemar 检验表:")
print(mcnemar_table)

mcnemar_test <- mcnemar.test(mcnemar_table)
print("McNemar 检验结果:")
print(mcnemar_test)

# 提取 McNemar 表中的值
a <- mcnemar_table[1, 2] # 不一致对：对照组有病，处理组无病
b <- mcnemar_table[2, 1] # 不一致对：处理组有病，对照组无病

# 计算 Rosenbaum 边界
gamma_values <- seq(1, 3, by = 0.1)
p_values <- sapply(gamma_values, function(gamma) {
  # 计算调整后的 McNemar 统计量
  adjusted_stat <- (abs(a - b * gamma) - 1)^2 / (a + b * gamma)
  1 - pchisq(adjusted_stat, 1)
})

# 创建数据框用于绘图
sensitivity_df <- data.frame(Gamma = gamma_values, P_value = p_values)

# 找出 Gamma 值使得 p 值刚好超过 0.05
critical_gamma <- sensitivity_df$Gamma[which(sensitivity_df$P_value > 0.05)[1]]

# 绘制敏感性分析图
p <- ggplot(sensitivity_df, aes(x = Gamma, y = P_value)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "red", size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Rosenbaum Sensitivity Analysis for CHG Effect on Cardiometabolic Diseases",
       subtitle = paste("McNemar test p-value:", round(mcnemar_test$p.value, 4)),
       x = "Gamma (Hidden Bias Strength)",
       y = "P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  scale_x_continuous(breaks = seq(1, max(gamma_values), by = 0.5)) +
  annotate("text", x = max(gamma_values) * 0.7, y = 0.07, 
           label = "Significance threshold (0.05)", color = "red", size = 5)

# 添加关键Gamma值的标注
if (!is.na(critical_gamma)) {
  p <- p + 
    geom_vline(xintercept = critical_gamma, linetype = "dotted", color = "green", size = 1) +
    annotate("text", x = critical_gamma + 0.1, y = max(p_values) * 0.8, 
             label = paste("Critical Gamma =", round(critical_gamma, 2)), 
             color = "green", size = 5, angle = 90)
}

# 显示图形
print(p)

# 输出结论
if (!is.na(critical_gamma)) {
  cat(sprintf("\n结论: 结果对隐藏偏倚稳健，直到 Gamma = %.2f\n", critical_gamma))
  cat("这意味着需要有一个强度至少为", round(critical_gamma, 2), 
      "倍的未观测混杂因素才能推翻我们的结论。\n")
} else {
  cat("\n结论: 即使在最大 Gamma 值下，结果仍然显著\n")
  cat("这表明我们的结论非常稳健，即使存在较强的未观测混杂因素。\n")
}

# 保存图形
ggsave("rosenbaum_sensitivity_analysis.png", plot = p, width = 10, height = 8, dpi = 300)

