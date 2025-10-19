rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/Scientific Reports")

# 加载必要包
library(haven)
library(tidyverse)
library(tableone)
library(nortest)  # 用于正态性检验
library(dplyr)
library(gtsummary)

# 读取数据
included_data1 <- read.csv("CHG-CMD_imputed.csv", stringsAsFactors = FALSE)

# 查看Diabetes_2020的编码
table(included_data1$Stroke_2020)

# 重新编码Diabetes_2020: "1" = 1 (有糖尿病), "2" = 0 (无糖尿病)
included_data1$Stroke_2020_numeric <- ifelse(included_data1$Stroke_2020 == "1", 1, 
                                               ifelse(included_data1$Stroke_2020 == "2", 0, NA))

# 检查编码结果
table(included_data1$Stroke_2020, included_data1$Stroke_2020_numeric, useNA = "always")

selected_vars <- c("Age", "Gender", "Education", "Smoking", "Drinking", 
                   "Residence", "Married_status", "Sleep",
                   "BMI", "BRI", "Glu", "TG", "HDL", "LDL", 
                   "TyG_BRI", "TyG","TyG_BMI","CTI","CRP","CHG","TC")

# 定义连续变量
continuous_vars <- c("Age", "Sleep", "BMI", "BRI", "Glu", "TG", "HDL", "LDL", 
                     "TyG", "TyG_BMI","TyG_BRI","CTI","CRP","CHG","TC")

# 进行正态性检验
nonnormal_vars <- character()

for (var in continuous_vars) {
  if (length(na.omit(included_data1[[var]])) > 3) {
    test_result <- ad.test(na.omit(included_data1[[var]]))
    if (test_result$p.value < 0.05) {
      nonnormal_vars <- c(nonnormal_vars, var)
      message(paste0("变量 ", var, " 是非正态分布的 (p = ", round(test_result$p.value, 4), ")"))
    } else {
      message(paste0("变量 ", var, " 是正态分布的 (p = ", round(test_result$p.value, 4), ")"))
    }
  } else {
    message(paste0("变量 ", var, " 没有足够数据进行分析"))
  }
}

# 创建TableOne表格 - 使用重新编码的结局变量
table1 <- CreateTableOne(
  vars = selected_vars,
  strata = "Stroke_2020_numeric",  # 按重新编码的结局分组
  data = included_data1,
  addOverall = TRUE,
  factorVars = c("Gender", "Education", "Smoking", "Drinking", 
                 "Residence", "Married_status")
)

# 打印结果
table_final1 <- print(table1,
                      nonnormal = nonnormal_vars,
                      showAllLevels = TRUE,
                      pDigits = 3,
                      catDigits = 1,
                      contDigits = 1,
                      smd = FALSE,
                      # 添加分组标签
                      groupLabels = c("No Stroke", "Stroke", "Overall"))

# 显示表格
print(table_final1)

# 创建包含分位数和边界值的详细报告
# 对CHG进行四分位数分组
quantiles_CHG <- quantile(included_data1$CHG, probs = seq(0, 1, 0.25), na.rm = TRUE)
included_data1$CHG_4 <- cut(included_data1$CHG,
                            breaks = quantiles_CHG,
                            labels = c('Q1', 'Q2', 'Q3', 'Q4'),
                            include.lowest = TRUE)

# 处理CHG的缺失值
included_data1$CHG_4[which(is.na(included_data1$CHG_4))] <- 'Q1'
table(included_data1$CHG_4, useNA = "always")

# 对TyG进行四分位数分组
quantiles_TyG <- quantile(included_data1$TyG, probs = seq(0, 1, 0.25), na.rm = TRUE)
included_data1$TyG_4 <- cut(included_data1$TyG,
                            breaks = quantiles_TyG,
                            labels = c('Q1', 'Q2', 'Q3', 'Q4'),
                            include.lowest = TRUE)

# 处理TyG的缺失值
included_data1$TyG_4[which(is.na(included_data1$TyG_4))] <- 'Q1'
table(included_data1$TyG_4, useNA = "always")

# 逻辑回归模型
# Model 1: 单变量模型
fit1 <- glm(
  Stroke_2020_numeric ~ factor(TyG_4),
  data = included_data1,
  family = quasibinomial
)
tbl1 <- tbl_regression(fit1, exponentiate = TRUE)
print(tbl1)

# Model 3: 调整人口学因素
fit3 <- glm(
  Stroke_2020_numeric ~ factor(TyG_4) + Age + Gender + Education + Smoking + Drinking,
  data = included_data1,
  family = quasibinomial
)
tbl3 <- tbl_regression(fit3, exponentiate = TRUE)
print(tbl3)

# Model 4: 全模型
fit4 <- glm(
  Stroke_2020_numeric ~ factor(TyG_4) +
    Age + Gender + Education + Smoking + Drinking + Residence + Married_status + 
    Sleep + BMI,
  data = included_data1,
  family = quasibinomial
)
tbl4 <- tbl_regression(fit4, exponentiate = TRUE)
print(tbl4)

# 可选：合并三个模型的结果
tbl_merged <- tbl_merge(
  list(tbl1, tbl3, tbl4),
  tab_spanner = c("**Model 1**", "**Model 3**", "**Model 4**")
)
print(tbl_merged)
