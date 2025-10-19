#
rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/")

# 加载必要包
library(haven)
library(tidyverse)
library(tableone)
library(nortest)  # 用于正态性检验
library(dplyr)
included_data1 <- read.csv("CHG-CMD_imputed.csv", stringsAsFactors = FALSE)

selected_vars <- c("Age", "Gender", "Education", "Smoking", "Drinking", 
                   "Residence", "Married_status", "Sleep",
                   "BMI", "BRI", "Glu", "TG", "HDL", "LDL", 
                   "TyG_BRI", "TyG","TyG_BMI","CTI","CRP","CHG","TC")

# 定义连续变量
continuous_vars <- c("Age", "Sleep", "BMI", "BRI", "Glu", "TG", "HDL", "LDL", "TyG", "TyG_BMI","TyG_BRI","CTI","CRP","CHG","TC")

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

# 创建TableOne表格
table1 <- CreateTableOne(
  vars = selected_vars,
  strata = "cardiometabolic_diseases",  # 按结局分组
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
                      groupLabels = c("No CMD", "CMD", "Overall"))

# 显示表格
print(table_final1)

library(gtsummary)
library(broom)
library(purrr)

included_data1$cardiometabolic_diseases<-ifelse(included_data1$cardiometabolic_diseases=="No",0,1)



# 创建包含分位数和边界值的详细报告
quantiles<-quantile(included_data1$CHG,probs = seq(0,1,0.25))
included_data1$CHG_4<-cut(included_data1$CHG,
                          breaks = quantiles,
                          labels=c('Q1','Q2','Q3','Q4'))
forcats::fct_count(included_data1$CHG_4)
included_data1$CHG_4[which(is.na(included_data1$CHG_4))]<-'Q1'
forcats::fct_count(included_data1$CHG_4) 




quantiles<-quantile(included_data1$TyG,probs = seq(0,1,0.25))
included_data1$TyG_4<-cut(included_data1$TyG,
                          breaks = quantiles,
                          labels=c('Q1','Q2','Q3','Q4'))
forcats::fct_count(included_data1$TyG_4)
included_data1$TyG_4[which(is.na(included_data1$TyG_4))]<-'Q1'
forcats::fct_count(included_data1$TyG_4) 


#model1
# Corrected Model 1
fit1 <- glm(
  cardiometabolic_diseases ~ factor(TyG_4) ,
  data = included_data1,
  family = quasibinomial
)
tbl_regression(fit1, exponentiate = TRUE)


fit3 <- glm(
  cardiometabolic_diseases ~ factor(TyG_4) + Age + Gender + Education + Smoking + Drinking,
  data = included_data1,
  family = quasibinomial
)
tbl_regression(fit3, exponentiate = TRUE)

# Corrected Model 4 (full model)
fit4 <- glm(
  cardiometabolic_diseases ~ factor(TyG_4) +
    Age + Gender + Education + Smoking + Drinking + Residence + Married_status + 
    Sleep+BMI,
  data = included_data1,
  family = quasibinomial
)
tbl_regression(fit4, exponentiate = TRUE)

