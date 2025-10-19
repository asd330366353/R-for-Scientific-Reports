# 安装并加载所需包
library(EValue)
library(dplyr)

# 读取数据
data <- read.csv("CHG-CMD_imputed.csv", stringsAsFactors = FALSE)

# 将 cardiometabolic_diseases 转为二分类变量（0/1）
data <- data %>%
  mutate(cardiometabolic_diseases = ifelse(cardiometabolic_diseases == "Yes", 1, 0))

# 检查结局变量是否罕见（患病率 < 10%）
prevalence <- mean(data$cardiometabolic_diseases, na.rm = TRUE)
is_rare <- prevalence < 0.1
print(paste("Outcome prevalence:", round(prevalence, 3)))
print(paste("Rare outcome assumption:", is_rare))

# 删除关键变量的缺失值
data_clean <- data %>%
  filter(!is.na(CHG), !is.na(TyG), !is.na(cardiometabolic_diseases))

# 拟合 logistic 回归模型（以 CHG 为例）
model_chg <- glm(cardiometabolic_diseases ~ CHG+Age + Gender + Education + Smoking + Drinking
                 + Residence + Married_status + Sleep, 
                 family = binomial(link = "logit"), 
                 data = data_clean)

# 提取 OR 和置信区间
or_chg <- exp(coef(model_chg)["CHG"])
ci_chg <- exp(confint(model_chg)["CHG", ])

# 计算 E-value - 明确指定 rare 参数
evalue_chg <- evalues.OR(est = or_chg, lo = ci_chg[1], hi = ci_chg[2], rare = is_rare)

# 打印结果
print("E-value for CHG:")
print(evalue_chg)

# 对 TyG 做同样分析
model_tyg <- glm(cardiometabolic_diseases ~ TyG+Age + Gender + Education + Smoking + Drinking
                 + Residence + Married_status + Sleep, 
                 family = binomial(link = "logit"), 
                 data = data_clean)

or_tyg <- exp(coef(model_tyg)["TyG"])
ci_tyg <- exp(confint(model_tyg)["TyG", ])

evalue_tyg <- evalues.OR(est = or_tyg, lo = ci_tyg[1], hi = ci_tyg[2], rare = is_rare)

print("E-value for TyG:")
print(evalue_tyg)

# 如果还想调整协变量，可以使用以下代码
# model_adjusted <- glm(cardiometabolic_diseases ~ CHG + Age + Gender + BMI,
#                      family = binomial(link = "logit"), 
#                      data = data_clean)
# 
# or_adjusted <- exp(coef(model_adjusted)["CHG"])
# ci_adjusted <- exp(confint(model_adjusted)["CHG", ])
# 
# evalue_adjusted <- evalues.OR(est = or_adjusted, 
#                              lo = ci_adjusted[1], 
#                              hi = ci_adjusted[2], 
#                              rare = is_rare)