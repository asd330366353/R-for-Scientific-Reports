# 读取数据并预处理
data <- read.csv("CHG-CMD_imputed.csv", header = TRUE)
data$cardiometabolic_diseases <- ifelse(data$cardiometabolic_diseases == "Yes", 1, 0)

# 进行多元逻辑回归分析
model <- glm(cardiometabolic_diseases ~ CHG + Age + Gender + Education + 
               Smoking + Drinking + Residence + Married_status + Sleep + BMI, 
             data = data, family = binomial)

# 输出回归结果摘要
summary(model)

# 计算并输出 OR 值和置信区间
or_ci <- exp(cbind(OR = coef(model), confint(model)))
print("Odds Ratios with 95% Confidence Intervals:")
print(or_ci)
