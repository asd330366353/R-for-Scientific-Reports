# 加载必要的包
library(mice)
library(dplyr)

# 读取数据
data <- read.csv("CHG-CMD.csv", na.strings = c("NA", ""))

# 查看数据结构及缺失值模式
str(data)
md.pattern(data)

# 将分类变量转换为因子
categorical_vars <- c("Gender", "Education", "Smoking", "Drinking", "Residence", 
                      "Married_status", "Diabetes", "Hypertension", "Stroke", 
                      "heart_problems", "Hypertension_2020", "Stroke_2020", 
                      "heart_problems_2020", "Diabetes_2020", "cardiometabolic_diseases")

data <- data %>%
  mutate(across(all_of(categorical_vars), as.factor))

# 执行多重插补（MICE）
set.seed(123)  # 确保结果可重现
# 创建自定义预测矩阵（默认全1矩阵）

imputed_data <- mice(
  data,
  m = 5,
  maxit = 10,
  method = "rf",  # 随机森林（推荐）或 "cart"
  print = FALSE
)

# 获取完整数据集（使用第一个插补结果）
complete_data <- complete(imputed_data, 1)

# 检查插补后的缺失值
sum(is.na(complete_data))

# 保存插补后的数据
write.csv(complete_data, "CHG-CMD_imputed.csv", row.names = FALSE)

# 查看插补效果示例
# 比较原始数据和插补数据
original_missing <- data %>% 
  select(height, weight, BMI, WC, BRI) %>% 
  summarise(across(everything(), ~sum(is.na(.))))

imputed_missing <- complete_data %>% 
  select(height, weight, BMI, WC, BRI) %>% 
  summarise(across(everything(), ~sum(is.na(.))))

print("原始数据缺失值数量:")
print(original_missing)
print("插补后缺失值数量:")
print(imputed_missing)