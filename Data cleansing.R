# 清空环境并设置工作目录
##
rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/CHARLS")

# 加载必要包
library(haven)
library(tidyverse)
library(tableone)
library(nortest)  # 用于正态性检验
library(dplyr)

# 读取2011年基线数据
household <- read_dta("2011/household_roster.dta")
family <- read_dta("2011/family_information.dta")
demographic <- read_dta("2011/demographic_background.dta")
health <- read_dta("2011/health_status_and_functioning.dta")
family_transfer <- read_dta("2011/family_transfer.dta")
health_care_and_insurance <- read_dta("2011/health_care_and_insurance.dta")
individual_income <- read_dta("2011/individual_income.dta") 
interviewer_observation <- read_dta("2011/interviewer_observation.dta") 
work_retirement_and_pension <- read_dta("2011/work_retirement_and_pension.dta")
biomarkers <- read_dta("2011/biomarkers.dta")
Blood <- read_dta("2011/Blood_20140429.dta")

# 合并2011年基线数据
Cleandata <- demographic %>%
  left_join(health, by = 'ID') %>%
  left_join(family, by = 'ID') %>%
  left_join(household, by = 'ID') %>%
  left_join(family_transfer, by = 'ID') %>%
  left_join(health_care_and_insurance, by = 'ID') %>%
  left_join(individual_income, by = 'ID') %>%
  left_join(interviewer_observation, by = 'ID') %>%
  left_join(work_retirement_and_pension, by = 'ID') %>%
  left_join(biomarkers, by = 'ID') %>%
  left_join(Blood, by = 'ID')

# 提取并创建所需变量
Cleandata2 <- Cleandata %>%
  mutate(ID = paste0(householdID.x, '0', str_sub(ID, -2, -1))) %>%
  transmute(
    ID,
    Age = coalesce(ifelse(!is.na(ba002_1), 2011 - ba002_1, NA_real_), ba004),
    Gender = ifelse(rgender == 1, "Male", "Female"),
    Education = case_when(
      bd001 %in% 1:7 ~ "Middle school and low",
      bd001 %in% 8:11 ~ "Middle school and above",
      TRUE ~ NA_character_
    ),
    Smoking = ifelse(da059 == 1, "Yes", "No"),
    Drinking = ifelse(da067 == 3, "No", "Yes"),
    Residence = ifelse(bc001 == 1, "Rural Village", "Urban Community"),
    Married_status = ifelse(be001 %in% c(1, 2), "Married", "Single"),
    Sleep = da049,
    height = if_else(qi002 > 220 | qi002 <= 50, NA_real_, qi002),
    weight = if_else(ql002 > 150 | ql002 <= 50, NA_real_, ql002),
    BMI = ifelse(height > 0 & weight > 0, weight / (height/100)^2, NA_real_),
    
    # 医学指标
    Diabetes = da007_3_,
    Hypertension = da007_1_,
    Stroke = da007_8_,
    heart_problems = da007_7_,
    Liver= da007_6_,
    digestive_disease=da007_10_,
    Glu = newglu,
    TG = newtg,
    HDL = newhdl,
    LDL = newldl,
    CRP=newcrp,
    TC=newcho,
    TyG = ifelse(TG > 0 & Glu > 0, log((TG * Glu)/2), NA_real_),
    CHG = ifelse(TC > 0 & Glu > 0 & HDL > 0, 
                 log((TC * Glu) / (2 * HDL)),  # 修正括号和函数调用
                 NA_real_))



# 读取并处理2020年随访数据
follow_up <- read_dta("2020/Health_Status_and_Functioning.dta")
follow_up_diseases <- follow_up %>%
  transmute(
    ID,
    Hypertension_2020 = da003_1_,
    Stroke_2020 = da003_8_,
    heart_problems_2020 = da003_7_,
    Diabetes_2020 = da003_3_
  )

# 合并随访数据并创建心血管疾病结局变量
Cleandata2 <- Cleandata2 %>%
  left_join(follow_up_diseases, by = "ID") %>%
  mutate(
    cardiometabolic_diseases = case_when(
      Diabetes_2020==2& Hypertension_2020 == 2 & Stroke_2020 == 2 & heart_problems_2020 == 2 ~ 0,  # 无疾病
      Diabetes_2020 == 1 | Hypertension_2020 == 1 | Stroke_2020 == 1 | heart_problems_2020 == 1 ~ 1,  # 有疾病
      TRUE ~ NA_real_
    )
  ) %>%
  # 转换为因子并添加标签
  mutate(
    cardiometabolic_diseases = factor(
      cardiometabolic_diseases ,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

# 添加变量标签
attr(Cleandata2$cardiometabolic_diseases, "label") <- 
  "2020 cardiometabolic diseases"

# 筛选分析人群
included_data <- Cleandata2 %>%
  filter(Age >= 45) %>%
  filter(Hypertension == 2) %>%     # 基线无高血压
  filter(Stroke == 2) %>%           # 基线无中风
  filter(heart_problems == 2) %>%   # 基线无心脏问题
  filter(Diabetes == 2) %>%         #基线无糖尿病
  filter(!is.na(cardiometabolic_diseases))


library(dplyr)

# 计算每个步骤的例数变化
step_counts <- list()
step_counts[["初始"]] <- nrow(Cleandata2)

included_data <- Cleandata2 %>%
  {step_counts[["Age>=45"]] <<- nrow(filter(., Age >= 45)); filter(., Age >= 45)} %>%
  {step_counts[["Hypertension==2"]] <<- nrow(filter(., Hypertension == 2)); filter(., Hypertension == 2)} %>%
  {step_counts[["Stroke==2"]] <<- nrow(filter(., Stroke == 2)); filter(., Stroke == 2)} %>%
  {step_counts[["heart_problems==2"]] <<- nrow(filter(., heart_problems == 2)); filter(., heart_problems == 2)} %>%
  {step_counts[["Diabetes==2"]] <<- nrow(filter(., Diabetes == 2)); filter(., Diabetes == 2)} %>%
  {step_counts[["!is.na(cardiometabolic_diseases)"]] <<- nrow(filter(., !is.na(cardiometabolic_diseases))); filter(., !is.na(cardiometabolic_diseases))}

# 计算减少的例数
reduction_df <- data.frame(
  Step = names(step_counts),
  Remaining = unlist(step_counts),
  Reduced = c(0, diff(-unlist(step_counts)))
)

print(reduction_df)


# 输出处理后的数据
write.csv(included_data, "CHG-CMD.csv", row.names = FALSE)

