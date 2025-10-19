# 清理环境并设置工作目录
rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/CHARLS")

# 1. 加载必要的包 ----
library(tidyverse)
library(jstable)
library(forestploter)
library(data.table)
library(grid)

# 2. 高效数据加载与预处理 ----
# 使用fread替代read_csv以提高读取速度
data <- fread("CHG-CMD_imputed.csv", na.strings = c("NA", "")) %>%
  mutate(
    cardiometabolic_diseases = as.integer(cardiometabolic_diseases == "Yes"),
    Age_group = factor(ifelse(Age < 60, "<60 years", "≥60 years"),
                       levels = c("<60 years", "≥60 years")),
    Sleep_group = factor(ifelse(Sleep < 7, "<7h/d", "≥7h/d"),
                         levels = c("<7h/d", "≥7h/d")),
    BMI_group = factor(
      ifelse(BMI <= 25, "≤25",
             ifelse(BMI > 25 & BMI <= 30, "25-30", ">30")),
      levels = c("≤25", "25-30", ">30")),
    across(c(Gender, Education, Smoking, Drinking, Residence, Married_status), as.factor)
  )

# 3. 亚组分析优化 ----
subgroups <- c("Age_group", "Gender", "Education", "Smoking", "Drinking", 
               "Residence", "Married_status", "Sleep_group","BMI_group")

# 使用更稳健的GLM分析
res <- TableSubgroupMultiGLM(
  formula = cardiometabolic_diseases ~CHG,
  var_subgroups = subgroups,
  data = data,
  family = "poisson"
)

# 4. 优化森林图数据处理 ----
# 检查结果结构并处理
if (!"Variable" %in% colnames(res)) {
  stop("结果中缺少'Variable'列，请检查TableSubgroupMultiGLM的输出结构")
}

# 确保结果列是数值型
dt <- res %>%
  slice(-1) %>%  # 删除第一行总结果
  rename(Subgroups = Variable) %>%
  mutate(
    across(c(RR, Lower, Upper, `P for interaction`), ~ as.numeric(as.character(.x))),
    # 格式化效应值和置信区间
    `β(95%CI)` = ifelse(
      is.na(RR) | is.na(Lower) | is.na(Upper), 
      "", 
      sprintf("%.2f (%.2f, %.2f)", RR, Lower, Upper)
    ),
    # 添加缩进以增强可读性
    Subgroups = ifelse(is.na(RR), Subgroups, paste0("   ", Subgroups)),
    # 创建空白列用于绘图
    ` ` = paste(rep(" ", 20), collapse = ""),
    # 安全地格式化 P 值
    `P for interaction` = case_when(
      is.na(`P for interaction`) ~ "",
      `P for interaction` < 0.001 ~ "<0.001",
      `P for interaction` < 0.01 ~ sprintf("%.3f", `P for interaction`),
      `P for interaction` < 0.05 ~ sprintf("%.3f", `P for interaction`),
      TRUE ~ sprintf("%.2f", `P for interaction`)
    ),
    # 添加显著性标记
    Significance = case_when(
      is.na(`P for interaction`) ~ "",
      `P for interaction` < 0.001 ~ "***",
      `P for interaction` < 0.01 ~ "**",
      `P for interaction` < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 5. 优化森林图可视化 ----
# 自定义主题 - 优化布局
tm <- forest_theme(
  base_size = 10,
  ci_pch = 18,
  ci_col = "#2E86AB",
  ci_fill = "#A23B72",
  ci_lty = 1,
  ci_lwd = 1.8,
  ci_alpha = 0.8,
  ci_Theight = 0.3,
  core = list(
    bg_params = list(fill = c("#F8F9FA", "#E9ECEF")),
    fg_params = list(hjust = 0, x = 0.05)
  ),
  refline_gp = gpar(col = "#F24236", lty = 1, lwd = 1.8),
  footnote_cex = 0.8,
  footnote_fontface = "italic",
  footnote_col = "#495057",
  vertline_gp = gpar(lty = 2, col = "gray70")
)

# 创建森林图 - 直接在forest函数中添加标题和脚注
forest_plot <- forest(
  dt[, c("Subgroups", "β(95%CI)", " ", "P for interaction")],
  est = dt$RR,
  lower = dt$Lower,
  upper = dt$Upper,
  ci_column = 3,
  xlim = c(0.5, 3.0),
  ticks_at = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
  ref_line = 1,
  arrow_lab = c("Lower Risk", "Higher Risk"),
  theme = tm,
  colgap = unit(6, "mm"),
  graphwidth = unit(10, "cm"),
  lineheight = "auto",
  clip = c(0.5, 3.0),
  col_names = c("Subgroup", "β (95% CI)", "", "P for interaction"),
  # 直接在forest函数中添加标题和脚注
  title = "Association between CHG and Cardiometabolic Diseases by Subgroups")

# 6. 动态计算高度并保存
n_subgroups <- nrow(dt)
plot_height <- max(8, n_subgroups * 0.5 + 3)

# 保存为高分辨率图像
ggsave("Forest_Plot_Full.png", 
       forest_plot,
       width = 18,
       height = plot_height,
       dpi = 300,
       bg = "white")

ggsave("Forest_Plot_Full.pdf", 
       forest_plot,
       width = 18,
       height = plot_height,
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

# 7. 创建简化版本用于演示
if(n_subgroups > 15) {
  # 只保留显著结果和几个代表性亚组
  dt_simple <- dt %>%
    filter(`P for interaction` < 0.1 | row_number() <= 5 | row_number() > n() - 5)
  
  simple_plot <- forest(
    dt_simple[, c("Subgroups", "β(95%CI)", " ", "P for interaction")],
    est = dt_simple$RR,
    lower = dt_simple$Lower,
    upper = dt_simple$Upper,
    ci_column = 3,
    xlim = c(0.5, 3.0),
    ref_line = 1,
    theme = tm,
    title = "Simplified Forest Plot: Association between CHG and Cardiometabolic Diseases"
  )
  
  ggsave("Forest_Plot_Simple.png", simple_plot, width = 16, height = 10, dpi = 300)
}

# 8. 添加结果摘要
cat("亚组分析结果摘要:\n")
cat(sprintf("总亚组数: %d\n", n_subgroups))
cat(sprintf("显著交互作用(P < 0.05): %d\n", sum(dt$`P for interaction` < 0.05, na.rm = TRUE)))
cat(sprintf("高度显著交互作用(P < 0.01): %d\n", sum(dt$`P for interaction` < 0.01, na.rm = TRUE)))

# 9. 显示图形
print(forest_plot)

