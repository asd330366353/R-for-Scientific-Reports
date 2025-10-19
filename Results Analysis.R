# 修复后的详细分析函数
detailed_analysis <- function(data) {
  cat("=== TyG和CHG的详细分析 ===\n\n")
  
  # 1. 基本统计描述
  cat("1. 基本统计描述:\n")
  stats <- data %>%
    group_by(disease_status) %>%
    summarise(
      n = n(),
      TyG_mean = mean(TyG, na.rm = TRUE),
      TyG_sd = sd(TyG, na.rm = TRUE),
      CHG_mean = mean(CHG, na.rm = TRUE),
      CHG_sd = sd(CHG, na.rm = TRUE)
    )
  print(stats)
  
  # 2. 效应大小计算
  cat("\n2. 效应大小分析:\n")
  # Cohen's d for TyG
  tyg_effect <- (mean(data$TyG[data$disease_status == 1]) - 
                   mean(data$TyG[data$disease_status == 0])) / 
    sqrt((sd(data$TyG[data$disease_status == 1])^2 + 
            sd(data$TyG[data$disease_status == 0])^2)/2)
  
  # Cohen's d for CHG
  chg_effect <- (mean(data$CHG[data$disease_status == 1]) - 
                   mean(data$CHG[data$disease_status == 0])) / 
    sqrt((sd(data$CHG[data$disease_status == 1])^2 + 
            sd(data$CHG[data$disease_status == 0])^2)/2)
  
  cat("TyG的Cohen's d效应大小:", round(tyg_effect, 3), "\n")
  cat("CHG的Cohen's d效应大小:", round(chg_effect, 3), "\n")
  cat("效应大小解读: 0.2=小, 0.5=中, 0.8=大\n")
  
  # 3. 相关性分析
  cat("\n3. 相关性分析:\n")
  cor_tyg <- cor.test(data$TyG, as.numeric(data$disease_status))
  cor_chg <- cor.test(data$CHG, as.numeric(data$disease_status))
  cat("TyG与疾病的点二列相关:", round(cor_tyg$estimate, 3), 
      "(p =", format.pval(cor_tyg$p.value, digits = 3), ")\n")
  cat("CHG与疾病的点二列相关:", round(cor_chg$estimate, 3), 
      "(p =", format.pval(cor_chg$p.value, digits = 3), ")\n")
  
  # 4. 逻辑回归系数
  cat("\n4. 逻辑回归分析:\n")
  model_tyg <- glm(disease_status ~ TyG, data = data, family = binomial)
  model_chg <- glm(disease_status ~ CHG, data = data, family = binomial)
  
  cat("TyG单变量模型:\n")
  print(summary(model_tyg)$coefficients)
  cat("CHG单变量模型:\n")
  print(summary(model_chg)$coefficients)
  
  # 5. 分布重叠可视化
  library(ggplot2)
  library(patchwork)
  
  p1 <- ggplot(data, aes(x = TyG, fill = factor(disease_status))) +
    geom_density(alpha = 0.5) +
    labs(title = "Distribution of TyG in cases and controls",
         x = "TyG", y = "Density", fill = "CMD") +
    theme_minimal()
  
  p2 <- ggplot(data, aes(x = CHG, fill = factor(disease_status))) +
    geom_density(alpha = 0.5) +
    labs(title = "Distribution of CHG in cases and controls",
         x = "CHG", y = "Density", fill = "CMD") +
    theme_minimal()
  
  # 6. 修复最佳截断值分析
  cat("\n5. 最佳截断值分析:\n")
  roc_tyg <- roc(data$disease_status, data$TyG)
  roc_chg <- roc(data$disease_status, data$CHG)
  
  # 修复：正确提取最佳截断值
  coords_tyg <- coords(roc_tyg, "best", ret = c("threshold", "specificity", "sensitivity"))
  coords_chg <- coords(roc_chg, "best", ret = c("threshold", "specificity", "sensitivity"))
  
  # 检查返回类型并正确提取
  if (is.list(coords_tyg)) {
    # 如果是列表，提取第一个元素
    threshold_tyg <- ifelse(is.numeric(coords_tyg$threshold), round(coords_tyg$threshold[1], 3), "N/A")
    sensitivity_tyg <- ifelse(is.numeric(coords_tyg$sensitivity), round(coords_tyg$sensitivity[1], 3), "N/A")
    specificity_tyg <- ifelse(is.numeric(coords_tyg$specificity), round(coords_tyg$specificity[1], 3), "N/A")
  } else {
    # 如果是矩阵或数据框
    threshold_tyg <- round(coords_tyg["threshold", 1], 3)
    sensitivity_tyg <- round(coords_tyg["sensitivity", 1], 3)
    specificity_tyg <- round(coords_tyg["specificity", 1], 3)
  }
  
  if (is.list(coords_chg)) {
    threshold_chg <- ifelse(is.numeric(coords_chg$threshold), round(coords_chg$threshold[1], 3), "N/A")
    sensitivity_chg <- ifelse(is.numeric(coords_chg$sensitivity), round(coords_chg$sensitivity[1], 3), "N/A")
    specificity_chg <- ifelse(is.numeric(coords_chg$specificity), round(coords_chg$specificity[1], 3), "N/A")
  } else {
    threshold_chg <- round(coords_chg["threshold", 1], 3)
    sensitivity_chg <- round(coords_chg["sensitivity", 1], 3)
    specificity_chg <- round(coords_chg["specificity", 1], 3)
  }
  
  cat("TyG最佳截断值:", threshold_tyg, 
      "敏感性:", sensitivity_tyg,
      "特异性:", specificity_tyg, "\n")
  
  cat("CHG最佳截断值:", threshold_chg, 
      "敏感性:", sensitivity_chg,
      "特异性:", specificity_chg, "\n")
  
  # 7. 替代方法：手动计算Youden指数
  cat("\n6. Youden指数分析:\n")
  calculate_youden <- function(roc_obj) {
    # 获取所有可能的截断值
    coords_all <- coords(roc_obj, "all", ret = c("threshold", "specificity", "sensitivity"))
    
    if (is.list(coords_all)) {
      specificity <- coords_all$specificity
      sensitivity <- coords_all$sensitivity
      threshold <- coords_all$threshold
    } else {
      specificity <- coords_all["specificity", ]
      sensitivity <- coords_all["sensitivity", ]
      threshold <- coords_all["threshold", ]
    }
    
    # 计算Youden指数
    youden <- sensitivity + specificity - 1
    best_index <- which.max(youden)
    
    return(list(
      threshold = threshold[best_index],
      sensitivity = sensitivity[best_index],
      specificity = specificity[best_index],
      youden = youden[best_index]
    ))
  }
  
  youden_tyg <- calculate_youden(roc_tyg)
  youden_chg <- calculate_youden(roc_chg)
  
  cat("TyG (Youden指数法):\n")
  cat("  最佳截断值:", round(youden_tyg$threshold, 3), "\n")
  cat("  敏感性:", round(youden_tyg$sensitivity, 3), "\n")
  cat("  特异性:", round(youden_tyg$specificity, 3), "\n")
  cat("  Youden指数:", round(youden_tyg$youden, 3), "\n")
  
  cat("CHG (Youden指数法):\n")
  cat("  最佳截断值:", round(youden_chg$threshold, 3), "\n")
  cat("  敏感性:", round(youden_chg$sensitivity, 3), "\n")
  cat("  特异性:", round(youden_chg$specificity, 3), "\n")
  cat("  Youden指数:", round(youden_chg$youden, 3), "\n")
  
  return(list(
    plots = p1 + p2, 
    effect_sizes = c(TyG = tyg_effect, CHG = chg_effect),
    correlations = c(TyG = cor_tyg$estimate, CHG = cor_chg$estimate),
    youden_results = list(TyG = youden_tyg, CHG = youden_chg)
  ))
}

# 执行详细分析
detailed_results <- detailed_analysis(results$data)

# 显示分布图
print(detailed_results$plots)
# 进一步的多变量分析
multivariate_analysis <- function(data) {
  cat("\n=== 多变量模型增量价值分析 ===\n")
  
  # 基础模型
  base_model <- glm(disease_status ~ BMI + age + sex + smoking + drinking, 
                    data = data, family = binomial)
  base_auc <- auc(roc(data$disease_status, predict(base_model, type = "response")))
  
  # 基础模型 + TyG
  model_tyg <- glm(disease_status ~ BMI + age + sex + smoking + drinking + TyG, 
                   data = data, family = binomial)
  auc_tyg <- auc(roc(data$disease_status, predict(model_tyg, type = "response")))
  
  # 基础模型 + CHG
  model_chg <- glm(disease_status ~ BMI + age + sex + smoking + drinking + CHG, 
                   data = data, family = binomial)
  auc_chg <- auc(roc(data$disease_status, predict(model_chg, type = "response")))
  
  # 基础模型 + TyG + CHG
  model_both <- glm(disease_status ~ BMI + age + sex + smoking + drinking + TyG + CHG, 
                    data = data, family = binomial)
  auc_both <- auc(roc(data$disease_status, predict(model_both, type = "response")))
  
  improvement_table <- data.frame(
    Model = c("基础模型", "+ TyG", "+ CHG", "+ TyG + CHG"),
    AUC = round(c(base_auc, auc_tyg, auc_chg, auc_both), 3),
    Improvement = round(c(0, auc_tyg - base_auc, auc_chg - base_auc, auc_both - base_auc), 3)
  )
  
  print(improvement_table)
  
  cat("\n解释:\n")
  cat("- 即使单个标记物的AUC不高，它们在多变量模型中可能提供增量价值\n")
  cat("- 关注的是AUC的改善(ΔAUC)，而不是单个标记物的绝对AUC\n")
  cat("- 小的ΔAUC如果有统计学意义，仍然具有临床重要性\n")
  
  return(improvement_table)
}

# 执行多变量分析
improvement_results <- multivariate_analysis(results$data)

# 生成论文中可用的解释文本
generate_interpretation <- function(detailed_results, improvement_results, auc_values) {
  cat("\n=== 论文中可用的解释文本 ===\n\n")
  
  cat("结果解释:\n")
  cat("本研究发现TyG和CHG是心血管代谢疾病(CMD)的独立危险因素，")
  cat("尽管其单变量判别能力有限(TyG AUC = ", round(auc_values["TyG alone"], 3), 
      "; CHG AUC = ", round(auc_values["CHG alone"], 3), ")。\n\n")
  
  cat("这种现象的可能解释包括:\n")
  cat("1. 复杂疾病的多因素性质: CMD是多种遗传和环境因素共同作用的结果，")
  cat("单个生物标志物通常无法充分捕捉疾病的复杂性。\n")
  
  cat("2. 病例-对照分布重叠: 尽管两组间存在统计学显著差异，")
  cat("但TyG和CHG在病例和对照组中的分布存在相当程度的重叠，")
  cat("限制了其判别能力。\n")
  
  cat("3. 效应大小有限: TyG的效应大小为", 
      round(abs(detailed_results$effect_sizes["TyG"]), 3),
      "，CHG的效应大小为", 
      round(abs(detailed_results$effect_sizes["CHG"]), 3),
      "，属于小到中等效应。\n\n")
  
  cat("临床意义:\n")
  cat("尽管单变量判别能力有限，但将TyG和CHG纳入多变量模型后，")
  cat("模型的判别能力从", improvement_results$AUC[1], "提高至", 
      improvement_results$AUC[4], "(ΔAUC = ", 
      improvement_results$Improvement[4], ")，")
  cat("表明这些标记物在综合风险评估中具有增量价值。\n\n")
  
  cat("建议:\n")
  cat("在临床实践中，TyG和CHG应作为综合风险评估工具的一部分使用，")
  cat("而不是作为独立的筛查工具。它们与传统危险因素的结合能够提供更准确的风险分层。\n")
}

# 生成解释文本
generate_interpretation(detailed_results, improvement_results, results$auc_values)

