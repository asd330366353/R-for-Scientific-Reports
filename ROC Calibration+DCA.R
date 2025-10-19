# 加载必要的包
library(pROC)
library(ggplot2)
library(randomForest)
library(dplyr)
library(dcurves)
library(purrr)
library(tidyr)

# 设置主题和颜色
theme_set(theme_bw())
model_colors <- c("Base Model" = "#1f77b4", "With CHG+TyG" = "#ff7f0e", 
                  "Complete Model" = "#2ca02c", "Random Forest" = "#d62728",
                  "CHG alone" = "#9467bd", "TyG alone" = "#8c564b")

# 读取数据并预处理
load_and_preprocess_data <- function(file_path) {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  data_clean <- data %>%
    mutate(
      disease_status = ifelse(cardiometabolic_diseases == "Yes", 1, 0),
      sex = ifelse(Gender == "Male", 1, 0),
      smoking = ifelse(Smoking == "Yes", 1, 0),
      drinking = ifelse(Drinking == "Yes", 1, 0),
      residence_urban = ifelse(Residence == "Urban Community", 1, 0),
      married = ifelse(Married_status == "Married", 1, 0),
      education_high = ifelse(Education == "Middle school and above", 1, 0),
      across(c(age = Age, BMI, CHG, TyG, sleep = Sleep), as.numeric)
    ) %>%
    select(disease_status, BMI, age, sex, education_high, smoking, drinking, 
           residence_urban, married, sleep, CHG, TyG) %>%
    na.omit()
  
  # 数据质量检查
  cat("数据概览:\n")
  cat("样本数:", nrow(data_clean), "\n")
  cat("疾病患病率:", round(mean(data_clean$disease_status) * 100, 1), "%\n")
  cat("缺失值: 无\n\n")
  
  return(data_clean)
}

# 构建多个模型
build_models <- function(data) {
  models <- list()
  
  # 定义模型公式
  formulas <- list(
    base = disease_status ~ BMI + age + sex + smoking + drinking,
    extended = disease_status ~ BMI + age + sex + smoking + drinking + CHG + TyG,
    complete = disease_status ~ BMI + age + sex + education_high + smoking + 
      drinking + residence_urban + married + sleep + CHG + TyG
  )
  
  # 训练所有模型
  models <- map(formulas, ~glm(.x, data = data, family = binomial))
  
  # 添加预测概率到数据中
  for (name in names(models)) {
    prob_col <- paste0(name, "_prob")
    data[[prob_col]] <- predict(models[[name]], type = "response")
  }
  
  return(list(models = models, data = data))
}

# 计算ROC性能
calculate_roc_performance <- function(data) {
  roc_curves <- list(
    "CHG alone" = roc(data$disease_status, data$CHG),
    "TyG alone" = roc(data$disease_status, data$TyG),
    "Base Model" = roc(data$disease_status, data$base_prob),
    "With CHG+TyG" = roc(data$disease_status, data$extended_prob),
    "Complete Model" = roc(data$disease_status, data$complete_prob)
  )
  
  # 计算AUC
  auc_values <- map_dbl(roc_curves, auc)
  
  return(list(roc_curves = roc_curves, auc_values = auc_values))
}

# 修复的NRI计算函数
calculate_nri <- function(prob1, prob2, outcome, cutpoints = c(0.2, 0.8)) {
  # 确保概率值在有效范围内
  prob1 <- pmin(pmax(prob1, 1e-10), 1-1e-10)
  prob2 <- pmin(pmax(prob2, 1e-10), 1-1e-10)
  
  categorize <- function(prob) {
    # 使用明确的边界条件
    category <- rep("Intermediate", length(prob))
    category[prob <= cutpoints[1]] <- "Low"
    category[prob > cutpoints[2]] <- "High"
    return(factor(category, levels = c("Low", "Intermediate", "High")))
  }
  
  cat1 <- categorize(prob1)
  cat2 <- categorize(prob2)
  
  events <- outcome == 1
  nonevents <- outcome == 0
  
  # 检查是否有足够的事件和非事件
  if (sum(events) == 0 || sum(nonevents) == 0) {
    warning("没有足够的事件或非事件来计算NRI")
    return(data.frame(
      NRI_events = NA,
      NRI_nonevents = NA,
      NRI_total = NA
    ))
  }
  
  event_up <- mean(cat2[events] > cat1[events], na.rm = TRUE)
  event_down <- mean(cat2[events] < cat1[events], na.rm = TRUE)
  nonevent_down <- mean(cat2[nonevents] < cat1[nonevents], na.rm = TRUE)
  nonevent_up <- mean(cat2[nonevents] > cat1[nonevents], na.rm = TRUE)
  
  nri_events <- event_up - event_down
  nri_nonevents <- nonevent_down - nonevent_up
  nri_total <- nri_events + nri_nonevents
  
  return(data.frame(
    NRI_events = nri_events,
    NRI_nonevents = nri_nonevents,
    NRI_total = nri_total
  ))
}

# 替代的简化NRI计算（如果上面的仍然有问题）
calculate_simple_nri <- function(prob1, prob2, outcome) {
  # 使用连续NRI（不需要分类）
  events <- outcome == 1
  nonevents <- outcome == 0
  
  # 事件组的NRI
  event_nri <- mean(prob2[events] > prob1[events]) - mean(prob2[events] < prob1[events])
  
  # 非事件组的NRI
  nonevent_nri <- mean(prob2[nonevents] < prob1[nonevents]) - mean(prob2[nonevents] > prob1[nonevents])
  
  total_nri <- event_nri + nonevent_nri
  
  return(data.frame(
    NRI_events = event_nri,
    NRI_nonevents = nonevent_nri,
    NRI_total = total_nri
  ))
}

# 使用dcurves进行决策曲线分析
perform_dca <- function(data) {
  # 创建DCA数据
  dca_data <- data %>%
    select(disease_status, base_prob, extended_prob, complete_prob)
  
  # 使用dcurves进行DCA
  dca_results <- dca(
    disease_status ~ base_prob + extended_prob + complete_prob,
    data = dca_data,
    thresholds = seq(0.01, 0.5, by = 0.01)
  )
  
  return(dca_results)
}

# 训练随机森林模型
train_random_forest <- function(data) {
  rf_data <- data %>%
    select(disease_status, BMI, age, sex, education_high, smoking, drinking, 
           residence_urban, married, sleep, CHG, TyG) %>%
    mutate(disease_status = as.factor(disease_status))
  
  set.seed(123)
  rf_model <- randomForest(disease_status ~ ., data = rf_data, 
                           ntree = 300,
                           importance = TRUE)
  
  # 获取预测概率
  rf_prob <- predict(rf_model, type = "prob")[,2]
  
  return(list(model = rf_model, probabilities = rf_prob))
}

# 修复的ROC曲线绘制函数
create_roc_plot <- function(roc_curves, auc_values) {
  # 提取ROC曲线数据
  roc_list <- list()
  for (i in 1:length(roc_curves)) {
    model_name <- names(roc_curves)[i]
    roc_obj <- roc_curves[[i]]
    roc_list[[i]] <- data.frame(
      sensitivity = roc_obj$sensitivities,
      specificity = roc_obj$specificities,
      model = model_name
    )
  }
  
  roc_data <- bind_rows(roc_list) %>%
    mutate(`1_specificity` = 1 - specificity)
  
  # 创建AUC标签
  auc_labels <- data.frame(
    model = names(auc_values),
    label = paste0(names(auc_values), ": AUC = ", round(auc_values, 3)),
    x = 0.55,
    y = seq(0.35, 0.35 - (length(auc_values)-1)*0.08, length.out = length(auc_values))
  )
  
  # 绘制
  p <- ggplot(roc_data, aes(x = `1_specificity`, y = sensitivity, color = model)) +
    geom_line(size = 1.1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
    geom_label(data = auc_labels, 
               aes(x = x, y = y, label = label, color = model),
               show.legend = FALSE, size = 3.5, hjust = 0) +
    labs(title = "Model ROC Curve Comparison",
         subtitle = "Univariate Models vs. Multivariate Models vs. Machine Learning Models",
         x = "1 - specificity",
         y = "Sensitivity") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "none") +
    scale_color_manual(values = model_colors) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1)
  
  return(p)
}

create_importance_plot <- function(model, title = "Random Forest Feature Importance") {
  imp <- importance(model)
  imp_df <- data.frame(
    Variable = rownames(imp),
    Importance = imp[, "MeanDecreaseGini"]
  ) %>%
    arrange(desc(Importance))
  
  #改进变量名的可读性
  imp_df$Variable <- recode(imp_df$Variable,
                            "age" = "Age",
                            "BMI" = "BMI",
                            "CHG" = "CHG",
                            "TyG" = "TyG",
                            "sex" = "Gender",
                            "smoking" = "Smoking",
                            "drinking" = "Drinking",
                            "education_high" = "High-education",
                            "residence_urban" = "Urban",
                            "married" = "Married",
                            "sleep" = "Sleep")
  
  ggplot(imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = title,
         x = "Variable",
         y = "Importance (Mean Decrease Gini)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# 创建改进的DCA图
create_dca_plot <- function(dca_results) {
  dca_plot <- plot(dca_results) +
    labs(title = "Decision Curve Analysis",
         subtitle = "Evaluating the clinical utility of the model under different decision thresholds",
         x = " Threshold Probability",
         y = "Net Benefit") +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12)) +
    scale_color_brewer(palette = "Set1")
  
  return(dca_plot)
}

# 决策曲线解读函数
interpret_dca <- function(dca_results, auc_values, data) {
  cat("\n=== 决策曲线分析解读 ===\n")
  cat("决策曲线分析评估模型在不同决策阈值下的临床净收益。\n\n")
  
  cat("关键解读要点:\n")
  cat("1. 横坐标(决策阈值): 表示临床医生或患者认为需要采取干预措施的概率阈值\n")
  cat("2. 纵坐标(净收益): 综合考虑真阳性和假阳性的临床价值\n")
  cat("3. 参考线:\n")
  cat("   - '全部治疗'(灰色实线): 假设所有患者都接受治疗的净收益\n")
  cat("   - '全部不治疗'(灰色虚线): 假设所有患者都不接受治疗的净收益\n\n")
  
  cat("解读指南:\n")
  cat("- 在特定阈值范围内，曲线位置越高表示临床净收益越大\n")
  cat("- 模型曲线高于'全部治疗'和'全部不治疗'线，说明在该阈值范围内使用该模型有临床价值\n")
  
  # 计算患病率
  prevalence <- mean(data$disease_status)
  cat("- 疾病患病率:", round(prevalence * 100, 1), "% - 通常关注患病率附近的阈值范围\n\n")
  
  # 基于AUC提供额外见解
  best_model <- names(which.max(auc_values))
  best_auc <- round(max(auc_values), 3)
  cat("基于分析结果:\n")
  cat("- 最佳性能模型:", best_model, "(AUC =", best_auc, ")\n")
  cat("- 决策曲线应显示该模型在较宽阈值范围内提供最高的净收益\n")
  
  # 提供具体的决策建议
  if (best_auc > 0.8) {
    cat("- 模型具有优秀的判别能力，在临床决策中具有重要价值\n")
  } else if (best_auc > 0.7) {
    cat("- 模型具有良好的判别能力，在临床决策中具有实用价值\n")
  } else {
    cat("- 模型的判别能力一般，需要结合其他临床信息进行决策\n")
  }
  
  # 检查CHG和TyG的增量价值
  base_auc <- auc_values["Base Model"]
  extended_auc <- auc_values["With CHG+TyG"]
  if (extended_auc > base_auc) {
    improvement <- round((extended_auc - base_auc) * 100, 1)
    cat("- 添加CHG和TyG标记物使AUC提高了", improvement, "%，表明这些标记物提供了增量价值\n")
  }
}

# 主分析流程
main_analysis <- function(file_path) {
  cat("开始数据分析...\n")
  
  # 1. 数据加载和预处理
  data_clean <- load_and_preprocess_data(file_path)
  
  # 2. 构建逻辑回归模型
  model_results <- build_models(data_clean)
  data_with_probs <- model_results$data
  models <- model_results$models
  
  # 3. 计算ROC性能
  roc_results <- calculate_roc_performance(data_with_probs)
  
  # 4. 计算NRI - 使用修复的函数
  cat("计算NRI...\n")
  
  # 先检查概率范围
  cat("概率范围检查:\n")
  cat("基础模型概率范围:", round(range(data_with_probs$base_prob), 4), "\n")
  cat("扩展模型概率范围:", round(range(data_with_probs$extended_prob), 4), "\n")
  cat("完整模型概率范围:", round(range(data_with_probs$complete_prob), 4), "\n")
  
  # 尝试使用修复的NRI函数
  nri_results <- list()
  
  tryCatch({
    nri_results$base_vs_extended <- calculate_nri(data_with_probs$base_prob, 
                                                  data_with_probs$extended_prob, 
                                                  data_with_probs$disease_status)
  }, error = function(e) {
    cat("标准NRI计算失败，使用简化版本\n")
    nri_results$base_vs_extended <- calculate_simple_nri(data_with_probs$base_prob, 
                                                         data_with_probs$extended_prob, 
                                                         data_with_probs$disease_status)
  })
  
  tryCatch({
    nri_results$base_vs_complete <- calculate_nri(data_with_probs$base_prob, 
                                                  data_with_probs$complete_prob, 
                                                  data_with_probs$disease_status)
  }, error = function(e) {
    nri_results$base_vs_complete <- calculate_simple_nri(data_with_probs$base_prob, 
                                                         data_with_probs$complete_prob, 
                                                         data_with_probs$disease_status)
  })
  
  tryCatch({
    nri_results$extended_vs_complete <- calculate_nri(data_with_probs$extended_prob, 
                                                      data_with_probs$complete_prob, 
                                                      data_with_probs$disease_status)
  }, error = function(e) {
    nri_results$extended_vs_complete <- calculate_simple_nri(data_with_probs$extended_prob, 
                                                             data_with_probs$complete_prob, 
                                                             data_with_probs$disease_status)
  })
  
  # 5. 决策曲线分析
  dca_results <- perform_dca(data_with_probs)
  
  # 6. 随机森林模型
  rf_results <- train_random_forest(data_with_probs)
  data_with_probs$rf_prob <- rf_results$probabilities
  
  # 7. 最终ROC比较（包含随机森林）
  final_roc_curves <- c(
    roc_results$roc_curves,
    list("Random Forest" = roc(data_with_probs$disease_status, data_with_probs$rf_prob))
  )
  final_auc <- map_dbl(final_roc_curves, auc)
  
  # 输出结果
  cat("\n=== 模型性能总结 ===\n")
  auc_table <- data.frame(
    Model = names(final_auc),
    AUC = round(final_auc, 3)
  )
  print(auc_table)
  
  cat("\n=== NRI 结果 ===\n")
  print(nri_results)
  
  cat("\n=== 最佳逻辑回归模型摘要 (Complete Model) ===\n")
  print(summary(models$complete))
  
  # 生成图形
  roc_plot <- create_roc_plot(final_roc_curves, final_auc)
  importance_plot <- create_importance_plot(rf_results$model)
  dca_plot <- create_dca_plot(dca_results)
  
  # 决策曲线解读
  interpret_dca(dca_results, final_auc, data_with_probs)
  
  return(list(
    data = data_with_probs,
    models = c(models, list(random_forest = rf_results$model)),
    roc_curves = final_roc_curves,
    auc_values = final_auc,
    nri_results = nri_results,
    plots = list(
      roc = roc_plot,
      importance = importance_plot,
      dca = dca_plot
    )
  ))
}

# 执行分析
results <- main_analysis("CHG-CMD_imputed.csv")

# 显示图形
print(results$plots$roc)
print(results$plots$importance)
print(results$plots$dca)

# 保存重要结果
cat("\n分析完成！主要结果:\n")
cat("最佳模型AUC:", round(max(results$auc_values), 3), "\n")
cat("样本数量:", nrow(results$data), "\n")

# 变量重要性比较
cat("\n变量重要性排名 (随机森林):\n")
imp <- importance(results$models$random_forest)
imp_sorted <- imp[order(-imp[, "MeanDecreaseGini"]), ]
print(imp_sorted)

# 保存图形（可选）
#ggsave("ROC_Curves.png", results$plots$roc, width = 10, height = 8, dpi = 300)
#ggsave("Variable_Importance.png", results$plots$importance, width = 10, height = 8, dpi = 300)
#ggsave("DCA_Plot.png", results$plots$dca, width = 10, height = 8, dpi = 300)
