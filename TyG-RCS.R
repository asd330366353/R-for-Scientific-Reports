# 加载必要的包
library(ggplot2)
library(rms)

# 创建只包含模型中实际使用变量的数据子集，避免常量变量警告
model_vars <- c("cardiometabolic_diseases", "TyG", "Age", "Gender", "Education", 
                "Smoking", "Drinking", "Residence", "Married_status", "Sleep", "BMI")
model_data <- included_data1[, model_vars]

# 设置数据分布 - 只使用模型中的变量
dd <- datadist(model_data)
options(datadist = "dd")

# 构建逻辑回归模型
fit <- lrm(cardiometabolic_diseases ~ rcs(TyG, 4) + 
             Age + Gender + Education + Smoking + 
             Drinking + Residence + Married_status + 
             Sleep + BMI,
           data = model_data, x = TRUE, y = TRUE)

# 模型检验
anova_results <- anova(fit)
print(anova_results)

# 提取P值
TyG_p_overall <- anova_results["TyG", "P"]
TyG_p_nonlinear <- anova_results[" Nonlinear", "P"]

# 格式化P值
format_p_value <- function(p) {
  if (p < 0.001) {
    return("<0.001")
  } else {
    return(sprintf("%.3f", p))
  }
}

p_overall_text <- format_p_value(TyG_p_overall)
p_nonlinear_text <- format_p_value(TyG_p_nonlinear)

# 获取RCS节点位置
knots_positions <- attr(rcs(model_data$TyG, 4), "parms")
cat("RCS节点位置 (knots):", knots_positions, "\n")

# 生成预测值
OR <- Predict(fit, TyG, ref.zero = TRUE, fun = exp)

# 确定参考值位置 (通常是最小值或中位数)
reference_value <- min(OR$TyG)
reference_index <- which.min(OR$TyG)

# 创建P值标签数据框
p_label_data <- data.frame(
  x = min(OR$TyG) + (max(OR$TyG) - min(OR$TyG)) * 0.05,
  y = max(OR$upper) * 0.95,
  label = paste0("P-overall = ", p_overall_text, "\nP-non-linear = ", p_nonlinear_text)
)

# 创建图形对象并赋值给变量p
p <- ggplot(OR, aes(x = TyG, y = yhat)) + 
  # 添加参考线
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 1.2) +
  # 添加置信区间阴影 - 使用更明显的颜色
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.4, fill = "#FF6B6B", color = NA) +
  # 添加趋势线 - 加粗
  geom_line(linewidth = 2, color = "#E41A1C") +
  # 标注参考值点
  geom_point(data = OR[reference_index, ], 
             aes(x = TyG, y = yhat), 
             color = "blue", size = 3, shape = 19) +
  # 添加参考值标注
  annotate("text", 
           x = reference_value, 
           y = OR$yhat[reference_index],
           label = "Reference\n(OR = 1.0)",
           vjust = -1, hjust = 0.5, size = 4.5,
           fontface = "bold", color = "blue") +
  # 添加节点位置标记（可选）
  geom_vline(xintercept = knots_positions, 
             linetype = "dotted", color = "gray30", alpha = 0.7) +
  # 添加P值标签 - 使用geom_text而不是annotate
  geom_text(data = p_label_data,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 6,
            fontface = "bold", color = "black") +
  # 添加节点信息标注（可选）
  annotate("text",
           x = min(OR$TyG) + (max(OR$TyG) - min(OR$TyG)) * 0.65,
           y = max(OR$upper) * 0.7,
           label = paste("Knots at:", paste(round(knots_positions, 2), collapse = ", ")),
           hjust = 0, vjust = 1, size = 4.5,
           fontface = "plain", color = "gray40") +
  labs(
    x = "TyG Index", 
    y = "Odds Ratio (95% CI)",
    title = "Association Between TyG Index and Cardiometabolic Diseases",
    subtitle = "Restricted Cubic Spline Analysis with 4 Knots"
  ) +
  theme_minimal(base_size = 16) +  # 基础字号从13增加到16
  theme(
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),  # 标题字号加大
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 14, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    # 增大坐标轴数字的字号并加粗
    axis.text = element_text(size = 14, face = "bold"),  # 字号加大并加粗
    axis.title = element_text(size = 18, face = "bold"),  # 坐标轴标题字号加大并加粗
    # 图例相关设置（如果有的话）
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 8),
    trans = "log",
    limits = c(min(OR$lower), max(OR$upper)),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x),
    breaks = scales::pretty_breaks(n = 8)  # 增加x轴刻度数量
  ) +
  coord_cartesian(clip = "off")

# 显示图形
print(p)

# 保存为高分辨率PDF
ggsave("TyG_RCS_Plot.pdf", plot = p, 
       width = 12, height = 12, device = cairo_pdf)  # 增加图像尺寸

# 保存为高分辨率TIFF
ggsave("TyG_RCS_Plot.tiff", plot = p, 
       width = 12, height = 12, dpi = 600,  # 提高DPI到600
       compression = "lzw")

# 输出节点位置到文本文件
writeLines(paste("RCS节点位置 (knots):", paste(round(knots_positions, 4), collapse = ", ")), 
           "TyG_RCS_knots.txt")

# 打印节点位置到控制台
cat("RCS节点位置 (knots):", paste(round(knots_positions, 4), collapse = ", "), "\n")

