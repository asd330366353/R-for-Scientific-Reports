# Install and load required packages
rm(list = ls())
gc()
setwd("C:/Users/Tang/Desktop/Scientific Reports")
library(PredictABEL)
library(pROC)
library(nricens)
library(survival)
data <- read.csv("CHG-CMD_imputed.csv", na.strings = c("NA", ""))

# Select complete data for model fitting
complete_data <- na.omit(data[, c("cardiometabolic_diseases", "Age", "Gender", "Education", 
                                  "Smoking", "Drinking", "Residence", "Married_status", 
                                  "Sleep", "Glu", "TG", "HDL", "LDL", "CRP", "BMI", "TyG_BMI","BRI","CHG","TyG","TyG_BRI")])
complete_data$cardiometabolic_diseases<-ifelse(complete_data$cardiometabolic_diseases=="No",0,1)

# ==================== Fit three models ====================
# Base model
base_model <- glm(
  cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
    Drinking + Residence + Married_status + Sleep + BMI,
  data = complete_data, 
  family = binomial
)

# Extended model 1 (Base model + TyG)
extended_model_TyG <- glm(
  cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
    Drinking + Residence + Married_status + Sleep + BMI + TyG,
  data = complete_data, 
  family = binomial
)

# Extended model 2 (Base model + CHG)
extended_model_CHG <- glm(
  cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
    Drinking + Residence + Married_status + Sleep + BMI + CHG,
  data = complete_data, 
  family = binomial
)

# Get predicted probabilities
base_predictions <- predict(base_model, type = "response")
extended_TyG_predictions <- predict(extended_model_TyG, type = "response")
extended_CHG_predictions <- predict(extended_model_CHG, type = "response")

# Set risk classification thresholds
cutpoints <- c(0.10, 0.20, 0.30)

# ==================== NRI Calculation ====================
cat("===== NRI Method Description =====\n")
cat("We use categorical NRI method based on the following risk categories:\n")
cat("• Category 1: <10% risk (Low risk)\n")
cat("• Category 2: 10%-20% risk (Low-medium risk)\n") 
cat("• Category 3: 20%-30% risk (Medium-high risk)\n")
cat("• Category 4: >30% risk (High risk)\n\n")

cat("NRI Definition: Net Reclassification Improvement index, evaluating the reclassification ability of new models across predefined risk categories\n")
cat("Calculation: NRI = (Correct upward movement in event group - Incorrect downward movement in event group) + (Correct downward movement in non-event group - Incorrect upward movement in non-event group)\n\n")

# Calculate NRI
nri_result_TyG <- nribin( 
  event = complete_data$cardiometabolic_diseases,
  p.std = base_predictions,
  p.new = extended_TyG_predictions,
  cut = cutpoints,
  niter = 1000
)

nri_result_CHG <- nribin( 
  event = complete_data$cardiometabolic_diseases,
  p.std = base_predictions,
  p.new = extended_CHG_predictions,
  cut = cutpoints,
  niter = 1000
)

# ==================== Safe NRI extraction function ====================
safe_extract_nri <- function(nri_result, model_name) {
  cat(paste("Processing NRI results for", model_name, "model:\n"))
  
  # Safely extract total NRI
  nri_est <- tryCatch({
    if (!is.null(nri_result$nri$est)) as.numeric(nri_result$nri$est) else NA
  }, error = function(e) NA)
  
  nri_lower <- tryCatch({
    if (!is.null(nri_result$nri$lower)) as.numeric(nri_result$nri$lower) else NA
  }, error = function(e) NA)
  
  nri_upper <- tryCatch({
    if (!is.null(nri_result$nri$upper)) as.numeric(nri_result$nri$upper) else NA
  }, error = function(e) NA)
  
  nri_p <- tryCatch({
    if (!is.null(nri_result$nri$p.value)) as.numeric(nri_result$nri$p.value) else NA
  }, error = function(e) NA)
  
  # Safely extract event group NRI
  nri_event_est <- tryCatch({
    if (!is.null(nri_result$nri_event$est)) as.numeric(nri_result$nri_event$est) else NA
  }, error = function(e) NA)
  
  nri_event_lower <- tryCatch({
    if (!is.null(nri_result$nri_event$lower)) as.numeric(nri_result$nri_event$lower) else NA
  }, error = function(e) NA)
  
  nri_event_upper <- tryCatch({
    if (!is.null(nri_result$nri_event$upper)) as.numeric(nri_result$nri_event$upper) else NA
  }, error = function(e) NA)
  
  nri_event_p <- tryCatch({
    if (!is.null(nri_result$nri_event$p.value)) as.numeric(nri_result$nri_event$p.value) else NA
  }, error = function(e) NA)
  
  # Safely extract non-event group NRI
  nri_nonevent_est <- tryCatch({
    if (!is.null(nri_result$nri_nonevent$est)) as.numeric(nri_result$nri_nonevent$est) else NA
  }, error = function(e) NA)
  
  nri_nonevent_lower <- tryCatch({
    if (!is.null(nri_result$nri_nonevent$lower)) as.numeric(nri_result$nri_nonevent$lower) else NA
  }, error = function(e) NA)
  
  nri_nonevent_upper <- tryCatch({
    if (!is.null(nri_result$nri_nonevent$upper)) as.numeric(nri_result$nri_nonevent$upper) else NA
  }, error = function(e) NA)
  
  nri_nonevent_p <- tryCatch({
    if (!is.null(nri_result$nri_nonevent$p.value)) as.numeric(nri_result$nri_nonevent$p.value) else NA
  }, error = function(e) NA)
  
  return(list(
    nri_est = nri_est,
    nri_lower = nri_lower, 
    nri_upper = nri_upper,
    nri_p = nri_p,
    nri_event_est = nri_event_est,
    nri_event_lower = nri_event_lower,
    nri_event_upper = nri_event_upper,
    nri_event_p = nri_event_p,
    nri_nonevent_est = nri_nonevent_est,
    nri_nonevent_lower = nri_nonevent_lower,
    nri_nonevent_upper = nri_nonevent_upper,
    nri_nonevent_p = nri_nonevent_p
  ))
}

# Extract NRI results
nri_TyG <- safe_extract_nri(nri_result_TyG, "Base model + TyG")
nri_CHG <- safe_extract_nri(nri_result_CHG, "Base model + CHG")

# ==================== Safe NRI output function ====================
safe_print_nri <- function(nri_values, model_name) {
  cat(paste(model_name, "\n"))
  
  # Safely output total NRI
  if (!is.na(nri_values$nri_est)) {
    cat("Total NRI:", round(nri_values$nri_est, 4))
    if (!is.na(nri_values$nri_lower) && !is.na(nri_values$nri_upper)) {
      cat(" (95% CI:", round(nri_values$nri_lower, 4), "-", 
          round(nri_values$nri_upper, 4), ")")
    }
    if (!is.na(nri_values$nri_p)) {
      cat(", p =", round(nri_values$nri_p, 4))
    }
    cat("\n")
  } else {
    cat("Total NRI: Cannot calculate\n")
  }
  
  # Safely output event group NRI
  if (!is.na(nri_values$nri_event_est)) {
    cat("Event group NRI:", round(nri_values$nri_event_est, 4))
    if (!is.na(nri_values$nri_event_lower) && !is.na(nri_values$nri_event_upper)) {
      cat(" (95% CI:", round(nri_values$nri_event_lower, 4), "-", 
          round(nri_values$nri_event_upper, 4), ")")
    }
    if (!is.na(nri_values$nri_event_p)) {
      cat(", p =", round(nri_values$nri_event_p, 4))
    }
    cat("\n")
  } else {
    cat("Event group NRI: Cannot calculate\n")
  }
  
  # Safely output non-event group NRI
  if (!is.na(nri_values$nri_nonevent_est)) {
    cat("Non-event group NRI:", round(nri_values$nri_nonevent_est, 4))
    if (!is.na(nri_values$nri_nonevent_lower) && !is.na(nri_values$nri_nonevent_upper)) {
      cat(" (95% CI:", round(nri_values$nri_nonevent_lower, 4), "-", 
          round(nri_values$nri_nonevent_upper, 4), ")")
    }
    if (!is.na(nri_values$nri_nonevent_p)) {
      cat(", p =", round(nri_values$nri_nonevent_p, 4))
    }
    cat("\n")
  } else {
    cat("Non-event group NRI: Cannot calculate\n")
  }
  
  cat("\n")
}

# ==================== IDI Calculation function ====================
calculate_idi <- function(y, p_old, p_new) {
  # Extract case and control groups
  cases <- which(y == 1)
  controls <- which(y == 0)
  
  # Calculate mean difference in predicted probabilities for cases and controls
  idi_case <- mean(p_new[cases]) - mean(p_old[cases])
  idi_control <- mean(p_old[controls]) - mean(p_new[controls])
  
  # IDI = Mean prediction probability increase in cases + Mean prediction probability decrease in controls
  idi <- idi_case + idi_control
  
  return(idi)
}

# Calculate IDI point estimates
idi_estimate_TyG <- calculate_idi(
  y = complete_data$cardiometabolic_diseases,
  p_old = base_predictions,
  p_new = extended_TyG_predictions
)

idi_estimate_CHG <- calculate_idi(
  y = complete_data$cardiometabolic_diseases,
  p_old = base_predictions,
  p_new = extended_CHG_predictions
)

# ==================== Bootstrap IDI Calculation ====================
cat("===== Bootstrap Method Description =====\n")
cat("We use 1000 bootstrap resamples to calculate confidence intervals and p-values:\n")
cat("• In each bootstrap, resample with replacement from original data with same sample size\n")
cat("• Refit base model and extended model on new sample\n")
cat("• Calculate NRI and IDI statistics\n")
cat("• 95% confidence intervals calculated using percentile method\n")
cat("• p-values obtained through two-sided test: p = 2 × min(Pr(estimate ≤ 0), Pr(estimate ≥ 0))\n\n")

bootstrap_idi <- function(y, data, base_formula, extended_formula, n_boot = 1000) {
  set.seed(123)
  idi_boot <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    # Resample data with replacement
    boot_idx <- sample(1:nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]
    
    # Fit models on new bootstrap sample
    base_model_boot <- glm(base_formula, data = boot_data, family = binomial)
    extended_model_boot <- glm(extended_formula, data = boot_data, family = binomial)
    
    # Get predicted probabilities
    base_pred_boot <- predict(base_model_boot, type = "response")
    extended_pred_boot <- predict(extended_model_boot, type = "response")
    
    # Calculate IDI for current bootstrap sample
    idi_boot[i] <- calculate_idi(
      y = boot_data$cardiometabolic_diseases,
      p_old = base_pred_boot,
      p_new = extended_pred_boot
    )
  }
  
  # Calculate 95% confidence interval
  idi_ci <- quantile(idi_boot, probs = c(0.025, 0.975))
  
  # Calculate p-value (two-sided test)
  idi_p_value <- 2 * min(
    mean(idi_boot <= 0),
    mean(idi_boot >= 0)
  )
  
  if (idi_p_value > 1) idi_p_value <- 1
  
  return(list(idi_boot = idi_boot, idi_ci = idi_ci, idi_p_value = idi_p_value))
}

# Define model formulas
base_formula <- cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
  Drinking + Residence + Married_status + Sleep + BMI

extended_formula_TyG <- cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
  Drinking + Residence + Married_status + Sleep + BMI + TyG

extended_formula_CHG <- cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
  Drinking + Residence + Married_status + Sleep + BMI + CHG

# Perform bootstrap for both comparisons
cat("Performing Bootstrap calculation for TyG model...\n")
bootstrap_TyG <- bootstrap_idi(
  y = complete_data$cardiometabolic_diseases,
  data = complete_data,
  base_formula = base_formula,
  extended_formula = extended_formula_TyG
)

cat("Performing Bootstrap calculation for CHG model...\n")
bootstrap_CHG <- bootstrap_idi(
  y = complete_data$cardiometabolic_diseases,
  data = complete_data,
  base_formula = base_formula,
  extended_formula = extended_formula_CHG
)

# ==================== Reclassification table function ====================
create_reclassification_table <- function(y, p_old, p_new, cutpoints, comparison_name) {
  # Classify predicted probabilities into risk categories
  categories_old <- cut(p_old, 
                        breaks = c(-Inf, cutpoints, Inf),
                        labels = c(paste0("<", cutpoints[1]*100, "%"),
                                   paste0(cutpoints[1]*100, "-", cutpoints[2]*100, "%"),
                                   paste0(cutpoints[2]*100, "-", cutpoints[3]*100, "%"),
                                   paste0(">", cutpoints[3]*100, "%")))
  
  categories_new <- cut(p_new, 
                        breaks = c(-Inf, cutpoints, Inf),
                        labels = c(paste0("<", cutpoints[1]*100, "%"),
                                   paste0(cutpoints[1]*100, "-", cutpoints[2]*100, "%"),
                                   paste0(cutpoints[2]*100, "-", cutpoints[3]*100, "%"),
                                   paste0(">", cutpoints[3]*100, "%")))
  
  # Create reclassification tables for event and non-event groups
  event_table <- table(Base_model = categories_old[y == 1], 
                       Extended_model = categories_new[y == 1])
  
  nonevent_table <- table(Base_model = categories_old[y == 0], 
                          Extended_model = categories_new[y == 0])
  
  return(list(Comparison = comparison_name, Event_group = event_table, Non_event_group = nonevent_table))
}

# Generate reclassification tables
reclass_tables_TyG <- create_reclassification_table(
  y = complete_data$cardiometabolic_diseases,
  p_old = base_predictions,
  p_new = extended_TyG_predictions,
  cutpoints = cutpoints,
  comparison_name = "Base model vs Base model + TyG"
)

reclass_tables_CHG <- create_reclassification_table(
  y = complete_data$cardiometabolic_diseases,
  p_old = base_predictions,
  p_new = extended_CHG_predictions,
  cutpoints = cutpoints,
  comparison_name = "Base model vs Base model + CHG"
)

# ==================== Results Output ====================
cat("===== NRI Results Summary =====\n\n")

cat("Comparison 1: Base model vs Base model + TyG\n")
safe_print_nri(nri_TyG, "Base model + TyG")

cat("Comparison 2: Base model vs Base model + CHG\n")
safe_print_nri(nri_CHG, "Base model + CHG")

# Display IDI results
cat("===== IDI (Integrated Discrimination Improvement) Results =====\n")
cat("IDI Definition: Integrated Discrimination Improvement index, evaluating overall improvement in model's ability to distinguish cases from controls\n")
cat("Calculation: IDI = (Mean prediction probability increase in cases) + (Mean prediction probability decrease in controls)\n\n")

cat("Comparison 1: Base model vs Base model + TyG\n")
cat("IDI point estimate:", round(idi_estimate_TyG, 4), 
    "(95% CI:", round(bootstrap_TyG$idi_ci[1], 4), "-", round(bootstrap_TyG$idi_ci[2], 4), 
    ", p =", round(bootstrap_TyG$idi_p_value, 4), ")\n\n")

cat("Comparison 2: Base model vs Base model + CHG\n")
cat("IDI point estimate:", round(idi_estimate_CHG, 4), 
    "(95% CI:", round(bootstrap_CHG$idi_ci[1], 4), "-", round(bootstrap_CHG$idi_ci[2], 4), 
    ", p =", round(bootstrap_CHG$idi_p_value, 4), ")\n\n")

# ==================== Reclassification Table Output ====================
cat("===== Reclassification Tables (Suggested for Supplementary Materials) =====\n")
cat("\nComparison 1: Base model vs Base model + TyG\n")
cat("Event group reclassification table:\n")
print(reclass_tables_TyG$Event_group)
cat("\nNon-event group reclassification table:\n")
print(reclass_tables_TyG$Non_event_group)

cat("\nComparison 2: Base model vs Base model + CHG\n")
cat("Event group reclassification table:\n")
print(reclass_tables_CHG$Event_group)
cat("\nNon-event group reclassification table:\n")
print(reclass_tables_CHG$Non_event_group)

# ==================== Model Fit Information and AUC Comparison ====================
cat("\n===== Model Fit Information and AUC Comparison =====\n")

# Calculate AUC
roc_base <- roc(complete_data$cardiometabolic_diseases ~ base_predictions)
roc_TyG <- roc(complete_data$cardiometabolic_diseases ~ extended_TyG_predictions)
roc_CHG <- roc(complete_data$cardiometabolic_diseases ~ extended_CHG_predictions)

# Calculate AUC confidence intervals
auc_ci_base <- ci.auc(roc_base)
auc_ci_TyG <- ci.auc(roc_TyG)
auc_ci_CHG <- ci.auc(roc_CHG)

cat("Base model:\n")
cat("  AIC:", round(AIC(base_model), 2), "\n")
cat("  AUC:", round(auc(roc_base), 4), 
    "(95% CI:", round(auc_ci_base[1], 4), "-", round(auc_ci_base[3], 4), ")\n\n")

cat("Base model + TyG:\n")
cat("  AIC:", round(AIC(extended_model_TyG), 2), "\n")
cat("  AUC:", round(auc(roc_TyG), 4), 
    "(95% CI:", round(auc_ci_TyG[1], 4), "-", round(auc_ci_TyG[3], 4), ")\n\n")

cat("Base model + CHG:\n")
cat("  AIC:", round(AIC(extended_model_CHG), 2), "\n")
cat("  AUC:", round(auc(roc_CHG), 4), 
    "(95% CI:", round(auc_ci_CHG[1], 4), "-", round(auc_ci_CHG[3], 4), ")\n\n")

# Statistical test for AUC comparison
cat("AUC comparison (DeLong test):\n")
cat("Base model vs Base model + TyG: p =", 
    round(roc.test(roc_base, roc_TyG)$p.value, 4), "\n")
cat("Base model vs Base model + CHG: p =", 
    round(roc.test(roc_base, roc_CHG)$p.value, 4), "\n")

# ==================== Results Interpretation ====================
cat("\n===== Results Interpretation =====\n")

# TyG model interpretation
cat("\n--- Base model + TyG Analysis ---\n")
if (!is.na(nri_TyG$nri_est) && nri_TyG$nri_est > 0) {
  cat(sprintf("Total NRI = %.4f > 0, indicates improved reclassification ability after adding TyG\n", nri_TyG$nri_est))
} else if (!is.na(nri_TyG$nri_est)) {
  cat(sprintf("Total NRI = %.4f < 0, indicates decreased reclassification ability after adding TyG\n", nri_TyG$nri_est))
}

if (!is.na(nri_TyG$nri_p) && nri_TyG$nri_p < 0.05) {
  cat(sprintf("NRI p-value = %.4f < 0.05, indicates statistically significant improvement\n", nri_TyG$nri_p))
} else if (!is.na(nri_TyG$nri_p)) {
  cat(sprintf("NRI p-value = %.4f > 0.05, indicates improvement is not statistically significant\n", nri_TyG$nri_p))
}

if (idi_estimate_TyG > 0) {
  improvement <- ifelse(idi_estimate_TyG > 0.01, "significant clinical improvement", "minor improvement")
  cat(sprintf("IDI = %.4f > 0, indicates %s in discrimination ability after adding TyG\n", idi_estimate_TyG, improvement))
} else {
  cat(sprintf("IDI = %.4f < 0, indicates decreased discrimination ability after adding TyG\n", idi_estimate_TyG))
}

if (bootstrap_TyG$idi_p_value < 0.05) {
  cat(sprintf("IDI p-value = %.4f < 0.05, indicates statistically significant improvement\n", bootstrap_TyG$idi_p_value))
} else {
  cat(sprintf("IDI p-value = %.4f > 0.05, indicates improvement is not statistically significant\n", bootstrap_TyG$idi_p_value))
}

# CHG model interpretation
cat("\n--- Base model + CHG Analysis ---\n")
if (!is.na(nri_CHG$nri_est) && nri_CHG$nri_est > 0) {
  cat(sprintf("Total NRI = %.4f > 0, indicates improved reclassification ability after adding CHG\n", nri_CHG$nri_est))
} else if (!is.na(nri_CHG$nri_est)) {
  cat(sprintf("Total NRI = %.4f < 0, indicates decreased reclassification ability after adding CHG\n", nri_CHG$nri_est))
}

if (!is.na(nri_CHG$nri_p) && nri_CHG$nri_p < 0.05) {
  cat(sprintf("NRI p-value = %.4f < 0.05, indicates statistically significant improvement\n", nri_CHG$nri_p))
} else if (!is.na(nri_CHG$nri_p)) {
  cat(sprintf("NRI p-value = %.4f > 0.05, indicates improvement is not statistically significant\n", nri_CHG$nri_p))
}

if (idi_estimate_CHG > 0) {
  improvement <- ifelse(idi_estimate_CHG > 0.01, "significant clinical improvement", "minor improvement")
  cat(sprintf("IDI = %.4f > 0, indicates %s in discrimination ability after adding CHG\n", idi_estimate_CHG, improvement))
} else {
  cat(sprintf("IDI = %.4f < 0, indicates decreased discrimination ability after adding CHG\n", idi_estimate_CHG))
}

if (bootstrap_CHG$idi_p_value < 0.05) {
  cat(sprintf("IDI p-value = %.4f < 0.05, indicates statistically significant improvement\n", bootstrap_CHG$idi_p_value))
} else {
  cat(sprintf("IDI p-value = %.4f > 0.05, indicates improvement is not statistically significant\n", bootstrap_CHG$idi_p_value))
}

# Clinical significance interpretation
cat("\n--- Clinical Significance Assessment ---\n")
if (idi_estimate_TyG > 0.01) {
  cat("TyG model IDI > 0.01, generally considered clinically meaningful\n")
} else if (idi_estimate_TyG > 0) {
  cat("TyG model IDI > 0 but < 0.01, indicates limited improvement\n")
}

if (idi_estimate_CHG > 0.01) {
  cat("CHG model IDI > 0.01, generally considered clinically meaningful\n")
} else if (idi_estimate_CHG > 0) {
  cat("CHG model IDI > 0 but < 0.01, indicates limited improvement\n")
}

# ==================== Visualization ====================
# ROC curve comparison
cat("\nGenerating visualization charts...\n")

# Set graphics parameters
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# 1. ROC curve comparison
plot(roc_base, col = "blue", main = "ROC Curve Comparison", 
     xlab = "1 - Specificity", ylab = "Sensitivity",
     legacy.axes = TRUE, print.auc = FALSE)
lines(roc_TyG, col = "red")
lines(roc_CHG, col = "green")
legend("bottomright", 
       legend = c(paste0("Base model (AUC = ", round(auc(roc_base), 3), ")"),
                  paste0("Base model + TyG (AUC = ", round(auc(roc_TyG), 3), ")"),
                  paste0("Base model + CHG (AUC = ", round(auc(roc_CHG), 3), ")")), 
       col = c("blue", "red", "green"), lwd = 2, cex = 0.8)

# 2. IDI distribution visualization - TyG
hist(bootstrap_TyG$idi_boot, breaks = 30, col = "#4DAF4A", border = "white",
     main = paste("IDI Bootstrap Distribution - TyG\nEstimate:", round(idi_estimate_TyG, 4)),
     xlab = paste("IDI value (95% CI: [", round(bootstrap_TyG$idi_ci[1], 4), ", ", 
                  round(bootstrap_TyG$idi_ci[2], 4), "])"),
     ylab = "Frequency")
abline(v = idi_estimate_TyG, col = "#E41A1C", lwd = 3)
abline(v = bootstrap_TyG$idi_ci, col = "#377EB8", lwd = 2, lty = 2)
legend("topright", bty = "n", cex = 0.7,
       legend = c(paste("Estimate:", round(idi_estimate_TyG, 4)),
                  paste("95% CI: [", round(bootstrap_TyG$idi_ci[1], 4), ", ", 
                        round(bootstrap_TyG$idi_ci[2], 4), "]")),
       col = c("#E41A1C", "#377EB8"), lwd = c(3, 2), lty = c(1, 2))
grid(col = "gray80", lty = 3)

# 3. IDI distribution visualization - CHG
hist(bootstrap_CHG$idi_boot, breaks = 30, col = "#FF7F00", border = "white",
     main = paste("IDI Bootstrap Distribution - CHG\nEstimate:", round(idi_estimate_CHG, 4)),
     xlab = paste("IDI value (95% CI: [", round(bootstrap_CHG$idi_ci[1], 4), ", ", 
                  round(bootstrap_CHG$idi_ci[2], 4), "])"),
     ylab = "Frequency")
abline(v = idi_estimate_CHG, col = "#E41A1C", lwd = 3)
abline(v = bootstrap_CHG$idi_ci, col = "#377EB8", lwd = 2, lty = 2)
legend("topright", bty = "n", cex = 0.7,
       legend = c(paste("Estimate:", round(idi_estimate_CHG, 4)),
                  paste("95% CI: [", round(bootstrap_CHG$idi_ci[1], 4), ", ", 
                        round(bootstrap_CHG$idi_ci[2], 4), "]")),
       col = c("#E41A1C", "#377EB8"), lwd = c(3, 2), lty = c(1, 2))
grid(col = "gray80", lty = 3)

# 4. Model performance comparison chart
performance_data <- data.frame(
  Model = c("Base", "Base+TyG", "Base+CHG"),
  AUC = c(auc(roc_base), auc(roc_TyG), auc(roc_CHG)),
  AIC = c(AIC(base_model), AIC(extended_model_TyG), AIC(extended_model_CHG))
)

# Standardize AIC for plotting (smaller is better)
performance_data$AIC_scaled <- 1 - (performance_data$AIC - min(performance_data$AIC)) / 
  (max(performance_data$AIC) - min(performance_data$AIC))

barplot(t(as.matrix(performance_data[, c("AUC", "AIC_scaled")])), 
        beside = TRUE, names.arg = performance_data$Model,
        col = c("#377EB8", "#4DAF4A"), 
        main = "Model Performance Comparison",
        ylab = "Standardized value", ylim = c(0, 1))
legend("topright", legend = c("AUC", "1 - Standardized AIC"), 
       fill = c("#377EB8", "#4DAF4A"), cex = 0.8)

# Reset graphics parameters
par(mfrow = c(1, 1))

# ==================== Save Results ====================
cat("\n===== Saving Analysis Results =====\n")

# Create results summary table
results_summary <- data.frame(
  Comparison = c("Base vs Base+TyG", "Base vs Base+CHG"),
  NRI_Total = c(ifelse(!is.na(nri_TyG$nri_est), round(nri_TyG$nri_est, 4), NA),
                ifelse(!is.na(nri_CHG$nri_est), round(nri_CHG$nri_est, 4), NA)),
  NRI_P_Value = c(ifelse(!is.na(nri_TyG$nri_p), round(nri_TyG$nri_p, 4), NA),
                  ifelse(!is.na(nri_CHG$nri_p), round(nri_CHG$nri_p, 4), NA)),
  IDI_Estimate = c(round(idi_estimate_TyG, 4), round(idi_estimate_CHG, 4)),
  IDI_P_Value = c(round(bootstrap_TyG$idi_p_value, 4), round(bootstrap_CHG$idi_p_value, 4)),
  AUC_Base = c(round(auc(roc_base), 4), round(auc(roc_base), 4)),
  AUC_Extended = c(round(auc(roc_TyG), 4), round(auc(roc_CHG), 4)),
  AUC_P_Value = c(round(roc.test(roc_base, roc_TyG)$p.value, 4), 
                  round(roc.test(roc_base, roc_CHG)$p.value, 4))
)

print(results_summary)

# Save results to CSV file
write.csv(results_summary, "Model_Comparison_Results_Summary.csv", row.names = FALSE)
cat("Results saved to: Model_Comparison_Results_Summary.csv\n")

cat("\n===== Analysis Complete =====\n")
cat("Please include reclassification tables in supplementary materials of the paper\n")
cat("Report NRI and IDI point estimates, confidence intervals, and p-values in the results section\n")
cat("Visualization charts have been generated and can be used as figures in the paper\n")


# ==================== Enhanced Reclassification Table Function ====================
create_enhanced_reclassification_table <- function(y, p_old, p_new, cutpoints, comparison_name) {
  # Classify predicted probabilities into risk categories
  categories_old <- cut(p_old, 
                        breaks = c(-Inf, cutpoints, Inf),
                        labels = c(paste0("<", cutpoints[1]*100, "%"),
                                   paste0(cutpoints[1]*100, "-", cutpoints[2]*100, "%"),
                                   paste0(cutpoints[2]*100, "-", cutpoints[3]*100, "%"),
                                   paste0(">", cutpoints[3]*100, "%")))
  
  categories_new <- cut(p_new, 
                        breaks = c(-Inf, cutpoints, Inf),
                        labels = c(paste0("<", cutpoints[1]*100, "%"),
                                   paste0(cutpoints[1]*100, "-", cutpoints[2]*100, "%"),
                                   paste0(cutpoints[2]*100, "-", cutpoints[3]*100, "%"),
                                   paste0(">", cutpoints[3]*100, "%")))
  
  # Create reclassification tables for events and non-events
  event_table <- table(Base_model = categories_old[y == 1], 
                       Extended_model = categories_new[y == 1])
  
  nonevent_table <- table(Base_model = categories_old[y == 0], 
                          Extended_model = categories_new[y == 0])
  
  # Calculate reclassification improvement statistics
  calculate_reclassification_stats <- function(event_table, nonevent_table) {
    n_event <- sum(event_table)
    n_nonevent <- sum(nonevent_table)
    
    # Proportion of upward movement in events
    event_up <- sum(event_table[lower.tri(event_table)]) / n_event
    
    # Proportion of downward movement in events  
    event_down <- sum(event_table[upper.tri(event_table)]) / n_event
    
    # Proportion of downward movement in non-events
    nonevent_down <- sum(nonevent_table[upper.tri(nonevent_table)]) / n_nonevent
    
    # Proportion of upward movement in non-events
    nonevent_up <- sum(nonevent_table[lower.tri(nonevent_table)]) / n_nonevent
    
    return(list(
      event_up = event_up,
      event_down = event_down,
      nonevent_down = nonevent_down,
      nonevent_up = nonevent_up
    ))
  }
  
  stats <- calculate_reclassification_stats(event_table, nonevent_table)
  
  return(list(
    Comparison = comparison_name,
    Events = event_table,
    Non_events = nonevent_table,
    Reclassification_stats = stats
  ))
}

# ==================== Generate Enhanced Reclassification Tables ====================
enhanced_reclass_TyG <- create_enhanced_reclassification_table(
  y = complete_data$cardiometabolic_diseases,
  p_old = base_predictions,
  p_new = extended_TyG_predictions,
  cutpoints = cutpoints,
  comparison_name = "Base model vs Base model + TyG"
)

enhanced_reclass_CHG <- create_enhanced_reclassification_table(
  y = complete_data$cardiometabolic_diseases,
  p_old = base_predictions,
  p_new = extended_CHG_predictions,
  cutpoints = cutpoints,
  comparison_name = "Base model vs Base model + CHG"
)

# ==================== Detailed Reclassification Table Output ====================
cat("\n===== Detailed Reclassification Tables (by Predicted Risk Categories) =====\n")

# TyG model reclassification table
cat("\n", enhanced_reclass_TyG$Comparison, "\n")
cat("Events Reclassification Table (showing observed number of events):\n")
cat("Rows: Base model predicted risk categories\n")
cat("Columns: Extended model (Base + TyG) predicted risk categories\n")
print(enhanced_reclass_TyG$Events)
cat("Total events:", sum(enhanced_reclass_TyG$Events), "events\n")

cat("\nNon-events Reclassification Table (showing observed number of non-events):\n")
cat("Rows: Base model predicted risk categories\n")  
cat("Columns: Extended model (Base + TyG) predicted risk categories\n")
print(enhanced_reclass_TyG$Non_events)
cat("Total non-events:", sum(enhanced_reclass_TyG$Non_events), "non-events\n")

# CHG model reclassification table
cat("\n", enhanced_reclass_CHG$Comparison, "\n")
cat("Events Reclassification Table (showing observed number of events):\n")
cat("Rows: Base model predicted risk categories\n")
cat("Columns: Extended model (Base + CHG) predicted risk categories\n")
print(enhanced_reclass_CHG$Events)
cat("Total events:", sum(enhanced_reclass_CHG$Events), "events\n")

cat("\nNon-events Reclassification Table (showing observed number of non-events):\n")
cat("Rows: Base model predicted risk categories\n")
cat("Columns: Extended model (Base + CHG) predicted risk categories\n")
print(enhanced_reclass_CHG$Non_events)
cat("Total non-events:", sum(enhanced_reclass_CHG$Non_events), "non-events\n")

# ==================== Reclassification Statistics Summary ====================
cat("\n===== Reclassification Statistics Summary =====\n")

print_reclassification_stats <- function(reclass_obj, model_name) {
  stats <- reclass_obj$Reclassification_stats
  cat(model_name, ":\n")
  cat(sprintf("  Event upward reclassification proportion: %.4f (%d/%d)\n", 
              stats$event_up,
              round(stats$event_up * sum(reclass_obj$Events)),
              sum(reclass_obj$Events)))
  cat(sprintf("  Event downward reclassification proportion: %.4f (%d/%d)\n", 
              stats$event_down,
              round(stats$event_down * sum(reclass_obj$Events)), 
              sum(reclass_obj$Events)))
  cat(sprintf("  Non-event downward reclassification proportion: %.4f (%d/%d)\n", 
              stats$nonevent_down,
              round(stats$nonevent_down * sum(reclass_obj$Non_events)),
              sum(reclass_obj$Non_events)))
  cat(sprintf("  Non-event upward reclassification proportion: %.4f (%d/%d)\n", 
              stats$nonevent_up,
              round(stats$nonevent_up * sum(reclass_obj$Non_events)),
              sum(reclass_obj$Non_events)))
  cat("\n")
}

print_reclassification_stats(enhanced_reclass_TyG, "Base model + TyG")
print_reclassification_stats(enhanced_reclass_CHG, "Base model + CHG")

# ==================== Reclassification Visualization ====================
cat("Generating reclassification visualizations...\n")

# Reclassification heatmap function
plot_reclassification_heatmap <- function(reclass_obj, title) {
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
  
  # Events heatmap
  event_matrix <- as.matrix(reclass_obj$Events)
  event_prop <- prop.table(event_matrix, 1)  # Standardize by row
  
  image(1:ncol(event_prop), 1:nrow(event_prop), t(event_prop), 
        col = heat.colors(12), axes = FALSE, 
        main = paste(title, "- Events"), 
        xlab = "Extended model risk categories", ylab = "Base model risk categories")
  axis(1, at = 1:ncol(event_prop), labels = colnames(event_matrix))
  axis(2, at = 1:nrow(event_prop), labels = rownames(event_matrix))
  
  # Add values
  for (i in 1:nrow(event_matrix)) {
    for (j in 1:ncol(event_matrix)) {
      text(j, i, event_matrix[i, j], cex = 0.8, font = 2)
    }
  }
  
  # Non-events heatmap
  nonevent_matrix <- as.matrix(reclass_obj$Non_events)
  nonevent_prop <- prop.table(nonevent_matrix, 1)
  
  image(1:ncol(nonevent_prop), 1:nrow(nonevent_prop), t(nonevent_prop), 
        col = heat.colors(12), axes = FALSE, 
        main = paste(title, "- Non-events"), 
        xlab = "Extended model risk categories", ylab = "Base model risk categories")
  axis(1, at = 1:ncol(nonevent_prop), labels = colnames(nonevent_matrix))
  axis(2, at = 1:nrow(nonevent_prop), labels = rownames(nonevent_matrix))
  
  # Add values
  for (i in 1:nrow(nonevent_matrix)) {
    for (j in 1:ncol(nonevent_matrix)) {
      text(j, i, nonevent_matrix[i, j], cex = 0.8, font = 2)
    }
  }
  
  par(mfrow = c(1, 1))
}

# Generate reclassification heatmaps
plot_reclassification_heatmap(enhanced_reclass_TyG, "Base + TyG")
plot_reclassification_heatmap(enhanced_reclass_CHG, "Base + CHG")

# ==================== Save Reclassification Tables to Files ====================
cat("Saving reclassification tables to files...\n")

# Create reclassification results summary
reclassification_summary <- data.frame(
  Model_comparison = c("Base vs Base+TyG", "Base vs Base+CHG"),
  Event_upward_reclassification = c(
    enhanced_reclass_TyG$Reclassification_stats$event_up,
    enhanced_reclass_CHG$Reclassification_stats$event_up
  ),
  Event_downward_reclassification = c(
    enhanced_reclass_TyG$Reclassification_stats$event_down, 
    enhanced_reclass_CHG$Reclassification_stats$event_down
  ),
  Non_event_downward_reclassification = c(
    enhanced_reclass_TyG$Reclassification_stats$nonevent_down,
    enhanced_reclass_CHG$Reclassification_stats$nonevent_down
  ),
  Non_event_upward_reclassification = c(
    enhanced_reclass_TyG$Reclassification_stats$nonevent_up,
    enhanced_reclass_CHG$Reclassification_stats$nonevent_up
  )
)

# Save to CSV
write.csv(reclassification_summary, "reclassification_statistics_summary1.csv", row.names = FALSE)

# Save detailed reclassification tables
sink("detailed_reclassification_tables.txt")
cat("===== Detailed Reclassification Tables =====\n\n")

cat(enhanced_reclass_TyG$Comparison, "\n\n")
cat("Events Reclassification Table:\n")
print(enhanced_reclass_TyG$Events)
cat("\nNon-events Reclassification Table:\n") 
print(enhanced_reclass_TyG$Non_events)

cat("\n", enhanced_reclass_CHG$Comparison, "\n\n")
cat("Events Reclassification Table:\n")
print(enhanced_reclass_CHG$Events)
cat("\nNon-events Reclassification Table:\n")
print(enhanced_reclass_CHG$Non_events)
sink()

cat("Reclassification tables saved to: detailed_reclassification_tables.txt\n")
cat("Reclassification statistics saved to: reclassification_statistics_summary.csv\n")


