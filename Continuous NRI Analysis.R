# ==================== Continuous NRI Analysis ====================
# Complete script for continuous Net Reclassification Improvement analysis

# Install and load required packages
rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/Scientific Reports")
library(PredictABEL)
library(pROC)
library(nricens)
library(boot)

# Load data
data <- read.csv("CHG-CMD_imputed.csv", na.strings = c("NA", ""))

# Select complete data for model fitting
complete_data <- na.omit(data[, c("cardiometabolic_diseases", "Age", "Gender", "Education", 
                                  "Smoking", "Drinking", "Residence", "Married_status", 
                                  "Sleep", "Glu", "TG", "HDL", "LDL", "CRP", "BMI", "TyG_BMI","BRI","CHG","TyG","TyG_BRI")])
complete_data$cardiometabolic_diseases <- ifelse(complete_data$cardiometabolic_diseases == "No", 0, 1)

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

# ==================== Continuous NRI Calculation ====================
cat("===== Continuous NRI (cNRI) Method Description =====\n")
cat("Continuous NRI Definition:\n")
cat("cNRI measures the proportion of individuals correctly reclassified without using predefined risk categories\n")
cat("Calculation: cNRI = [P(up|event) - P(down|event)] + [P(down|nonevent) - P(up|nonevent)]\n")
cat("Where:\n")
cat("  P(up|event): Proportion of events with increased predicted probability in new model\n")
cat("  P(down|event): Proportion of events with decreased predicted probability in new model\n")
cat("  P(down|nonevent): Proportion of non-events with decreased predicted probability in new model\n")
cat("  P(up|nonevent): Proportion of non-events with increased predicted probability in new model\n\n")

# Function to calculate continuous NRI
calculate_continuous_nri <- function(y, p_old, p_new) {
  # Identify events and non-events
  events <- which(y == 1)
  nonevents <- which(y == 0)
  
  # Calculate proportions for events
  events_up <- mean(p_new[events] > p_old[events])    # Events moving up
  events_down <- mean(p_new[events] < p_old[events])  # Events moving down
  events_same <- mean(p_new[events] == p_old[events]) # Events staying same
  
  # Calculate proportions for non-events
  nonevents_up <- mean(p_new[nonevents] > p_old[nonevents])    # Non-events moving up
  nonevents_down <- mean(p_new[nonevents] < p_old[nonevents])  # Non-events moving down
  nonevents_same <- mean(p_new[nonevents] == p_old[nonevents]) # Non-events staying same
  
  # Calculate continuous NRI components
  cNRI_events <- events_up - events_down
  cNRI_nonevents <- nonevents_down - nonevents_up
  cNRI_total <- cNRI_events + cNRI_nonevents
  
  return(list(
    cNRI_total = cNRI_total,
    cNRI_events = cNRI_events,
    cNRI_nonevents = cNRI_nonevents,
    events_up = events_up,
    events_down = events_down,
    events_same = events_same,
    nonevents_up = nonevents_up,
    nonevents_down = nonevents_down,
    nonevents_same = nonevents_same,
    n_events = length(events),
    n_nonevents = length(nonevents)
  ))
}

# Calculate point estimates for continuous NRI
cat("Calculating continuous NRI point estimates...\n")
cNRI_TyG <- calculate_continuous_nri(
  complete_data$cardiometabolic_diseases,
  base_predictions,
  extended_TyG_predictions
)

cNRI_CHG <- calculate_continuous_nri(
  complete_data$cardiometabolic_diseases,
  base_predictions,
  extended_CHG_predictions
)

# ==================== Bootstrap for Confidence Intervals ====================
cat("Performing bootstrap for confidence intervals...\n")

# Bootstrap function for continuous NRI
bootstrap_continuous_nri <- function(y, data, base_formula, extended_formula, n_boot = 1000) {
  set.seed(123)
  cNRI_total_boot <- numeric(n_boot)
  cNRI_events_boot <- numeric(n_boot)
  cNRI_nonevents_boot <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    # Bootstrap sample with replacement
    boot_idx <- sample(1:nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]
    
    tryCatch({
      # Fit models on bootstrap sample
      base_model_boot <- glm(base_formula, data = boot_data, family = binomial)
      extended_model_boot <- glm(extended_formula, data = boot_data, family = binomial)
      
      # Get predictions
      base_pred_boot <- predict(base_model_boot, type = "response")
      extended_pred_boot <- predict(extended_model_boot, type = "response")
      
      # Calculate continuous NRI for bootstrap sample
      cNRI_boot <- calculate_continuous_nri(
        boot_data$cardiometabolic_diseases,
        base_pred_boot,
        extended_pred_boot
      )
      
      cNRI_total_boot[i] <- cNRI_boot$cNRI_total
      cNRI_events_boot[i] <- cNRI_boot$cNRI_events
      cNRI_nonevents_boot[i] <- cNRI_boot$cNRI_nonevents
    }, error = function(e) {
      cNRI_total_boot[i] <- NA
      cNRI_events_boot[i] <- NA
      cNRI_nonevents_boot[i] <- NA
    })
  }
  
  # Remove NA values
  cNRI_total_boot <- na.omit(cNRI_total_boot)
  cNRI_events_boot <- na.omit(cNRI_events_boot)
  cNRI_nonevents_boot <- na.omit(cNRI_nonevents_boot)
  
  # Calculate 95% confidence intervals
  ci_total <- quantile(cNRI_total_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  ci_events <- quantile(cNRI_events_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  ci_nonevents <- quantile(cNRI_nonevents_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Calculate p-values (two-sided)
  p_total <- 2 * min(mean(cNRI_total_boot <= 0), mean(cNRI_total_boot >= 0))
  p_events <- 2 * min(mean(cNRI_events_boot <= 0), mean(cNRI_events_boot >= 0))
  p_nonevents <- 2 * min(mean(cNRI_nonevents_boot <= 0), mean(cNRI_nonevents_boot >= 0))
  
  # Ensure p-values don't exceed 1
  p_total <- min(p_total, 1)
  p_events <- min(p_events, 1)
  p_nonevents <- min(p_nonevents, 1)
  
  return(list(
    cNRI_total_boot = cNRI_total_boot,
    cNRI_events_boot = cNRI_events_boot,
    cNRI_nonevents_boot = cNRI_nonevents_boot,
    ci_total = ci_total,
    ci_events = ci_events,
    ci_nonevents = ci_nonevents,
    p_total = p_total,
    p_events = p_events,
    p_nonevents = p_nonevents
  ))
}

# Define model formulas
base_formula <- cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
  Drinking + Residence + Married_status + Sleep + BMI

extended_formula_TyG <- cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
  Drinking + Residence + Married_status + Sleep + BMI + TyG

extended_formula_CHG <- cardiometabolic_diseases ~ Age + Gender + Education + Smoking + 
  Drinking + Residence + Married_status + Sleep + BMI + CHG

# Perform bootstrap for both comparisons
cat("Bootstrapping for Base vs Base+TyG...\n")
bootstrap_TyG <- bootstrap_continuous_nri(
  complete_data$cardiometabolic_diseases,
  complete_data,
  base_formula,
  extended_formula_TyG
)

cat("Bootstrapping for Base vs Base+CHG...\n")
bootstrap_CHG <- bootstrap_continuous_nri(
  complete_data$cardiometabolic_diseases,
  complete_data,
  base_formula,
  extended_formula_CHG
)

# ==================== Detailed Reclassification Analysis ====================
cat("Performing detailed reclassification analysis...\n")

# Function for detailed reclassification analysis
detailed_reclassification_analysis <- function(y, p_old, p_new, model_name) {
  # Calculate probability differences
  prob_diff <- p_new - p_old
  
  # Create analysis dataframe
  analysis_df <- data.frame(
    y = y,
    p_old = p_old,
    p_new = p_new,
    prob_diff = prob_diff,
    reclass_direction = ifelse(prob_diff > 0, "Up", 
                               ifelse(prob_diff < 0, "Down", "Same")),
    reclass_magnitude = abs(prob_diff)
  )
  
  # Separate events and non-events
  events_df <- analysis_df[y == 1, ]
  nonevents_df <- analysis_df[y == 0, ]
  
  # Summary statistics for events
  events_summary <- list(
    n = nrow(events_df),
    mean_diff = mean(events_df$prob_diff),
    median_diff = median(events_df$prob_diff),
    sd_diff = sd(events_df$prob_diff),
    up_prop = mean(events_df$reclass_direction == "Up"),
    down_prop = mean(events_df$reclass_direction == "Down"),
    same_prop = mean(events_df$reclass_direction == "Same"),
    mean_magnitude = mean(events_df$reclass_magnitude)
  )
  
  # Summary statistics for non-events
  nonevents_summary <- list(
    n = nrow(nonevents_df),
    mean_diff = mean(nonevents_df$prob_diff),
    median_diff = median(nonevents_df$prob_diff),
    sd_diff = sd(nonevents_df$prob_diff),
    up_prop = mean(nonevents_df$reclass_direction == "Up"),
    down_prop = mean(nonevents_df$reclass_direction == "Down"),
    same_prop = mean(nonevents_df$reclass_direction == "Same"),
    mean_magnitude = mean(nonevents_df$reclass_magnitude)
  )
  
  return(list(
    model_name = model_name,
    events = events_summary,
    nonevents = nonevents_summary,
    full_data = analysis_df
  ))
}

# Perform detailed analysis for both models
detailed_TyG <- detailed_reclassification_analysis(
  complete_data$cardiometabolic_diseases,
  base_predictions,
  extended_TyG_predictions,
  "Base + TyG"
)

detailed_CHG <- detailed_reclassification_analysis(
  complete_data$cardiometabolic_diseases,
  base_predictions,
  extended_CHG_predictions,
  "Base + CHG"
)

# ==================== Results Output ====================
cat("\n===== CONTINUOUS NRI RESULTS =====\n\n")

# Function to print continuous NRI results
print_continuous_nri_results <- function(cNRI_point, bootstrap_results, model_name) {
  cat("MODEL:", model_name, "\n")
  cat("=", rep("=", nchar(model_name) + 7), "\n", sep = "")
  
  # Total cNRI
  cat(sprintf("Total cNRI: %.4f (95%% CI: %.4f to %.4f), p = %.4f\n",
              cNRI_point$cNRI_total,
              bootstrap_results$ci_total[1],
              bootstrap_results$ci_total[2],
              bootstrap_results$p_total))
  
  # Events cNRI
  cat(sprintf("Events cNRI: %.4f (95%% CI: %.4f to %.4f), p = %.4f\n",
              cNRI_point$cNRI_events,
              bootstrap_results$ci_events[1],
              bootstrap_results$ci_events[2],
              bootstrap_results$p_events))
  
  # Non-events cNRI
  cat(sprintf("Non-events cNRI: %.4f (95%% CI: %.4f to %.4f), p = %.4f\n",
              cNRI_point$cNRI_nonevents,
              bootstrap_results$ci_nonevents[1],
              bootstrap_results$ci_nonevents[2],
              bootstrap_results$p_nonevents))
  
  cat("\nDetailed reclassification proportions:\n")
  cat(sprintf("  Events: Up=%.3f, Down=%.3f, Same=%.3f\n",
              cNRI_point$events_up, cNRI_point$events_down, cNRI_point$events_same))
  cat(sprintf("  Non-events: Up=%.3f, Down=%.3f, Same=%.3f\n",
              cNRI_point$nonevents_up, cNRI_point$nonevents_down, cNRI_point$nonevents_same))
  cat(sprintf("  Sample sizes: Events=%d, Non-events=%d\n",
              cNRI_point$n_events, cNRI_point$n_nonevents))
  cat("\n")
}

# Print results for both models
print_continuous_nri_results(cNRI_TyG, bootstrap_TyG, "Base model vs Base + TyG")
print_continuous_nri_results(cNRI_CHG, bootstrap_CHG, "Base model vs Base + CHG")

# ==================== Detailed Analysis Output ====================
cat("\n===== DETAILED RECLASSIFICATION ANALYSIS =====\n\n")

print_detailed_analysis <- function(detailed_analysis, model_name) {
  cat("MODEL:", model_name, "\n")
  cat("=", rep("=", nchar(model_name) + 7), "\n", sep = "")
  
  cat("EVENTS (n =", detailed_analysis$events$n, "):\n")
  cat(sprintf("  Mean probability change: %.4f (SD=%.4f)\n", 
              detailed_analysis$events$mean_diff, detailed_analysis$events$sd_diff))
  cat(sprintf("  Median probability change: %.4f\n", detailed_analysis$events$median_diff))
  cat(sprintf("  Mean reclassification magnitude: %.4f\n", detailed_analysis$events$mean_magnitude))
  
  cat("\nNON-EVENTS (n =", detailed_analysis$nonevents$n, "):\n")
  cat(sprintf("  Mean probability change: %.4f (SD=%.4f)\n", 
              detailed_analysis$nonevents$mean_diff, detailed_analysis$nonevents$sd_diff))
  cat(sprintf("  Median probability change: %.4f\n", detailed_analysis$nonevents$median_diff))
  cat(sprintf("  Mean reclassification magnitude: %.4f\n", detailed_analysis$nonevents$mean_magnitude))
  cat("\n")
}

print_detailed_analysis(detailed_TyG, "Base + TyG")
print_detailed_analysis(detailed_CHG, "Base + CHG")

# ==================== Statistical Significance Interpretation ====================
cat("\n===== STATISTICAL SIGNIFICANCE INTERPRETATION =====\n\n")

interpret_cNRI <- function(cNRI_point, bootstrap_results, model_name) {
  cat(model_name, ":\n")
  
  # Total cNRI interpretation
  if (bootstrap_results$p_total < 0.05) {
    if (cNRI_point$cNRI_total > 0) {
      cat(sprintf("  Total cNRI: Significant improvement (p=%.4f)\n", bootstrap_results$p_total))
    } else {
      cat(sprintf("  Total cNRI: Significant deterioration (p=%.4f)\n", bootstrap_results$p_total))
    }
  } else {
    cat(sprintf("  Total cNRI: No significant change (p=%.4f)\n", bootstrap_results$p_total))
  }
  
  # Events cNRI interpretation
  if (bootstrap_results$p_events < 0.05) {
    if (cNRI_point$cNRI_events > 0) {
      cat(sprintf("  Events cNRI: Significant improvement (p=%.4f)\n", bootstrap_results$p_events))
    } else {
      cat(sprintf("  Events cNRI: Significant deterioration (p=%.4f)\n", bootstrap_results$p_events))
    }
  } else {
    cat(sprintf("  Events cNRI: No significant change (p=%.4f)\n", bootstrap_results$p_events))
  }
  
  # Non-events cNRI interpretation
  if (bootstrap_results$p_nonevents < 0.05) {
    if (cNRI_point$cNRI_nonevents > 0) {
      cat(sprintf("  Non-events cNRI: Significant improvement (p=%.4f)\n", bootstrap_results$p_nonevents))
    } else {
      cat(sprintf("  Non-events cNRI: Significant deterioration (p=%.4f)\n", bootstrap_results$p_nonevents))
    }
  } else {
    cat(sprintf("  Non-events cNRI: No significant change (p=%.4f)\n", bootstrap_results$p_nonevents))
  }
  cat("\n")
}

interpret_cNRI(cNRI_TyG, bootstrap_TyG, "Base + TyG")
interpret_cNRI(cNRI_CHG, bootstrap_CHG, "Base + CHG")

# ==================== Visualization ====================
cat("Generating visualizations...\n")

# Set up plotting parameters
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# 1. Bootstrap distribution for Total cNRI - TyG
hist(bootstrap_TyG$cNRI_total_boot, breaks = 30, 
     col = "#4DAF4A", border = "white",
     main = paste("Bootstrap Distribution - TyG\nTotal cNRI:", 
                  round(cNRI_TyG$cNRI_total, 4)),
     xlab = paste("Total cNRI (95% CI: [", 
                  round(bootstrap_TyG$ci_total[1], 4), ", ",
                  round(bootstrap_TyG$ci_total[2], 4), "])"),
     ylab = "Frequency")
abline(v = cNRI_TyG$cNRI_total, col = "#E41A1C", lwd = 3)
abline(v = bootstrap_TyG$ci_total, col = "#377EB8", lwd = 2, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 3)
legend("topright", bty = "n", cex = 0.7,
       legend = c(paste("Estimate:", round(cNRI_TyG$cNRI_total, 4)),
                  paste("95% CI: [", round(bootstrap_TyG$ci_total[1], 4), ", ",
                        round(bootstrap_TyG$ci_total[2], 4), "]")),
       col = c("#E41A1C", "#377EB8"), lwd = c(3, 2), lty = c(1, 2))

# 2. Bootstrap distribution for Total cNRI - CHG
hist(bootstrap_CHG$cNRI_total_boot, breaks = 30, 
     col = "#FF7F00", border = "white",
     main = paste("Bootstrap Distribution - CHG\nTotal cNRI:", 
                  round(cNRI_CHG$cNRI_total, 4)),
     xlab = paste("Total cNRI (95% CI: [", 
                  round(bootstrap_CHG$ci_total[1], 4), ", ",
                  round(bootstrap_CHG$ci_total[2], 4), "])"),
     ylab = "Frequency")
abline(v = cNRI_CHG$cNRI_total, col = "#E41A1C", lwd = 3)
abline(v = bootstrap_CHG$ci_total, col = "#377EB8", lwd = 2, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 3)
legend("topright", bty = "n", cex = 0.7,
       legend = c(paste("Estimate:", round(cNRI_CHG$cNRI_total, 4)),
                  paste("95% CI: [", round(bootstrap_CHG$ci_total[1], 4), ", ",
                        round(bootstrap_CHG$ci_total[2], 4), "]")),
       col = c("#E41A1C", "#377EB8"), lwd = c(3, 2), lty = c(1, 2))

# 3. Reclassification proportions comparison
cNRI_components <- data.frame(
  Model = rep(c("TyG", "CHG"), each = 3),
  Component = rep(c("Events Up", "Events Down", "Non-events Down"), 2),
  Value = c(cNRI_TyG$events_up, cNRI_TyG$events_down, cNRI_TyG$nonevents_down,
            cNRI_CHG$events_up, cNRI_CHG$events_down, cNRI_CHG$nonevents_down)
)

barplot(matrix(c(cNRI_TyG$events_up, cNRI_TyG$events_down, cNRI_TyG$nonevents_down,
                 cNRI_CHG$events_up, cNRI_CHG$events_down, cNRI_CHG$nonevents_down),
               nrow = 3, ncol = 2),
        beside = TRUE, names.arg = c("TyG", "CHG"),
        col = c("#1B9E77", "#D95F02", "#7570B3"),
        main = "cNRI Components Comparison",
        ylab = "Proportion", ylim = c(0, max(c(cNRI_TyG$events_up, cNRI_TyG$events_down, 
                                               cNRI_TyG$nonevents_down, cNRI_CHG$events_up,
                                               cNRI_CHG$events_down, cNRI_CHG$nonevents_down)) * 1.1))
legend("topright", legend = c("Events Up", "Events Down", "Non-events Down"),
       fill = c("#1B9E77", "#D95F02", "#7570B3"), cex = 0.8)

# 4. Probability change distributions
plot(density(detailed_TyG$full_data$prob_diff[detailed_TyG$full_data$y == 1]),
     col = "#E41A1C", lwd = 2, main = "Probability Change Distributions",
     xlab = "Probability Change (New - Old)", ylab = "Density",
     xlim = range(c(detailed_TyG$full_data$prob_diff, detailed_CHG$full_data$prob_diff)))
lines(density(detailed_TyG$full_data$prob_diff[detailed_TyG$full_data$y == 0]),
      col = "#377EB8", lwd = 2)
lines(density(detailed_CHG$full_data$prob_diff[detailed_CHG$full_data$y == 1]),
      col = "#E41A1C", lwd = 2, lty = 2)
lines(density(detailed_CHG$full_data$prob_diff[detailed_CHG$full_data$y == 0]),
      col = "#377EB8", lwd = 2, lty = 2)
abline(v = 0, col = "black", lty = 3)
legend("topright", 
       legend = c("TyG Events", "TyG Non-events", "CHG Events", "CHG Non-events"),
       col = c("#E41A1C", "#377EB8", "#E41A1C", "#377EB8"),
       lwd = 2, lty = c(1, 1, 2, 2), cex = 0.7)

# Reset plotting parameters
par(mfrow = c(1, 1))

# ==================== Save Results ====================
cat("Saving results to files...\n")

# Create comprehensive results summary
results_summary <- data.frame(
  Model_comparison = c("Base vs Base+TyG", "Base vs Base+CHG"),
  cNRI_total = c(cNRI_TyG$cNRI_total, cNRI_CHG$cNRI_total),
  cNRI_total_CI_lower = c(bootstrap_TyG$ci_total[1], bootstrap_CHG$ci_total[1]),
  cNRI_total_CI_upper = c(bootstrap_TyG$ci_total[2], bootstrap_CHG$ci_total[2]),
  cNRI_total_p = c(bootstrap_TyG$p_total, bootstrap_CHG$p_total),
  cNRI_events = c(cNRI_TyG$cNRI_events, cNRI_CHG$cNRI_events),
  cNRI_events_CI_lower = c(bootstrap_TyG$ci_events[1], bootstrap_CHG$ci_events[1]),
  cNRI_events_CI_upper = c(bootstrap_TyG$ci_events[2], bootstrap_CHG$ci_events[2]),
  cNRI_events_p = c(bootstrap_TyG$p_events, bootstrap_CHG$p_events),
  cNRI_nonevents = c(cNRI_TyG$cNRI_nonevents, cNRI_CHG$cNRI_nonevents),
  cNRI_nonevents_CI_lower = c(bootstrap_TyG$ci_nonevents[1], bootstrap_CHG$ci_nonevents[1]),
  cNRI_nonevents_CI_upper = c(bootstrap_TyG$ci_nonevents[2], bootstrap_CHG$ci_nonevents[2]),
  cNRI_nonevents_p = c(bootstrap_TyG$p_nonevents, bootstrap_CHG$p_nonevents),
  events_up_prop = c(cNRI_TyG$events_up, cNRI_CHG$events_up),
  events_down_prop = c(cNRI_TyG$events_down, cNRI_CHG$events_down),
  nonevents_up_prop = c(cNRI_TyG$nonevents_up, cNRI_CHG$nonevents_up),
  nonevents_down_prop = c(cNRI_TyG$nonevents_down, cNRI_CHG$nonevents_down)
)

# Save to CSV
write.csv(results_summary, "continuous_NRI_results_summary.csv", row.names = FALSE)

# Save detailed analysis
detailed_results <- list(
  TyG_detailed = detailed_TyG,
  CHG_detailed = detailed_CHG,
  bootstrap_TyG = bootstrap_TyG,
  bootstrap_CHG = bootstrap_CHG
)

saveRDS(detailed_results, "continuous_NRI_detailed_results.rds")

# Save session info for reproducibility
sink("session_info.txt")
print(sessionInfo())
sink()

cat("\n===== ANALYSIS COMPLETE =====\n")
cat("Results saved to:\n")
cat("  - continuous_NRI_results_summary.csv\n")
cat("  - continuous_NRI_detailed_results.rds\n")
cat("  - session_info.txt\n")
cat("Visualizations have been generated and displayed.\n")
cat("Use continuous NRI results for your manuscript reporting.\n")

