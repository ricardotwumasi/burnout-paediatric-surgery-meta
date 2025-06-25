# ============================================================================
# BURNOUT AMONG PEDIATRIC SURGEONS: ADVANCED META-ANALYSIS (PART 2)
# Version: 1.0.1
# Authors: Ricardo Twumasi; Sebastian Kirdar-Smith
# Repository: https://github.com/ricardotwumasi/burnout-paediatric-surgery-meta
# ============================================================================

# ============================================================================
# PART 2: SUBGROUP ANALYSIS, PUBLICATION BIAS, AND META-REGRESSION
# ============================================================================

cat("🔥 PEDIATRIC SURGEON BURNOUT META-ANALYSIS - PART 2\n")
cat("===================================================\n")
cat("Advanced Analysis: Subgroups, Bias Assessment, Meta-Regression\n")
cat("Version 1.0.0 - Initial Release\n\n")

# Load required libraries
required_packages <- c("metafor", "tidyverse", "meta", "dplyr", "ggplot2", "viridis")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
    cat("📦 Installed package:", pkg, "\n")
  }
  library(pkg, character.only = TRUE)
}

# Load workspace from Part 1
if (file.exists("pediatric_burnout_part1_workspace.RData")) {
  load("pediatric_burnout_part1_workspace.RData")
  cat("✅ Loaded workspace from Part 1\n")
} else if (file.exists("pediatric_burnout_prepared_data.RData")) {
  load("pediatric_burnout_prepared_data.RData")
  cat("✅ Loaded prepared data from Part 1\n")
  # Re-run main analysis if needed
  res_overall <- rma(yi, vi, data = data_final, method = "REML", test = "knha")
} else {
  stop("❌ Part 1 data not found. Please run Part 1 first.")
}

cat("📊 Dataset loaded:", nrow(data_final), "studies with", sum(data_final$number_of_participants), "participants\n\n")

# ============================================================================
# INFLUENCE ANALYSIS
# ============================================================================

cat("=== INFLUENCE ANALYSIS ===\n")

# Comprehensive influence diagnostics
influence_results <- influence(res_overall)

# Create influence diagnostics table
tryCatch({
  influence_table <- data.frame(
    Study = data_final$study_label_clean,
    Cook_D = round(as.numeric(influence_results$cook.d), 3),
    DFFITS = round(as.numeric(influence_results$dffits), 3),
    Hat = round(as.numeric(influence_results$hat), 3),
    Weight = round(weights(res_overall), 1),
    stringsAsFactors = FALSE
  )
  
  # Sort by Cook's D to identify most influential
  influence_table <- influence_table[order(-influence_table$Cook_D), ]
  cat("📊 Top 5 most influential studies (by Cook's D):\n")
  print(head(influence_table, 5))
  
  # Identify most influential study
  most_influential_idx <- which.max(as.numeric(influence_results$cook.d))
  most_influential_study <- data_final$study_label_clean[most_influential_idx]
  cat("\n🎯 Most influential study:", most_influential_study, "\n")
  cat("Cook's D =", round(as.numeric(influence_results$cook.d)[most_influential_idx], 3), "\n")
  
  # Save influence results
  write.csv(influence_table, "pediatric_burnout_influence_analysis.csv", row.names = FALSE)
  
}, error = function(e) {
  cat("⚠️ Error in influence analysis. Using leave-one-out analysis instead.\n")
  loo_results <- leave1out(res_overall)
  
  original_est <- res_overall$b[1]
  est_changes <- abs(loo_results$estimate - original_est)
  
  influence_table <- data.frame(
    Study = data_final$study_label_clean,
    Est_Change = round(est_changes, 3),
    I2_Change = round(abs(loo_results$I2 - res_overall$I2), 1),
    stringsAsFactors = FALSE
  )
  
  influence_table <- influence_table[order(-influence_table$Est_Change), ]
  cat("📊 Top 5 most influential studies (by estimate change):\n")
  print(head(influence_table, 5))
  
  most_influential_idx <- which.max(est_changes)
  most_influential_study <- data_final$study_label_clean[most_influential_idx]
  cat("\n🎯 Most influential study:", most_influential_study, "\n")
  cat("Estimate change when removed:", round(est_changes[most_influential_idx], 3), "\n")
})

# ============================================================================
# SUBGROUP ANALYSES
# ============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SUBGROUP ANALYSES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Function to perform and summarize subgroup analysis
perform_subgroup_analysis <- function(data, grouping_var, var_name) {
  cat("\n--- Subgroup Analysis by", var_name, "---\n")
  
  # Create formula for mixed-effects model
  formula_str <- paste("~ ", grouping_var)
  
  # Fit mixed-effects model
  subgroup_model <- rma(yi, vi, mods = as.formula(formula_str), data = data, method = "REML")
  print(subgroup_model)
  
  # Calculate individual subgroup estimates
  grouping_levels <- unique(data[[grouping_var]])
  subgroup_results <- list()
  
  for(level in grouping_levels) {
    subset_data <- data[data[[grouping_var]] == level, ]
    if(nrow(subset_data) > 1) {
      subgroup_res <- rma(yi, vi, data = subset_data, method = "REML")
      subgroup_est <- predict(subgroup_res, transf = transf.ilogit)
      
      result_text <- sprintf("%s (n=%d): %.2f%% (95%% CI: %.2f%% to %.2f%%), I² = %.1f%%", 
                            level, nrow(subset_data),
                            subgroup_est$pred * 100, 
                            subgroup_est$ci.lb * 100, 
                            subgroup_est$ci.ub * 100,
                            subgroup_res$I2)
      cat(result_text, "\n")
      
      subgroup_results[[level]] <- list(
        n_studies = nrow(subset_data),
        estimate = subgroup_est$pred * 100,
        ci_lb = subgroup_est$ci.lb * 100,
        ci_ub = subgroup_est$ci.ub * 100,
        i2 = subgroup_res$I2
      )
    } else {
      result_text <- sprintf("%s (n=%d): %.2f%% (single study)", 
                            level, nrow(subset_data),
                            subset_data$burnout_proportion * 100)
      cat(result_text, "\n")
      
      subgroup_results[[level]] <- list(
        n_studies = nrow(subset_data),
        estimate = subset_data$burnout_proportion * 100,
        ci_lb = NA,
        ci_ub = NA,
        i2 = NA
      )
    }
  }
  
  return(list(model = subgroup_model, results = subgroup_results))
}

# 1. Subgroup by Study Quality
quality_analysis <- perform_subgroup_analysis(data_final, "nos_quality", "Study Quality")

# 2. Subgroup by Sample Size (median split)
median_n <- median(data_final$number_of_participants)
data_final$size_group <- ifelse(data_final$number_of_participants >= median_n, 
                                paste0("Large (≥", median_n, ")"), 
                                paste0("Small (<", median_n, ")"))

size_analysis <- perform_subgroup_analysis(data_final, "size_group", "Sample Size")

# 3. Subgroup by Publication Year (median split)
median_year <- median(data_final$publication_year)
data_final$year_group <- ifelse(data_final$publication_year >= median_year, 
                                paste0("Recent (≥", median_year, ")"), 
                                paste0("Older (<", median_year, ")"))

year_analysis <- perform_subgroup_analysis(data_final, "year_group", "Publication Year")

# 4. Subgroup by Burnout Measurement Tool
cat("\n--- Subgroup Analysis by Burnout Measurement Tool ---\n")

# Display MBI classification summary
mbi_summary <- table(data_final$mbi_used)
cat("📊 BURNOUT MEASUREMENT TOOL DISTRIBUTION:\n")
print(mbi_summary)
cat("MBI studies:", sum(data_final$mbi_used == "MBI"), "\n")
cat("Non-MBI studies:", sum(data_final$mbi_used == "Non-MBI"), "\n\n")

# Perform MBI subgroup analysis
mbi_analysis <- perform_subgroup_analysis(data_final, "mbi_used", "Burnout Measurement Tool")

# Calculate difference between MBI and Non-MBI
mbi_data <- data_final[data_final$mbi_used == "MBI", ]
non_mbi_data <- data_final[data_final$mbi_used == "Non-MBI", ]

if(nrow(mbi_data) > 1 && nrow(non_mbi_data) > 1) {
  mbi_res <- rma(yi, vi, data = mbi_data, method = "REML")
  non_mbi_res <- rma(yi, vi, data = non_mbi_data, method = "REML")
  
  mbi_est <- predict(mbi_res, transf = transf.ilogit)
  non_mbi_est <- predict(non_mbi_res, transf = transf.ilogit)
  
  diff_estimate <- (mbi_est$pred - non_mbi_est$pred) * 100
  cat("🔍 Difference (MBI - Non-MBI):", sprintf("%.2f%%", diff_estimate), "\n")
  
  # Test for subgroup differences
  cat("Test for subgroup differences:\n")
  cat("QM =", sprintf("%.3f", mbi_analysis$model$QM), 
      ", df =", mbi_analysis$model$m-1, 
      ", p =", sprintf("%.3f", mbi_analysis$model$QMp), "\n")
  cat("Variance explained by burnout tool: R² =", sprintf("%.1f%%", mbi_analysis$model$R2), "\n")
}

# ============================================================================
# FOREST PLOTS BY SUBGROUPS
# ============================================================================

cat("\n=== CREATING SUBGROUP FOREST PLOTS ===\n")

# Create forest plot by study quality
pdf("forest_plot_by_quality.pdf", width = 12, height = 8)
data_ordered <- data_final[order(data_final$nos_quality, -data_final$burnout_proportion), ]
res_ordered <- rma(yi, vi, data = data_ordered, method = "REML")

par(mar = c(10, 0, 4, 2))
forest(res_ordered,
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-3, 2),
       alim = c(0, 1),
       xlab = "Proportion of Burnout",
       header = c("Study", "Proportion [95% CI]"),
       mlab = "Overall Random Effects Model",
       ilab = cbind(data_ordered$nos_quality, data_ordered$number_of_participants),
       ilab.xpos = c(-2.2, -1.6),
       ilab.pos = 2,
       cex = 0.75,
       slab = data_ordered$study_label_clean)

text(c(-2.2, -1.6), res_ordered$k + 2, 
     c("Quality", "N"), pos = 2, cex = 0.75, font = 2)

dev.off()

# Create forest plot by MBI status
pdf("forest_plot_by_mbi.pdf", width = 12, height = 8)
data_mbi_ordered <- data_final[order(data_final$mbi_used, -data_final$burnout_proportion), ]
res_mbi_ordered <- rma(yi, vi, data = data_mbi_ordered, method = "REML")

par(mar = c(10, 0, 4, 2))
forest(res_mbi_ordered,
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-3, 2),
       alim = c(0, 1),
       xlab = "Proportion of Burnout",
       header = c("Study", "Proportion [95% CI]"),
       mlab = "Overall Random Effects Model",
       ilab = cbind(data_mbi_ordered$mbi_used, data_mbi_ordered$number_of_participants),
       ilab.xpos = c(-2.2, -1.6),
       ilab.pos = 2,
       cex = 0.75,
       slab = data_mbi_ordered$study_label_clean)

text(c(-2.2, -1.6), res_mbi_ordered$k + 2, 
     c("Tool", "N"), pos = 2, cex = 0.75, font = 2)

dev.off()
par(mar = c(5, 4, 4, 2))  # Reset margins

cat("📊 Subgroup forest plots saved\n")

# ============================================================================
# PUBLICATION BIAS ASSESSMENT
# ============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PUBLICATION BIAS ASSESSMENT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Create comprehensive publication bias plots
pdf("publication_bias_assessment.pdf", width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

# 1. Standard funnel plot
funnel(res_overall, xlab = "Logit Proportion", 
       main = "Funnel Plot", cex = 0.8, pch = 19)

# 2. Contour-enhanced funnel plot
funnel(res_overall, level = c(90, 95, 99), shade = c("white", "gray75", "gray50"),
       xlab = "Logit Proportion", main = "Contour-Enhanced Funnel Plot")

# 3. Egger's Regression Test
egger_test <- regtest(res_overall, model = "lm")
cat("\n--- Egger's Regression Test ---\n")
print(egger_test)
egger_interpretation <- ifelse(egger_test$pval < 0.05, 
                              "Significant asymmetry detected (p < 0.05)", 
                              "No significant asymmetry detected (p ≥ 0.05)")
cat("🔍 Interpretation:", egger_interpretation, "\n")

# 4. Trim-and-Fill Analysis
trim_fill <- trimfill(res_overall)
cat("\n--- Trim-and-Fill Analysis ---\n")
print(trim_fill)

# Plot trim-and-fill results
funnel(trim_fill, xlab = "Logit Proportion", 
       main = "Trim-and-Fill Analysis")

# Calculate adjusted estimate
tf_est <- predict(trim_fill, transf = transf.ilogit)
cat("📊 Adjusted estimate after trim-and-fill:", 
    sprintf("%.2f%% (95%% CI: %.2f%% to %.2f%%)", 
            tf_est$pred * 100, tf_est$ci.lb * 100, tf_est$ci.ub * 100), "\n")

# Check if studies were trimmed
studies_trimmed <- max(0, trim_fill$k - nrow(data_final))
cat("Studies estimated to be missing:", studies_trimmed, "\n")

# 5. Rank Correlation Test (Begg's test)
rank_test <- ranktest(res_overall)
cat("\n--- Rank Correlation Test (Begg's Test) ---\n")
print(rank_test)
rank_interpretation <- ifelse(rank_test$pval < 0.05, 
                              "Significant rank correlation detected (p < 0.05)", 
                              "No significant rank correlation detected (p ≥ 0.05)")
cat("🔍 Interpretation:", rank_interpretation, "\n")

dev.off()
par(mfrow = c(1, 1))

# Create publication bias summary
bias_summary <- data.frame(
  Test = c("Egger's Regression", "Begg's Rank", "Trim-and-Fill"),
  P_Value = c(egger_test$pval, rank_test$pval, NA),
  Result = c(egger_interpretation, rank_interpretation, 
             paste("Estimated missing studies:", studies_trimmed)),
  stringsAsFactors = FALSE
)

write.csv(bias_summary, "pediatric_burnout_bias_assessment.csv", row.names = FALSE)

# ============================================================================
# SENSITIVITY ANALYSIS
# ============================================================================

cat("\n=== SENSITIVITY ANALYSIS ===\n")

# Leave-one-out analysis
loo_results <- leave1out(res_overall)

# Create sensitivity analysis summary
loo_summary <- data.frame(
  Study_Removed = data_final$study_label_clean,
  Estimate_Without = sprintf("%.2f%%", transf.ilogit(loo_results$estimate) * 100),
  CI_Lower = sprintf("%.2f%%", transf.ilogit(loo_results$ci.lb) * 100),
  CI_Upper = sprintf("%.2f%%", transf.ilogit(loo_results$ci.ub) * 100),
  I2_Without = sprintf("%.1f%%", loo_results$I2),
  stringsAsFactors = FALSE
)

# Find studies that substantially change the overall estimate
original_est <- transf.ilogit(res_overall$b[1]) * 100
loo_ests <- transf.ilogit(loo_results$estimate) * 100
changes <- abs(loo_ests - original_est)

substantial_change_threshold <- 2  # 2% threshold
substantial_change_idx <- which(changes > substantial_change_threshold)

if(length(substantial_change_idx) > 0) {
  cat("📊 Studies causing >", substantial_change_threshold, "% change in overall estimate:\n")
  print(loo_summary[substantial_change_idx, ])
} else {
  cat("✅ No studies cause >", substantial_change_threshold, "% change in overall estimate\n")
}

# Most impactful study
max_change_idx <- which.max(changes)
cat("\n🎯 Most impactful study when removed:", data_final$study_label_clean[max_change_idx], "\n")
cat("Change in estimate:", sprintf("%.2f%%", changes[max_change_idx]), "\n")
cat("Original:", sprintf("%.2f%%", original_est), "→ Without:", 
    sprintf("%.2f%%", loo_ests[max_change_idx]), "\n")

# Save sensitivity analysis results
write.csv(loo_summary, "pediatric_burnout_sensitivity_analysis.csv", row.names = FALSE)

# ============================================================================
# GENDER DATA PREPARATION AND ANALYSIS
# ============================================================================

cat("\n=== GENDER DATA ANALYSIS ===\n")

# Convert female_percentage to numeric, handling missing/non-numeric values
data_final$female_percentage_numeric <- suppressWarnings({
  ifelse(is.na(data_final$female_percentage) | 
         data_final$female_percentage == "" | 
         data_final$female_percentage == "n/a", 
         NA, 
         as.numeric(data_final$female_percentage))
})

# Check gender data availability
gender_available <- !is.na(data_final$female_percentage_numeric)
cat("📊 Studies with gender data:", sum(gender_available), "out of", nrow(data_final), "\n")

if(sum(gender_available) > 0) {
  cat("Female percentage range:", 
      range(data_final$female_percentage_numeric, na.rm = TRUE), "%\n")
  
  # Gender data summary
  gender_summary <- data.frame(
    Study = data_final$study_label_clean,
    Female_Percent = ifelse(is.na(data_final$female_percentage_numeric), 
                            "Missing", 
                            paste0(data_final$female_percentage_numeric, "%")),
    Has_Gender_Data = !is.na(data_final$female_percentage_numeric),
    stringsAsFactors = FALSE
  )
  
  print(gender_summary)
  write.csv(gender_summary, "pediatric_burnout_gender_data.csv", row.names = FALSE)
} else {
  cat("⚠️ No usable gender data available\n")
}

# ============================================================================
# META-REGRESSION ANALYSIS
# ============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("META-REGRESSION ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Function to perform and report meta-regression
perform_meta_regression <- function(data, predictor, predictor_name) {
  cat("\n--- Meta-Regression:", predictor_name, "---\n")
  
  formula_str <- paste("~", predictor)
  meta_reg <- rma(yi, vi, mods = as.formula(formula_str), data = data, method = "REML")
  
  print(meta_reg)
  cat("R² =", sprintf("%.1f%%", meta_reg$R2), "\n")
  
  significance <- ifelse(meta_reg$pval[2] < 0.05, 
                        paste("Significant predictor (p =", sprintf("%.3f", meta_reg$pval[2]), ")"),
                        paste("Not significant (p =", sprintf("%.3f", meta_reg$pval[2]), ")"))
  cat("🔍 Interpretation:", significance, "\n")
  
  return(meta_reg)
}

# 1. Univariate meta-regressions
cat("\n--- Univariate Meta-Regressions ---\n")

meta_reg_year <- perform_meta_regression(data_final, "publication_year", "Publication Year")
meta_reg_size <- perform_meta_regression(data_final, "number_of_participants", "Sample Size")
meta_reg_nos <- perform_meta_regression(data_final, "newcastle_ottawa_scale_nos_total", "NOS Total Score")

# Individual NOS components
meta_reg_selection <- perform_meta_regression(data_final, "nos_selection", "NOS Selection")
meta_reg_comparability <- perform_meta_regression(data_final, "nos_comparability", "NOS Comparability")
meta_reg_outcome <- perform_meta_regression(data_final, "nos_outcome", "NOS Outcome")

# Gender meta-regression (if data available)
if(sum(!is.na(data_final$female_percentage_numeric)) >= 3) {
  meta_reg_gender <- perform_meta_regression(data_final, "female_percentage_numeric", 
                                           paste("Female Percentage (", sum(!is.na(data_final$female_percentage_numeric)), "studies)"))
} else {
  cat("\n⚠️ Insufficient gender data for meta-regression\n")
  meta_reg_gender <- NULL
}

# MBI as moderator
meta_reg_mbi <- perform_meta_regression(data_final, "mbi_used", "MBI Usage")

# ============================================================================
# MULTIVARIABLE META-REGRESSION
# ============================================================================

cat("\n--- Multivariable Meta-Regression ---\n")

# Model without gender
meta_reg_multi_no_gender <- rma(yi, vi, 
                                mods = ~ publication_year + number_of_participants + newcastle_ottawa_scale_nos_total,
                                data = data_final, method = "REML")

cat("\nMultivariable model (Year + Sample Size + NOS Total):\n")
print(meta_reg_multi_no_gender)
cat("R² =", sprintf("%.1f%%", meta_reg_multi_no_gender$R2), "\n")

# Quality + MBI model
meta_reg_quality_mbi <- rma(yi, vi, 
                            mods = ~ nos_quality + mbi_used, 
                            data = data_final, method = "REML")
cat("\n--- Quality + MBI Model ---\n")
print(meta_reg_quality_mbi)
cat("R² =", sprintf("%.1f%%", meta_reg_quality_mbi$R2), "\n")

# Full model with all predictors
meta_reg_full <- rma(yi, vi, 
                     mods = ~ mbi_used + newcastle_ottawa_scale_nos_total + 
                       number_of_participants + publication_year,
                     data = data_final, method = "REML")
cat("\n--- Full Multivariable Model ---\n")
print(meta_reg_full)
cat("R² =", sprintf("%.1f%%", meta_reg_full$R2), "\n")

# Model comparison
models_list <- list(
  "Null model" = res_overall,
  "Publication year" = meta_reg_year,
  "Sample size" = meta_reg_size,
  "NOS total" = meta_reg_nos,
  "MBI only" = meta_reg_mbi,
  "Quality + MBI" = meta_reg_quality_mbi,
  "Multivariable" = meta_reg_multi_no_gender,
  "Full model" = meta_reg_full
)

# Add gender model if available
if(!is.null(meta_reg_gender)) {
  models_list[["Female percentage"]] <- meta_reg_gender
}

# AIC comparison
aic_comparison <- data.frame(
  Model = names(models_list),
  AIC = round(sapply(models_list, AIC), 2),
  R_squared = round(sapply(models_list, function(x) ifelse(is.null(x$R2), 0, x$R2)), 1),
  Parameters = sapply(models_list, function(x) x$p),
  stringsAsFactors = FALSE
)
aic_comparison$Delta_AIC <- round(aic_comparison$AIC - min(aic_comparison$AIC), 2)
aic_comparison <- aic_comparison[order(aic_comparison$AIC), ]

cat("\n--- Model Comparison (AIC) ---\n")
print(aic_comparison)
cat("\n🏆 Best model (lowest AIC):", aic_comparison$Model[1], "\n")
cat("📊 Models with ΔAICc < 2 are considered equally supported\n")

# Save model comparison
write.csv(aic_comparison, "pediatric_burnout_model_comparison.csv", row.names = FALSE)

# ============================================================================
# REGRESSION PLOTS
# ============================================================================

cat("\n=== Creating Meta-Regression Plots ===\n")

pdf("meta_regression_plots.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(5, 4, 4, 2))

# Publication year
regplot(meta_reg_year, xlab = "Publication Year", ylab = "Logit Proportion",
        main = "Meta-Regression: Publication Year", 
        pch = 19, col = "darkblue", shade = TRUE)

# Sample size
regplot(meta_reg_size, xlab = "Sample Size", ylab = "Logit Proportion",
        main = "Meta-Regression: Sample Size",
        pch = 19, col = "darkred", shade = TRUE)

# NOS total score
regplot(meta_reg_nos, xlab = "Newcastle-Ottawa Scale Total", ylab = "Logit Proportion",
        main = "Meta-Regression: Study Quality",
        pch = 19, col = "darkgreen", shade = TRUE)

# Female percentage (if available)
if(!is.null(meta_reg_gender)) {
  regplot(meta_reg_gender, xlab = "Female Percentage (%)", ylab = "Logit Proportion",
          main = paste("Meta-Regression: Female %"),
          pch = 19, col = "purple", shade = TRUE)
}

# Bubble plot: Sample size vs Year
plot(data_final$publication_year, data_final$number_of_participants,
     cex = 1/sqrt(data_final$vi), 
     pch = 19, col = rgb(0.5, 0.5, 0.5, 0.7),
     xlab = "Publication Year", ylab = "Sample Size",
     main = "Study Characteristics\n(Bubble size = Precision)")

# Gender vs Burnout scatter (if available)
if(sum(!is.na(data_final$female_percentage_numeric)) > 2) {
  gender_data <- data_final[!is.na(data_final$female_percentage_numeric), ]
  plot(gender_data$female_percentage_numeric, gender_data$burnout_proportion * 100,
       pch = 19, col = "orange", cex = 1.2,
       xlab = "Female Percentage (%)", ylab = "Burnout Rate (%)",
       main = "Gender vs Burnout Rate")
  # Add regression line
  if(nrow(gender_data) > 2) {
    abline(lm(burnout_proportion * 100 ~ female_percentage_numeric, data = gender_data), 
           col = "red", lwd = 2)
  }
}

dev.off()
par(mfrow = c(1, 1))

cat("📊 Meta-regression plots saved\n")

# ============================================================================
# CUMULATIVE META-ANALYSIS
# ============================================================================

cat("\n=== CUMULATIVE META-ANALYSIS BY PUBLICATION YEAR ===\n")

# Order studies by publication year
data_by_year <- data_final[order(data_final$publication_year), ]

# Perform cumulative meta-analysis
cum_ma <- cumul(rma(yi, vi, data = data_by_year, method = "REML"))

# Cumulative results table
cum_results <- data.frame(
  Year = data_by_year$publication_year,
  Study = data_by_year$study_label_clean,
  Cumulative_Proportion = sprintf("%.2f%%", transf.ilogit(cum_ma$estimate) * 100),
  CI_Lower = sprintf("%.2f%%", transf.ilogit(cum_ma$ci.lb) * 100),
  CI_Upper = sprintf("%.2f%%", transf.ilogit(cum_ma$ci.ub) * 100),
  I2 = sprintf("%.1f%%", cum_ma$I2),
  stringsAsFactors = FALSE
)

cat("📊 Cumulative results by publication year:\n")
print(cum_results)

# Save cumulative results
write.csv(cum_results, "pediatric_burnout_cumulative_analysis.csv", row.names = FALSE)

# Create cumulative forest plot
pdf("cumulative_forest_plot.pdf", width = 10, height = 8)
par(mar = c(6, 0, 4, 2))
forest(cum_ma,
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-2.5, 2),
       alim = c(0, 1),
       xlab = "Cumulative Proportion of Burnout",
       header = c("Study Added", "Cumulative Proportion [95% CI]"),
       mlab = "Cumulative Random Effects Estimate",
       cex = 0.75,
       slab = paste(data_by_year$publication_year, data_by_year$clean_authors, sep = ": "))

dev.off()
par(mar = c(5, 4, 4, 2))

# ============================================================================
# COMPREHENSIVE SUMMARY REPORT
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("COMPREHENSIVE SUMMARY REPORT\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Overall findings
overall_results <- predict(res_overall, transf = transf.ilogit)
cat("\n🔥 MAIN FINDINGS:\n")
cat("• Pooled burnout prevalence:", sprintf("%.1f%% (95%% CI: %.1f%% to %.1f%%)", 
                                            overall_results$pred * 100,
                                            overall_results$ci.lb * 100,
                                            overall_results$ci.ub * 100), "\n")
cat("• Prediction interval:", sprintf("%.1f%% to %.1f%%", 
                                      overall_results$pi.lb * 100,
                                      overall_results$pi.ub * 100), "\n")
cat("• Heterogeneity: I² =", sprintf("%.1f%%", res_overall$I2), "\n")

# Key moderators
cat("\n📊 KEY MODERATORS:\n")
significant_mods <- c()
if(quality_analysis$model$QMp < 0.05) significant_mods <- c(significant_mods, "Study Quality")
if(mbi_analysis$model$QMp < 0.05) significant_mods <- c(significant_mods, "Burnout Measurement")
if(meta_reg_nos$pval[2] < 0.05) significant_mods <- c(significant_mods, "NOS Score")

if(length(significant_mods) > 0) {
  cat("• Significant moderators:", paste(significant_mods, collapse = ", "), "\n")
} else {
  cat("• No significant moderators identified\n")
}

# Publication bias
cat("\n📈 PUBLICATION BIAS:\n")
cat("• Egger's test: p =", sprintf("%.3f", egger_test$pval), 
    ifelse(egger_test$pval < 0.05, "(significant)", "(not significant)"), "\n")
cat("• Begg's test: p =", sprintf("%.3f", rank_test$pval), 
    ifelse(rank_test$pval < 0.05, "(significant)", "(not significant)"), "\n")
cat("• Trim-and-fill estimated missing studies:", studies_trimmed, "\n")

# Model selection
cat("\n🏆 BEST MODEL:", aic_comparison$Model[1], "\n")
cat("• R² =", sprintf("%.1f%%", aic_comparison$R_squared[1]), "\n")

# Clinical implications
cat("\n🏥 CLINICAL IMPLICATIONS:\n")
cat("• Burnout affects approximately", round(overall_results$pred * 100), "% of pediatric healthcare workers\n")
cat("• Wide variation between studies suggests multiple contributing factors\n")
cat("• Study methodology influences findings - standardization needed\n")

# Recommendations
cat("\n💡 RECOMMENDATIONS:\n")
cat("• Standardize burnout measurement tools (preferably MBI)\n")
cat("• Investigate subspecialty and training level differences\n")
cat("• Consider workplace and demographic factors\n")
cat("• Develop targeted intervention studies\n")

# Files created
cat("\n📁 FILES CREATED:\n")
output_files <- c(
  "forest_plot_by_quality.pdf",
  "forest_plot_by_mbi.pdf", 
  "publication_bias_assessment.pdf",
  "meta_regression_plots.pdf",
  "cumulative_forest_plot.pdf",
  "pediatric_burnout_influence_analysis.csv",
  "pediatric_burnout_bias_assessment.csv",
  "pediatric_burnout_sensitivity_analysis.csv",
  "pediatric_burnout_model_comparison.csv",
  "pediatric_burnout_cumulative_analysis.csv"
)

for(file in output_files) {
  cat("  •", file, "\n")
}

# Save final workspace
save.image("pediatric_burnout_complete_analysis.RData")

cat("\n🎉 PART 2: ADVANCED META-ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("🚀 All files saved and ready for publication/presentation\n")
cat("📊 Repository ready for Git version control\n")
