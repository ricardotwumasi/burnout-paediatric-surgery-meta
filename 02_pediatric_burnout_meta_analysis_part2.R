# ============================================================================
# BURNOUT AMONG PEDIATRIC SURGEONS: ADVANCED META-ANALYSIS (PART 2)
# Version: 1.1.0 - UPDATED WITH COVID AND INCOME PLOTS
# Authors: Ricardo Twumasi; Sebastian Kirdar-Smith
# Repository: https://github.com/ricardotwumasi/burnout-paediatric-surgery-meta
# ============================================================================

# ============================================================================
# PART 2: SUBGROUP ANALYSIS, PUBLICATION BIAS, AND META-REGRESSION
# ============================================================================

cat("🔥 PEDIATRIC SURGEON BURNOUT META-ANALYSIS - PART 2\n")
cat("===================================================\n")
cat("Advanced Analysis: Subgroups, Bias Assessment, Meta-Regression\n")
cat("Version 1.1.0 - Added COVID and Income Plots\n\n")

# Load required libraries
required_packages <- c("metafor", "tidyverse", "ggplot2")
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
} else if (file.exists("meta_analysis_workspace.RData")) {
  load("meta_analysis_workspace.RData")
  cat("✅ Loaded workspace from COVID/Income analysis\n")
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
  write.csv(influence_table, "outputs/pediatric_burnout_influence_analysis.csv", row.names = FALSE)
  
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
  cat("\n📊 MBI vs Non-MBI difference:", sprintf("%.2f%%", diff_estimate), "\n")
}

# 5. NEW: COVID Timing Subgroup Analysis
cat("\n--- Subgroup Analysis by COVID Timing ---\n")

if("COVID_status_combined" %in% names(data_final)) {
  covid_analysis <- perform_subgroup_analysis(data_final, "COVID_status_combined", "COVID Timing")
} else {
  cat("⚠️ COVID timing variable not found in dataset\n")
  covid_analysis <- NULL
}

# 6. NEW: Country Income Level Subgroup Analysis
cat("\n--- Subgroup Analysis by Country Income Level ---\n")

if("country_income_level_clean" %in% names(data_final)) {
  # Filter out studies with missing income data
  data_income_filtered <- data_final %>% filter(!is.na(country_income_level_clean))
  
  if(nrow(data_income_filtered) > 0) {
    cat("Studies with income data:", nrow(data_income_filtered), "\n")
    cat("Studies excluded (multi-country/missing):", nrow(data_final) - nrow(data_income_filtered), "\n\n")
    
    income_analysis <- perform_subgroup_analysis(data_income_filtered, 
                                                  "country_income_level_clean", 
                                                  "Country Income Level")
  } else {
    cat("⚠️ No studies with valid income level data\n")
    income_analysis <- NULL
  }
} else {
  cat("⚠️ Income level variable not found in dataset\n")
  income_analysis <- NULL
}

# ============================================================================
# FOREST PLOTS FOR SUBGROUPS
# ============================================================================

cat("\n=== Creating Subgroup Forest Plots ===\n")

# Forest plot by Quality
pdf("outputs/forest_plot_by_quality.pdf", width = 10, height = 8)
par(mar = c(5, 4, 4, 2))
data_by_quality <- data_final[order(data_final$nos_quality), ]
forest(rma(yi, vi, data = data_by_quality, method = "REML"),
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-2, 3),
       xlab = "Proportion of Burnout",
       slab = paste(data_by_quality$study_label_clean, 
                    paste0("(", data_by_quality$nos_quality, ")")),
       header = c("Study (Quality)", "Prevalence [95% CI]"),
       cex = 0.75)
dev.off()

# Forest plot by MBI
pdf("outputs/forest_plot_by_mbi.pdf", width = 10, height = 8)
par(mar = c(5, 4, 4, 2))
data_by_mbi <- data_final[order(data_final$mbi_used), ]
forest(rma(yi, vi, data = data_by_mbi, method = "REML"),
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-2, 3),
       xlab = "Proportion of Burnout",
       slab = paste(data_by_mbi$study_label_clean, 
                    paste0("(", data_by_mbi$mbi_used, ")")),
       header = c("Study (Tool)", "Prevalence [95% CI]"),
       cex = 0.75)
dev.off()

cat("✅ Subgroup forest plots created\n")

# ============================================================================
# PUBLICATION BIAS ASSESSMENT
# ============================================================================

cat("\n=== PUBLICATION BIAS ASSESSMENT ===\n")

# Funnel plot
pdf("outputs/publication_bias_assessment.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

# Standard funnel plot
funnel(res_overall, main = "Funnel Plot", xlab = "Logit Proportion")

# Contour-enhanced funnel plot
funnel(res_overall, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"),
       refline = 0, main = "Contour-Enhanced Funnel Plot")

# Egger's test
egger_test <- regtest(res_overall, model = "lm")
cat("Egger's regression test for funnel plot asymmetry:\n")
print(egger_test)

# Begg's rank correlation test
rank_test <- ranktest(res_overall)
cat("\nBegg's rank correlation test:\n")
print(rank_test)

# Trim-and-fill analysis
taf <- trimfill(res_overall)
funnel(taf, main = "Trim-and-Fill Analysis")
cat("\nTrim-and-fill analysis:\n")
print(taf)
studies_trimmed <- taf$k0

# Influence plot
plot(influence_results)

dev.off()
par(mfrow = c(1, 1))

# Save bias assessment results
bias_results <- data.frame(
  Test = c("Egger's test", "Begg's test", "Trim-and-fill"),
  Statistic = c(egger_test$zval, rank_test$tau, NA),
  P_value = c(egger_test$pval, rank_test$pval, NA),
  Result = c(
    ifelse(egger_test$pval < 0.05, "Significant asymmetry", "No significant asymmetry"),
    ifelse(rank_test$pval < 0.05, "Significant correlation", "No significant correlation"),
    paste(studies_trimmed, "studies trimmed")
  ),
  stringsAsFactors = FALSE
)
write.csv(bias_results, "outputs/pediatric_burnout_bias_assessment.csv", row.names = FALSE)

cat("📊 Publication bias assessment complete\n")

# ============================================================================
# SENSITIVITY ANALYSIS
# ============================================================================

cat("\n=== SENSITIVITY ANALYSIS ===\n")

# Leave-one-out analysis
loo_results <- leave1out(res_overall, transf = transf.ilogit)
loo_data <- data.frame(
  Study = data_final$study_label_clean,
  Estimate = sprintf("%.2f%%", loo_results$estimate * 100),
  CI_Lower = sprintf("%.2f%%", loo_results$ci.lb * 100),
  CI_Upper = sprintf("%.2f%%", loo_results$ci.ub * 100),
  I2 = sprintf("%.1f%%", loo_results$I2),
  stringsAsFactors = FALSE
)

cat("📊 Leave-one-out sensitivity analysis:\n")
cat("Range of estimates:", 
    sprintf("%.2f%% to %.2f%%", min(loo_results$estimate * 100), max(loo_results$estimate * 100)), "\n")

# Find most influential
original_est <- transf.ilogit(res_overall$b[1]) * 100
est_changes <- abs(loo_results$estimate * 100 - original_est)
most_influential <- which.max(est_changes)

cat("Most influential study:", data_final$study_label_clean[most_influential], "\n")
cat("Change when removed:", sprintf("%.2f%%", est_changes[most_influential]), "\n")

# Save sensitivity results
write.csv(loo_data, "outputs/pediatric_burnout_sensitivity_analysis.csv", row.names = FALSE)

# ============================================================================
# META-REGRESSION ANALYSES
# ============================================================================

cat("\n=== META-REGRESSION ANALYSES ===\n")

# 1. Publication year
cat("\n--- Meta-Regression: Publication Year ---\n")
meta_reg_year <- rma(yi, vi, mods = ~ publication_year, data = data_final, method = "REML")
print(meta_reg_year)

# 2. Sample size
cat("\n--- Meta-Regression: Sample Size ---\n")
meta_reg_size <- rma(yi, vi, mods = ~ number_of_participants, data = data_final, method = "REML")
print(meta_reg_size)

# 3. Study quality (NOS total)
cat("\n--- Meta-Regression: Study Quality (NOS) ---\n")
meta_reg_nos <- rma(yi, vi, mods = ~ newcastle_ottawa_scale_nos_total, data = data_final, method = "REML")
print(meta_reg_nos)

# 4. Gender percentage (if available)
if(sum(!is.na(data_final$female_percentage_numeric)) > 5) {
  cat("\n--- Meta-Regression: Female Percentage ---\n")
  gender_data <- data_final[!is.na(data_final$female_percentage_numeric), ]
  meta_reg_gender <- rma(yi, vi, mods = ~ female_percentage_numeric, 
                         data = gender_data, method = "REML")
  print(meta_reg_gender)
} else {
  meta_reg_gender <- NULL
  cat("\n⚠️ Insufficient data for gender meta-regression\n")
}

# 5. NEW: COVID Timing Meta-Regression
if("COVID_status_combined" %in% names(data_final)) {
  cat("\n--- Meta-Regression: COVID Timing ---\n")
  meta_reg_covid <- rma(yi, vi, mods = ~ COVID_status_combined, 
                        data = data_final, method = "REML", test = "knha")
  print(meta_reg_covid)
  cat(sprintf("Variance explained by COVID timing: R² = %.1f%%\n", meta_reg_covid$R2))
} else {
  meta_reg_covid <- NULL
}

# 6. NEW: Country Income Level Meta-Regression
if("country_income_level_clean" %in% names(data_final)) {
  cat("\n--- Meta-Regression: Country Income Level ---\n")
  
  # Use only studies with valid income data
  data_income_mr <- data_final %>% 
    filter(!is.na(country_income_level_clean)) %>%
    filter(country_income_level_clean %in% c("high", "medium"))  # Exclude low if only 1 study
  
  if(nrow(data_income_mr) > 2) {
    meta_reg_income <- rma(yi, vi, mods = ~ country_income_level_clean, 
                          data = data_income_mr, method = "REML", test = "knha")
    print(meta_reg_income)
    cat(sprintf("Variance explained by income level: R² = %.1f%%\n", meta_reg_income$R2))
  } else {
    meta_reg_income <- NULL
    cat("⚠️ Insufficient studies for income level meta-regression\n")
  }
} else {
  meta_reg_income <- NULL
}

# Multivariable models
cat("\n--- Multivariable Meta-Regression Models ---\n")

# Model 1: Year + Quality
model_1 <- rma(yi, vi, mods = ~ publication_year + newcastle_ottawa_scale_nos_total, 
               data = data_final, method = "REML")
cat("\nModel 1: Year + Quality\n")
print(model_1)

# Model 2: Size + Quality
model_2 <- rma(yi, vi, mods = ~ number_of_participants + newcastle_ottawa_scale_nos_total, 
               data = data_final, method = "REML")
cat("\nModel 2: Sample Size + Quality\n")
print(model_2)

# Model 3: MBI + Quality
model_3 <- rma(yi, vi, mods = ~ mbi_used + newcastle_ottawa_scale_nos_total, 
               data = data_final, method = "REML")
cat("\nModel 3: MBI + Quality\n")
print(model_3)

# Model 4: NEW - COVID + MBI (if COVID available)
if(!is.null(meta_reg_covid)) {
  model_4 <- rma(yi, vi, mods = ~ COVID_status_combined + mbi_used, 
                 data = data_final, method = "REML")
  cat("\nModel 4: COVID + MBI\n")
  print(model_4)
} else {
  model_4 <- NULL
}

# Model 5: NEW - Income + MBI (if income available)
if(!is.null(meta_reg_income)) {
  model_5 <- rma(yi, vi, mods = ~ country_income_level_clean + mbi_used, 
                 data = data_income_mr, method = "REML")
  cat("\nModel 5: Income Level + MBI\n")
  print(model_5)
} else {
  model_5 <- NULL
}

# Model comparison
models_list <- list(
  "Year only" = meta_reg_year,
  "Size only" = meta_reg_size,
  "NOS only" = meta_reg_nos,
  "Year + NOS" = model_1,
  "Size + NOS" = model_2,
  "MBI + NOS" = model_3
)

# Add COVID and Income models if available
if(!is.null(meta_reg_covid)) {
  models_list[["COVID only"]] <- meta_reg_covid
  if(!is.null(model_4)) models_list[["COVID + MBI"]] <- model_4
}
if(!is.null(meta_reg_income)) {
  models_list[["Income only"]] <- meta_reg_income
  if(!is.null(model_5)) models_list[["Income + MBI"]] <- model_5
}

aic_comparison <- data.frame(
  Model = names(models_list),
  AIC = sapply(models_list, AIC),
  R_squared = sapply(models_list, function(x) ifelse(is.null(x$R2), 0, x$R2)),
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
write.csv(aic_comparison, "outputs/pediatric_burnout_model_comparison.csv", row.names = FALSE)

# ========================= FDR helper (BH) =========================
suppressPackageStartupMessages({ library(dplyr) })

.collect_tests <- function() {
  tests <- list()
  add <- function(name, p, family, note = NA_character_) {
    if (!is.null(p) && is.finite(p)) {
      tests[[length(tests) + 1]] <<- data.frame(
        Test   = name,
        p_raw  = as.numeric(p),
        Family = family,
        Note   = note,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ## ===== PRIMARY (confirmatory) subgroup tests =====
  if (exists("subgroup_covid_2way"))   add("COVID Timing (Pre vs During/Post)",  subgroup_covid_2way$QMp,  "Primary")
  if (exists("subgroup_income_2way"))  add("Income Level (High vs Medium)",      subgroup_income_2way$QMp, "Primary")
  
  ## ===== EXPLORATORY moderators / meta-regressions =====
  if (exists("meta_reg_year"))         add("Meta-reg: Publication year",         meta_reg_year$QMp,        "Moderators")
  if (exists("meta_reg_size"))         add("Meta-reg: Sample size",              meta_reg_size$QMp,        "Moderators")
  if (exists("meta_reg_nos"))          add("Meta-reg: Study quality (NOS)",      meta_reg_nos$QMp,         "Moderators")
  if (exists("meta_reg_gender") && !is.null(meta_reg_gender))
    add("Meta-reg: Female %",                 meta_reg_gender$QMp,      "Moderators")
  if (exists("meta_reg_covid")  && !is.null(meta_reg_covid))
    add("Meta-reg: COVID timing",             meta_reg_covid$QMp,       "Moderators")
  if (exists("meta_reg_income") && !is.null(meta_reg_income))
    add("Meta-reg: Country income level",     meta_reg_income$QMp,      "Moderators")
  
  ## ===== SMALL-STUDY/PUBLICATION BIAS =====
  if (exists("egger_test"))            add("Egger’s regression test",            egger_test$pval,          "Bias")
  if (exists("rank_test"))             add("Begg’s rank correlation test",       rank_test$pval,           "Bias")
  
  ## Bind & adjust (BH) within each family
  if (length(tests) == 0) return(invisible(NULL))
  out <- dplyr::bind_rows(tests) %>%
    dplyr::group_by(Family) %>%
    dplyr::mutate(
      p_FDR    = p.adjust(p_raw, method = "BH"),
      q        = p_FDR,
      sig_unadj= p_raw < 0.05,
      sig_FDR  = p_FDR < 0.05
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(match(Family, c("Primary","Moderators","Bias")), p_FDR, p_raw)
  out
}

fdr_table <- .collect_tests()

if (!is.null(fdr_table)) {
  cat("\n=== FDR (BH) summary by family ===\n")
  print(dplyr::mutate(fdr_table,
                      p_raw = sprintf('%.4f', p_raw),
                      p_FDR = sprintf('%.4f', p_FDR),
                      q     = sprintf('%.4f', q)))
  if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)
  utils::write.csv(fdr_table, "outputs/fdr_corrected_pvalues.csv", row.names = FALSE)
  cat("\nSaved: outputs/fdr_corrected_pvalues.csv\n")
} else {
  cat("\n[Note] No tests collected for FDR adjustment.\n")
}
# ======================= end FDR helper (BH) =======================

# ============================================================================
# META-REGRESSION PLOTS (ENHANCED WITH COVID AND INCOME)
# ============================================================================

cat("\n=== Creating Enhanced Meta-Regression Plots ===\n")

pdf("outputs/meta_regression_plots.pdf", width = 14, height = 10)
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
          main = "Meta-Regression: Female %",
          pch = 19, col = "purple", shade = TRUE)
}

# NEW: COVID Timing (if available)
if(!is.null(meta_reg_covid)) {
  # Create boxplot for categorical variable
  covid_prev_data <- data_final %>%
    mutate(prev_pct = burnout_proportion * 100)
  
  boxplot(prev_pct ~ COVID_status_combined, data = covid_prev_data,
          col = c("lightblue", "lightcoral"),
          main = "Burnout by COVID Timing",
          xlab = "COVID Period",
          ylab = "Burnout Prevalence (%)",
          las = 1)
  
  # Add points
  points(jitter(as.numeric(as.factor(covid_prev_data$COVID_status_combined)), 0.3), 
         covid_prev_data$prev_pct,
         pch = 19, col = rgb(0, 0, 0, 0.5))
  
  # Add p-value
  text(x = 1.5, y = max(covid_prev_data$prev_pct, na.rm = TRUE) * 0.95,
       labels = sprintf("p = %.3f", meta_reg_covid$QMp),
       pos = 3, cex = 1.2, font = 2)
}

# NEW: Income Level (if available)
if(!is.null(meta_reg_income)) {
  # Create boxplot for categorical variable
  income_prev_data <- data_income_mr %>%
    mutate(prev_pct = burnout_proportion * 100,
           income_label = tools::toTitleCase(country_income_level_clean))
  
  boxplot(prev_pct ~ income_label, data = income_prev_data,
          col = c("lightgreen", "lightyellow"),
          main = "Burnout by Income Level",
          xlab = "Country Income Level",
          ylab = "Burnout Prevalence (%)",
          las = 1)
  
  # Add points
  points(jitter(as.numeric(as.factor(income_prev_data$income_label)), 0.3), 
         income_prev_data$prev_pct,
         pch = 19, col = rgb(0, 0, 0, 0.5))
  
  # Add p-value
  text(x = 1.5, y = max(income_prev_data$prev_pct, na.rm = TRUE) * 0.95,
       labels = sprintf("p = %.3f", meta_reg_income$QMp),
       pos = 3, cex = 1.2, font = 2)
}

dev.off()
par(mfrow = c(1, 1))

cat("📊 Enhanced meta-regression plots saved\n")

# ============================================================================
# ADDITIONAL VISUALIZATION: COVID AND INCOME FOREST PLOTS
# ============================================================================

cat("\n=== Creating COVID and Income Specific Visualizations ===\n")

# COVID timing forest plot with subgroups
if(!is.null(meta_reg_covid)) {
  pdf("outputs/forest_plot_covid_detailed.pdf", width = 12, height = 10)
  par(mar = c(6, 4, 4, 2))
  
  data_covid_ordered <- data_final %>%
    arrange(COVID_status_combined, publication_year)
  
  metafor::forest(rma(yi, vi, data = data_covid_ordered, method = "REML"),
         transf = transf.ilogit,
         at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
         xlim = c(-2, 3),
         alim = c(0, 1),
         xlab = "Burnout Prevalence",
         slab = paste(data_covid_ordered$study_label_clean,
                      paste0("(", data_covid_ordered$COVID_status_combined, ")")),
         header = c("Study (COVID Period)", "Prevalence [95% CI]"),
         cex = 0.75,
         main = "Burnout Prevalence by COVID Timing")
  
  dev.off()
  cat("✅ COVID forest plot created\n")
}

# Income level forest plot with subgroups
if(!is.null(meta_reg_income)) {
  pdf("outputs/forest_plot_income_detailed.pdf", width = 12, height = 10)
  par(mar = c(6, 4, 4, 2))
  
  data_income_ordered <- data_income_mr %>%
    arrange(country_income_level_clean, publication_year) %>%
    mutate(income_label = tools::toTitleCase(country_income_level_clean))
  
  metafor::forest(rma(yi, vi, data = data_income_ordered, method = "REML"),
         transf = transf.ilogit,
         at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
         xlim = c(-2, 3),
         alim = c(0, 1),
         xlab = "Burnout Prevalence",
         slab = paste(data_income_ordered$study_label_clean,
                      paste0("(", data_income_ordered$income_label, ")")),
         header = c("Study (Income Level)", "Prevalence [95% CI]"),
         cex = 0.75,
         main = "Burnout Prevalence by Country Income Level")
  
  dev.off()
  cat("✅ Income level forest plot created\n")
}

# ============================================================================
# ADVANCED VISUALIZATION: ggplot2 PLOTS FOR COVID AND INCOME
# ============================================================================

cat("\n=== Creating ggplot2 Visualizations ===\n")

# COVID timing comparison plot
if(!is.null(meta_reg_covid)) {
  covid_plot_data <- data_final %>%
    mutate(
      prev_pct = burnout_proportion * 100,
      ci_lower = ci_lb * 100,
      ci_upper = ci_ub * 100,
      weight = 1/sqrt(vi)
    )
  
  p_covid <- ggplot(covid_plot_data, aes(x = COVID_status_combined, y = prev_pct)) +
    geom_jitter(aes(size = weight), width = 0.2, alpha = 0.6, color = "steelblue") +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "lightblue") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, alpha = 0.4) +
    labs(
      title = "Burnout Prevalence by COVID Timing",
      subtitle = sprintf("QM = %.2f, p = %.3f | R² = %.1f%%", 
                        meta_reg_covid$QM, meta_reg_covid$QMp, meta_reg_covid$R2),
      x = "COVID Period",
      y = "Burnout Prevalence (%)",
      size = "Study Weight"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      legend.position = "right"
    )
  
  ggsave("outputs/covid_comparison_ggplot.pdf", p_covid, 
         width = 10, height = 7)
  cat("✅ COVID ggplot2 visualization created\n")
}

# Income level comparison plot
if(!is.null(meta_reg_income)) {
  income_plot_data <- data_income_mr %>%
    mutate(
      prev_pct = burnout_proportion * 100,
      ci_lower = ci_lb * 100,
      ci_upper = ci_ub * 100,
      weight = 1/sqrt(vi),
      income_label = tools::toTitleCase(country_income_level_clean)
    )
  
  p_income <- ggplot(income_plot_data, aes(x = income_label, y = prev_pct)) +
    geom_jitter(aes(size = weight), width = 0.2, alpha = 0.6, color = "darkgreen") +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "lightgreen") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, alpha = 0.4) +
    labs(
      title = "Burnout Prevalence by Country Income Level",
      subtitle = sprintf("QM = %.2f, p = %.3f | R² = %.1f%%", 
                        meta_reg_income$QM, meta_reg_income$QMp, meta_reg_income$R2),
      x = "Country Income Level",
      y = "Burnout Prevalence (%)",
      size = "Study Weight"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      legend.position = "right"
    )
  
  ggsave("outputs/income_comparison_ggplot.pdf", p_income, 
         width = 10, height = 7)
  cat("✅ Income level ggplot2 visualization created\n")
}

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
write.csv(cum_results, "outputs/pediatric_burnout_cumulative_analysis.csv", row.names = FALSE)

# Create cumulative forest plot
pdf("outputs/cumulative_forest_plot.pdf", width = 10, height = 8)
par(mar = c(6, 0, 4, 2))
metafor::forest(cum_ma,
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
if(!is.null(meta_reg_covid) && meta_reg_covid$QMp < 0.05) significant_mods <- c(significant_mods, "COVID Timing")
if(!is.null(meta_reg_income) && meta_reg_income$QMp < 0.05) significant_mods <- c(significant_mods, "Income Level")

if(length(significant_mods) > 0) {
  cat("• Significant moderators:", paste(significant_mods, collapse = ", "), "\n")
} else {
  cat("• No significant moderators identified\n")
}

# COVID findings (if available)
if(!is.null(covid_analysis)) {
  cat("\n🦠 COVID TIMING:\n")
  for(period in names(covid_analysis$results)) {
    result <- covid_analysis$results[[period]]
    cat(sprintf("• %s: %.1f%% (k=%d)\n", period, result$estimate, result$n_studies))
  }
  cat(sprintf("• Test for difference: p = %.3f\n", meta_reg_covid$QMp))
}

# Income findings (if available)
if(!is.null(income_analysis)) {
  cat("\n💰 COUNTRY INCOME LEVEL:\n")
  for(level in names(income_analysis$results)) {
    result <- income_analysis$results[[level]]
    cat(sprintf("• %s: %.1f%% (k=%d)\n", 
                tools::toTitleCase(level), result$estimate, result$n_studies))
  }
  if(!is.null(meta_reg_income)) {
    cat(sprintf("• Test for difference: p = %.3f\n", meta_reg_income$QMp))
  }
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
cat("• Burnout affects approximately", round(overall_results$pred * 100), "% of pediatric surgeons\n")
cat("• Wide variation between studies suggests multiple contributing factors\n")
if(!is.null(meta_reg_covid)) {
  cat("• COVID pandemic timing shows", 
      ifelse(meta_reg_covid$QMp < 0.05, "significant", "no significant"),
      "association with burnout rates\n")
}
if(!is.null(meta_reg_income)) {
  cat("• Country income level shows",
      ifelse(meta_reg_income$QMp < 0.05, "significant", "no significant"),
      "differences in burnout prevalence\n")
}

# Recommendations
cat("\n💡 RECOMMENDATIONS:\n")
cat("• Standardize burnout measurement tools (preferably MBI)\n")
cat("• Investigate subspecialty and training level differences\n")
cat("• Consider workplace and demographic factors\n")
if(!is.null(meta_reg_covid)) {
  cat("• Monitor long-term effects of COVID-19 pandemic on burnout\n")
}
if(!is.null(meta_reg_income)) {
  cat("• Examine system-level factors across different income settings\n")
}
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

# Add COVID and income files if created
if(!is.null(meta_reg_covid)) {
  output_files <- c(output_files, 
                   "forest_plot_covid_detailed.pdf",
                   "covid_comparison_ggplot.pdf")
}
if(!is.null(meta_reg_income)) {
  output_files <- c(output_files,
                   "forest_plot_income_detailed.pdf",
                   "income_comparison_ggplot.pdf")
}

for(file in output_files) {
  cat("  •", file, "\n")
}

# Save final workspace
save.image("outputs/pediatric_burnout_complete_analysis.RData")

cat("\n🎉 PART 2: ADVANCED META-ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("🚀 All files saved and ready for publication/presentation\n")
cat("📊 Repository ready for Git version control\n")
cat("\n✨ NEW: COVID and Income Level analyses included with visualizations!\n")
