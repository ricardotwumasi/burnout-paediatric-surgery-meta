# ============================================================================
# BURNOUT AMONG PEDIATRIC SURGEONS: COMPREHENSIVE META-ANALYSIS
# Version: 1.0.1
# Authors: Ricardo Twumasi; Sebastian Kirdar-Smith
# Repository: https://github.com/ricardotwumasi/burnout-paediatric-surgery-meta
# ============================================================================

# ============================================================================
# PART 1: DATA PREPARATION AND EFFECT SIZE CALCULATION
# ============================================================================

cat("🔥 PEDIATRIC SURGEON BURNOUT META-ANALYSIS\n")
cat("==========================================\n")
cat("Version 1.0.0 - Initial Release\n")
cat("Repository: https://github.com/ricardotwumasi/burnout-paediatric-surgery-meta\n\n")

# Load required libraries
required_packages <- c("metafor", "tidyverse", "meta", "dplyr", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
    cat("📦 Installed package:", pkg, "\n")
  }
  library(pkg, character.only = TRUE)
}

cat("✅ All required packages loaded successfully!\n\n")

# Read and prepare the data
cat("📊 Loading pediatric burnout meta-analysis data...\n")
tryCatch({
  # Try multiple potential file paths
  if (file.exists("Burnout Peads Meta Data.csv")) {
    data <- read.csv("Burnout Peads Meta Data.csv", stringsAsFactors = FALSE)
  } else if (file.exists("BurnoutPeads/Burnout Peads Meta Data.csv")) {
    data <- read.csv("BurnoutPeads/Burnout Peads Meta Data.csv", stringsAsFactors = FALSE)
  } else if (file.exists("data/Burnout Peads Meta Data.csv")) {
    data <- read.csv("data/Burnout Peads Meta Data.csv", stringsAsFactors = FALSE)
  } else {
    stop("Data file not found. Please ensure 'Burnout Peads Meta Data.csv' is in the working directory.")
  }
  
  cat("✅ Data loaded successfully!\n")
  cat("📈 Dataset contains:", nrow(data), "studies with", sum(data$number_of_participants), "total participants\n\n")
  
}, error = function(e) {
  cat("❌ Error loading data:", e$message, "\n")
  stop("Please check file path and ensure data file exists.")
})

# Examine the data structure
cat("=== DATA STRUCTURE OVERVIEW ===\n")
str(data)
cat("\n=== FIRST FEW ROWS ===\n")
print(head(data))
cat("\n=== SUMMARY STATISTICS ===\n")
summary(data)

# ============================================================================
# DATA CLEANING AND PREPARATION
# ============================================================================

cat("\n=== DATA CLEANING AND PREPARATION ===\n")

data_clean <- data %>%
  mutate(
    # Convert burnout prevalence from percentage to proportion
    burnout_proportion = burnout_prevalence_percent / 100,
    
    # Calculate number of events (burnout cases)
    events = round(burnout_proportion * number_of_participants),
    
    # Create study labels for forest plots
    study_label = paste(authors, publication_year, sep = ", "),
    
    # Clean author names (remove et al. if present for cleaner labels)
    clean_authors = str_replace(authors, " et al\\.", ""),
    study_label_clean = paste(clean_authors, publication_year, sep = ", "),
    
    # Categorize NOS scores for quality assessment
    nos_quality = case_when(
      newcastle_ottawa_scale_nos_total >= 7 ~ "High Quality",
      newcastle_ottawa_scale_nos_total >= 4 ~ "Moderate Quality", 
      newcastle_ottawa_scale_nos_total < 4 ~ "Low Quality"
    ),
    
    # Create risk of bias categories based on NOS components
    selection_bias = case_when(
      nos_selection >= 3 ~ "Low Risk",
      nos_selection == 2 ~ "Moderate Risk",
      nos_selection < 2 ~ "High Risk"
    ),
    
    comparability_bias = case_when(
      nos_comparability >= 1 ~ "Low Risk",
      nos_comparability == 0 ~ "High Risk"
    ),
    
    outcome_bias = case_when(
      nos_outcome >= 2 ~ "Low Risk",
      nos_outcome == 1 ~ "Moderate Risk",
      nos_outcome < 1 ~ "High Risk"
    )
  )

# Data validation
cat("\n=== DATA VALIDATION ===\n")
cat("Missing data check:\n")
missing_summary <- colSums(is.na(data_clean))
print(missing_summary[missing_summary > 0])

cat("\nRange validation:\n")
cat("• Burnout proportions range:", range(data_clean$burnout_proportion), "\n")
cat("• Sample sizes range:", range(data_clean$number_of_participants), "\n")
cat("• Events range:", range(data_clean$events), "\n")
cat("• Publication years range:", range(data_clean$publication_year), "\n")

# ============================================================================
# EFFECT SIZE CALCULATION
# ============================================================================

cat("\n=== EFFECT SIZE CALCULATION ===\n")

# Calculate effect sizes using logit transformation for proportions
data_es <- escalc(
  measure = "PLO",           # Proportion logit transformation
  xi = data_clean$events,    # Number of events (burnout cases)
  ni = data_clean$number_of_participants,  # Sample size
  data = data_clean,
  slab = data_clean$study_label_clean,     # Study labels
  digits = 4
)

# Add the effect size data back to our main dataframe
data_final <- data_clean %>%
  bind_cols(
    yi = data_es$yi,    # Log odds (logit proportions)
    vi = data_es$vi,    # Sampling variances
    sei = sqrt(data_es$vi)  # Standard errors
  ) %>%
  mutate(
    # Back-transform logit to proportion for interpretation
    prop_est = transf.ilogit(yi),
    # Calculate 95% CI bounds
    ci_lb = transf.ilogit(yi - 1.96 * sei),
    ci_ub = transf.ilogit(yi + 1.96 * sei)
  )

# Display effect size calculations
cat("Effect size calculations completed!\n")
effect_size_table <- data.frame(
  Study = data_final$study_label_clean,
  Events = data_final$events,
  N = data_final$number_of_participants,
  Proportion = round(data_final$burnout_proportion, 3),
  LogitProp = round(data_final$yi, 3),
  SE = round(data_final$sei, 3),
  NOS_Total = data_final$newcastle_ottawa_scale_nos_total,
  Quality = data_final$nos_quality
)

print(effect_size_table)

cat("\n=== STUDY-LEVEL RESULTS WITH 95% CIs ===\n")
ci_table <- data.frame(
  Study = data_final$study_label_clean,
  Proportion = sprintf("%.2f%%", data_final$prop_est * 100),
  CI_95 = sprintf("[%.2f%%, %.2f%%]", data_final$ci_lb * 100, data_final$ci_ub * 100),
  Quality = data_final$nos_quality
)
print(ci_table)

# Save prepared data for analysis continuation
cat("\n=== DATA PREPARATION SUMMARY ===\n")
cat("✅ Data preparation complete. Ready for meta-analysis.\n")
cat("📊 Total studies:", nrow(data_final), "\n")
cat("👥 Total participants:", sum(data_final$number_of_participants), "\n")
cat("🔥 Total burnout cases:", sum(data_final$events), "\n")
cat("📅 Study period:", min(data_final$publication_year), "to", max(data_final$publication_year), "\n")

# Save intermediate results for part 2
save(data_final, file = "pediatric_burnout_prepared_data.RData")
cat("💾 Prepared data saved to 'pediatric_burnout_prepared_data.RData'\n")

# ============================================================================
# MAIN RANDOM EFFECTS META-ANALYSIS
# ============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("MAIN RANDOM EFFECTS META-ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Fit random-effects model using REML estimator (recommended for small samples)
res_overall <- rma(yi, vi, data = data_final, method = "REML", test = "knha")

# Display detailed results
cat("\n=== OVERALL META-ANALYSIS RESULTS ===\n")
print(res_overall)

# Calculate overall proportion with confidence interval (back-transformed)
overall_results <- predict(res_overall, transf = transf.ilogit, digits = 3)
cat("\n=== OVERALL BURNOUT PREVALENCE ===\n")
cat("🔥 Pooled Proportion:", sprintf("%.2f%% (95%% CI: %.2f%% to %.2f%%)", 
                                     overall_results$pred * 100, 
                                     overall_results$ci.lb * 100, 
                                     overall_results$ci.ub * 100), "\n")
cat("📊 Prediction Interval:", sprintf("%.2f%% to %.2f%%", 
                                       overall_results$pi.lb * 100, 
                                       overall_results$pi.ub * 100), "\n")

# Heterogeneity statistics
cat("\n=== HETEROGENEITY ASSESSMENT ===\n")
cat("I² =", sprintf("%.1f%%", res_overall$I2), "\n")
cat("τ² =", sprintf("%.2f", res_overall$tau2), "\n") 
cat("Q =", sprintf("%.2f (p = %.3f)", res_overall$QE, res_overall$QEp), "\n")

# Interpret heterogeneity
heterogeneity_interpretation <- case_when(
  res_overall$I2 >= 75 ~ "Very high heterogeneity",
  res_overall$I2 >= 50 ~ "Substantial heterogeneity",
  res_overall$I2 >= 25 ~ "Moderate heterogeneity",
  TRUE ~ "Low heterogeneity"
)
cat("📈 Interpretation:", heterogeneity_interpretation, "\n")

# ============================================================================
# FOREST PLOT CREATION
# ============================================================================

cat("\n=== CREATING FOREST PLOT ===\n")

# Prepare risk of bias data for forest plot
rob_symbols <- data_final %>%
  mutate(
    rob_selection = case_when(
      nos_selection >= 3 ~ "+",
      nos_selection == 2 ~ "?", 
      nos_selection < 2 ~ "-"
    ),
    rob_comparability = case_when(
      nos_comparability >= 1 ~ "+",
      nos_comparability == 0 ~ "-"
    ),
    rob_outcome = case_when(
      nos_outcome >= 2 ~ "+",
      nos_outcome == 1 ~ "?",
      nos_outcome < 1 ~ "-"
    )
  )

# Calculate weights for display
weights_pct <- paste0(formatC(weights(res_overall), format="f", digits=1), "%")

# Create publication-quality forest plot
pdf("forest_plot_pediatric_burnout.pdf", width = 14, height = 10)

# Set up plotting parameters for RevMan-style forest plot
par(mar = c(12, 0, 4, 2), mgp = c(3, 0.3, 0), tcl = -0.2)

# Create enhanced forest plot with risk of bias
forest_plot <- forest(res_overall,
                      transf = transf.ilogit,
                      at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                      xlim = c(-3.5, 2.5),
                      alim = c(0, 1),
                      xlab = "Proportion of Burnout in Pediatric Healthcare Workers",
                      header = c("Study", "Proportion [95% CI]"),
                      mlab = "Random Effects Model (REML)",
                      ilab = cbind(data_final$events, 
                                   data_final$number_of_participants,
                                   weights_pct),
                      ilab.xpos = c(-2.8, -2.4, -2.0),
                      ilab.pos = 2,
                      cex = 0.75,
                      efac = c(0, 1),
                      refline = transf.ilogit(res_overall$b[1]))

# Add column headers
text(c(-2.8, -2.4, -2.0), 
     res_overall$k + 2, 
     c("Events", "Total", "Weight"), 
     pos = 2, cex = 0.75, font = 2)

# Add risk of bias visualization
rob_positions <- seq(1.2, 1.8, length.out = 3)
text(rob_positions, res_overall$k + 2, c("S", "C", "O"), cex = 0.75, font = 2)
text(mean(rob_positions), res_overall$k + 3, "Risk of Bias", cex = 0.75, font = 2)

# Add risk of bias symbols
colors_rob <- c("+" = "darkgreen", "?" = "orange", "-" = "red")
for(i in 1:nrow(data_final)) {
  # Selection bias
  points(rob_positions[1], res_overall$k + 1 - i, 
         pch = 19, col = colors_rob[rob_symbols$rob_selection[i]], cex = 1.2)
  text(rob_positions[1], res_overall$k + 1 - i, 
       rob_symbols$rob_selection[i], cex = 0.6, font = 2, col = "white")
  
  # Comparability
  points(rob_positions[2], res_overall$k + 1 - i, 
         pch = 19, col = colors_rob[rob_symbols$rob_comparability[i]], cex = 1.2)
  text(rob_positions[2], res_overall$k + 1 - i, 
       rob_symbols$rob_comparability[i], cex = 0.6, font = 2, col = "white")
  
  # Outcome assessment
  points(rob_positions[3], res_overall$k + 1 - i, 
         pch = 19, col = colors_rob[rob_symbols$rob_outcome[i]], cex = 1.2)
  text(rob_positions[3], res_overall$k + 1 - i, 
       rob_symbols$rob_outcome[i], cex = 0.6, font = 2, col = "white")
}

# Add heterogeneity statistics
par(xpd = TRUE)
text(-2.5, -2, pos = 4, cex = 0.75,
     bquote(paste("Heterogeneity: ", tau^2, " = ", .(sprintf("%.2f", res_overall$tau2)),
                  "; ", chi^2, " = ", .(sprintf("%.2f", res_overall$QE)),
                  ", df = ", .(res_overall$k - res_overall$p),
                  " (P = ", .(sprintf("%.3f", res_overall$QEp)), "); ",
                  I^2, " = ", .(sprintf("%.0f%%", res_overall$I2)))))

# Test for overall effect
text(-2.5, -3, pos = 4, cex = 0.75,
     bquote(paste("Test for overall effect: Z = ", .(sprintf("%.2f", res_overall$zval)),
                  " (P = ", .(sprintf("%.3f", res_overall$pval)), ")")))

# Add risk of bias legend
text(-2.5, -5, pos = 4, cex = 0.7, font = 2, "Risk of bias legend:")
text(-2.5, -6, pos = 4, cex = 0.65, "(S) Selection of participants")
text(-2.5, -7, pos = 4, cex = 0.65, "(C) Comparability of groups") 
text(-2.5, -8, pos = 4, cex = 0.65, "(O) Outcome assessment")
text(-2.5, -9, pos = 4, cex = 0.65, "+ Low risk  ? Unclear risk  - High risk")

par(xpd = FALSE)
dev.off()

cat("📊 Forest plot saved to 'forest_plot_pediatric_burnout.pdf'\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== STUDY CHARACTERISTICS SUMMARY ===\n")
summary_stats <- data_final %>%
  summarise(
    n_studies = n(),
    total_participants = sum(number_of_participants),
    median_sample_size = median(number_of_participants),
    range_sample_size = paste(min(number_of_participants), "to", max(number_of_participants)),
    median_burnout_rate = sprintf("%.1f%%", median(burnout_proportion * 100)),
    range_burnout_rate = paste0(sprintf("%.1f%%", min(burnout_proportion * 100)), 
                                " to ", 
                                sprintf("%.1f%%", max(burnout_proportion * 100))),
    high_quality_studies = sum(nos_quality == "High Quality"),
    moderate_quality_studies = sum(nos_quality == "Moderate Quality"),
    low_quality_studies = sum(nos_quality == "Low Quality")
  )

print(as.data.frame(summary_stats))

# Export summary for reporting
write.csv(summary_stats, "pediatric_burnout_summary_stats.csv", row.names = FALSE)
write.csv(data_final, "pediatric_burnout_final_dataset.csv", row.names = FALSE)

cat("\n🎯 Part 1 Meta-analysis completed successfully!\n")
cat("📁 Files created:\n")
cat("  • forest_plot_pediatric_burnout.pdf\n")
cat("  • pediatric_burnout_prepared_data.RData\n")
cat("  • pediatric_burnout_summary_stats.csv\n")
cat("  • pediatric_burnout_final_dataset.csv\n")

heterogeneity_level <- if(res_overall$I2 >= 50) "substantial" else "moderate"
cat(paste0("\n📈 Results show ", heterogeneity_level, " heterogeneity between studies.\n"))
cat("🚀 Ready for Part 2: Advanced analysis (subgroups, bias assessment, meta-regression)\n")

# Save workspace for part 2
save.image("pediatric_burnout_part1_workspace.RData")
cat("💾 Workspace saved for Part 2 analysis\n")
