# ============================================================================
# BURNOUT AMONG PEDIATRIC SURGEONS: META-ANALYSIS (PART 1)
# Version: 1.1.0
# Authors: Ricardo Twumasi; Sebastian Kirdar-Smith
# Repository: https://github.com/ricardotwumasi/burnout-paediatric-surgery-meta
# ============================================================================

# Load required libraries
required_packages <- c("metafor", "tidyverse", "grid")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

cat("=============================================================================\n")
cat("META-ANALYSIS: BURNOUT AMONG PEDIATRIC SURGEONS\n")
cat("With COVID Timing and Country Income Level Subgroup Analyses\n")
cat("=============================================================================\n\n")

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

cat("PART 1: DATA PREPARATION\n")
cat("-------------------------\n")

# Read the data (handling BOM if present)
data <- read.csv("Burnout Peads Meta Data.csv", fileEncoding = "UTF-8-BOM")

# Clean the data and prepare variables
data_clean <- data %>%
  mutate(
    # Convert burnout prevalence from percentage to proportion
    burnout_proportion = burnout_prevalence_percent / 100,
    
    # Calculate number of events (burnout cases)
    events = round(burnout_proportion * number_of_participants),
    
    # Create study labels for forest plots
    study_label = paste(authors, publication_year, sep = ", "),
    
    # Clean author names
    clean_authors = str_replace(authors, " et al\\.", ""),
    study_label_clean = paste(clean_authors, publication_year, sep = ", "),
    
    # Categorize NOS scores for quality assessment
    nos_quality = case_when(
      newcastle_ottawa_scale_nos_total >= 7 ~ "High Quality",
      newcastle_ottawa_scale_nos_total >= 4 ~ "Moderate Quality", 
      newcastle_ottawa_scale_nos_total < 4 ~ "Low Quality"
    ),
    
    # Clean COVID status variable
    COVID_status = str_trim(COVID_status),
    COVID_status = str_replace(COVID_status, "Post COVID", "Post-COVID"),
    COVID_status = str_replace(COVID_status, "Pre-COVID ", "Pre-COVID"),
    
    # Create combined COVID category (Intra + Post)
    COVID_status_combined = case_when(
      COVID_status == "Pre-COVID" ~ "Pre-COVID",
      COVID_status %in% c("Intra-COVID", "Post-COVID") ~ "During/Post-COVID"
    ),
    
    # Clean country income level
    country_income_level = str_trim(country_income_level),
    country_income_level_clean = case_when(
      country_income_level == "na" ~ NA_character_,
      TRUE ~ country_income_level
    ),
    
    # Create risk of bias categories
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
    ),
    
    # Convert female percentage to numeric
    female_percentage_numeric = as.numeric(female_percentage)
  )

# Check COVID status distribution
cat("\nCOVID STATUS DISTRIBUTION:\n")
cat("Original categories:\n")
print(table(data_clean$COVID_status, useNA = "always"))
cat("\nCombined categories:\n")
print(table(data_clean$COVID_status_combined, useNA = "always"))

# Check income level distribution
cat("\nCOUNTRY INCOME LEVEL DISTRIBUTION:\n")
print(table(data_clean$country_income_level_clean, useNA = "always"))

# Calculate effect sizes using logit transformation
data_es <- escalc(
  measure = "PLO",
  xi = data_clean$events,
  ni = data_clean$number_of_participants,
  data = data_clean,
  slab = data_clean$study_label_clean,
  digits = 4
)

# Combine effect sizes with main dataframe
data_final <- data_clean %>%
  bind_cols(
    yi = data_es$yi,
    vi = data_es$vi,
    sei = sqrt(data_es$vi)
  ) %>%
  mutate(
    prop_est = transf.ilogit(yi),
    ci_lb = transf.ilogit(yi - 1.96 * sei),
    ci_ub = transf.ilogit(yi + 1.96 * sei)
  )

cat("\nData preparation complete.\n")
cat("Total studies:", nrow(data_final), "\n")
cat("Total participants:", sum(data_final$number_of_participants), "\n\n")

# ============================================================================
# PART 2: OVERALL META-ANALYSIS
# ============================================================================

cat("PART 2: OVERALL META-ANALYSIS\n")
cat("-------------------------------\n")

# Fit overall random-effects model
res_overall <- rma(yi, vi, data = data_final, method = "REML", test = "knha")

# Display results
overall_results <- predict(res_overall, transf = transf.ilogit, digits = 3)
cat("\nOVERALL BURNOUT PREVALENCE:\n")
cat("Pooled Proportion:", sprintf("%.1f%% (95%% CI: %.1f%% to %.1f%%)\n", 
                                  overall_results$pred * 100, 
                                  overall_results$ci.lb * 100, 
                                  overall_results$ci.ub * 100))
cat("Prediction Interval:", sprintf("%.1f%% to %.1f%%\n", 
                                    overall_results$pi.lb * 100, 
                                    overall_results$pi.ub * 100))

cat("\nHETEROGENEITY:\n")
cat("I² =", sprintf("%.1f%%\n", res_overall$I2))
cat("τ² =", sprintf("%.2f\n", res_overall$tau2))
cat("Q =", sprintf("%.2f (p < %.3f)\n\n", res_overall$QE, res_overall$QEp))

# Store results for FDR correction
covid_combined_results_list <- list()
income_results_list <- list()

# ============================================================================
# PART 3: COVID TIMING SUBGROUP ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("PART 3: COVID TIMING SUBGROUP ANALYSIS\n")
cat("============================================================================\n\n")

# -------------------------
# Analysis 1: Three-way comparison (Pre/Intra/Post)
# -------------------------

cat("ANALYSIS 3A: THREE-WAY COVID COMPARISON\n")
cat("----------------------------------------\n")

cat("\nSample sizes by COVID timing:\n")
covid_table <- data_final %>%
  group_by(COVID_status) %>%
  summarise(
    n_studies = n(),
    n_participants = sum(number_of_participants),
    median_burnout = median(burnout_proportion * 100)
  )
print(as.data.frame(covid_table))

# Mixed-effects model for three-way comparison
cat("\nMixed-effects model testing for differences between COVID periods:\n")
subgroup_covid_3way <- rma(yi, vi, mods = ~ COVID_status, 
                            data = data_final, method = "REML", test = "knha")
print(subgroup_covid_3way)

cat("\nVariance explained by COVID timing (3-way): R² =", 
    sprintf("%.1f%%\n", subgroup_covid_3way$R2))

# Calculate pooled estimates for each COVID period
cat("\nPooled estimates by COVID timing:\n\n")

covid_results_list <- list()

for(period in unique(data_final$COVID_status)) {
  subset_data <- data_final[data_final$COVID_status == period, ]
  n_studies <- nrow(subset_data)
  
  if(n_studies > 0) {
    cat(sprintf("%s (k = %d):\n", period, n_studies))
    
    if(n_studies > 1) {
      period_res <- rma(yi, vi, data = subset_data, method = "REML", test = "knha")
      period_est <- predict(period_res, transf = transf.ilogit)
      
      cat(sprintf("  Prevalence: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", 
                  period_est$pred * 100, 
                  period_est$ci.lb * 100, 
                  period_est$ci.ub * 100))
      cat(sprintf("  I² = %.1f%%, τ² = %.2f\n", period_res$I2, period_res$tau2))
      
      covid_results_list[[period]] <- list(
        n = n_studies,
        prev = period_est$pred * 100,
        ci_lb = period_est$ci.lb * 100,
        ci_ub = period_est$ci.ub * 100,
        i2 = period_res$I2,
        estimate = period_res$b[1],
        se = period_res$se
      )
    } else {
      # Single study - just report the values
      cat(sprintf("  Prevalence: %.1f%% (single study)\n", 
                  subset_data$burnout_proportion[1] * 100))
      
      covid_results_list[[period]] <- list(
        n = 1,
        prev = subset_data$burnout_proportion[1] * 100,
        ci_lb = NA,
        ci_ub = NA,
        i2 = NA,
        estimate = subset_data$yi[1],
        se = subset_data$sei[1]
      )
    }
  }
  cat("\n")
}

# -------------------------
# Analysis 2: Two-way comparison (Pre vs During/Post)
# -------------------------

cat("\nANALYSIS 3B: TWO-WAY COVID COMPARISON (RECOMMENDED)\n")
cat("----------------------------------------------------\n")
cat("Combining Intra-COVID and Post-COVID for more robust analysis\n\n")

cat("Sample sizes by COVID timing (combined):\n")
covid_combined_table <- data_final %>%
  group_by(COVID_status_combined) %>%
  summarise(
    n_studies = n(),
    n_participants = sum(number_of_participants),
    median_burnout = median(burnout_proportion * 100)
  )
print(as.data.frame(covid_combined_table))

# Mixed-effects model for two-way comparison
cat("\nMixed-effects model testing Pre-COVID vs During/Post-COVID:\n")
subgroup_covid_2way <- rma(yi, vi, mods = ~ COVID_status_combined, 
                            data = data_final, method = "REML", test = "knha")
print(subgroup_covid_2way)

cat("\nVariance explained by COVID timing (2-way): R² =", 
    sprintf("%.1f%%\n", subgroup_covid_2way$R2))

# Calculate pooled estimates for combined categories
cat("\nPooled estimates by COVID timing (combined categories):\n\n")

for(period in c("Pre-COVID", "During/Post-COVID")) {
  subset_data <- data_final[data_final$COVID_status_combined == period, ]
  n_studies <- nrow(subset_data)
  
  cat(sprintf("%s (k = %d):\n", period, n_studies))
  
  period_res <- rma(yi, vi, data = subset_data, method = "REML", test = "knha")
  period_est <- predict(period_res, transf = transf.ilogit)
  
  cat(sprintf("  Prevalence: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", 
              period_est$pred * 100, 
              period_est$ci.lb * 100, 
              period_est$ci.ub * 100))
  cat(sprintf("  I² = %.1f%%, τ² = %.2f\n", period_res$I2, period_res$tau2))
  cat(sprintf("  Studies included: %s\n", 
              paste(subset_data$clean_authors, collapse = ", ")))
  
  covid_combined_results_list[[period]] <- list(
    n = n_studies,
    prev = period_est$pred * 100,
    ci_lb = period_est$ci.lb * 100,
    ci_ub = period_est$ci.ub * 100,
    i2 = period_res$I2,
    estimate = period_res$b[1],
    se = period_res$se
  )
  cat("\n")
}

# Calculate difference between pre and during/post COVID
diff_covid <- covid_combined_results_list[["Pre-COVID"]]$prev - 
              covid_combined_results_list[["During/Post-COVID"]]$prev
cat("Difference in prevalence (Pre-COVID minus During/Post-COVID):", 
    sprintf("%.1f%%\n", diff_covid))

# Statistical test interpretation
if(subgroup_covid_2way$QMp < 0.05) {
  cat("Interpretation: SIGNIFICANT difference in burnout rates between COVID periods (p < 0.05)\n")
} else {
  cat("Interpretation: No significant difference in burnout rates between COVID periods (p ≥ 0.05)\n")
}

# -------------------------
# Forest plot for COVID subgroups (FIXED VERSION - NO OVERLAPPING LABELS)
# -------------------------

cat("\nCreating forest plot for COVID timing subgroups...\n")

# Create subgroup forest plot for two-way COVID comparison
pdf("outputs/forest_plot_covid_subgroups.pdf", width = 11, height = 10)

# Set margins with extra space at bottom
par(mar = c(7, 4, 4, 2))

# Order data by COVID status for plotting
data_plot_ordered <- data_final %>%
  arrange(desc(COVID_status_combined), publication_year)

# Create simple forest plot with all studies
forest(rma(yi, vi, data = data_plot_ordered, method = "REML"),
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-2, 3),
       alim = c(0, 1),
       xlab = "Proportion of Burnout",
       slab = paste(data_plot_ordered$study_label_clean, 
                    paste0("(", data_plot_ordered$COVID_status_combined, ")")),
       header = c("Study (COVID Period)", "Prevalence [95% CI]"),
       cex = 0.75,
       refline = transf.ilogit(res_overall$b[1]))

# Add summary statistics below plot with adequate spacing
par(xpd = TRUE)
text(-2, -2, pos = 4, cex = 0.7,
     sprintf("Subgroup estimates - Pre-COVID: %.1f%% (%.1f%%-%.1f%%), k=%d", 
             covid_combined_results_list[["Pre-COVID"]]$prev,
             covid_combined_results_list[["Pre-COVID"]]$ci_lb,
             covid_combined_results_list[["Pre-COVID"]]$ci_ub,
             covid_combined_results_list[["Pre-COVID"]]$n))
text(-2, -3, pos = 4, cex = 0.7,
     sprintf("During/Post-COVID: %.1f%% (%.1f%%-%.1f%%), k=%d", 
             covid_combined_results_list[["During/Post-COVID"]]$prev,
             covid_combined_results_list[["During/Post-COVID"]]$ci_lb,
             covid_combined_results_list[["During/Post-COVID"]]$ci_ub,
             covid_combined_results_list[["During/Post-COVID"]]$n))
text(-2, -4.5, pos = 4, cex = 0.7,
     sprintf("Test for subgroup differences: QM = %.2f, p = %.3f (unadjusted)", 
             subgroup_covid_2way$QM, subgroup_covid_2way$QMp))
par(xpd = FALSE)

dev.off()
cat("Forest plot saved: forest_plot_covid_subgroups.pdf\n\n")

# ============================================================================
# PART 4: COUNTRY INCOME LEVEL SUBGROUP ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("PART 4: COUNTRY INCOME LEVEL SUBGROUP ANALYSIS\n")
cat("============================================================================\n\n")

# Filter out studies with multiple countries (marked as "na")
data_income <- data_final %>%
  filter(!is.na(country_income_level_clean))

cat("Studies included in income analysis:", nrow(data_income), "\n")
cat("Studies excluded (multi-country or missing data):", nrow(data_final) - nrow(data_income), "\n\n")

# -------------------------
# Analysis 1: Three-way comparison (High/Medium/Low)
# -------------------------

cat("ANALYSIS 4A: THREE-WAY INCOME COMPARISON\n")
cat("-----------------------------------------\n")

cat("\nSample sizes by income level:\n")
income_table <- data_income %>%
  group_by(country_income_level_clean) %>%
  summarise(
    n_studies = n(),
    n_participants = sum(number_of_participants),
    median_burnout = median(burnout_proportion * 100)
  )
print(as.data.frame(income_table))

# Check if we have enough studies for 3-way comparison
if(all(table(data_income$country_income_level_clean) > 1)) {
  # Mixed-effects model for three-way comparison
  cat("\nMixed-effects model testing for differences between income levels:\n")
  subgroup_income <- rma(yi, vi, mods = ~ country_income_level_clean, 
                        data = data_income, method = "REML", test = "knha")
  print(subgroup_income)
  
  cat("\nVariance explained by income level (3-way): R² =", 
      sprintf("%.1f%%\n", subgroup_income$R2))
} else {
  cat("\nNote: Three-way comparison not conducted due to small sample sizes in some categories.\n")
  cat("Proceeding with two-way comparison (High vs Medium) only.\n\n")
}

# -------------------------
# Analysis 2: Two-way comparison (High vs Medium)
# -------------------------

cat("\nANALYSIS 4B: TWO-WAY INCOME COMPARISON (RECOMMENDED)\n")
cat("-----------------------------------------------------\n")
cat("Comparing high vs medium income countries (excluding low income, k=1)\n\n")

# Filter to high and medium income only
data_income_2way <- data_income %>%
  filter(country_income_level_clean %in% c("high", "medium"))

cat("Sample sizes by income level (high vs medium):\n")
income_2way_table <- data_income_2way %>%
  group_by(country_income_level_clean) %>%
  summarise(
    n_studies = n(),
    n_participants = sum(number_of_participants),
    median_burnout = median(burnout_proportion * 100)
  )
print(as.data.frame(income_2way_table))

# Mixed-effects model for two-way comparison
cat("\nMixed-effects model testing High vs Medium income:\n")
subgroup_income_2way <- rma(yi, vi, mods = ~ country_income_level_clean, 
                             data = data_income_2way, method = "REML", test = "knha")
print(subgroup_income_2way)

cat("\nVariance explained by income level (2-way): R² =", 
    sprintf("%.1f%%\n", subgroup_income_2way$R2))

# Calculate pooled estimates for each income level
cat("\nPooled estimates by income level:\n\n")

high_data <- data_income_2way[data_income_2way$country_income_level_clean == "high", ]
medium_data <- data_income_2way[data_income_2way$country_income_level_clean == "medium", ]

res_high <- rma(yi, vi, data = high_data, method = "REML", test = "knha")
res_medium <- rma(yi, vi, data = medium_data, method = "REML", test = "knha")

high_est <- predict(res_high, transf = transf.ilogit)
medium_est <- predict(res_medium, transf = transf.ilogit)

cat(sprintf("High income (k = %d):\n", nrow(high_data)))
cat(sprintf("  Prevalence: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", 
            high_est$pred * 100, high_est$ci.lb * 100, high_est$ci.ub * 100))
cat(sprintf("  I² = %.1f%%, τ² = %.2f\n\n", res_high$I2, res_high$tau2))

cat(sprintf("Medium income (k = %d):\n", nrow(medium_data)))
cat(sprintf("  Prevalence: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", 
            medium_est$pred * 100, medium_est$ci.lb * 100, medium_est$ci.ub * 100))
cat(sprintf("  I² = %.1f%%, τ² = %.2f\n\n", res_medium$I2, res_medium$tau2))

# Store results for FDR correction
income_results_list[["high"]] <- list(
  n = nrow(high_data),
  prev = high_est$pred * 100,
  ci_lb = high_est$ci.lb * 100,
  ci_ub = high_est$ci.ub * 100,
  i2 = res_high$I2
)

income_results_list[["medium"]] <- list(
  n = nrow(medium_data),
  prev = medium_est$pred * 100,
  ci_lb = medium_est$ci.lb * 100,
  ci_ub = medium_est$ci.ub * 100,
  i2 = res_medium$I2
)

diff_income <- high_est$pred * 100 - medium_est$pred * 100
cat(sprintf("Difference (high - medium): %.1f%%\n", diff_income))

if(subgroup_income_2way$QMp < 0.05) {
  cat("Statistical significance: p < 0.05 (SIGNIFICANT)\n\n")
} else {
  cat(sprintf("Statistical significance: p = %.3f (NOT significant)\n\n", subgroup_income_2way$QMp))
}

# -------------------------
# Forest plot for income level (FIXED VERSION)
# -------------------------

cat("Creating forest plot for income level subgroups...\n")

pdf("outputs/forest_plot_income_subgroups.pdf", width = 11, height = 10)

# Set up plotting area with extra bottom margin
par(mar = c(8, 4, 4, 2))

# Order data by income level for plotting
data_income_plot_ordered <- data_income_2way %>%
  arrange(desc(country_income_level_clean), publication_year)

# Create simple forest plot with all studies
res_income_plot <- rma(yi, vi, data = data_income_plot_ordered, method = "REML")

forest(res_income_plot,
       transf = transf.ilogit,
       at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
       xlim = c(-2, 3),
       alim = c(0, 1),
       xlab = "Proportion of Burnout",
       slab = paste(data_income_plot_ordered$study_label_clean, 
                    paste0("(", tools::toTitleCase(data_income_plot_ordered$country_income_level_clean), " income)")),
       header = c("Study (Income Level)", "Prevalence [95% CI]"),
       cex = 0.75,
       refline = transf.ilogit(res_overall$b[1]))

# Add summary statistics below plot with adequate spacing
par(xpd = TRUE)
text(-2, -2, pos = 4, cex = 0.7,
     sprintf("Subgroup estimates - High income: %.1f%% (%.1f%%-%.1f%%), k=%d", 
             high_est$pred * 100,
             high_est$ci.lb * 100,
             high_est$ci.ub * 100,
             nrow(high_data)))
text(-2, -3, pos = 4, cex = 0.7,
     sprintf("Medium income: %.1f%% (%.1f%%-%.1f%%), k=%d", 
             medium_est$pred * 100,
             medium_est$ci.lb * 100,
             medium_est$ci.ub * 100,
             nrow(medium_data)))
text(-2, -4.5, pos = 4, cex = 0.7,
     sprintf("Test for subgroup differences: QM = %.2f, p = %.3f (unadjusted)", 
             subgroup_income_2way$QM, subgroup_income_2way$QMp))
text(-2, -6, pos = 4, cex = 0.65, font = 3,
     "Note: 4 multi-country studies and 1 low-income study excluded from this analysis")
par(xpd = FALSE)

dev.off()
cat("Forest plot saved: forest_plot_income_subgroups.pdf\n\n")

# ============================================================================
# PART 5: FALSE DISCOVERY RATE CORRECTION FOR MULTIPLE COMPARISONS
# ============================================================================

cat("============================================================================\n")
cat("PART 5: FALSE DISCOVERY RATE (FDR) CORRECTION\n")
cat("============================================================================\n\n")

cat("BENJAMINI-HOCHBERG PROCEDURE\n")
cat("-----------------------------\n\n")

cat("Correcting for multiple comparisons across primary subgroup analyses.\n")
cat("This accounts for the fact that we are conducting multiple statistical tests,\n")
cat("which increases the risk of Type I errors (false positives).\n\n")

# Collect all primary p-values for correction
# Primary comparisons: COVID 2-way and Income 2-way
primary_tests <- data.frame(
  Comparison = c("COVID Timing (Pre vs During/Post)",
                 "Income Level (High vs Medium)"),
  Unadjusted_p = c(subgroup_covid_2way$QMp, 
                   subgroup_income_2way$QMp),
  stringsAsFactors = FALSE
)

# Apply Benjamini-Hochberg FDR correction
primary_tests$FDR_adjusted_p <- p.adjust(primary_tests$Unadjusted_p, method = "BH")

# Calculate q-values (proportion of false discoveries)
primary_tests$q_value <- primary_tests$FDR_adjusted_p

# Add significance indicators
primary_tests$Unadjusted_sig <- ifelse(primary_tests$Unadjusted_p < 0.05, "Yes", "No")
primary_tests$FDR_sig <- ifelse(primary_tests$FDR_adjusted_p < 0.05, "Yes", "No")

cat("PRIMARY SUBGROUP COMPARISONS (n = 2):\n")
cat("--------------------------------------\n")
print(primary_tests, row.names = FALSE)
cat("\n")

# Interpretation
cat("INTERPRETATION:\n")
cat("---------------\n")
for(i in 1:nrow(primary_tests)) {
  test_name <- primary_tests$Comparison[i]
  unadj_p <- primary_tests$Unadjusted_p[i]
  adj_p <- primary_tests$FDR_adjusted_p[i]
  
  cat(sprintf("%s:\n", test_name))
  cat(sprintf("  Unadjusted p-value: %.4f\n", unadj_p))
  cat(sprintf("  FDR-adjusted p-value: %.4f\n", adj_p))
  
  if(adj_p < 0.05) {
    cat("  FDR-adjusted result: SIGNIFICANT (q < 0.05)\n")
    cat("  Interpretation: This finding survives FDR correction at alpha=0.05\n")
  } else if(adj_p < 0.10) {
    cat("  FDR-adjusted result: MARGINALLY SIGNIFICANT (q < 0.10)\n")
    cat("  Interpretation: Suggestive evidence, but does not meet strict FDR threshold\n")
  } else {
    cat("  FDR-adjusted result: NOT SIGNIFICANT (q >= 0.10)\n")
    cat("  Interpretation: This finding does not survive FDR correction\n")
  }
  cat("\n")
}

# Summary
cat("SUMMARY OF FDR-CORRECTED RESULTS:\n")
cat("----------------------------------\n")
significant_after_fdr <- sum(primary_tests$FDR_sig == "Yes")
cat(sprintf("Number of tests: %d\n", nrow(primary_tests)))
cat(sprintf("Tests significant before FDR correction: %d\n", sum(primary_tests$Unadjusted_sig == "Yes")))
cat(sprintf("Tests significant after FDR correction: %d\n", significant_after_fdr))
cat(sprintf("Expected proportion of false discoveries: %.1f%%\n\n", 
            ifelse(significant_after_fdr > 0, max(primary_tests$FDR_adjusted_p[primary_tests$FDR_sig == "Yes"]) * 100, 0)))

# Create detailed results table for export
fdr_results_table <- primary_tests
fdr_results_table$Unadjusted_p <- sprintf("%.4f", fdr_results_table$Unadjusted_p)
fdr_results_table$FDR_adjusted_p <- sprintf("%.4f", fdr_results_table$FDR_adjusted_p)
fdr_results_table$q_value <- sprintf("%.4f", fdr_results_table$q_value)

cat("DETAILED RESULTS TABLE (formatted for export):\n")
print(fdr_results_table, row.names = FALSE)
cat("\n")

# Additional sensitivity: Also check exploratory 3-way comparisons
cat("EXPLORATORY: INCLUDING 3-WAY COMPARISONS\n")
cat("-----------------------------------------\n")
cat("Note: 3-way comparisons are exploratory due to small subgroup sizes\n\n")

all_tests <- data.frame(
  Comparison = c("COVID Timing (3-way: Pre/Intra/Post)",
                 "COVID Timing (2-way: Pre vs During/Post)",
                 "Income Level (2-way: High vs Medium)"),
  Unadjusted_p = c(subgroup_covid_3way$QMp,
                   subgroup_covid_2way$QMp,
                   subgroup_income_2way$QMp),
  Type = c("Exploratory", "Primary", "Primary"),
  stringsAsFactors = FALSE
)

all_tests$FDR_adjusted_p <- p.adjust(all_tests$Unadjusted_p, method = "BH")
all_tests$Unadjusted_sig <- ifelse(all_tests$Unadjusted_p < 0.05, "Yes", "No")
all_tests$FDR_sig <- ifelse(all_tests$FDR_adjusted_p < 0.05, "Yes", "No")

cat("ALL COMPARISONS (n = 3):\n")
print(all_tests, row.names = FALSE)
cat("\n")

cat("Note: When including exploratory 3-way comparisons, FDR correction becomes\n")
cat("more conservative. Primary 2-way comparisons are recommended for reporting.\n\n")

# ============================================================================
# PART 6: SENSITIVITY ANALYSES
# ============================================================================

cat("============================================================================\n")
cat("PART 6: SENSITIVITY ANALYSES\n")
cat("============================================================================\n\n")

# Leave-one-out analyses for COVID subgroups
cat("LEAVE-ONE-OUT SENSITIVITY ANALYSIS: COVID SUBGROUPS\n")
cat("----------------------------------------------------\n\n")

cat("Testing stability of COVID subgroup estimates:\n\n")

# Pre-COVID leave-one-out
cat("Pre-COVID studies:\n")
pre_covid_data <- data_final[data_final$COVID_status_combined == "Pre-COVID", ]
pre_estimates <- numeric(nrow(pre_covid_data))

for(i in 1:nrow(pre_covid_data)) {
  loo_data <- pre_covid_data[-i, ]
  loo_res <- rma(yi, vi, data = loo_data, method = "REML", test = "knha")
  loo_est <- predict(loo_res, transf = transf.ilogit)
  pre_estimates[i] <- loo_est$pred * 100
  cat(sprintf("  Removing %s: %.1f%%\n", 
              pre_covid_data$clean_authors[i], 
              loo_est$pred * 100))
}

pre_range <- range(pre_estimates)
cat(sprintf("Range: %.1f%% to %.1f%% (difference = %.1f%%)\n\n", 
            pre_range[1], pre_range[2], diff(pre_range)))

# During/Post-COVID leave-one-out
cat("During/Post-COVID studies:\n")
during_post_covid_data <- data_final[data_final$COVID_status_combined == "During/Post-COVID", ]
during_post_estimates <- numeric(nrow(during_post_covid_data))

for(i in 1:nrow(during_post_covid_data)) {
  loo_data <- during_post_covid_data[-i, ]
  loo_res <- rma(yi, vi, data = loo_data, method = "REML", test = "knha")
  loo_est <- predict(loo_res, transf = transf.ilogit)
  during_post_estimates[i] <- loo_est$pred * 100
  cat(sprintf("  Removing %s: %.1f%%\n", 
              during_post_covid_data$clean_authors[i], 
              loo_est$pred * 100))
}

during_post_range <- range(during_post_estimates)
cat(sprintf("Range: %.1f%% to %.1f%% (difference = %.1f%%)\n\n", 
            during_post_range[1], during_post_range[2], diff(during_post_range)))

# Leave-one-out for income subgroups
cat("LEAVE-ONE-OUT SENSITIVITY ANALYSIS: INCOME SUBGROUPS\n")
cat("-----------------------------------------------------\n\n")

cat("Testing stability of income level subgroup estimates:\n\n")

# High income leave-one-out
cat("High income studies:\n")
high_estimates <- numeric(nrow(high_data))

for(i in 1:nrow(high_data)) {
  loo_data <- high_data[-i, ]
  loo_res <- rma(yi, vi, data = loo_data, method = "REML", test = "knha")
  loo_est <- predict(loo_res, transf = transf.ilogit)
  high_estimates[i] <- loo_est$pred * 100
  cat(sprintf("  Removing %s: %.1f%%\n", 
              high_data$clean_authors[i], 
              loo_est$pred * 100))
}

high_range <- range(high_estimates)
cat(sprintf("Range: %.1f%% to %.1f%% (difference = %.1f%%)\n\n", 
            high_range[1], high_range[2], diff(high_range)))

# Medium income leave-one-out
cat("Medium income studies:\n")
medium_estimates <- numeric(nrow(medium_data))

for(i in 1:nrow(medium_data)) {
  loo_data <- medium_data[-i, ]
  loo_res <- rma(yi, vi, data = loo_data, method = "REML", test = "knha")
  loo_est <- predict(loo_res, transf = transf.ilogit)
  medium_estimates[i] <- loo_est$pred * 100
  cat(sprintf("  Removing %s: %.1f%%\n", 
              medium_data$clean_authors[i], 
              loo_est$pred * 100))
}

medium_range <- range(medium_estimates)
cat(sprintf("Range: %.1f%% to %.1f%% (difference = %.1f%%)\n\n", 
            medium_range[1], medium_range[2], diff(medium_range)))

cat("SUMMARY OF SENSITIVITY ANALYSES:\n")
cat("--------------------------------\n")
cat("Range of estimates when removing one study:\n")
cat(sprintf("Pre-COVID: %.1f%% to %.1f%% (range = %.1f%%)\n", 
            pre_range[1], pre_range[2], diff(pre_range)))
cat(sprintf("During/Post-COVID: %.1f%% to %.1f%% (range = %.1f%%)\n", 
            during_post_range[1], during_post_range[2], diff(during_post_range)))
cat(sprintf("High income: %.1f%% to %.1f%% (range = %.1f%%)\n", 
            high_range[1], high_range[2], diff(high_range)))
cat(sprintf("Medium income: %.1f%% to %.1f%% (range = %.1f%%)\n\n", 
            medium_range[1], medium_range[2], diff(medium_range)))

# ============================================================================
# PART 7: MULTIVARIABLE META-REGRESSION (EXPLORATORY)
# ============================================================================

cat("============================================================================\n")
cat("PART 7: MULTIVARIABLE META-REGRESSION\n")
cat("============================================================================\n\n")

cat("EXPLORATORY ANALYSIS: Multiple predictors simultaneously\n")
cat("---------------------------------------------------------\n\n")

# Model 1: COVID + MBI
cat("Model 1: COVID timing + Burnout measurement tool\n")
multi_reg_1 <- rma(yi, vi, mods = ~ COVID_status_combined + mbi_used, 
                   data = data_final, method = "REML", test = "knha")
print(multi_reg_1)
cat("R² =", sprintf("%.1f%%\n\n", multi_reg_1$R2))

# Model 2: Income + MBI (for income-available subset)
cat("Model 2: Income level + Burnout measurement tool\n")
multi_reg_2 <- rma(yi, vi, mods = ~ country_income_level_clean + mbi_used, 
                   data = data_income, method = "REML", test = "knha")
print(multi_reg_2)
cat("R² =", sprintf("%.1f%%\n\n", multi_reg_2$R2))

# Model 3: COVID + Study Quality
cat("Model 3: COVID timing + Study quality\n")
multi_reg_3 <- rma(yi, vi, mods = ~ COVID_status_combined + newcastle_ottawa_scale_nos_total, 
                   data = data_final, method = "REML", test = "knha")
print(multi_reg_3)
cat("R² =", sprintf("%.1f%%\n\n", multi_reg_3$R2))

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
# PART 8: COMPREHENSIVE SUMMARY
# ============================================================================

cat("============================================================================\n")
cat("PART 8: COMPREHENSIVE SUMMARY OF FINDINGS\n")
cat("============================================================================\n\n")

cat("OVERALL META-ANALYSIS:\n")
cat("----------------------\n")
cat(sprintf("• Pooled burnout prevalence: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n",
            overall_results$pred * 100, 
            overall_results$ci.lb * 100, 
            overall_results$ci.ub * 100))
cat(sprintf("• Total studies: %d (n = %d participants)\n", 
            nrow(data_final), 
            sum(data_final$number_of_participants)))
cat(sprintf("• Heterogeneity: I² = %.1f%%, τ² = %.2f\n", 
            res_overall$I2, res_overall$tau2))
cat(sprintf("• Prediction interval: %.1f%% to %.1f%%\n\n", 
            overall_results$pi.lb * 100, 
            overall_results$pi.ub * 100))

cat("COVID TIMING ANALYSIS:\n")
cat("----------------------\n")
cat("Two-way comparison (Pre-COVID vs During/Post-COVID):\n")
cat(sprintf("• Pre-COVID (k=%d): %.1f%% (95%% CI: %.1f%% to %.1f%%)\n",
            covid_combined_results_list[["Pre-COVID"]]$n,
            covid_combined_results_list[["Pre-COVID"]]$prev,
            covid_combined_results_list[["Pre-COVID"]]$ci_lb,
            covid_combined_results_list[["Pre-COVID"]]$ci_ub))
cat(sprintf("• During/Post-COVID (k=%d): %.1f%% (95%% CI: %.1f%% to %.1f%%)\n",
            covid_combined_results_list[["During/Post-COVID"]]$n,
            covid_combined_results_list[["During/Post-COVID"]]$prev,
            covid_combined_results_list[["During/Post-COVID"]]$ci_lb,
            covid_combined_results_list[["During/Post-COVID"]]$ci_ub))
cat(sprintf("• Difference: %.1f%% (%s)\n",
            abs(diff_covid),
            ifelse(diff_covid > 0, "Pre-COVID higher", "During/Post-COVID higher")))
cat(sprintf("• Test for subgroup difference: QM = %.2f, p = %.3f (unadjusted)\n",
            subgroup_covid_2way$QM, subgroup_covid_2way$QMp))
cat(sprintf("• Variance explained: R² = %.1f%%\n",
            subgroup_covid_2way$R2))
cat(sprintf("• Statistical significance (unadjusted): %s\n",
            ifelse(subgroup_covid_2way$QMp < 0.05, 
                   "SIGNIFICANT (p < 0.05)", 
                   "NOT significant (p ≥ 0.05)")))
cat(sprintf("• Statistical significance (FDR-adjusted): %s\n\n",
            ifelse(primary_tests$FDR_sig[1] == "Yes", 
                   "SIGNIFICANT (q < 0.05)", 
                   "NOT significant (q ≥ 0.05)")))

cat("COUNTRY INCOME LEVEL ANALYSIS:\n")
cat("------------------------------\n")
cat("High vs Medium income comparison (n = 10 studies):\n")
cat(sprintf("• High income (k=%d): %.1f%% (95%% CI: %.1f%% to %.1f%%)\n",
            income_results_list[["high"]]$n,
            income_results_list[["high"]]$prev,
            income_results_list[["high"]]$ci_lb,
            income_results_list[["high"]]$ci_ub))
cat(sprintf("• Medium income (k=%d): %.1f%% (95%% CI: %.1f%% to %.1f%%)\n",
            income_results_list[["medium"]]$n,
            income_results_list[["medium"]]$prev,
            income_results_list[["medium"]]$ci_lb,
            income_results_list[["medium"]]$ci_ub))
cat(sprintf("• Difference: %.1f%%\n", abs(diff_income)))
cat(sprintf("• Test for subgroup difference: QM = %.2f, p = %.3f (unadjusted)\n",
            subgroup_income_2way$QM, subgroup_income_2way$QMp))
cat(sprintf("• Variance explained: R² = %.1f%%\n",
            subgroup_income_2way$R2))
cat(sprintf("• Statistical significance (unadjusted): %s\n",
            ifelse(subgroup_income_2way$QMp < 0.05, 
                   "SIGNIFICANT (p < 0.05)", 
                   "NOT significant (p ≥ 0.05)")))
cat(sprintf("• Statistical significance (FDR-adjusted): %s\n\n",
            ifelse(primary_tests$FDR_sig[2] == "Yes", 
                   "SIGNIFICANT (q < 0.05)", 
                   "NOT significant (q ≥ 0.05)")))

cat("Note: Low income category (k=1, Pakistan) and multi-country studies (k=4)\n")
cat("      excluded from formal income comparison.\n\n")

cat("KEY INTERPRETATION POINTS:\n")
cat("--------------------------\n")

# COVID interpretation
if(primary_tests$FDR_sig[1] == "Yes") {
  if(diff_covid > 0) {
    cat("• COVID TIMING: Burnout rates were SIGNIFICANTLY HIGHER before COVID (survives FDR correction)\n")
    cat("  This unexpected finding may reflect:\n")
    cat("  - Selection bias in post-pandemic studies\n")
    cat("  - Adaptation and resilience building during pandemic\n")
    cat("  - Changes in work practices or reporting\n")
    cat("  - Limited post-COVID data (k=3 studies)\n")
  } else {
    cat("• COVID TIMING: Burnout rates SIGNIFICANTLY INCREASED during/after COVID (survives FDR correction)\n")
    cat("  This aligns with expected pandemic-related stressors\n")
  }
} else if(subgroup_covid_2way$QMp < 0.05) {
  cat("• COVID TIMING: Nominally significant difference (p < 0.05) but does NOT survive FDR correction\n")
  cat("  This suggests the finding should be interpreted cautiously and may reflect:\n")
  cat("  - Type I error (false positive) due to multiple comparisons\n")
  cat("  - Insufficient statistical power (limited post-COVID studies, k=3)\n")
  cat("  - High heterogeneity masking true effects\n")
} else {
  cat("• COVID TIMING: No significant difference detected between pre and during/post COVID\n")
  cat("  This suggests either:\n")
  cat("  - Burnout rates remained stable across COVID periods\n")
  cat("  - Insufficient statistical power (limited post-COVID studies, k=3)\n")
  cat("  - High heterogeneity masking true differences\n")
}
cat("\n")

# Income interpretation
if(primary_tests$FDR_sig[2] == "Yes") {
  if(diff_income > 0) {
    cat("• INCOME LEVEL: Burnout SIGNIFICANTLY HIGHER in high-income countries (survives FDR correction)\n")
    cat("  Possible explanations:\n")
    cat("  - Higher work demands and expectations\n")
    cat("  - Greater administrative burden\n")
    cat("  - Different work-life balance norms\n")
  } else {
    cat("• INCOME LEVEL: Burnout SIGNIFICANTLY HIGHER in medium-income countries (survives FDR correction)\n")
    cat("  Possible explanations:\n")
    cat("  - Resource limitations\n")
    cat("  - Workforce shortages\n")
    cat("  - Healthcare system strain\n")
  }
} else if(subgroup_income_2way$QMp < 0.05) {
  cat("• INCOME LEVEL: Nominally significant difference (p < 0.05) but does NOT survive FDR correction\n")
  cat("  This finding should be interpreted cautiously as it may reflect Type I error\n")
  cat("  Limited sample (k=10) reduces power to detect true differences\n")
} else {
  cat("• INCOME LEVEL: No significant difference between high and medium income countries\n")
  cat("  This suggests burnout is a universal challenge across healthcare systems\n")
  cat("  Limited sample (k=10) may reduce power to detect differences\n")
}
cat("\n")

cat("METHODOLOGICAL CONSIDERATIONS:\n")
cat("-------------------------------\n")
cat("• COVID analysis limited by small post-pandemic sample (k=3)\n")
cat("• Income analysis excludes 5 studies (1 low-income, 4 multi-country)\n")
cat("• High residual heterogeneity indicates unmeasured moderators\n")
cat("• Leave-one-out analyses show generally stable estimates\n")
cat("• FDR correction accounts for multiple comparisons\n\n")

cat("RECOMMENDATIONS FOR REPORTING:\n")
cat("------------------------------\n")
cat("1. Report both unadjusted and FDR-adjusted p-values for transparency\n")
cat("2. Report COVID findings with caution given limited post-pandemic data\n")
cat("3. Acknowledge excluded studies in income analysis\n")
cat("4. Emphasize need for more post-COVID longitudinal research\n")
cat("5. Consider income analysis exploratory given small subgroups\n")
cat("6. Explore other potential moderators (specialty, training level, etc.)\n\n")

cat("FILES CREATED:\n")
cat("--------------\n")
cat("1. forest_plot_covid_subgroups.pdf - COVID timing forest plot (NO overlapping labels)\n")
cat("2. forest_plot_income_subgroups.pdf - Income level forest plot (NO overlapping labels)\n")
cat("3. meta_analysis_workspace.RData - Complete workspace for further analysis\n\n")

cat("=============================================================================\n")
cat("META-ANALYSIS COMPLETE\n")
cat("=============================================================================\n")

# Save workspace for further analysis
save.image(file = "outputs/meta_analysis_workspace.RData")
cat("\nWorkspace saved: meta_analysis_workspace.RData\n")
cat("Load with: load('meta_analysis_workspace.RData')\n")
