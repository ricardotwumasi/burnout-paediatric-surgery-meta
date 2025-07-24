# Load required libraries with checks
required_packages <- c("tidyverse")
lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

# Load required libraries
library(tidyverse)

# Read data
burnout_data <- read.csv("surgical_specialties_comparison.csv")

# Data preprocessing
burnout_data <- burnout_data %>%
  mutate(
    # Convert percentage to proportion for plotting
    burnout_prevalence = burnout_prevalence_percent / 100,
    ci_lower = CI_Lower / 100,
    ci_upper = CI_Upper / 100,
    
    # Create indicator for studies with CIs
    has_ci = !is.na(CI_Lower) & !is.na(CI_Upper),
    
    # Create combined specialty and author label
    specialty_label = paste0(specialty, "\n(", authors, ", ", publication_year, ")"),
    
    # Create a shorter label for better readability
    short_label = paste0(specialty, "\n(", authors, ")")
  )

# Create forest plot
forest_plot <- ggplot(burnout_data, aes(y = reorder(short_label, burnout_prevalence))) +
  # Add points and error bars
  geom_point(aes(x = burnout_prevalence, shape = has_ci, color = has_ci), size = 3) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper, x = burnout_prevalence),
    height = 0.2,
    color = "black",
    data = burnout_data %>% filter(has_ci)
  ) +
  
  # Add value labels on the right
  geom_text(aes(x = burnout_prevalence, 
                label = paste0(round(burnout_prevalence_percent, 1), "%")),
            hjust = -0.3, size = 3, fontface = "bold") +
  
  # Labels and formatting
  labs(
    title = "Burnout Prevalence Across Surgical Specialties",
    subtitle = "Comparison of Systematic Reviews and Meta-Analyses",
    x = "Burnout Prevalence",
    y = NULL,
    caption = "Note: Studies without confidence intervals shown as triangles"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(size = 9, lineheight = 0.9),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 10),
    plot.margin = margin(5.5, 30, 5.5, 5.5, "pt")  # Extra right margin for labels
  ) +
  
  scale_x_continuous(
    limits = c(0, 0.8), 
    breaks = seq(0, 0.8, 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  scale_shape_manual(
    values = c(17, 19),  # Triangle for no CI, circle for CI
    labels = c("No confidence interval", "Confidence interval reported"),
    name = "Study precision"
  ) +
  
  scale_color_manual(
    values = c("red", "blue"),  # Red for no CI, blue for CI
    labels = c("No confidence interval", "Confidence interval reported"),
    name = "Study precision"
  ) +
  
  guides(
    shape = guide_legend(override.aes = list(color = c("red", "blue"))),
    color = "none"  # Hide the color legend since shape legend shows both
  )

# Print the plot
print(forest_plot)

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Number of studies:", nrow(burnout_data), "\n")
cat("Studies with confidence intervals:", sum(burnout_data$has_ci), "\n")
cat("Studies without confidence intervals:", sum(!burnout_data$has_ci), "\n\n")

# Show which studies lack confidence intervals
studies_no_ci <- burnout_data %>% 
  filter(!has_ci) %>% 
  select(authors, specialty, burnout_prevalence_percent)

if(nrow(studies_no_ci) > 0) {
  cat("Studies without confidence intervals:\n")
  print(studies_no_ci)
}

cat("\nBurnout prevalence range:", 
    round(min(burnout_data$burnout_prevalence_percent), 1), "% -", 
    round(max(burnout_data$burnout_prevalence_percent), 1), "%\n")

cat("Mean burnout prevalence:", 
    round(mean(burnout_data$burnout_prevalence_percent), 1), "%\n")

cat("Median burnout prevalence:", 
    round(median(burnout_data$burnout_prevalence_percent), 1), "%\n")

# Save the plot (uncomment to use)
# ggsave("surgical_specialties_burnout_comparison.pdf", forest_plot, 
#        width = 12, height = 8, units = "inches", dpi = 300)
# ggsave("surgical_specialties_burnout_comparison.png", forest_plot, 
#        width = 12, height = 8, units = "inches", dpi = 300)