library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)

# ------------------------------------------------------------------------------

# Read the data 
raw_df <-fread( "~/Desktop/Projects/EnvPredict/EnvPredictSubmission/Supplementary_Table_Simulation_Scenario_Summary.tsv" )

# ------------------------------------------------------------------------------
# Reshape and Standardize Data
# ------------------------------------------------------------------------------
plot_data <- bind_rows(
  raw_df %>% select(conf_level, genetic_arch, FPR = med_FPR_RDA, FNR = med_FNR_RDA) %>% mutate(Method = "RDA"),
  raw_df %>% select(conf_level, genetic_arch, FPR = med_FPR_pRDA, FNR = med_FNR_pRDA) %>% mutate(Method = "pRDA"),
  raw_df %>% select(conf_level, genetic_arch, FPR = med_FPR_LFMM, FNR = med_FNR_LFMM) %>% mutate(Method = "LFMM")
) %>%
  mutate(
    # Clean up any potential string truncation/case differences
    conf_clean = case_when(
      str_detect(conf_level, "(?i)extreme")  ~ "Extreme Confounding",
      str_detect(conf_level, "(?i)high")     ~ "High Confounding",
      str_detect(conf_level, "(?i)moderate") ~ "Moderate Confounding",
      str_detect(conf_level, "(?i)low")      ~ "Low Confounding",
      TRUE ~ as.character(conf_level)
    ),
    arch_clean = case_when(
      str_detect(genetic_arch, "(?i)oligogenic") ~ "Oligogenic",
      str_detect(genetic_arch, "(?i)moderate")   ~ "Moderately Polygenic",
      str_detect(genetic_arch, "(?i)polygenic")  ~ "Polygenic",
      TRUE ~ as.character(genetic_arch)
    ),
    
    # Set factor levels for consistent legend ordering
    conf_clean = factor(
      conf_clean, 
      levels = c("Extreme Confounding", "High Confounding", "Moderate Confounding", "Low Confounding")
    ),
    arch_clean = factor(
      arch_clean, 
      levels = c("Oligogenic", "Moderately Polygenic", "Polygenic")
    ),
    
    # Calculate True Positive Rate (Sensitivity)
    TPR = 1 - FNR,
    
    # Combine Method and Confounding level for color mapping
    Method_Conf = factor(
      paste(Method, conf_clean, sep = " - "),
      levels = c(
        "RDA - Low Confounding", "RDA - Moderate Confounding", "RDA - High Confounding", "RDA - Extreme Confounding",
        "pRDA - Low Confounding", "pRDA - Moderate Confounding", "pRDA - High Confounding", "pRDA - Extreme Confounding",
        "LFMM - Low Confounding", "LFMM - Moderate Confounding", "LFMM - High Confounding", "LFMM - Extreme Confounding"
      )
    )
  )

# ------------------------------------------------------------------------------
# 3. Custom Color Palette Mapping
# ------------------------------------------------------------------------------
color_palette <- c(
  
  # RDA (Sea Blue -> Navy)
  "RDA - Low Confounding"      = "#5DA5DA",
  "RDA - Moderate Confounding" = "#2B83BA",
  "RDA - High Confounding"     = "#1D5F8A",
  "RDA - Extreme Confounding"  = "#08306B",
  
  # pRDA (Gold -> Burnished Gold)
  "pRDA - Low Confounding"      = "#E6B800",
  "pRDA - Moderate Confounding" = "#D39B00",
  "pRDA - High Confounding"     = "#B8860B",
  "pRDA - Extreme Confounding"  = "#8C6500",
  
  # LFMM (Red -> Deep Crimson)
  "LFMM - Low Confounding"      = "#EF3B2C",
  "LFMM - Moderate Confounding" = "#D7301F",
  "LFMM - High Confounding"     = "#B30000",
  "LFMM - Extreme Confounding"  = "#67000D"
)
# ------------------------------------------------------------------------------
# 4. Generate Plot
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 4. Generate Plot
# ------------------------------------------------------------------------------
ggplot(plot_data, aes(x = FNR, y = FPR, color = Method_Conf, shape = arch_clean)) +
  geom_point(
    position = position_jitter(width = 0.003, height = 0.003, seed = 123),
    size = 6,
    alpha = 0.95
  ) +
  scale_color_manual(values = color_palette, name = "Method & Confounding") +
  scale_shape_manual(
    values = c(16, 17, 15),
    name = "Genetic Architecture"
  ) +
  coord_cartesian(
    xlim = c(0, 1.02),
    ylim = c(0, 1.02)
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Median False Negative Rate (FNR)",
    y = "Median False Positive Rate (FPR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  ) +
  theme_classic()

