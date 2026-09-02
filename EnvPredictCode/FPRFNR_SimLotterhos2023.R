#setwd
setwd("~/Desktop/Projects/EnvPredict/Lotterhos2023Simulation/")
library(tidyverse)
library(ggridges)
library(broom.mixed)
library(dplyr)

#read in simulations
sim <- read_csv("summary_20220428_20220726(1).csv")

#make an index of potential for confounding based on clinality of structure and env
sim2 <- sim %>%
  mutate(confounding_index = abs(cor_PC1_temp)) %>%
  rowid_to_column("sim_id")

#make a polygenicity index based on their median temperature correlation and number of causal loci
sim2 <- sim2 %>%
  mutate(
    # total number of causal loci
    n_causal = num_causal_temp,
    
    # typical effect size proxy
    median_effect = abs(median_causal_temp_cor),
    
    polygenicity_index =
      n_causal / (median_effect + 1e-6)
  )


# Calculate thresholds based on  Polygenicity

poly_hi <- quantile(sim2$polygenicity_index, 0.66, na.rm = TRUE)
poly_lo <- quantile(sim2$polygenicity_index, 0.33, na.rm = TRUE)

# absolute correlation/confounding thresholds
value_extreme = 0.80  # Extreme confounding
value_hi      = 0.60  # High confounding (0.6 to 0.8)
value_lo      = 0.30  # Low confounding

# Classify with confounding level
sim_classified <- sim2 %>%
  mutate(
    # 1. Define Confounding Levels (Rows)
    conf_level = case_when(
      confounding_index >= 0.80 ~ "Extreme Confounding",
      confounding_index >= 0.60 ~ "High Confounding",
      confounding_index <= 0.30 ~ "Low Confounding",
      TRUE                      ~ "Moderate Confounding"
    ),
    # Set levels for Top-to-Bottom order
    conf_level = factor(conf_level, levels = c("Extreme Confounding", "High Confounding", "Moderate Confounding", "Low Confounding")),
    
    # 2. Define Genetic Architecture (Columns)
    genetic_arch = case_when(
      polygenicity_index >= poly_hi ~ "Polygenic",
      polygenicity_index <= poly_lo ~ "Oligogenic",
      TRUE                          ~ "Moderately polygenic"
    ),
    # Set levels for Left-to-Right order
    genetic_arch = factor(genetic_arch, levels = c("Oligogenic", "Moderately polygenic", "Polygenic"))
  )

#estimate false negative rates
sim_classified <- sim_classified %>%
  mutate(
    # 1. LEA temp false negative rate (FNR)
    LEA3.2_lfmm2_FNR_temp = 1 - LEA3_2_lfmm2_TPR_temp,
    
    # 2. RDA with structure correction false negative rate
    RDA_FNR_corr = 1 - RDA_TPR_corr,
    
    # 3. RDA without structure correction false negative rate
    RDA_FNR = 1 - RDA_TPR
  )



#pivot long for plotting FNR per scenario 
sim_long <- sim_classified %>%
  select(sim_id, conf_level, genetic_arch, RDA_FNR, RDA_FNR_corr, LEA3.2_lfmm2_FNR_temp) %>%
  pivot_longer(
    cols = starts_with(c("RDA", "LEA")),
    names_to = "Method",
    values_to = "FNR"
  )



#now with FDR
sim_long_FDR <- sim_classified %>%
  select(sim_id, conf_level, genetic_arch, RDA_FDR_allSNPs, RDA_FDR_allSNPs_corr, LEA3_2_lfmm2_FDR_allSNPs_temp) %>%
  pivot_longer(
    cols = starts_with(c("RDA", "LEA")),
    names_to = "Method",
    values_to = "FPR"
  )


ggplot(sim_long, aes(x = FNR, fill = Method)) +
  geom_density(alpha = 0.5) +
  # facet_wrap allows 'free' scales for EVERY individual panel
  facet_wrap(conf_level ~ genetic_arch, 
             ncol = 3, 
             scales = "free_y") + 
  theme_minimal() +
  labs(
    title = "Shift in FNR across scenarios",
    x = "False Negative Rate (1.0 = total failure)",
    y = "Density"
  ) +
  scale_fill_manual(values = c(
    "RDA_FNR" = "#A6CEE3", 
    "RDA_FNR_corr" = "#FDBF6F", 
    "LEA3.2_lfmm2_FNR_temp" = "#E31A1C"
  )) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(1, "lines") # Adds a bit of breathing room between plots
  ) + theme_classic()



ggplot(sim_long_FDR, aes(x = FPR, fill = Method)) +
  geom_density(alpha = 0.5) +
  # facet_wrap allows 'free' scales for EVERY individual panel
  facet_wrap(conf_level ~ genetic_arch, 
             ncol = 3, 
             scales = "free_y") + 
  theme_minimal() +
  labs(
    title = "Shift in FNR across scenarios",
    x = "False Positive Rate",
    y = "Density"
  ) +
  scale_fill_manual(values = c(
    "RDA_FDR_allSNPs" = "#A6CEE3", 
    "RDA_FDR_allSNPs_corr" = "#FDBF6F", 
    "LEA3_2_lfmm2_FDR_allSNPs_temp" = "#E31A1C"
  )) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(1, "lines") # Adds a bit of breathing room between plots
  ) + theme_classic()


# Assess separation of parameters across scenarios
scenario_summary <- sim_classified %>%
  group_by(conf_level, genetic_arch) %>%
  summarize(
    n_simulations = n(),
    
    # --- Raw Median FNRs per Method ---
    med_FNR_RDA      = median(RDA_FNR, na.rm = TRUE),
    med_FNR_pRDA     = median(RDA_FNR_corr, na.rm = TRUE),
    med_FNR_LFMM     = median(LEA3.2_lfmm2_FNR_temp, na.rm = TRUE),
    
    # --- Scenario Parameters ---
    med_n_causal     = median(n_causal, na.rm = TRUE),
    med_effect       = median(median_effect, na.rm = TRUE),
    med_confounding  = median(confounding_index, na.rm = TRUE),
    
    .groups = "drop"
  )


scenario_summary_FDR <- sim_classified_FDR %>%
  group_by(conf_level, genetic_arch) %>%
  summarize(
    n_simulations = n(),
    
    # --- Raw Median FNRs per Method ---
    med_FPR_RDA      = median(RDA_FDR_allSNPs, na.rm = TRUE),
    med_FPR_pRDA     = median(RDA_FDR_allSNPs_corr, na.rm = TRUE),
    med_FPR_LFMM     = median(LEA3_2_lfmm2_FDR_allSNPs_temp, na.rm = TRUE),
    
    # --- Scenario Parameters ---
    med_n_causal     = median(n_causal, na.rm = TRUE),
    med_effect       = median(median_effect, na.rm = TRUE),
    med_confounding  = median(confounding_index, na.rm = TRUE),
    
    .groups = "drop"
  )
# Assess separation of parameters across scenarios


# View the table
print(scenario_summary)
write.table(scenario_summary, "Table_Scenario_Summary.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(scenario_summary_FDR, "Table_Scenario_Summary_FPR.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
      