library(tidyverse)
getwd()
#Load random forest env prediction data
data <- read_tsv("AllSpeciesRandomForestTable.tsv") %>%  filter(Variable %in% c("Comp.1", "Comp.2"))

#clean and factorize
# We need to create a specific order for the legend to match the blue/orange gradient
plot_data <- data %>%
  filter(Variable %in% c("Comp.1", "Comp.2")) %>%
  mutate(
    Species = str_replace_all(Dataset, "_", " "),
    # Ensure the factor order matches the visual gradient: 
    # Corrected 500 down to 100, then Uncorrected 100 up to 500
    Panel = factor(Outlier_set, levels = c(
      "PCcor_OL_500", "PCcor_OL_400", "PCcor_OL_300", "PCcor_OL_200", "PCcor_OL_100",
      "OL_100", "OL_200", "OL_300", "OL_400", "OL_500"
    ))
  )
blue_ramp <- colorRampPalette(c("lightblue", "blue"))(5)
orange_ramp <- colorRampPalette(c("bisque", "orange"))(5)

# Map them to the levels (Order: 500, 400, 300, 200, 100 for Blue; 100...500 for Orange)
# Note: blue_ramp[5] is the "bluest" (500), blue_ramp[1] is "lightblue" (100)
custom_colors <- c(
  "PCcor_OL_500" = blue_ramp[5], "PCcor_OL_400" = blue_ramp[4], 
  "PCcor_OL_300" = blue_ramp[3], "PCcor_OL_200" = blue_ramp[2], 
  "PCcor_OL_100" = blue_ramp[1],
  "OL_100" = orange_ramp[1], "OL_200" = orange_ramp[2], 
  "OL_300" = orange_ramp[3], "OL_400" = orange_ramp[4], 
  "OL_500" = orange_ramp[5]
)
# create the Plot
ggplot(plot_data, aes(x = R2, y = Species, color = Panel)) +
  # Use geom_jitter to prevent dots from overlapping perfectly
  geom_jitter(size = 1.5, height = 0.2) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Predictive Accuracy across Species and SNP Panels",
    subtitle = "Points represent individual simulations (Runs/Variables)",
    x = expression(Predictive~Accuracy~(R^2)),
    y = NULL,
    color = "SNP Panel"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90"),
    axis.text.y = element_text(face = "bold")
  ) + facet_wrap(~Variable) + theme_classic()
