setwd("~/Desktop/Projects/EnvPredict/")
library(superheat)
library(wesanderson)
cortable <- fread("Cortabletopivot.tsv")

cortable2 <- cortable %>%
  mutate(PC_Species = paste0(Var1, "_", str_replace(Species, " ", "_"))) %>%
  dplyr::select(PC_Species, Var2, R_Square)

mat <- cortable2 %>%
  pivot_wider(names_from = Var2, values_from = R_Square) %>%
  column_to_rownames("PC_Species")

cormat_clean <- mat %>% mutate(todrop = rownames(mat)) %>%  mutate(todrop = str_replace(todrop, "PC4.*", "PC4")) %>%  filter(!todrop %in% "PC4") %>%  dplyr::select(-todrop)

superheat(X = cormat_clean, pretty.order.rows = T, pretty.order.cols = F, row.dendrogram = T,col.dendrogram = T, scale = F, heat.lim = c(0, 1), heat.pal.values = c(0, 0.25, 0.5, 1), force.left.label = TRUE, heat.pal = wes_palette("Zissou1"))
