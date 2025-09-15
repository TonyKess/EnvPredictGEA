setwd("~/Desktop/Projects/EnvPredict/")
dataset <- fread("Supplementary_Table_1_Allspecies_cortests.csv")
PCs <- paste0("PC", rep(1:4))
Comps <- paste0("Comp.", rep(1:4))
PCcomp <- c(PCs, Comps)

Envassoc <- dataset %>%
  group_by(Species) %>%
  filter(Var1 %in% PCs  | Var2 %in% PCs ) %>% 
    ungroup() 

Envassoc_pruned_bysig <- Envassoc %>%
  group_by(Species) %>%
  filter(P_Value < (0.05/160)) %>% 
  ungroup() 

Envassoc_pruned_bysig %>%
  arrange(desc(R_Square)) %>% group_by(Species) %>%  summarise(medianR2perspec = max(R_Square))

Envassoc_pruned_byassoc <- Envassoc_pruned_bysig %>%
  group_by(Species) %>%
  slice_max(order_by = R_Square, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(R_Square))


AllSpeciesR2RF <- fread("AllSpeciesResnFigTable.tsv")

AllSpeciesR2RF %>%  group_by(Dataset) %>%  slice_max(R2)


PCCor <- paste0("PCcor_OL_", c(100, 200, 300, 400, 500))
OL <-  paste0("OL_", c(100, 200, 300, 400, 500))

table(AllSpeciesR2RF$Dataset)

Specieslist <-  c("AtlanticSalmon", "BlackSpruce","CohoSalmon","KessCod2020","Layton2021Charr", "Mahony2020_Pinus",   "MVWScallop2024", "Schweizer_Wolves")

Getwilcox <- function(Species){
  Singlespecies <- AllSpeciesR2RF %>%  filter(Dataset %in% Species)
  OLset <- Singlespecies %>%  filter(Outlier_set %in% OL)
  PCcorOLset <- Singlespecies %>%  filter(Outlier_set %in% PCCor)
  WCtest <- wilcox.test(OLset$R2, PCcorOLset$R2)
  return(WCtest)
}

WilcoxResults <- map_dfr(Specieslist, function(sp) {
  test <- Getwilcox(sp)
  tibble(
    Species = sp,
    W = test$statistic,
    p_value = test$p.value
  )
})


WilcoxResults <- WilcoxResults %>%  mutate(Dataset = Species) %>%  select(-Species)

MeanR2OL <- AllSpeciesR2RF %>%  filter(Outlier_set %in% OL) %>%  group_by(Dataset) %>%  summarize(MeanR2OL = mean(R2))
SDR2OL <- AllSpeciesR2RF %>%  filter(Outlier_set %in% OL) %>%  group_by(Dataset) %>%  summarize(SDR2OL = sd(R2))
MeanR2PCcor <- AllSpeciesR2RF %>%  filter(Outlier_set %in% PCCor ) %>%  group_by(Dataset) %>%  summarize(MeanR2PCcor = mean(R2))
SDR2PCcor  <- AllSpeciesR2RF %>%  filter(Outlier_set %in% PCCor) %>%  group_by(Dataset) %>%  summarize(SDR2PCcor = sd(R2))

Table_2 <-  inner_join(MeanR2OL, MeanR2PCcor)  %>%  inner_join(WilcoxResults) %>%  arrange(desc(MeanR2OL))
Table_3 <- AllSpeciesR2RF %>%  group_by(Dataset)  %>% slice_max(R2) %>%  arrange(desc(R2)) %>%  select(-Run)


write.table(Table_2, "Table_2.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(Table_3, "Table_3.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
