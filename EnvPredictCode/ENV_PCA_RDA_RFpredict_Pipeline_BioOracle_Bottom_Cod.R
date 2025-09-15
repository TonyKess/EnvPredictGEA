#load packages
library(tidyverse)
library(data.table)
library(vegan)
library(pcadapt)
library(caret)
library(robust)
library(qvalue)
library(ranger)
library(patchwork)

#set environment
wd = "~/Desktop/Projects/EnvPredict/"
subproj = "KessCod2020"
setwd(paste0(wd, subproj))

system("mkdir resfigntable")

#bring in geo data
Env<- data.frame(fread(paste0(subproj, "_Env.tsv")))

#scale bioclim data for PCA
BioOracle2.bottom.var<- c("BO2_tempmax_bdmean", 
                           "BO2_tempmean_bdmean", 
                           "BO2_tempmin_bdmean",
                           "BO2_dissoxmax_bdmean", 
                           "BO2_dissoxmean_bdmean", 
                           "BO2_dissoxmin_bdmean", 
                           "BO2_salinitymax_bdmean", 
                           "BO2_salinitymean_bdmean", 
                           "BO2_salinitymin_bdmean") 


BioOracle2.bottom.sd <- paste0(BioOracle2.bottom.var, ".scale")

Env<- Env %>%
  mutate(across(all_of(BioOracle2.bottom.var), ~ as.numeric(scale(.)), .names = "{.col}.scale"))


#PCA bioclim variables

Env_SDmatrix <-Env  %>% select(all_of(BioOracle2.bottom.sd))

Env_SDmatrix <- as.matrix(Env_SDmatrix)

Env_PCs <- princomp(Env_SDmatrix) 

EnvPCvar <- summary(Env_PCs)

#get loadings of bioclim variables and identify top variable for plotting

loadings  <- data.frame(unclass(Env_PCs$loadings)) %>% 
  mutate(BIOvar = rownames(.), Study = subproj)

topBIOvar <- loadings %>% 
  slice_max(abs(Comp.1)) %>%  
  select(BIOvar)

topBIOvarname <- topBIOvar$BIOvar

#combine env and PC env data
EnvPC_Biovars <- bind_cols(Env,Env_PCs$scores[,1:5])

#plot
ENVPC_C1C2_topbiovar <- ggplot() + 
  geom_point(data = EnvPC_Biovars, aes(x = Comp.1, y = Comp.2, colour = .data[[topBIOvarname]])) + 
  scale_color_gradient(low = "blue", high = "red") + 
  labs(color = topBIOvarname) + 
  theme_classic()
#save plot
ggsave(plot = ENVPC_C1C2_topbiovar, filename = paste0("resfigntable/", subproj, "ENVPC_C1C2_topbiovar.pdf"), 
       device = pdf, 
       height = 8, 
       width = 13)

 
#add kmeans clusters for splitting test and training data
kmeans_addclusters_func <- function(n_clusters, PCA){
  set.seed(123)
  kmeanrun <-  kmeans(PCA, centers = n_clusters, nstart = 1000)
  plinkPCAdataclustered <- bind_cols(PCA, cluster= kmeanrun$cluster)
  return(plinkPCAdataclustered)}


kmeans_evalfunc <- function(n_clusters,PCA){
  
  # Initialize total within sum of squares error: wss
  wss <- numeric(n_clusters)
  
  set.seed(123)
  
  # Look over 1 to n possible clusters
  for (i in 1:n_clusters) {
    # Fit the model: km.out
    km.out <- kmeans(PCA, centers = i, nstart = 1000,iter.max = 100, algorithm = "MacQueen")
    # Save the within cluster sum of squares
    wss[i] <- km.out$tot.withinss
  }
  
  # Produce a scree plot
  wss_df <- tibble(clusters = 1:n_clusters, wss = wss)
  
  scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
    geom_point(size = 4)+
    geom_line() +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
    xlab('Number of clusters')
  scree_plot
}

kmeans_evalfunc(10, Env_PCs$scores) + theme_classic()
Env_PCs_clusters <-  kmeans_addclusters_func(n_clusters = 5, PCA = Env_PCs$scores)
Env_PCs_clusters  <- bind_cols(EnvPC_Biovars, cluster = Env_PCs_clusters$cluster)

#plot
ENVPC_C1C2_clusterk5 <- ggplot() + 
  geom_point(data = Env_PCs_clusters, aes(x = Comp.1, y = Comp.2, colour = as.factor(cluster))) + 
  theme_classic()


#now pcadapt
pcadapt.object <-  read.pcadapt(paste0(subproj, ".bed"), type = "bed")
Fam  <- fread(paste0(subproj, ".fam")) %>%  
  mutate(FID = V1, IID = V2) %>%
  select(FID, IID)
Env_PCs_clusters <- inner_join(Fam, Env_PCs_clusters)
PCA <- pcadapt(input = pcadapt.object , K = 5)
PCAScores <- PCA$scores
colnames(PCAScores) <- paste0("PC", rep(1:5))
PCA$singular.values^2
plot(PCA, option = "screeplot") + theme_classic()

GenoPCA_EnvPCA<- data.frame(bind_cols(Env_PCs_clusters, PCAScores))


GenoPCA_12 <- ggplot() + 
  geom_point(data = GenoPCA_EnvPCA, aes(x = PC1, y = PC2, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")



GenoPCA_23 <- ggplot() + 
  geom_point(data = GenoPCA_EnvPCA, aes(x = PC2, y = PC3, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")


GenoPCA_34 <- ggplot() + geom_point(data = GenoPCA_EnvPCA, aes(x = PC3, y = PC4, colour = Comp.1)) + theme_classic() +
  scale_color_gradient(low = "blue", high = "red") + theme_classic()


# Assuming GenoPCA_12, GenoPCA_23, GenoPCA_34 are ggplot objects
GenoPCAplots <- GenoPCA_12 + GenoPCA_23 + GenoPCA_34 + plot_layout(ncol = 1)


#save plot
ggsave(plot = GenoPCAplots, filename = paste0("resfigntable/", subproj, "GenoPCAplots.pdf"), 
       device = pdf, 
       height = 8, 
       width = 13)


#Check correlations between PCA and environment
# Create an empty dataframe to store the results
results <- data.frame()

# Loop over the PC variables
for(i in 1:4){
  # Loop over the Comp variables
  for(j in 1:3){
    # Get the PC and Comp column names
    pc_name <- paste0("PC", i)
    comp_name <- paste0("Comp.", j)
    
    # Run the correlation test
    test <- cor.test(GenoPCA_EnvPCA[[pc_name]], GenoPCA_EnvPCA[[comp_name]])
    
    # Store the results in the dataframe
    results <- rbind(results, data.frame(Var1 = pc_name, Var2 = comp_name, R_Square = test$estimate^2, P_Value = test$p.value))
  }
} 

# Loop over the Lat/Lon variables
for(latlon in c("Lat", "Lon")){
  # Loop over the PC variables
  for(i in 1:4){
    # Get the PC column names
    pc_name <- paste0("PC", i)
    
    # Run the correlation test
    test <- cor.test(GenoPCA_EnvPCA[[latlon]], GenoPCA_EnvPCA[[pc_name]])
    
    # Store the results in the dataframe
    results <- rbind(results, data.frame(Var1 = latlon, Var2 = pc_name, R_Square = test$estimate^2, P_Value = test$p.value))
  }
  
  # Loop over the Comp variables
  for(j in 1:3){
    # Get the Comp column names
    comp_name <- paste0("Comp.", j)
    
    # Run the correlation test
    test <- cor.test(GenoPCA_EnvPCA[[latlon]], GenoPCA_EnvPCA[[comp_name]])
    
    # Store the results in the dataframe
    results <- rbind(results, data.frame(Var1 = latlon, Var2 = comp_name, R_Square = test$estimate^2, P_Value = test$p.value))
  }
}

# Print the results
Cortests <- results 

write.table(Cortests, paste0("resfigntable/", subproj, "_cortests.tsv"), col.names = T, 
            row.names = F, quote = F, sep = "\t")


#RDA for general pop structure
system(paste0("source ~/.zshrc; conda activate plink; plink --bfile ", subproj, 
              " --recodeA --allow-extra-chr --maf 0.01 --out ", subproj))

genos.dose <- fread(paste0(subproj, ".raw"), sep = " ") %>%  select (-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
genos.dose.imp<- apply(genos.dose, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

#Run RDAs
Full.rda <- rda(genos.dose.imp ~ GenoPCA_EnvPCA$Comp.1 + GenoPCA_EnvPCA$Comp.2  + GenoPCA_EnvPCA$Comp.3  ,scale=T )
Full.PCcor.rda <- rda(genos.dose.imp ~  GenoPCA_EnvPCA$Comp.1 + GenoPCA_EnvPCA$Comp.2  + GenoPCA_EnvPCA$Comp.3  + Condition(GenoPCA_EnvPCA$PC1 + GenoPCA_EnvPCA$PC2+ GenoPCA_EnvPCA$PC3 ) ,scale=T )
Full.Geocor.rda <- rda(genos.dose.imp ~  GenoPCA_EnvPCA$Comp.1 + GenoPCA_EnvPCA$Comp.2  + GenoPCA_EnvPCA$Comp.3  + Condition(GenoPCA_EnvPCA$Lat + GenoPCA_EnvPCA$Lon) ,scale=T )

#Get 
Full.rda.r2 <- bind_cols(adjr2 = unclass(RsquareAdj(Full.rda))$adj.r.squared, analysis = "fullrda")
Full.PCcor.rda.r2 <- bind_cols(adjr2 = unclass(RsquareAdj(Full.PCcor.rda))$adj.r.squared, analysis = "fullrda")
Full.Geocor.rda.r2 <- bind_cols(adjr2 = unclass(RsquareAdj(Full.Geocor.rda))$adj.r.squared, analysis = "fullrda")
AllRDA.r2 <- bind_rows(Full.rda.r2, Full.PCcor.rda.r2, Full.Geocor.rda.r2) %>%  mutate(dataset = subproj)
write.table(AllRDA.r2, paste0("resfigntable/", subproj, ".RDA.r2.tsv"), col.names = F, row.names = F, sep = "\t", quote = F)

Full.rda.site_scores <- as.data.frame(unclass(scores(Full.rda, display = "sites")))
Full.rda.PCcor.site_scores <- as.data.frame(unclass(scores(Full.PCcor.rda, display = "sites")))
Full.rda.Geocor.site_scores <- as.data.frame(unclass(scores(Full.Geocor.rda, display = "sites")))

colnames(Full.rda.PCcor.site_scores) <- c("RDA1.PCcor", "RDA2.PCcor")
colnames(Full.rda.Geocor.site_scores) <- c("RDA1.Geocor", "RDA2.Geocor")



GenoPCA_EnvPCA_RDA <- bind_cols(GenoPCA_EnvPCA, Full.rda.site_scores, Full.rda.PCcor.site_scores)

#new RDA plot bit nicer
GenoRDA_12 <- ggplot() + 
  geom_point(data =GenoPCA_EnvPCA_RDA , aes(x = RDA1, y = RDA2, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")
GenoRDAPCcor_12 <- ggplot() + 
  geom_point(data =GenoPCA_EnvPCA_RDA , aes(x = RDA1.PCcor, y = RDA2.PCcor, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")


# Create the plot - Full RDA
pdf(paste0("resfigntable/", subproj, ".Fullrda.pdf"), height = 8, width = 13)
plot(Full.rda, type="n", scaling=3)
color.gradient <- colorRampPalette(c("blue", "red"))(length(unique(GenoPCA_EnvPCA$Comp.1)))
color.match <- match(GenoPCA_EnvPCA$Comp.1, sort(unique(GenoPCA_EnvPCA$Comp.1)))
points(Full.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color.gradient[color.match])
text(Full.rda, scaling=3, display="bp", col="#0868ac", cex=1)

# Close the PDF device
dev.off()

# Create the plot - PCcor RDA
pdf(paste0("resfigntable/", subproj, ".PCcor.Fullrda.pdf"), height = 8, width = 13)
plot(Full.PCcor.rda, type="n", scaling=3)
color.gradient <- colorRampPalette(c("blue", "red"))(length(unique(GenoPCA_EnvPCA$Comp.1)))
color.match <- match(GenoPCA_EnvPCA$Comp.1, sort(unique(GenoPCA_EnvPCA$Comp.1)))
points(Full.PCcor.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color.gradient[color.match])
text(Full.PCcor.rda, scaling=3, display="bp", col="#0868ac", cex=1)
# Close the PDF device
dev.off()


pdf(paste0("resfigntable/", subproj, "Geocor.Fullrda.pdf"), height = 8, width = 13)
plot(Full.Geocor.rda, type="n", scaling=3)
color.gradient <- colorRampPalette(c("blue", "red"))(length(unique(GenoPCA_EnvPCA$Comp.1)))
color.match <- match(GenoPCA_EnvPCA$Comp.1, sort(unique(GenoPCA_EnvPCA$Comp.1)))
points(Full.Geocor.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color.gradient[color.match])
text(Full.Geocor.rda, scaling=3, display="bp", col="#0868ac", cex=1)
# Close the PDF device
dev.off()


#Now we do environmental predictions
# Create the data partitions for each cluster
set.seed(123)  # for reproducibility

# Split the data into a training set
GenoPCA_EnvPCA <- GenoPCA_EnvPCA

train_set <-GenoPCA_EnvPCA %>%
  group_by(cluster) %>%
  sample_frac(size = 0.5) %>%  ungroup() %>%  select(FID, IID)

# Get the remaining data as a test set
test_set <- anti_join(GenoPCA_EnvPCA, train_set)  %>%  select(FID, IID)
write.table(test_set, "testset.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(train_set, "trainset.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

system(paste0("source ~/.zshrc; conda activate plink; plink --bfile ", subproj, 
              " --keep trainset.tsv --allow-extra-chr --recodeA --make-bed --maf 0.01 --out ", subproj, ".train"))

##RDA and PCA

genos.dose.train <- fread(paste0(subproj, ".train.raw"), sep = " ") %>%  select (-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
genos.dose.imp.train<- apply(genos.dose.train, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

train.fams <- fread(paste0(subproj, ".train.fam")) %>% 
  mutate(FID = as.character(V1),  IID = V2) %>% 
  select(FID, IID)
train.pcadapt.object <-  read.pcadapt(paste0(subproj, ".train.bed"), type = "bed")
train_ENV <- inner_join(train.fams, Env_PCs_clusters)
train.PCA <- pcadapt(input = train.pcadapt.object , K = 5)
train.PCAScores <- train.PCA$scores
colnames(train.PCAScores) <- paste0("PC", rep(1:5))

train_ENV <- bind_cols(train_ENV, train.PCAScores)


##RDA
train.rda <- rda(genos.dose.imp.train ~ train_ENV$Comp.1 + train_ENV$Comp.2  + train_ENV$Comp.3  ,scale=T )
train.PCcor.rda <- rda(genos.dose.imp.train ~ train_ENV$Comp.1 + train_ENV$Comp.2  + train_ENV$Comp.3 + Condition(train_ENV$PC1 + train_ENV$PC2+ train_ENV$PC3), scale=T )

#outlier detection
snpmap <- fread(paste0(subproj, '.train.bim'))

rda_snploadings <- bind_cols(snpmap, train.rda$CCA$v)
OL_RDA1_100 <- rda_snploadings %>%   slice_max(abs(RDA1), n = 100)
OL_RDA1_200 <- rda_snploadings %>%   slice_max(abs(RDA1), n = 200)
OL_RDA1_300 <- rda_snploadings %>%   slice_max(abs(RDA1), n = 300)
OL_RDA1_400 <- rda_snploadings %>%   slice_max(abs(RDA1), n = 400)
OL_RDA1_500 <- rda_snploadings %>%   slice_max(abs(RDA1), n = 500)

OL_RDA2_100 <- rda_snploadings %>%   slice_max(abs(RDA2), n = 100)
OL_RDA2_200 <- rda_snploadings %>%   slice_max(abs(RDA2), n = 200)
OL_RDA2_300 <- rda_snploadings %>%   slice_max(abs(RDA2), n = 300)
OL_RDA2_400 <- rda_snploadings %>%   slice_max(abs(RDA2), n = 400)
OL_RDA2_500 <- rda_snploadings %>%   slice_max(abs(RDA2), n = 500)

OL_100 <- unique(c(OL_RDA1_100$V2, OL_RDA2_100$V2))
OL_200 <- unique(c(OL_RDA1_200$V2, OL_RDA2_200$V2))
OL_300 <- unique(c(OL_RDA1_300$V2, OL_RDA2_300$V2))
OL_400 <- unique(c(OL_RDA1_400$V2, OL_RDA2_400$V2))
OL_500 <- unique(c(OL_RDA1_500$V2, OL_RDA2_500$V2))


rda_snploadings
write.table(OL_100, "OLSNPs100.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(OL_200, "OLSNPs200.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(OL_300, "OLSNPs300.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(OL_400, "OLSNPs400.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(OL_500, "OLSNPs500.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

system(paste0("source ~/.zshrc; conda activate plink; for i in 100 200 300 400 500 ; do plink --bfile ", subproj, 
              " --keep testset.tsv --recodeA --allow-extra-chr --extract OLSNPs$i.tsv --make-bed  --out ", subproj, ".test.OL.$i ; done"))


rda_PCcor_snploadings <- bind_cols(snpmap, train.PCcor.rda$CCA$v)
PC_cor_OL_RDA1_100 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA1), n = 100)
PC_cor_OL_RDA1_200 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA1), n = 200)
PC_cor_OL_RDA1_300 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA1), n = 300)
PC_cor_OL_RDA1_400 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA1), n = 400)
PC_cor_OL_RDA1_500 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA1), n = 500)

PC_cor_OL_RDA2_100 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA2), n = 100)
PC_cor_OL_RDA2_200 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA2), n = 200)
PC_cor_OL_RDA2_300 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA2), n = 300)
PC_cor_OL_RDA2_400 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA2), n = 400)
PC_cor_OL_RDA2_500 <- rda_PCcor_snploadings %>%   slice_max(abs(RDA2), n = 500)

PC_cor_OL_100 <- unique(c(PC_cor_OL_RDA1_100$V2, PC_cor_OL_RDA2_100$V2))
PC_cor_OL_200 <- unique(c(PC_cor_OL_RDA1_200$V2, PC_cor_OL_RDA2_200$V2))
PC_cor_OL_300 <- unique(c(PC_cor_OL_RDA1_300$V2, PC_cor_OL_RDA2_300$V2))
PC_cor_OL_400 <- unique(c(PC_cor_OL_RDA1_400$V2, PC_cor_OL_RDA2_400$V2))
PC_cor_OL_500 <- unique(c(PC_cor_OL_RDA1_500$V2, PC_cor_OL_RDA2_500$V2))


write.table(PC_cor_OL_100, "PCcor.OLSNPs100.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(PC_cor_OL_200, "PCcor.OLSNPs200.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(PC_cor_OL_300, "PCcor.OLSNPs300.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(PC_cor_OL_400, "PCcor.OLSNPs400.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(PC_cor_OL_500, "PCcor.OLSNPs500.tsv", col.names = F, row.names = F, sep = "\t", quote = F)


system(paste0("source ~/.zshrc; conda activate plink; for i in 100 200 300 400 500 ; do plink --bfile ", subproj, 
              " --keep testset.tsv --recodeA --allow-extra-chr --extract PCcor.OLSNPs$i.tsv --make-bed  --out ", subproj, ".test.PCcor.OL.$i ; done"))


#import data

test.fams <- fread(paste0(subproj,  ".test.PCcor.OL.100.fam")) %>% 
  mutate(FID = as.character(V1),  IID = V2) %>% 
  select(FID, IID) 
test_ENV <- inner_join(test.fams, GenoPCA_EnvPCA)



# Function to import datasets
import_datasets <- function(suffixes, base_name) {
  
  # Use purrr::map to iterate over suffixes
  map(suffixes, ~ {
    # Create file name
    file_name <- paste0(subproj, base_name, ".", .x, ".raw")
    
    # Read data
    data <- fread(file_name, sep = " ") %>%
      select(-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
    
    # Replace NA values
    data.imp <- apply(data, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
    
    # Assign data frame to new variable in global environment
    assign(paste0("genos", base_name, "_", .x), data.imp, envir = .GlobalEnv)
  })
}

# List of suffixes
suffixes <- c(100, 200, 300, 400, 500)

# Import datasets
import_datasets(suffixes, ".test.PCcor.OL")
import_datasets(suffixes, ".test.OL")


#random foret function

# Define the variables to be used as y in the model
y_vars <- c(paste0("Comp.", 1:3), BioOracle2.bottom.var)

# Initialize a data frame to store the results
results <- data.frame()


# Define the function
run_rf_model <- function(y_vars, x_data, x_name, mtry) {
  
  # Initialize a data frame to store the results
  results <- data.frame()
  
  # Loop over the y variables
  for (y_var in y_vars) {
    
    # Run the model 5 times
    for (i in 1:5) {
      
      # Fit the model
      rf_model <- ranger(y = test_ENV[[y_var]], 
                         x = x_data,
                         mtry = mtry, num.trees = 1000, num.threads = 30)
      
      # Add the results to the data frame
      results <- rbind(results, data.frame(Variable = y_var, R2 = rf_model$r.squared, Outlier_set = x_name, Run = i))
    }
  }
  
  # Return the results
  return(results)
}


# Call the function
results_OL_100 <- run_rf_model(y_vars = y_vars, x_data = genos.test.OL_100, x_name = "OL_100", mtry = 160)
results_OL_200 <- run_rf_model(y_vars = y_vars, x_data = genos.test.OL_200, x_name = "OL_200", mtry = 320)
results_OL_300 <- run_rf_model(y_vars = y_vars, x_data = genos.test.OL_300, x_name = "OL_300", mtry = 480)
results_OL_400 <- run_rf_model(y_vars = y_vars, x_data = genos.test.OL_400, x_name = "OL_400", mtry = 640)
results_OL_500 <- run_rf_model(y_vars = y_vars, x_data = genos.test.OL_500, x_name = "OL_500", mtry = 800)

results_PC_cor100 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_100, x_name = "PCcor_OL_100", mtry = 160)
results_PC_cor200 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_200, x_name = "PCcor_OL_200", mtry = 320)
results_PC_cor300 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_300, x_name = "PCcor_OL_300", mtry = 480)
results_PC_cor400 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_400, x_name = "PCcor_OL_400", mtry = 640)
results_PC_cor500 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_500, x_name = "PCcor_OL_500", mtry = 800)


RF_env_res <- bind_rows(results_OL_100,
                        results_OL_200,
                        results_OL_300,
                        results_OL_400,
                        results_OL_500,
                        results_PC_cor100,
                        results_PC_cor200,
                        results_PC_cor300,
                        results_PC_cor400,
                        results_PC_cor500) %>%  mutate(Dataset = subproj)

color_palette <- c(colorRampPalette(c("lightblue", "blue"))(5),
                   colorRampPalette(c("bisque", "orange"))(5))

RF_result_plot <- ggplot() + 
  geom_point(data = RF_env_res, aes(x = Variable, y = R2, colour = Outlier_set, size = 3)) + 
  theme_classic() + 
  ylim(-0.25,1)  +
  scale_color_manual(values = color_palette)




write.table(RF_env_res, paste0("resfigntable/", subproj, "RF_env_res.tsv"), col.names = T,
            row.names = F, 
            sep = "\t", 
            quote = F)

#save plot
ggsave(plot = RF_result_plot, filename = paste0("resfigntable/", subproj, "RF_results_plot.pdf"), 
       device = pdf, 
       height = 8, 
       width = 13)

#separate and test:
RF_env_res <-  RF_env_res %>%  mutate(PCcor = str_replace(Outlier_set, "_.*", ""))
RF_env_res <-  RF_env_res %>%  mutate(Sampsize = str_replace(Outlier_set, ".*_", ""))

PCcor <- RF_env_res %>% filter(PCcor %in% "PCcor")
OL <- RF_env_res %>% filter(PCcor %in% "OL")

wilcox.test(PCcor$R2, OL$R2)



