This workflow goes through show code from EnvPredict performs:
1. Extraction of environmental variables
2. Scaling and PCA of envrionmental variables
3. PCA of genetic variation in each species,and correlation of genetic and environmental PCs per species.
4. RDA of genetic and environmental variation in each species, with and without structure correction
5. Splitting of datasets into detect and test sets
6. Outlier detection using RDA in detect set, with and without structure correction
7. Random forest comparison of outlier sets with and without structure correction to test environmental prediction capacity. 

# 1. Environmental variables:
Bring in environments and libraries 
```
library(tidyverse)
library(data.table)
library(geodata)
library(raster)
#set environment
wd = "<yourpath>/EnvPredict/"
subproj = "AtlanticSalmon"
setwd(paste0(wd, subproj))
data_path = "~/EnvPredict/Envdataset"
```
Next bring in geo data as Lats/Longs
```
Geo <- data.frame(fread(paste0(subproj, "_Geo.tsv"))) 


# res = Valid resolutions are 10, 5, 2.5 (minutes of a degree)
res = 2.5
```
Download Bioclim data:
```
bioclim_data <- worldclim_global(var = "bio",
                                 res = 2.5,
                                 path = data_path)

Env <- raster::extract(bioclim_data, Geo %>% dplyr::select(Lon, Lat)) 
colnames(Env)<- c("id", (paste0("BIO", 1:19)))
Env <- Env %>% dplyr::select(-id)

Env <- bind_cols(Geo, Env)
```
Write out the table of Environmental data for analysis:
```
write.table(Env, 
            paste0(subproj, "_Env.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
```

# 2. Scaling and PCA of envrionmental variables - 
[This code marks the begining of the Env_PCA_RDA_RFpredict_pipeline workflow](https://github.com/TonyKess/EnvPredictGEA/blob/main/EnvPredictCode/ENV_PCA_RDA_RFpredict_Pipeline_BioClim_AtlanticSalmon.R)

First, load packages
```

library(tidyverse)
library(data.table)
library(vegan)
library(pcadapt)
library(caret)
library(robust)
library(qvalue)
library(ranger)
library(patchwork)
```

Set environments
```
wd = "<yourpath>/EnvPredict/"
subproj = "AtlanticSalmon"
setwd(paste0(wd, subproj))
system("mkdir resfigntable") # save figures here
```

Bring in geographic and environmental data
```
Env<- data.frame(fread(paste0(subproj, "_Env.tsv")))
```
Scale bioclim data for PCA - reduces variation from variables on different scales
```
BIOvar <- paste0("BIO", 1:19)

BIOvar_sd <- paste0(BIOvar, ".scale")

Env<- Env %>%
  mutate(across(all_of(BIOvar), ~ as.numeric(scale(.)), .names = "{.col}.scale"))
```

PCA bioclim variables
```
Env_SDmatrix <-Env  %>% select(all_of(BIOvar_sd))

Env_SDmatrix <- as.matrix(Env_SDmatrix)

Env_PCs <- princomp(Env_SDmatrix) 
summary(Env_PCs)


#combine env and PC env data
EnvPC_Biovars <- bind_cols(Env,Env_PCs$scores[,1:5])
```
Now we build a function to evaluate the number of  kmeans clusters in environmental data for splitting test and training data
```
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
```
Apply to our environmental data
```
kmeans_evalfunc(10,Env_PCs$scores) + theme_classic()
```
Cool, K = 5 looks fine, wss just keeps going down but marginally.
<img width="1061" height="653" alt="image" src="https://github.com/user-attachments/assets/f06ef435-d64a-4297-b8cf-009d219634f9" />

Now we make a function to add clusters:
```
kmeans_addclusters_func <- function(n_clusters, PCA){
  set.seed(123)
  kmeanrun <-  kmeans(PCA, centers = n_clusters, nstart = 1000)
  plinkPCAdataclustered <- bind_cols(PCA, cluster= kmeanrun$cluster)
  return(plinkPCAdataclustered)}

Env_PCs_clusters <-  kmeans_addclusters_func(n_clusters = 5, PCA = Env_PCs$scores)
```
Add our clusters to our environmental PC data
```
Env_PCs_clusters  <- bind_cols(EnvPC_Biovars, cluster = Env_PCs_clusters$cluster)
```
# 3. PCA of genetic variation in each species,and correlation of genetic and environmental PCs per species.
Now we use pcadapt to carry out PCA of genomic data. First read data into pcadapt, and match to the environmental data:
```
pcadapt.object <-  read.pcadapt(paste0(subproj, ".bed"), type = "bed")
Fam  <- fread(paste0(subproj, ".fam")) %>%  
  mutate(FID = V1, IID = V2) %>%
  select(FID, IID)
Env_PCs_clusters <- inner_join(Fam, Env_PCs_clusters)
```
Now carry out PCA and get PCA scores across 5 axes:
```
PCA <- pcadapt(input = pcadapt.object , K = 5)
PCAScores <- PCA$scores
colnames(PCAScores) <- paste0("PC", rep(1:5))
```
Evaluate the Scree plot:
```
plot(PCA, option = "screeplot") + theme_classic()

```

<img width="895" height="542" alt="image" src="https://github.com/user-attachments/assets/3829b021-77eb-4473-80d7-d7fb30c7e828" />
Cattel method says K ~ 3

Now combine Genetic and Environmental PCAs and plot:
```
GenoPCA_EnvPCA<- data.frame(bind_cols(Env_PCs_clusters, PCAScores))

GenoPCA_12 <- ggplot() + 
  geom_point(data = GenoPCA_EnvPCA, aes(x = PC1, y = PC2, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")
```
<img width="804" height="475" alt="image" src="https://github.com/user-attachments/assets/58e530b0-4c4b-4a8b-8946-b2928e87235e" />
Lots of genome/environment correlation!

# 4. RDA of genetic and environmental variation in each species, with and without structure correction
First use plink to make a .raw dosage file, RDA likes these. You can do this right from the R console:
```
system(paste0("source ~/.zshrc; conda activate plink; plink --bfile ", subproj, 
              " --recodeA --allow-extra-chr --chr-set 29 --maf 0.01 --out ", subproj))
```
Read in the dosages and impute by common value using Forester script:
```
genos.dose <- fread(paste0(subproj, ".raw"), sep = " ") %>%  select (-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
genos.dose.imp<- apply(genos.dose, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
```

Run RDA, without structure correction:
```
Full.rda <- rda(genos.dose.imp ~ GenoPCA_EnvPCA$Comp.1 + GenoPCA_EnvPCA$Comp.2  + GenoPCA_EnvPCA$Comp.3  ,scale=T )
```
And with the first 3 genomic PCs:
```
Full.PCcor.rda <- rda(genos.dose.imp ~  GenoPCA_EnvPCA$Comp.1 + GenoPCA_EnvPCA$Comp.2  + GenoPCA_EnvPCA$Comp.3  + Condition(GenoPCA_EnvPCA$PC1 + GenoPCA_EnvPCA$PC2 + GenoPCA_EnvPCA$PC3) ,scale=T )
```
Get RDA scores per individual in the ordinations, combine them into data frames with environmental data:
```
Full.rda.site_scores <- as.data.frame(unclass(scores(Full.rda, display = "sites")))
Full.rda.PCcor.site_scores <- as.data.frame(unclass(scores(Full.PCcor.rda, display = "sites")))
colnames(Full.rda.PCcor.site_scores) <- c("RDA1.PCcor", "RDA2.PCcor")
colnames(Full.rda.Geocor.site_scores) <- c("RDA1.Geocor", "RDA2.Geocor")


GenoPCA_EnvPCA_RDA <- bind_cols(GenoPCA_EnvPCA, Full.rda.site_scores, Full.rda.PCcor.site_scores)
```
Now plot RDA:
```
GenoRDA_12 <- ggplot() + 
  geom_point(data =GenoPCA_EnvPCA_RDA , aes(x = RDA1, y = RDA2, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")
```
<img width="856" height="549" alt="image" src="https://github.com/user-attachments/assets/4edb4f54-f4e6-4469-86df-59c136e9c2ef" />
Wow, looks great! The environmental gradient and genetic clustering is clear across the first two RDA axes

Now lets plot the PC-corrected RDA:

```
GenoRDAPCcor_12 <- ggplot() + 
  geom_point(data =GenoPCA_EnvPCA_RDA , aes(x = RDA1.PCcor, y = RDA2.PCcor, colour = Comp.1), size = 3) + 
  theme_classic() +
  scale_color_gradient(low = "blue", high = "red")
```
<img width="841" height="522" alt="image" src="https://github.com/user-attachments/assets/e9d445ac-cc49-46c9-9881-0fa26c401675" />

This looks worse.
Note these are the "spruced up" figures from the manuscript but mostly I just made the font bigger and bolder.

# 5. Splitting of datasets into detect and test sets
We separate these for outlier detection and then testing whether those outliers are any good at predicting environments. We need to split the data, starting with a training set for detection of outliers:

```
set.seed(123)  # for reproducibility

GenoPCA_EnvPCA <- GenoPCA_EnvPCA %>%  mutate(FID = str_replace(FID, "pop_", ""))

train_set <-GenoPCA_EnvPCA %>%
  group_by(cluster) %>%
  sample_frac(size = 0.5) %>%  ungroup() %>%  select(FID, IID)
```

And then we get the remaining data as a test set
```
test_set <- anti_join(GenoPCA_EnvPCA, train_set)  %>%  select(FID, IID)
write.table(test_set, "testset.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(train_set, "trainset.tsv", col.names = F, row.names = F, sep = "\t", quote = F)
```
Subset the genotype file for the detection/train set:
```
system(paste0("source ~/.zshrc; conda activate plink; plink --bfile ", subproj, 
              " --keep trainset.tsv --allow-extra-chr --recodeA --make-bed --chr-set 29 --maf 0.01 --out ", subproj, ".train"))
```
# 6. Outlier detection using RDA in detect set, with and without structure correction

In the detection/train set, we do RDA to detect outliers. First import, impute:
```
genos.dose.train <- fread(paste0(subproj, ".train.raw"), sep = " ") %>%  select (-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
genos.dose.imp.train<- apply(genos.dose.train, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

train.fams <- fread(paste0(subproj, ".train.fam")) %>% 
  mutate(FID = as.character(V1),  IID = V2) %>% 
  select(FID, IID)
```
PCA of the test set for including a correction in pRDA:
```
train.pcadapt.object <-  read.pcadapt(paste0(subproj, ".train.bed"), type = "bed")
train.PCA <- pcadapt(input = train.pcadapt.object , K = 5)
train.PCAScores <- train.PCA$scores
```
Match to environmental data:
```
train_ENV <- inner_join(train.fams, Env_PCs_clusters)
colnames(train.PCAScores) <- paste0("PC", rep(1:5))
train_ENV <- bind_cols(train_ENV, train.PCAScores)
```
RDA:
```
train.rda <- rda(genos.dose.imp.train ~ train_ENV$Comp.1 + train_ENV$Comp.2  + train_ENV$Comp.3  ,scale=T )
train.PCcor.rda <- rda(genos.dose.imp.train ~ train_ENV$Comp.1 + train_ENV$Comp.2  + train_ENV$Comp.3 + Condition(train_ENV$PC1 + train_ENV$PC2 ,scale=T ))
```
Outlier detection using absolute loadings on RDA axis 1 and 2:
```
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
```

Subset the test individuals for these outliers:
```
system(paste0("source ~/.zshrc; conda activate plink; for i in 100 200 300 400 500 ; do plink --bfile ", subproj, 
              " --keep testset.tsv --recodeA --allow-extra-chr --extract OLSNPs$i.tsv --make-bed --chr-set 29  --out ", subproj, ".test.OL.$i ; done"))
```
Now do the same for RDA with population structure correction:
```
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
              " --keep testset.tsv --recodeA --allow-extra-chr --extract PCcor.OLSNPs$i.tsv --make-bed --chr-set 29  --out ", subproj, ".test.PCcor.OL.$i ; done"))
```
# 7. Random forest comparison of outlier sets with and without structure correction to test environmental prediction capacity. 
Import the IDs and match to metadata:
```
test.fams <- fread(paste0(subproj,  ".test.PCcor.OL.100.fam")) %>% 
  mutate(FID = V1,  IID = V2) %>% 
  select(FID, IID) 
test_ENV <- inner_join(test.fams, GenoPCA_EnvPCA)
```

We use this function to import all of the SNP dosage datasets for each set RDA outliers:
```
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
```

Then we train separate random forest models on genomic data to predict environmental variation:
```
# Define the variables to be used as y in the model
y_vars <- c(paste0("Comp.", 1:3), paste0("BIO", 1:19))

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
                         mtry = mtry, num.trees = 1000)
      
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

results_PC_cor100 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_100, x_name = "PCcor_OL_100", mtry = 100)
results_PC_cor200 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_200, x_name = "PCcor_OL_200", mtry = 320)
results_PC_cor300 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_300, x_name = "PCcor_OL_300", mtry = 480)
results_PC_cor400 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_400, x_name = "PCcor_OL_400", mtry = 640)
results_PC_cor500 <- run_rf_model(y_vars = y_vars, x_data = genos.test.PCcor.OL_500, x_name = "PCcor_OL_500", mtry = 800)

```
Combine and plot the model R2:
```
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
```
<img width="776" height="506" alt="image" src="https://github.com/user-attachments/assets/bbd6249b-7f6b-48a8-965a-7716a377d1d7" />

Now we can also test the model performance between corrected and uncorrected outlier sets:
```
RF_env_res <-  RF_env_res %>%  mutate(PCcor = str_replace(Outlier_set, "_.*", ""))
RF_env_res <-  RF_env_res %>%  mutate(Sampsize = str_replace(Outlier_set, ".*_", ""))

PCcor <- RF_env_res %>% filter(PCcor %in% "PCcor")
OL <- RF_env_res %>% filter(PCcor %in% "OL")

wilcox.test(PCcor$R2, OL$R2)
```
