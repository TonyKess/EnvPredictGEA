This workflow goes through show code from EnvPredict performs:
1. Extraction of environmental variables
2. Scaling and PCA of envrionmental variables
3. PCA of genetic variation in each species,and correlation of genetic and environmental PCs per species.
4. RDA of genetic and environmental variation in each species, with and without structure correction
5. Splitting of datasets to detect and test sets
6. Outlier detection using RDA in detect set, with and without structure correction
7. Random forest comparison of outlier sets with and without structure correction to test environmental prediction capacity. 

#1. Environmental variables:
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

#2. Scaling and PCA of envrionmental variables - 
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
