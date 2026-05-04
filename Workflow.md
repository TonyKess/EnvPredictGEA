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
