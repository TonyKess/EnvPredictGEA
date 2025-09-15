library(tidyverse)
library(data.table)
library(geodata)
library(raster)
#set environment
wd = "~/Desktop/Projects/EnvPredict/"
subproj = "AtlanticSalmon"
setwd(paste0(wd, subproj))
data_path = "~/Desktop/Projects/EnvPredict/Envdataset"

#bring in geo data
Geo <- data.frame(fread(paste0(subproj, "_Geo.tsv"))) 


# res = Valid resolutions are 10, 5, 2.5 (minutes of a degree)
res = 2.5

bioclim_data <- worldclim_global(var = "bio",
                                 res = 2.5,
                                 path = data_path)

Env <- raster::extract(bioclim_data, Geo %>% dplyr::select(Lon, Lat)) 
colnames(Env)<- c("id", (paste0("BIO", 1:19)))
Env <- Env %>% dplyr::select(-id)

Env <- bind_cols(Geo, Env)

write.table(Env, 
            paste0(subproj, "_Env.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
