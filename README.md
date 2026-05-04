# EnvPredictGEA
Code for testing the impact of population structure correction on genome environment association in eight species. This code carries out the following analyses of genomic and environmental data:
1. Extraction of environmental variables
2. Scaling and PCA of envrionmental variables
3. PCA of genetic variation in each species,and correlation of genetic and environmental PCs per species.
4. RDA of genetic and environmental variation in each species, with and without structure correction
5. Splitting of datasets to detect and test sets
6. Outlier detection using RDA in detect set, with and without structure correction
7. Random forest comparison of outlier sets with and without structure correction to test environmental prediction capacity. 

Scripts are structured such that all analyses begin at the top of the directory (here EnvPredict) and then work from each species project subdirectory.

##Data access
Data for each analysis is found at the following dois:
Atlantic salmon:
Kess et al. 2024; https://doi.org/10.5061/dryad.g1jwstqz3 

Coho salmon :
Xuereb et al. 2022; https://doi.org/10.5061/dryad.r4xgxd2gx 

Arctic charr :
Layton et al. 2021; https://datadryad.org/dataset/doi:10.5061/dryad.8sf7m0ckd 

Atlantic cod 
Kess et al. 2020; https://datadryad.org/dataset/doi:10.5061/dryad.f1vhhmgsg 

North American grey wolf 
Schweizer et al. 2016; https://datadryad.org/dataset/doi:10.5061/dryad.c9b25 

Sea scallop 
Van Wyngaarden et al. 2018; https://datadryad.org/dataset/doi:10.5061/dryad.c15v5 

Lodgepole pine 
Mahony et al. 2020; https://datadryad.org/dataset/doi:10.5061/dryad.56j8vq8 

Black spruce 
https://doi.org/10.5281/zenodo.19961100 

#Workflow
## 1. Extraction of environmental variables
Using geographic coordinates, after setting project information, do:
```
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
``` 
