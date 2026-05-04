# EnvPredictGEA
Code for testing the impact of population structure correction on genome environment association in eight species. This code carries out the following analyses of genomic and environmental data:
1. Extraction of environmental variables
2. Scaling and PCA of envrionmental variables
3. PCA of genetic variation in each species,and correlation of genetic and environmental PCs per species.
4. RDA of genetic and environmental variation in each species, with and without structure correction
5. Splitting of datasets to detect and test sets
6. Outlier detection using RDA in detect set, with and without structure correction
7. Random forest comparison of outlier sets with and without structure correction to test environmental prediction capacity. 

Scripts are structured such that all analyses begin at the top of the directory (here EnvPredictCode) and then work from each species project subdirectory. 
