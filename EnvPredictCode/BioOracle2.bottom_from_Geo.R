library(tidyverse)
library(data.table)
library(sdmpredictors)
library(raster)

wd = "~/Desktop/Projects/EnvPredict/"
subproj = "KessCod2020"
setwd(paste0(wd, subproj))
data_path = "~/Desktop/Projects/EnvPredict/Envdataset"

Geo <- data.frame(fread(paste0(subproj, "_Geo.tsv"))) 

#set this so downloads don't time out
options(timeout = max(1000, getOption("timeout")))

#get bottom layers
environment.bottom <- load_layers(c("BO2_tempmax_bdmean", 
                                    "BO2_tempmean_bdmean", 
                                    "BO2_tempmin_bdmean",
                                    "BO2_dissoxmax_bdmean", 
                                    "BO2_dissoxmean_bdmean", 
                                    "BO2_dissoxmin_bdmean", 
                                    "BO2_salinitymax_bdmean", 
                                    "BO2_salinitymean_bdmean", 
                                    "BO2_salinitymin_bdmean")) 
#get bathymetry
bathymetry <- load_layers("BO_bathymean")

my.sites.environment <- data.frame(IID=Geo$IID, depth=extract(bathymetry,Geo[,3:4]), extract(environment.bottom, Geo[,3:4]))

Env <- inner_join(Geo, my.sites.environment)

write.table(Env, 
            paste0(subproj, "_Env.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

