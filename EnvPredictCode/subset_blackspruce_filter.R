rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/23_filter")
dir()


library(tidyverse)
library(dartR)
library(cowplot)
library(readxl)
library(stringr)
library(adegenet)


############
### DATA ###
############

### This part is only to select a subset of samples from the whole dataset based 
### on previous admixture analysis and overall quality statistics
### Below is the code to filter SNPs starting with thr raw dataset


# Table with ok samples (excluding potential mixed samples)
df <- read.csv("01_remove_mixed_admixture_metadata_clean.tsv", sep="\t", header=T)
dim(df)
head(df)

# Removing red spruce, and subclusters (other than West,Central, or East)
df_sub <- df %>% filter(BestK6 %in% c("K6.West", "K6.Central", "K6.East")) %>% filter(BestK4 != "K4.RedSpruce")
table(df_sub$BestK6)
head(df_sub)

dg <- df_sub %>% group_by(POP) %>% dplyr::summarise(n=n(), med_cr = median(call_rate), strata = first(STRATA)) %>% arrange(strata, -med_cr)
head(dg,20)

# Selecting provenances with best median call rates
pops <- c("4420","7000","2209","6988","6983","6999","6971","6965","6927","6916","6914","333","6901","352","355","1530","6804")
pops.sel <- dg %>% filter(POP %in% pops)
pops.sel

df_pops <- df_sub %>% filter(POP %in% pops.sel$POP)
df_sel <- df_pops %>% arrange(-call_rate) %>% group_by(POP) %>% filter(row_number() %in% c(1:15))
dim(df_sel)
df_sel %>% group_by(POP) %>% summarize(n=n())


ggplot(df_sel) + aes(y = lat, x = lon) +
  geom_point()

ggplot(df_sel) + aes(x = POP, y =call_rate) +
  geom_boxplot()



#################
### FUNCTIONS ###
#################


# A function to obtain per sample call rate
ind_call_rate = function(row, output){
  x = length(row[!is.na(row)])/length(row)
  return(x)
}


# Function for calculating heterozygosity per sample
ind_freq_het = function(row, output){
  hets <- row[row == 2]
  hets_complete <- complete.cases(hets)
  hets_nonan <- hets[hets_complete]
  x = length(hets_nonan)/length(row[!is.na(row)])
  return(x)
}




##########################
#  Filtering genotypes   #
##########################

# Raw variants dataset
gel <- gl.load("../DATA_intermediate/04_merging.Rdata")
gel

samples.to.keep <- df_sel$id
gel.clean <- gl.keep.ind(gel, ind.list = samples.to.keep, recalc = TRUE, mono.rm = TRUE)
gel.clean



### Finding loci names with poor call rate

lociPoorCallRate <- function(gel, minCR) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_minCR <- gel_meta %>% filter(CallRate < minCR) %>% pull(locus)
  return(gel_minCR)
}

poor.loci <- lociPoorCallRate(gel.clean, 0.5)
length(poor.loci)




### Finding loci names with excess heterozygosity

head(gel.clean@other$loc.metrics)
lociHighHet <- function(gel, maxHet) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_maxHet <- gel_meta %>% filter(FreqHets > maxHet) %>% pull(locus)
  return(gel_maxHet)
}

highHet.loci.EPN <- lociHighHet(gel.clean, 0.5)
length(highHet.loci.EPN)

ggplot(gel.clean@other$loc.metrics) + aes(x = FreqHets) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  labs(y = "N SNPs") +
  geom_vline(xintercept = 0.5, colour = "skyblue") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))


### Finding loci names with low MAF (<0.01)

lociLowMAF <- function(gel, minMAF) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_minHet <- gel_meta %>% filter(maf < minMAF) %>% pull(locus)
  return(gel_minHet)
}

lowMAF.loci.EPN <- lociLowMAF(gel.clean, 0.01)
length(lowMAF.loci.EPN)

ggplot(gel.clean@other$loc.metrics) + aes(x = maf) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  labs(y = "N SNPs") +
  geom_vline(xintercept = 0.01, colour = "skyblue") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))



### Filtering poor SNPs (low call rate, high het, low MAF)

loci.to.drop <- unique(c(poor.loci, highHet.loci.EPN, lowMAF.loci.EPN))
length(loci.to.drop)


gel.CR <- gl.drop.loc(x = gel.clean, loci.to.drop)
gel.CR <- gl.filter.monomorphs(gel.CR)
gel.CR <- gl.recalc.metrics(gel.CR)
gel.CR





####################################
### Overall SNP call rate filter ###
####################################

# 50%
gel.CR.CR <- gl.filter.callrate(x = gel.CR,
                                method = "loc",
                                threshold = 0.5,
                                mono.rm = TRUE,
                                recalc = TRUE,
                                recursive = TRUE)
gel.CR.CR
# 2,322 genotypes,  36,430 SNPs


######################
### No secondaries ###
######################

#gl.report.secondaries(gel.CR.CR)
#gel.CR.CR.nosec <- gl.filter.secondaries(gel.CR.CR, method = "best", verbose=3)
#gel.CR.CR.nosec <- gl.filter.monomorphs(gel.CR.CR.nosec)
#gel.CR.CR.nosec <- gl.recalc.metrics(gel.CR.CR.nosec)

removeSecondaries <- function(gel) {
  df <- data.frame("locus" = locNames(gel))
  df$CallRate <- gel@other$loc.metrics$CallRate
  df[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(df$locus, '-', 3)
  df_grouped <- df %>% group_by(CloneID) %>% mutate(n = n()) %>% ungroup()
  high.CR.loci <- df_grouped %>% group_by(CloneID) %>% filter(CallRate == max(CallRate)) %>% pull(locus)
  
  gel.nosec <- gl.keep.loc(x = gel, high.CR.loci)
  gel.nosec <- gl.filter.monomorphs(gel.nosec)
  gel.nosec <- gl.recalc.metrics(gel.nosec)
  return(gel.nosec)
}

gel.CR.CR.nosec <- removeSecondaries(gel.CR.CR)


# Number of loci before
dim(gel.CR.CR)[2]
# 36430
# Number of loci after filtering
dim(gel.CR.CR.nosec)[2]
# 31294


### Random 5000

set.seed(18102024)
loc.sample <- sample(gel.CR.CR.nosec@loc.names, 5000)
loc.sample

gel.subset <- gl.keep.loc(x = gel.CR.CR.nosec, loc.sample)
gel.subset <- gl.filter.monomorphs(gel.subset)
gel.subset <- gl.recalc.metrics(gel.subset)
gel.subset



##############
### Saving ###
##############

gl.save(gel.subset, file="./03_filter_snps_TK.Rdata")





###########################
#    Converting to VCF    #
###########################


dmeta <- gel.subset@other$loc.metrics
head(dmeta)
dmeta$locus <- row.names(dmeta)
dmeta[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(dmeta$locus, '-', 3)
dmeta[c('REF','ALT')] <- str_split_fixed(dmeta$Alleles, '/', 2)
head(dmeta)


gel.subset$chromosome <- as.factor(paste0("Tag_",dmeta$CloneID))
gel.subset$chromosome 

gel.subset$position <- as.factor(dmeta$TagPos)
gel.subset$position
gel.subset



gl2plink(
  gel.subset,
  plink_path = "//wsl.localhost/Ubuntu/home/BlackSpruce/23_filter",
  bed_file = FALSE,
  outfile = "03_filter_snps_TK",
  outpath = getwd(),
  chr_format = "character",
  pos_cM = "0",
  ID_dad = "0",
  ID_mom = "0",
  sex_code = "unknown",
  phen_value = "0",
  verbose = NULL
)



#gl.save(genl, file="BlackSpruce_1431_genotypes_genlight.Rdata")



meta <- gel.subset@other$ind.metrics
head(meta)

meta <- meta %>% dplyr::select(id, POP, STRATA, lat, lon)
dim(meta)
meta$vcfsample <- paste0("pop1_",meta$id)
head(meta)

write.table(meta, file="03_filter_snps_TK_metadata.tsv", sep="\t", col.names = T, row.names = F, append=F, quote=F)

###########


