# This script performs normalization and filtering of CpG probes on EPIC methylation
# microarray data from various Treg adoptive cell products.
# * First compared different normalization methods. 
# * Filtered out CpGs that failed in more than one sample, overlap with SNPs, and are
#   predicted to be crossreactive. 
# * Did NOT remove CpGs on sex chromosomes! Do this during differential analysis.  

# SETUP ----
library("tidyverse")
library("minfi")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("IlluminaHumanMethylationEPICmanifest")
library("cowplot")
library("factoextra")
library("ggrepel")


# IMPORT DATA ----
# sample sheet in same directory as idat files
dir = "S:/AG/AG-Polansky-Biskup/Daten/KO+AS/EPIC_Array/idat"
targets = read.metharray.sheet(base = dir, pattern ="ko_all.csv", recursive = T)
rgSet = read.metharray.exp(targets = targets, force = T)

# QUALITY CONTROL ----
## Name individual elements/samples of rgSet list
sampleNames(rgSet) = targets$ID

## Generate default report 
qcReport(rgSet, sampNames = targets$ID, sampGroups = targets$Cell, pdf = "./210211_reshape/figs/all.qcReport.pdf")

## Remove samples with average detP > 0.05
### NOTE: detP dimensions: rows = cg, cols = sample
detP = detectionP(rgSet)
keep = colMeans(detP) < 0.05

# columns = samples
rgSet = rgSet[,keep]

# rows = samples
targets = targets[keep,]

# NORMALIZATION ----
## want to test whether quantile or functional normalization or noob normalization is more appropriate
mSetNorm1 = preprocessQuantile(rgSet)

mSetNorm2 = preprocessFunnorm(rgSet)

mSetNorm3 = preprocessNoob(rgSet, dyeMethod = "single")

## visualizing before and after normalization

# create function to plot density
plot_density = function(norm.obj){
  plot = getBeta(norm.obj) %>%
    as.data.frame() %>%
    pivot_longer(cols = 1:ncol(norm.obj), values_to = "beta", names_to = "ID") %>%
    left_join(., targets[,c("ID","Cell")], by = "ID") %>%
    ggplot(., aes(beta, group = ID, color = Cell))+
    geom_density()+
    labs(title = "Norm Beta Values", x = "Beta", y = "Density")+
    theme_classic()
  
  return(plot)
}

# save density plots
density.quantile = plot_density(mSetNorm1) + labs(subtitle = "Quantile Normalization")
density.funnorm = plot_density(mSetNorm2) + labs(subtitle = "Functional Normalization")
density.noob = plot_density(mSetNorm3) + labs(subtitle = "Noob Normalization")

plot_grid(density.quantile, density.funnorm, density.noob, align = "v", nrow = 3)

# create function to plot pca
plot_pca = function(norm.obj){
  analysis = t(getBeta(norm.obj)) %>% prcomp()
  
  variance.perc = t(get_eigenvalue(analysis)) %>%
    as.data.frame(rownames = "Feature") %>%
    rename_all(~str_replace(.,"Dim.","PC")) %>%
    round(digits = 2)
  
  plot = as_tibble(analysis$x, rownames = NA) %>% 
    rownames_to_column(var = "ID") %>% 
    left_join(targets, by = "ID") %>%
    ggplot(data = ., aes(PC1, PC2, group = ID, color = Cell))+
    geom_jitter(size = 3)+
    theme_bw() +
    xlab(paste("PC1 ",toString(variance.perc[2, "PC1"]),"%")) +
    ylab(paste("PC2 ",toString(variance.perc[2, "PC2"]), "%"))
  
  return(plot)
}

pca.quantile = plot_pca(mSetNorm1) + labs(subtitle = "Quantile Normalization")
pca.funnorm = plot_pca(mSetNorm2) + labs(subtitle = "Functional Normalization") 
pca.noob = plot_pca(mSetNorm3) + labs(subtitle = "Noob Normalization")

plot_grid(pca.quantile, pca.funnorm, pca.noob, align = "v", nrow = 3)

# chosen normalization method = single sample Noob, based on publication Fortin, et al (2017)

# FILTERING ----
## Probes are in the same order in the mSetNorm and detP objects 
detP = detP[match(featureNames(mSetNorm3), rownames(detP)),]

## Remove any probes that have failed in one or more samples, cutoff at pval = 0.01
keep = rowSums(detP < 0.01) == ncol(mSetNorm3)

## Create filtered dataset 
mSetFlt = mSetNorm3[keep,]

## Remove probes which overlap with SNPs 
snps = getSnpInfo(mSetFlt)

# check to see all rownames match
all(rownames(mSetFlt) == rownames(getBeta(mSetFlt)))

# filter out mSetFlt containing snps
mSetFlt = mSetFlt[is.na(snps$Probe_rs),]

## Remove cross reactive
cross.reactive.sites = read.csv("./crossreactive.csv")
keep = !(featureNames(mSetFlt) %in% cross.reactive.sites$X)

mSetFlt = mSetFlt[keep,]


# SAVE RDS for later processing and analysis
saveRDS(mSetFlt, file = "./210211_reshape/rds/mSetFlt.rds")
saveRDS(targets, file = "./210211_reshape/rds/targets.rds")
