# Exploratory analysis of Reshape GMP-cultured Treg cell lines. 
# Aim is to identify and understand unpredicted patterns in sample subgroups. 

# SETUP ----
library("tidyverse")
library("minfi")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("IlluminaHumanMethylationEPICmanifest")
library("cowplot")
library("factoextra")
library("ggrepel")
library("clusterProfiler")
library("org.Hs.eg.db")

# DATA IMPORT ----
mSetFlt = readRDS("./210211_reshape/rds/mSetFlt.rds")
targets = readRDS("./210211_reshape/rds/targets.rds")

# Import annotation
anno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Reorder levels
targets$Cell = str_remove(targets$Cell, pattern = " ") %>%
  factor(., levels = c("Neg_Fraction","CD4.na","CD4.mem","CD4.Treg","CD8.na","CD8.mem","CD8.Treg"))

# Create new identifier column that is more descriptive
targets$Desc = str_c(targets$Cell, targets$Expansion_Days, sep = "_")

# DATA FILTERING ----
# Get Beta vals for exploratory analysis
beta = getBeta(mSetFlt)

# Remove cpgs in sex chromosomes
xy.cpgs = rownames(anno[anno$chr %in% c("chrX","chrY"),])
  
beta.noXY = beta[!(rownames(beta) %in% xy.cpgs),]

# Remove samples that were not isolated/expanded under GMP conditions (so KO samples)
beta.GMP = beta.noXY[,str_detect(colnames(beta.noXY), pattern = "KO", negate = T)]
targets.GMP = targets[str_detect(targets$ID, pattern = "KO", negate = T),]

# PCA ANALYSIS ----
pca = prcomp(t(beta.GMP))

# get variance contribution of each PC from facto_extra
variance.perc = t(get_eigenvalue(pca)) %>%
  as.data.frame(rownames = "Feature") %>%
  rename_all(~str_replace(.,"Dim.","PC")) %>%
  round(digits = 2)

# plot pca
as_tibble(pca$x, rownames = NA) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(targets.GMP, by = "ID") %>%
  ggplot(data = ., aes(PC1, PC2, group = ID, color = Cell))+
  geom_point(aes(shape = Purification), size = 3) +
  theme_classic() + 
  #geom_text_repel(aes(label=ID), size = 4, segment.alpha = 0.5) +
  xlab(paste("PC1 ",toString(variance.perc[2, "PC1"]),"%")) +
  ylab(paste("PC2 ",toString(variance.perc[2, "PC2"]), "%")) +
  labs(subtitle = "chrXY removed")

# GO TERM ANALYSIS FOR TOP PCs ----
# extract contribution weights for first 5 PC
cpg.contrib = get_pca_var(pca)$contrib[,1:5]

# subset annotation file in order to get list of genes
anno.sub = anno[anno$Name %in% rownames(cpg.contrib), c("Name", "UCSC_RefGene_Name", "chr","pos")]

# check that cpg order matches before cbind
all(anno.sub$Name == rownames(cpg.contrib))

cpg.contrib = as.data.frame(cbind(cpg.contrib, anno.sub))

# extract genes for top 0.1% contributing cpgs for PC1 and PC2
topPC1.genes = dplyr::filter(cpg.contrib, Dim.1 >= quantile(cpg.contrib$Dim.1, 0.999),
                             UCSC_RefGene_Name != "")
  
topPC1.genes = stringr::word(topPC1.genes$UCSC_RefGene_Name, start = 1, sep = ";") %>%
  unique()

topPC2.genes = dplyr::filter(cpg.contrib, Dim.2 >= quantile(cpg.contrib$Dim.2, 0.999),
                             UCSC_RefGene_Name != "")

topPC2.genes = stringr::word(topPC2.genes$UCSC_RefGene_Name, start = 1, sep = ";") %>%
  unique()

# run GO term analysis on top PC1 and PC2 genes
pc1.go = enrichGO(gene = topPC1.genes, 
                  OrgDb = org.Hs.eg.db,
                  readable = F, 
                  keyType = "SYMBOL",
                  ont = "BP")

pc2.go = enrichGO(gene = topPC2.genes, 
                  OrgDb = org.Hs.eg.db,
                  readable = F, 
                  keyType = "SYMBOL",
                  ont = "BP")

# simplify terms
pc1.simplified = simplify(pc1.go, cutoff = 0.7, by = "p.adjust", select_fun=min)
pc2.simplified = simplify(pc2.go, cutoff = 0.7, by = "p.adjust", select_fun=min)

# plot GO term results
dotplot(pc1.simplified, showCategory = 5) +
  labs(title = "Top contributors to PC1")

dotplot(pc2.simplified, showCategory = 5) +
  labs(title = "Top contributors to PC2")

# CORRELATION HEATMAP ----
# get top 10% most variable CpGs
beta.most.variable = beta.GMP[rowSds(beta.GMP) >= quantile(rowSds(beta.GMP), 0.90),]

# perform pairwise correlation
cor.mat = cor(beta.most.variable, method = "pearson")

# define heatmap annotation rows
annotation = HeatmapAnnotation(Cell_Type = targets.GMP$Cell,
                               Disease = targets.GMP$Disease,
                               Treatment = targets.GMP$Treatment,
                               Purification = targets.GMP$Purification,
                               Sex = targets.GMP$pred_sex,
                               col = list(Cell_Type = c("Neg_Fraction" = "Dark Gray",
                                                        "CD4.Treg" = "Lime Green",
                                                        "CD8.Treg" = "Purple"),
                                          Purification = c("Bead" = "Orange",
                                                           "Sorted" = "Black")))

# plot heatmap
set.seed(56)
Heatmap(cor.mat, name = "r", top_annotation = annotation,
        column_title = "Pearson Correlation of Top 10% Most Variable CpGs",
        show_row_names = T, show_column_names = T,
        row_labels = targets.GMP$User_ID,
        column_labels = targets.GMP$User_ID,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))