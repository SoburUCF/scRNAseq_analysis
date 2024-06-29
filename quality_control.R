
## Install packages from Bioconductor using the the BiocManager::install() function
install.packages("BiocManager")

BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
BiocManager::install("multtest")
BiocManager::install("glmGamPoi")
BiocManager::install("Matrix")

## Install packages from CRAN

install.packages("Matrix")
install.packages("RCurl")
install.packages("scales")
install.packages("cowplot")
install.packages("Seurat")
install.packages("metap")
install.packages("tidyverse")

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(metap)
library(AnnotationHub)
library(ensembldb)
library(multtest)
library(glmGamPoi)


sessionInfo()

# Import data

rds_obj <- readRDS("data/GSE182275_seurat_object.rds")

str(rds_obj)

# View the metadata for gene counts
View(rds_obj@meta.data)

# Update Seurat Object

rds_obj <- UpdateSeuratObject(rds_obj)

# Quality Control
# Check the totoal number of RNA counts,unique gene counts and percent of mt RNA in the quality control step
VlnPlot(rds_obj, features = "nFeature_RNA")


VlnPlot(rds_obj, features = "nCount_RNA")

VlnPlot(rds_obj, features = "percent.mt")

# Another quality control step is to look by a scatter plot

FeatureScatter(rds_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')


#ElbowPlot is used to help determine the number of principal components (PCs)
#to use for downstream analyses such as clustering or visualization (e.g., UMAP or t-SNE).
#The goal of the ElbowPlot is to identify the "elbow point," which is the point where the rate of decrease sharply changes. 
#This point indicates where adding more components contributes marginally to explaining additional variance.

ElbowPlot(rds_obj)

#Normalization

rds_obj <- NormalizeData(rds_obj, normalization.method = "LogNormalize")

str(rds_obj)


#Clustering
#constructs a KNN graph based on the specified principal components.
#dims = 1:15 indicates that the first 15 PCs are used to define the cell neighborhoods.

rds_obj <- FindNeighbors(rds_obj, dims = 1:15)


# Setting Resolution
rds_obj <- FindClusters(rds_obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

View(rds_obj@meta.data)

# Clustering
my_cols <- c("#1F77B4","#A4DFF2","#FFBB78", "#FF7F0E", "#2CA02C", "#D62728", "#AC8F14", "#9467BD", 
             "#8C564B", "#BCBD22", "#17BECF", "blue")
cluster <- DimPlot(rds_obj, group.by = "RNA_snn_res.0.35", label = TRUE, cols = my_cols)


#Clustering by tissue
tissue <- DimPlot(rds_obj, group.by = "tissue", label = TRUE)

#Combine plost
fig2ab <- tissue|cluster

install.packages('devtools')
devtools::install_github('immunogenomics/presto')

# find markers for every cluster compared to all remaining cells, report only the positive ones
rds.markers <- FindAllMarkers(rds_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Top 5 genes

top5 <- rds.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

top5

# Taking only first 12 cluster
top5 <- top5 %>% 
  dplyr::filter(cluster %in% 0:11)

clusters_to_keep <- unique(top5$cluster)
rds_obj_subset <- subset(rds_obj, idents = clusters_to_keep)

DoHeatmap(rds_obj, features = top5$gene, group.bar.height = 0.01,size=3,combine = TRUE) 

p2 <- DoHeatmap(rds_obj_subset, features = top5$gene, group.bar.height = 0.01,size=3,combine = FALSE) 

p2 <- lapply(X = p2, FUN = function(x) x + 
               theme(plot.title = element_text(size = 8)) +
               theme(axis.title.y = element_text(size = 5)) +
               theme(axis.title.x = element_text(size = 5)) +
               theme(axis.text.y = element_text(size = 8)) +
               theme(legend.position = "none")  )

library(patchwork)

# Step 7: Combine the list of ggplot objects if necessary
if (is.list(p2)) {
  p2_combined <- wrap_plots(p2) + 
    scale_fill_gradientn(colors = c("lightblue", "white", "red4"))
} else {
  p2_combined <- p2 + 
    scale_fill_gradientn(colors = c("lightblue", "white", "red4"))
}

ggsave("figures/heatmap2.pdf", plot = p2_combined, dpi = 600)

