
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

View(rds_obj@meta.data)
rds_obj <- UpdateSeuratObject(rds_obj)

VlnPlot(rds_obj, features = "nFeature_RNA")


VlnPlot(rds_obj, features = "nCount_RNA")

VlnPlot(rds_obj, features = "percent.mt")



ElbowPlot(rds_obj)




#Clustering

rds_obj <- FindNeighbors(rds_obj, dims = 1:15)


# Setting Resolution
rds_obj <- FindClusters(rds_obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

View(rds_obj@meta.data)

# Clustering
cluster <- DimPlot(rds_obj, group.by = "RNA_snn_res.0.35", label = TRUE)


#Clustering by tissue
tissue <- DimPlot(rds_obj, group.by = "tissue", label = TRUE)

#Combine plost
fig2ab <- tissue|cluster

#Save the figure
ggsave("figures/Figure2ab.png", fig2ab, width = 12, height = 5.65)
