
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
# Reference Paper: https://www.nature.com/articles/s41590-022-01229-8#Sec2
# Seurat_object RDS data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182275 


cd8.obj <- readRDS("data/GSE182275_seurat_object.rds")

# Update Seurat Object
cd8.obj <- UpdateSeuratObject(cd8.obj)

str(cd8.obj)

# View the metadata for gene counts
View(cd8.obj@meta.data)


####### 1. Quality Control of the Data ###########

#Calculate the % of mitochondrial genes
#In this data, already calculated
#cd8.obj[["percent.mt"]] <- PercentageFeatureSet(cd8.obj, pattern = "^MT-")

# Check the total number of RNA counts,unique gene counts and percent of mt RNA in the quality control step
VlnPlot(cd8.obj, features = "nFeature_RNA")
VlnPlot(cd8.obj, features = "nCount_RNA")
VlnPlot(cd8.obj, features = "percent.mt")

# Another quality control step is to look by a scatter plot
FeatureScatter(cd8.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

#Filter out the low quality data
#Cells with <500 or >2,500 detected genes or a mitochondrial read percentage >10 were discarded
cd8.obj <- subset(cd8.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)



############## 2.Normalization################

cd8.obj <- NormalizeData(cd8.obj, normalization.method = "LogNormalize")



############### 3. Feature Selection ####################


# Identification of highly variable features (feature selection)
#Calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).

cd8.obj <- FindVariableFeatures(cd8.obj, selection.method = "vst", nfeatures = 2000)




############## 4. Scaling Data ##############

all.genescd8 <- rownames(cd8.obj)
cd8.obj <- ScaleData(cd8.obj, features = all.genescd8)


############# 5. Perform PCA- a linear dimensional reduction##

cd8.obj <- RunPCA(cd8.obj, features = VariableFeatures(object = cd8.obj))

# Visualize PCA
DimPlot(cd8.obj, reduction = "pca") + NoLegend()


############# 6. Clustering #############

#Determine the ‘dimensionality’ of the dataset

#ElbowPlot is used to help determine the number of principal components (PCs)
#to use for downstream analyses such as clustering or visualization (e.g., UMAP or t-SNE).
#The goal of the ElbowPlot is to identify the "elbow point," which is the point where the rate of decrease sharply changes. 
#This point indicates where adding more components contributes marginally to explaining additional variance.

ElbowPlot(cd8.obj)

#constructs a KNN graph based on the specified principal components.
#dims = 1:15 indicates that the first 15 PCs are used to define the cell neighborhoods.
cd8.obj <- FindNeighbors(cd8.obj, dims = 1:15)


# Setting Clustering Resolution 
cd8.obj <- FindClusters(cd8.obj, resolution = 0.35)

View(cd8.obj@meta.data)


#Run non-linear dimensional reduction (UMAP)
cd8.obj <- RunUMAP(cd8.obj, dims = 1:15)



# Set clustering colors
my_cols <- c("#1F77B4","#A4DFF2","#FFBB78", "#FF7F0E", "#2CA02C", "#D62728", "#AC8F14", "#9467BD", 
             "#8C564B", "#BCBD22", "#17BECF", "blue")

DimPlot(cd8.obj, reduction = "umap")

DimPlot(cd8.obj, group.by = "RNA_snn_res.0.35", label = TRUE, cols = my_cols)


#Clustering by tissue
DimPlot(cd8.obj, group.by = "tissue", label = TRUE)

#Combine plost
fig2ab <- tissue|cluster


#Create feature Plot

FeaturePlot(cd8.obj, features = c("Cd69", "Itgae", "Klrg1"), cols = c("blue", "yellow", "red"))

install.packages('devtools')
devtools::install_github('immunogenomics/presto')

# find markers for every cluster compared to all remaining cells, report only the positive ones
cd8.markers <- FindAllMarkers(cd8.obj, only.pos = TRUE)


# Top 5 genes

top5 <- rds.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

top5

# Taking only first 12 cluster
top5 <- top5 %>% 
  dplyr::filter(cluster %in% 0:11)

clusters_to_keep <- unique(top5$cluster)
cd8.obj_subset <- subset(cd8.obj, idents = clusters_to_keep)

DoHeatmap(cd8.obj, features = top5$gene, group.bar.height = 0.01,size=3,combine = TRUE) 

p2 <- DoHeatmap(cd8.obj_subset, features = top5$gene, group.bar.height = 0.01,size=3,combine = FALSE) 

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

