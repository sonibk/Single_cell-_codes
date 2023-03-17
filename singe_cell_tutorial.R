#### Install required packages
install.packages("tidyverse")
BiocManager::install("SeuratObject")
install.packages("Seurat")
install.packages("hdf5r")


#### Load required packages
library(tidyverse)
library(SeuratObject)
library(Seurat)
library(hdf5r)
library(grid)


###read in the data in a count matrix
nsclc.sparse.m <- Read10X_h5(filename = "../20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5")


####Identifying the different modalities present
str(nsclc.sparse.m)
 
#### Getting the gene expression matrix
cts <- nsclc.sparse.m$`Gene Expression`


### initialize the seurat object with the raw(non-normalized data)
nsclc.Seurat.obj <- CreateSeuratObject(counts = cts, 
                                       project = "NSCLC", 
                                       min.cells = 3,
                                       min.features = 200)



####check the structure of the seurat object
str(nsclc.Seurat.obj)

#### View the metadata
View(nsclc.Seurat.obj@meta.data)


#1 QC filtering----
#### percentage of mitochondrial genes 
nsclc.Seurat.obj[["percent.mt"]] <- PercentageFeatureSet(object = nsclc.Seurat.obj,
                                                         pattern = "^MT-")
View(nsclc.Seurat.obj@meta.data)


###visualize with a violin plot - features are from metadata
VlnPlot(object = nsclc.Seurat.obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)


###plot a scatter plot of no.of genes(features) vs counts of the genes
FeatureScatter(object = nsclc.Seurat.obj,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               pt.size = 1.5 )+
  geom_smooth(method = "lm")


#2 Filtering----
####subset the Seurat object: filter out low #of genes and mT genes
nsclc.Seurat.obj <- subset(nsclc.Seurat.obj,subset = nFeature_RNA >200 & 
                             nFeature_RNA < 2500 & percent.mt > 5)


#3 Normalization----
## gene counts in a given cell/total gene counts for that cell x 
#scaling factor then log transform
nsclc.Seurat.obj <- NormalizeData(nsclc.Seurat.obj, 
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)



#4 finding variable features----
## Identifies features that are outliers on a 'mean variability plot
nsclc.Seurat.obj <- FindVariableFeatures(nsclc.Seurat.obj,
                                         selection.method= "vst",
                                         nfeatures = 2000)


### Identify the top 10 variable features
Top10 <- head(VariableFeatures(nsclc.Seurat.obj), 10)


###plot variable features with or without lables
plot1 <- VariableFeaturePlot(nsclc.Seurat.obj)


#### label your plot with the top 10 variable genes
### xnudge and ynudge help in making the plot fit in the limits
### Always load the grid package
LabelPoints(plot = plot1, 
            points = Top10,
            repel = TRUE,
            xnudge = 0,
            ynudge = 0)



###Scale data-----
all.genes <- row.names(nsclc.Seurat.obj)
nsclc.Seurat.obj <- ScaleData(nsclc.Seurat.obj, features = all.genes)



##perform linear dimensionality reduction------
###PCA
nsclc.Seurat.obj <- RunPCA(object = nsclc.Seurat.obj,
                           features = VariableFeatures( object = nsclc.Seurat.obj))


###Visualize the first 5 principle components
print(nsclc.Seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)


####visualize if the first priciple component is showing heterogeity

DimHeatmap(object = nsclc.Seurat.obj,
           dims = 1, 
           cells = 500,
           balanced = TRUE)



###Determine the dimentionality
ElbowPlot(nsclc.Seurat.obj) #15



###clustering-----
##finding neighbours
nsclc.Seurat.obj<- FindNeighbors(object = nsclc.Seurat.obj,
                                 dims = 1:15)

nsclc.Seurat.obj<- FindClusters(object = nsclc.Seurat.obj,
                                resolution = c(0.1, 0.5, 0.8, 1 ))

#check metadata
View(nsclc.Seurat.obj@meta.data)


##visualize the clusters
DimPlot(object = nsclc.Seurat.obj,
        group.by = "RNA_snn_res.0.1",
        label = TRUE)


##  too much granularity
DimPlot(object = nsclc.Seurat.obj,
        group.by = "RNA_snn_res.0.5",
        label = TRUE)


## setting identities of cells
Idents(nsclc.Seurat.obj)<- "RNA_snn_res.0.1"

Idents(nsclc.Seurat.obj)


### non-linear dimentionality reduction-----
###UMAP
##install UMAP
reticulate::py_install(packages = "umap_learn")


####do umap
nsclc.Seurat.obj<- RunUMAP(object = nsclc.Seurat.obj,
                           dims = 1:15)


####visualize
DimPlot(object = nsclc.Seurat.obj,
        reduction = "umap")
