### load required packages
library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(tidyverse)
library(grid)
library(hdf5r)

### gzip the files----

# write the name of the file you want to zip
my_file <- "AT-01-Br/filtered_feature_bc_matrix/barcodes.tsv"

#write the intended name of the gz
gz_file_name <- "AT-01-Br/filtered_feature_bc_matrix/barcodes.tsv.gz"

# Open a connection to the compressed file using gzfile
gz_file <- gzfile(gz_file_name, "w")

# Read the contents of the original file and write them to the compressed file
writeLines(readLines(my_file), gz_file)

# Close the compressed file connection
close(gz_file)

### Set working directory
setwd("C:/Users/Yp/Desktop/Singlecell transcriptomics/Brain")

### Load the dataset----
brain.data <- Read10X(data.dir = "AT-01-Br/filtered_feature_bc_matrix/")


### seuratobject----
# Initialize the Seurat object with the raw (non-normalized data).
brain <- CreateSeuratObject(counts = brain.data, project = "brain_2020",
                            min.cells = 3, 
                            min.features = 200)
brain


#QC----
#### find the % mitochondrial reads
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")


# Visualize QC metrics as a violin plot
VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


brain <- subset(brain, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

brain

###normalizing data----
brain <- NormalizeData(brain,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

## Access normalised data
normalized_data <- brain[["RNA"]]@data


####head it
head(normalized_data)



#highly variable features-----
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(brain), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(brain)
plot2 <- LabelPoints(plot = plot1, 
                     points = top10,
                     repel = TRUE, 
                     xnudge = 0,
                     ynudge = 0,
                     max.overlaps = 20)

plot1 + plot2


## Scaling data----
#all.genes <- rownames(brain)
#brain <- ScaleData(brain, features = all.genes)

brain <- ScaleData(brain, vars.to.regress = "percent.mt")


### do linear dimensionality reduction using PCA----
brain <- RunPCA(brain, features = VariableFeatures(object = brain))


#### Examine and visualize PCA results a few different ways
print(brain[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(brain, dims = 1:2, reduction = "pca")


#### plot the PCA
DimPlot(brain, reduction = "pca")

DimHeatmap(brain, dims = 1:15, cells = 200, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
brain <- JackStraw(brain, num.replicate = 100)
brain <- ScoreJackStraw(brain, dims = 1:20)

### determine dimensionality
JackStrawPlot(brain, dims = 1:15)

##elbow plot
ElbowPlot(brain)#13 PCs

####clustering ----
brain <- FindNeighbors(brain, dims = 1:13)
brain <- FindClusters(brain, resolution = c(0.1, 0.3, 0.5, 0.8, 1))


#check metadata
View(brain@meta.data)


##visualize the clusters----
 #0.3
DimPlot(object = brain,
        group.by = "RNA_snn_res.0.3",
        label = TRUE)

#0.5

DimPlot(object = brain,
        group.by = "RNA_snn_res.0.5",
        label = TRUE)


#0.8
DimPlot(object = brain,
        group.by = "RNA_snn_res.0.8",
        label = TRUE)

## setting identities of cells
Idents(brain)<- "RNA_snn_res.0.3"


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
brain <- RunUMAP(brain, dims = 1:13)



# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(brain, reduction = "umap", label = T)

###save the project
saveRDS(brain, file = "C:/Users/Yp/Desktop/Singlecell transcriptomics/output")


##### Find the deferentially expressed genes
###one cluster versus 
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
brain.markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
brain.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

View(top10)



### Visualize using heatmap
brain.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(brain, features = top10$gene) + NoLegend()


### visualize on the UMAP plot

FeaturePlot(brain, features = c("SLC1A3", "SSH2", "PRKG1", "CLDN5", "SKAP1", "SNTG1", 
                                "MS4A4E", "BMPER"))

saveRDS(brain, file = "brain_final.rds")











