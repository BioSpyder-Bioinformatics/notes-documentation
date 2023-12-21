# Seurat analysis

Seurat is a R tool used to
- Get expression values in cells
- Scale their expression data
- Make PCA analysis
- Clusters cells based on PCA with T-SNE/U-MAP techniques
	+ Similar to PCA but it does not assume linear relation
- Allows to check overexpressed genes in clusters/shows markers between clusters
- Differential expression
	+ Check variability within groups against group variability

Steps:
- Load data
- Normalise
- Run PCA
	+ Determine the # of PCs based on ElbowPlot
- Find neighbours-clusters and run UMAP based on dimensionality decided above



### Preparation
Need to run the script tbl2mtx to obtain the 3 files required by seurat eg barcodes.tsv genes.tsv matrix.mtx


### From tutorial

*Load data in*
Start reading the data in as 10X data, then create the seurat object with the counts
```R
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./matrix_out/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 10, min.features = 2000)
pbmc
```
_Parameters_
- min.cells -> Keep probes being epressed at least in # cells (reduces noise later) -> preferred 10
- min.features -> min # of genes expressed in a cell -> dependent on the type of cells being sequenced; for EPA which was massive they used 5000



*QC*
_THIS IS NOT REQUIRED FOR OUR DATA, HERE FOR REFERENCE!!!!!!!_

There are then a few QC steps, in the tutorial they wiolin plot the rna features, count and the percentage of mitochondrial (^MT-) genes (in our data they shouldn't be as prevalent as there are no probes on them)
It also plots feature scatter plots.
```R
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


### !!!!!!!!!!!!! THIS IS THEM FILTERING BASED ON THE SHAPE OF THE DATA, excluding few-too many features and cells with >5% mtRNA
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```





*Normalisation of the data*
By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data
```R
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Which is pretty much the same as this, given that above are parameters set by default
pbmc <- NormalizeData(pbmc)
```



*Feature selection (identification of highly variable features)*
Procedure described in (https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue)
```R
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)  # nfeatures is the max # of genes that are checked for variability, the method is also standard

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
#plot with labels
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```




*Data Scaling*
Use the function ScaleData to scale the expression so that the mean is 0 and the variance across cells is 1. This is done to give equal weight to all genes in downstream analysis

```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```



*Linear dimensionality reduction*
First we calculate a PCA, we then have a look at the features that dictate each PC, we can also visualise it in plots
```R
# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways

# Shows first 5 genes of first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


#Plot the first 2 PCs dimensionality reduction
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")


# Plot actual PCA
DimPlot(pbmc, reduction = "pca")

# Look at heatmaps of PCs variations
# Only 1
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# All 15
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```



*Determine dimensionality of dataset*
We basically want to determine how many PCs we want to include in the analysis based on their p-value. This is done to avoid introducing excessive background noise.
In this snippet there are 2 approaches, one with the JackStraw procedure which is a bootstrap-like approach, the other is the classic elbowplot to STD which is much faster.

```R
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

# Easier
ElbowPlot(pbmc)
```


*Cluster the cells*
```R
# The dims were determined before
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# Resolution might increase with larger datasets
pbmc <- FindClusters(pbmc, resolution = 0.5) # We want to tune this to get 10-20 clusters

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```


*Run UMAP*
Use the dimensionality reduction you used before!
```R
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
```


*SAVE A BACKUP OF SESSION FILE*
```R
# Exporting the UMAP coordinates 
Write.csv(Embeddings(epa, reduction='umap'), "./umap_coordinates.csv") 
 
# Exporting cluster information 
Write.csv(epa@meta.data, "./metadata.csv") 

# Save session file
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
```


*Find differentially expressed features (cluster biomarkers)*
You can do it in a more controlled manner with FindMarkers or just run everything vs everything with FindAllMarkers

min.pct flag requires a feature to be detected at a minimum percentage in eitherof the 2 groups of cells, while the thresh.test argument requires a feature to be differentiall expressed by some amount between the two groups. 

```R
### ------ Ignore - specific only

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
### -------------------------------

# THE ONE OF INTEREST!!
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)




# Check 'classification power' (eg how determining the marker is) for any individual marker

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# Could also use RidgePlot(), CellScatter(), and DotPlot()


# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)



# COLOUR UMAP PLOTS BASED ON FEATURES (need to determine your ones though)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
# FeaturePlot(pbmc, features = c('DDX55-1828', 'GCLC-22866', 'TRBV2-88126', 'BRMS1-24014', 'CARHSP1-26965', 'SCPEP1-90417', 'SLFN11-15461', 'ATXN7L1-26874', 'ATXN7L1-26874', 'MMP9-15866', 'CTSH-1644'))


#DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 10 markers (or all markers if less than 20) for each cluster.
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

```



*Assign cell type identity to clusters*
```R
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#To save session
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
```











# Streamlined script

```R
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_cd3_cd19_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#------------

# Check amount of features
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#------------

# Subset the data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#--------------

# Normalise the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find 2000 most variable features and plot them
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#---------------

# Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) #npcs=50
# Check 5 most variable features per PC
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Get PCA
DimPlot(pbmc, reduction = "pca")

#---------------

# Check elbowplot
ElbowPlot(pbmc)


#--------------
# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Feed PCs to the UMAP Calculation
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# Get 'deterministic' markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

```


FeaturePlot(pbmc, features = c("IGHM-89466", "CD79A-87745", "MS4A1-14568", "TRAC-89080", "LAPTM5-23860", "CST7-26199", "CD1C-18380", "IL7R-3347", "CD3D-1127"))






kawaki_markers = ['CD34', 'PPBP', 'HBB', 'MKI67', 'MS4A7', 'FCGR3A', 'CD14', 'LYZ', 'CST3', 'XCR1', 'CLEC9A', 'FCER1A', 'CD160', 'KLRC1', 'KLRD1', 'NCR1', 'NKG7', 'GNLY', 'FOXP3', 'CCR7', 'IL7R', 'CD8A', 'CD3E', 'CD38', 'MZB1', 'CD79B', 'CD79A', 'CD19', 'MS4A1', 'MS4A1']

c("CD34-14176", "PPBP-24138", "HBB-33622", "MKI67-28355", "MS4A7-12103", "FCGR3A-13772", "CD14-18823", "LYZ-23572", "CST3-27071", "XCR1-10671", "CLEC9A-20989", "FCER1A-88227", "CD160-89972", "KLRC1-28394", "KLRD1-89433", "NCR1-4503", "NCR1-26566", "NKG7-23420", "GNLY-88176", "FOXP3-87523", "CCR7-89746", "IL7R-3347", "CD8A-1152", "CD3E-16113", "CD38-19506", "MZB1-11368", "CD79B-90964", "CD79A-87745", "CD19-28252", "MS4A18-33827", "MS4A1-14568")


# you can plot raw counts as well
VlnPlot(pbmc, features = c("CD34-14176", "PPBP-24138", "HBB-33622", "MKI67-28355", "MS4A7-12103", "FCGR3A-13772", "CD14-18823", "LYZ-23572", "CST3-27071", "XCR1-10671", "CLEC9A-20989", "FCER1A-88227", "CD160-89972", "KLRC1-28394", "KLRD1-89433", "NCR1-4503", "NCR1-26566", "NKG7-23420", "GNLY-88176", "FOXP3-87523", "CCR7-89746", "IL7R-3347", "CD8A-1152", "CD3E-16113", "CD38-19506", "MZB1-11368", "CD79B-90964", "CD79A-87745", "CD19-28252", "MS4A18-33827", "MS4A1-14568"), slot = "counts", log = TRUE)


# To save it
jpeg('umap_violins.jpg', height=3000, width=3000)
# Plot what you need
VlnPlot(pbmc, features = c("CD34-14176", "PPBP-24138", "HBB-33622", "MKI67-28355", "MS4A7-12103", "FCGR3A-13772", "CD14-18823", "LYZ-23572", "CST3-27071", "XCR1-10671", "CLEC9A-20989", "FCER1A-88227", "CD160-89972", "KLRC1-28394", "KLRD1-89433", "NCR1-4503", "NCR1-26566", "NKG7-23420", "GNLY-88176", "FOXP3-87523", "CCR7-89746", "IL7R-3347", "CD8A-1152", "CD3E-16113", "CD38-19506", "MZB1-11368", "CD79B-90964", "CD79A-87745", "CD19-28252", "MS4A18-33827", "MS4A1-14568"), slot = "counts", log = TRUE)
dev.off()



VlnPlot(pbmc, features = c('MS4A1-14568'), slot = "counts", log = TRUE)






Standard salvo
Rscript seurat_script.R ../../filtered_cd3_matrix 10 200 1500 5 0.1 50 0

standard 10x
Rscript ../seurat_script_modified.R ../../filtered_pbmc_no_cult_matrix 3 200 2000 10 0.5 30 0.3
# not working!!




