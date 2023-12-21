argument = commandArgs(trailingOnly = TRUE)

input_folder = argument[1]
min_cells = argument[2]
min_features = argument[3]
num_variable_features = argument[4]
ndim = argument[5]
resolution = argument[6]
num_neighbours = argument[7]
mindist  = argument[8]

output_folder <- paste("mincell_",min_cells, "_minfeat_",min_features,"_ngenes_",num_variable_features,"_ndim_",ndim,"_res_",resolution,"_nneigh_",num_neighbours,"_mindist_",mindist,sep="")
dir.create(output_folder)

library(dplyr)
library(Seurat)
library(patchwork)

# load and filter the data
pbmc.data <- Read10X(data.dir = input_folder)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


# normalize the data
pbmc <- NormalizeData(pbmc)

#find most variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Run PCA and produce plots
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
jpeg(paste(output_folder,"/PCA_plot.jpg",sep=""))
DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

jpeg(paste(output_folder,"/Elbow_plot.jpg",sep=""))
ElbowPlot(pbmc)
dev.off()

# cluster cells
pbmc <- FindNeighbors(pbmc, dims=1:ndim)
pbmc <- FindClusters(pbmc, resolution = as.double(resolution))


# run umap analysis and print plot
pbmc <- RunUMAP(pbmc, dims = 1:ndim, n_neighbohrs = num_neighbours, min_dist = mindist )

jpeg(paste(output_folder,"/UMAP_plot.jpg",sep=""))
DimPlot(pbmc, reduction = "umap")
dev.off()

# Find markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
    
write.csv(pbmc.markers,paste(output_folder,"/markers.csv",sep=""))


# To save it
jpeg(paste(output_folder,'/umap_violins.jpg',sep=""), height=3500, width=3500)
# Plot what you need
VlnPlot(pbmc, features = c("CD34-14176", "PPBP-24138", "HBB-33622", "MKI67-28355", "MS4A7-12103", "FCGR3A-13772", "CD14-18823", "LYZ-23572", "CST3-27071", "XCR1-10671", "CLEC9A-20989", "FCER1A-88227", "CD160-89972", "KLRC1-28394", "KLRD1-89433", "NCR1-4503", "NCR1-26566", "NKG7-23420", "GNLY-88176", "FOXP3-87523", "CCR7-89746", "IL7R-3347", "CD8A-1152", "CD3E-16113", "CD38-19506", "MZB1-11368", "CD79B-90964", "CD79A-87745", "CD19-28252", "MS4A18-33827", "MS4A1-14568"), slot = "counts", log = TRUE)
dev.off()



#seurat_tests/mincell_10_minfeat_200_ngenes_1500_ndim_5_res_0.1_nneigh_50_mindist_0/num_used_cells.txt
