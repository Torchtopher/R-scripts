library(dplyr)
library(Seurat)
library(patchwork)

# Load the data
pbmc.data <- Read10X(data.dir = "/home/chris/code/paul_data_anayalsis/data")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "hESCsRNASEQ", min.cells = 3, min.features = 200)

# Calculate the proportion of transcripts mapping to mitochondrial genes, probably not what we want here
# Wasn't sure what pattern to search for so just used mitochondrial
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# plot different QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# actually do the QC, really arbitrary here, just chose values that looked "good"
pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 5500 & percent.mt < 60)

# log normalize
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# finds features that are highly expressed in some cells and lowly expressed in others
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes from the above calculation
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scale data so that mean expression is 0, and variance of expression across cells is 1
# can't defend why we do this, but it is considered "standard"
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# run PCA compute to reduce dimensionality of data, lots of complicated math here
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# see Macosko et al if you actually care about the math behind this
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# closer to dashed line is better
# graph looks really bad, almost random, need to look into that or figure out if data has any real features
JackStrawPlot(pbmc, dims = 1:15)

# not really like an elbow, more like exponetial decay, not really what we are looking for here
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)

head(Idents(pbmc), 5)
# using PCs 1-10
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
pbmc <- FindClusters(pbmc, resolution = 0.5)

DimPlot(pbmc, reduction = "umap", label = TRUE)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# not sure why p=0?
#           p_val avg_log2FC pct.1 pct.2 p_val_adj
#PPAP2B      0  0.9186659 0.577 0.148         0
#IFI44L      0  0.8453088 0.528 0.124         0
#CNN3        0  1.2008168 0.911 0.425         0
#C1orf61     0  2.0903425 0.918 0.367         0
#NES         0  1.2651935 0.866 0.360         0


