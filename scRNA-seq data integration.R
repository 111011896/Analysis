# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = "/Users/uni1896/easy_r/data", recursive = F, full.names = F)

for (x in dirs) {
  name <- gsub("MGI0369_1_", "", x)
  
  pbmc <- ReadMtx(mtx = paste0("/Users/uni1896/easy_r/data/", x, "/matrix.mtx.gz"),
                  features = paste0("/Users/uni1896/easy_r/data/", x, "/features.tsv.gz"),
                  cells = paste0("/Users/uni1896/easy_r/data/", x, "/barcodes.tsv.gz"))
  
  assign(name, CreateSeuratObject(counts = pbmc))
}

# merge datasets
merged_seurat <- merge(`SLAB-C1_C1`, y = c(`SLAB-145_0`, `SLAB-145_7`),
      add.cell.ids = ls()[4:6],
      project = "SLAB")

# QC & filtering
View(merged_seurat@meta.data)
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data,
                                    col = 'sample',
                                    into = c("Patient", "Type", "Barcode"), sep = "_")

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                          nFeature_RNA < 5000 &
                          percent.mt < 15)

# perform integration to correct for batch effects
obj.list <- SplitObject(merged_seurat, split.by = "Patient")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
}

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

seurat.integrated <- IntegrateData(anchorset = anchors)

# scale data, run PCA and UMAP and visualize integrated data
DefaultAssay(seurat.integrated) <- "integrated"
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)

p1 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Type")
p2 <- DimPlot(seurat.integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2
