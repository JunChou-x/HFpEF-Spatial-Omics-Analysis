library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(plyr)
library(patchwork)

control1 <- Read10X("./Control1")
control2 <- Read10X("./Control2")
hfpef1 <- Read10X("./HFpEF1")
hfpef2 <- Read10X("./HFpEF2")

# Create a Seurat object
control1_obj <- CreateSeuratObject(counts = control1, project = "Control", min.cells = 3, min.features = 200)
control2_obj <- CreateSeuratObject(counts = control2, project = "Control", min.cells = 3, min.features = 200)
hfpef1_obj <- CreateSeuratObject(counts = hfpef1, project = "HFpEF", min.cells = 3, min.features = 200)
hfpef2_obj <- CreateSeuratObject(counts = hfpef2, project = "HFpEF", min.cells = 3, min.features = 200)
control1_obj$condition <- "Control"
control2_obj$condition <- "Control"
hfpef1_obj$condition <- "HFpEF"
hfpef2_obj$condition <- "HFpEF"
control1_obj$sample <- "Control1"
control2_obj$sample <- "Control2"
hfpef1_obj$sample <- "HFpEF1"
hfpef2_obj$sample <- "HFpEF2"

Merge <- merge(control1_obj, y = list(control2_obj, hfpef1_obj, hfpef2_obj), 
               add.cell.ids = c("Control1", "Control2", "HFpEF1", "HFpEF2"), 
               project = "HFpEF_Study")
Merge <- readRDS("./Merge.rds")
Merge[["percent.mt"]] <- PercentageFeatureSet(Merge, pattern = "^mt-")
VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
summary(Merge@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")])

# QC
Merge <- subset(Merge, subset = percent.mt < 10 & 
                  nFeature_RNA > 500 & nFeature_RNA < 6000 & 
                  nCount_RNA > 1000 & nCount_RNA < 50000)
Merge <- NormalizeData(Merge, normalization.method = "LogNormalize", scale.factor = 10000)
Merge <- FindVariableFeatures(Merge, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Merge), 10)
plot5 <- LabelPoints(VariableFeaturePlot(Merge), points = top10, repel = TRUE)
plot5

all.genes <- rownames(Merge)
Merge <- ScaleData(Merge, features = all.genes)
Merge <- RunPCA(Merge, features = VariableFeatures(object = Merge))
Merge <- RunHarmony(Merge, group.by.vars = "sample")  
VizDimLoadings(Merge, dims = 1:2, reduction = "pca")
DimPlot(Merge, reduction = "pca") + NoLegend()
DimHeatmap(Merge, dims = 1:10, cells = 500, balanced = TRUE)

ElbowPlot(Merge)

Merge <- FindNeighbors(Merge, dims = 1:20)
Merge <- FindClusters(Merge, resolution = 0.5)
Merge <- RunUMAP(Merge, reduction = "harmony", dims = 1:20)
color=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
        "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
        "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
        "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
DimPlot(Merge, reduction = "umap",cols = color)


umap_coords <- Embeddings(Merge, reduction = "umap")
cluster_centers <- aggregate(umap_coords, by = list(Merge$seurat_clusters), FUN = mean)

DimPlot(Merge, reduction = "umap", cols = color) +
  geom_text(data = cluster_centers, aes(x = umap_1, y = umap_2, label = Group.1), 
            size = 5, color = "black") +
  theme(legend.position = "none") 

sweep.res.list <- paramSweep(Merge, PCs = 1:20, sct = FALSE)
sweep.stats_Merge <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn_Merge <- find.pK(sweep.stats_Merge)
bcmvn_Merge
pK_Merge <- bcmvn_Merge %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK_Merge <- as.numeric(as.character(pK_Merge[[1]]))
ggplot(bcmvn_Merge, aes(pK, Bcmetric, group = 1)) +
  geom_point() +
  geom_line()
annotations_Merge <- Merge@meta.data$seurat_clusters
homotypic.prop_Merge <- modelHomotypic(annotations_Merge)
nExp_poi_Merge <- round(0.076 * nrow(Merge@meta.data)) 
nExp_poi.adj_Merge <- round(nExp_poi_Merge * (1 - homotypic.prop_Merge))

Merge <- doubletFinder(Merge,
                       PCs = 1:20,
                       pN = 0.25,
                       pK = pK_Merge,
                       nExp = nExp_poi.adj_Merge,
                       reuse.pANN = FALSE, sct = FALSE
)


Merge$celltype <- plyr::mapvalues(
  Merge$cluster_0.5,
  from = 0:23,
  to = c(
    "Cardiomyocyte",   # 0
    "Endothelial",     # 1
    "Fibroblast",      # 2
    "Cardiomyocyte",   # 3
    "Cardiomyocyte",   # 4
    "Cardiomyocyte",   # 5
    "Cardiomyocyte",   # 6
    "Endothelial",     # 7
    "Macrophage",      # 8
    "Pericyte",        # 9
    "Endothelial",     # 10
    "Fibroblast",      # 11
    "Macrophage",      # 12
    "Cardiomyocyte",   # 13
    "Endothelial",     # 14
    "Cardiomyocyte",   # 15
    "Endothelial",     # 16
    "Macrophage",      # 17
    "T cell",          # 18
    "Cardiomyocyte",   # 19
    "Pericyte",        # 20
    "Smooth Muscle",   # 21
    "Pericyte",        # 22
    "Schwann Cell"     # 23
  )
)

saveRDS(Merge,"scObj.rds")

