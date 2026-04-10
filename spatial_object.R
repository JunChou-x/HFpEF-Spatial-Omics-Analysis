library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(harmony)
library(BPCells)
library(plyr)
library(forcats)
library(paletteer)
library(arrow)
library(rjson)
library(future)

# Set global variables
options(future.globals.maxSize = 10 * 1024^3)
source("AuxFunctions.R")


# Sample configuration
base_path <- "./data/spatial_data"
output_path <- "./results"

SampleInfo <- data.frame(
  Sample = c("Control_Lower", "Control_Upper", "HFpEF_Lower", "HFpEF_Upper"),
  H5 = file.path(base_path, c("./filtered_feature_bc_matrix.h5", 
                              "./filtered_feature_bc_matrix.h5", 
                              "./filtered_feature_bc_matrix.h5", 
                              "./filtered_feature_bc_matrix.h5"))
)


# Data loading
Matrices <- vector("list", length = nrow(SampleInfo))
names(Matrices) <- SampleInfo$Sample

for (jj in 1:nrow(SampleInfo)) {
  
  MatData <- open_matrix_10x_hdf5(path = SampleInfo$H5[jj])
  
  temp_dir <- file.path("./temp_bpcells", SampleInfo$Sample[jj])
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  
  write_matrix_dir(mat = MatData, dir = temp_dir, overwrite = TRUE)
  mat <- open_matrix_dir(dir = temp_dir)
  
  colnames(mat) <- paste0(colnames(mat), "_", SampleInfo$Sample[jj])
  
  features_path <- file.path(dirname(SampleInfo$H5[jj]), "features.tsv.gz")
  if (!file.exists(features_path)) stop(paste("Missing features file:", features_path))
  
  Genes <- read.delim(features_path, sep = "\t", header = FALSE)
  rownames(mat) <- make.unique(Genes$V2[match(rownames(mat), Genes$V1)])
  
  Matrices[[jj]] <- mat
}

merged.object <- CreateSeuratObject(counts = Matrices)
merged.object <- JoinLayers(merged.object)
merged.object$sample <- gsub(".*_([^_]+_[^_]+)$", "\\1", colnames(merged.object))

# QC
merged.object[["percent.mt"]] <- PercentageFeatureSet(merged.object, pattern = "^mt-")
merged.object <- subset(merged.object, subset = percent.mt < 30)

# SketchData 
merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object)

merged.object <- SketchData(
  object = merged.object, 
  ncells = 300000, 
  method = "LeverageScore", 
  sketched.assay = "sketch"
)

DefaultAssay(merged.object) <- "sketch"

merged.object <- FindVariableFeatures(merged.object) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by = "sample") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 2) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, return.model = TRUE)

print(table(Idents(merged.object)))

# Cell type annotation
cluster_to_l1 <- c(
  "0"="Cardiomyocyte", "1"="Cardiomyocyte", "2"="Cardiomyocyte", "3"="Cardiomyocyte",
  "4"="Fibroblast", "5"="Cardiomyocyte", "6"="Cardiomyocyte", "7"="Cardiomyocyte",
  "8"="Cardiomyocyte", "9"="Endothelial Cell", "10"="Cardiomyocyte", "11"="Cardiomyocyte",
  "12"="Macrophage", "13"="Cardiomyocyte", "14"="Endothelial Cell", "15"="Smooth Muscle Cell",
  "16"="Macrophage", "17"="Cardiomyocyte", "18"="Macrophage", "19"="Fibroblast",
  "20"="Smooth Muscle Cell", "21"="Mesothelial Cell", "22"="Pericyte", "23"="Endothelial Cell",
  "24"="Cardiomyocyte", "25"="Cardiomyocyte", "26"="Fibroblast", "27"="Fibroblast",
  "28"="Smooth Muscle Cell", "29"="Doubled", "30"="Endothelial Cell", "31"="Adipocyte",
  "32"="Cardiomyocyte", "33"="Fibroblast", "34"="Adipocyte", "35"="Fibroblast",
  "36"="Doubled", "37"="Myofibroblast", "38"="Doubled", "39"="Macrophage",
  "40"="Endothelial Cell", "41"="Doubled", "42"="Low-quality cells", "43"="Low-quality cells"
)

merged.object$Level1 <- plyr::mapvalues(
  merged.object$sketch_snn_res.2, 
  from = names(cluster_to_l1), 
  to = cluster_to_l1
)


Idents(merged.object) <- "Level1"
all_levels <- levels(factor(merged.object$Level1))
MetaDataSubClusters <- list()

for(ClusterID in all_levels) {
  message("Processing: ", ClusterID)
  Subset <- subset(merged.object, idents = ClusterID) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:25) %>%
    FindClusters(resolution = 1)
  
  clean_id <- gsub(" ", "", ClusterID)
  Subset$Level2 <- paste0(clean_id, "_", Subset$seurat_clusters)
  MetaDataSubClusters[[ClusterID]] <- Subset@meta.data
}

# Merge metadata
combined_meta <- do.call(rbind, MetaDataSubClusters)
combined_meta$Barcode <- sapply(strsplit(rownames(combined_meta), "[.]"), `[`, 2)
merged.object$Level2 <- combined_meta$Level2[match(colnames(merged.object), combined_meta$Barcode)]


# Label cleaning
merged.object$Level2 <- as.character(merged.object$Level2)
merged.object$Level2[merged.object$Level2 %in% c("Doubled_0", "Doubled_1")] <- "Doubled"
merged.object$Level2[merged.object$Level2 %in% c("Low-qualitycells_0", "Low-qualitycells_1")] <- "Low-qualitycells"
merged.object$Level2[merged.object$Level2 == "MesothelialCell_0"] <- "MesothelialCell"
merged.object$Level2[merged.object$Level2 == "Myofibroblast_0"] <- "Myofibroblast"


ordered_levels <- sort(unique(merged.object$Level2))
merged.object$Level2 <- factor(merged.object$Level2, levels = ordered_levels)

print(table(merged.object$Level2))


# ProjectData
merged.object <- ProjectData(
  object = merged.object,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "harmony",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(L1 = "Level1", L2 = "Level2")
)




# Subgroup subdivision and annotation (using Cardiomyocyte as an example)
sub_obj <- subset(merged.object, subset = Level1 == "Cardiomyocyte")
DefaultAssay(sub_obj) <- "sketch"
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
sub_obj <- ScaleData(sub_obj)
sub_obj <- RunPCA(sub_obj)
sub_obj <- RunHarmony(sub_obj, group.by = "sample")
sub_obj <- FindNeighbors(sub_obj, dims = 1:25)
sub_obj <- FindClusters(sub_obj, resolution = 1)
sub_markers <- FindAllMarkers(sub_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cell_type_map <- c(
  "0" = "vCMs",
  "1" = "aCMs",
  "2" = "Endothelial", 
  "7" = "Pericytes",
  "8" = "Stress CMs"
  # ... 
)
sub_obj$L3 <- unname(cell_type_map[as.character(sub_obj$sketch_snn_res.0.5)])
if(!"Level3" %in% colnames(merged.object@meta.data)) {
  merged.object$Level3 <- NA
}
source_annotations <- sub_obj$L3
names(source_annotations) <- Cells(sub_obj)
new_annotations_aligned <- source_annotations[Cells(merged.object)]
merged.object$Level3 <- ifelse(
  is.na(merged.object$Level3), 
  new_annotations_aligned, 
  merged.object$Level3
)
names_to_standardize <- c(
  "Atrial Cardiomyocyte", "Ventricular Cardiomyocyte",
  "Mixed (Endo-Vwf + aCM + Fib)", "Endothelial (CM-mixed)",
  "Low-quality cells"
)
new_standard_names <- c(
  "Atrial Cardiomyocytes", "Ventricular Cardiomyocytes",
  "Endothelial Cells", "Endothelial Cells",
  "Low-Quality Cells"
)
merged.object$Level5 <- mapvalues(
  merged.object$Level4, 
  from = names_to_standardize, 
  to = new_standard_names,
  warn_missing = FALSE
)

# Eliminate low-quality and contaminated cells
Idents(merged.object) <- "L5"
merged.object <- subset(
  merged.object, 
  idents = c("Low-Quality Cells", "RBC Contamination"), 
  invert = TRUE
)
my_levels_order <- c(
  "Ventricular Cardiomyocytes",
  "Atrial Cardiomyocytes",
  "Stressed Cardiomyocytes",
  "Conduction System Cardiomyocytes",
  "Ventricular-associated Fibroblasts",
  "Perivascular Fibroblasts",
  "Myofibroblasts",
  "Stress-Response Fibroblasts",
  "Atrial-associated Fibroblasts",
  "Resident Fibroblasts",
  "Vascular Endothelial Cells",
  "Endocardial Endothelial Cells",
  "Lymphatic Endothelial Cells",
  "Pericytes",
  "Ventricular-associated SMCs",
  "Arterial SMCs",
  "Atrial-associated SMCs",
  "B Cells",
  "Brown Adipocytes",
  "White Adipocytes",
  "Schwann Cells",
  "Neurons",
  "Fibro-Immune Interface",
  "Proliferating Cells",
  "Mesothelial Cells"
)
merged.object$L5 <- factor(merged.object$L5, levels = my_levels_order)
