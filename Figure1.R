library(Seurat)
library(scattermore)
library(tidyverse)
library(data.table)
library(wesanderson)
library(patchwork)
library(RColorBrewer)
library(furrr)
library(paletteer)
library(arrow)
library(pheatmap)
library(RColorBrewer)
library(distances)
library(BPCells)
library(harmony)
library(plotly)
library(htmlwidgets)
library(monocle3)
library(SeuratWrappers)
library(scales)

# Set global variables
options(future.globals.maxSize = 10 * 1024^3)
source("AuxFunctions.R")

# Organize information ----------------------------------------------------------------

MetaData.merged <- merged.object@meta.data
umap_coords <- Embeddings(merged.object, reduction = "full.umap")
MetaData.merged <- cbind(MetaData.merged, umap_coords)
MetaData.merged$L2<-factor(MetaData.merged$L2,levels=sort(unique(MetaData.merged$L2)))
LabelsL2<-MetaData.merged %>% group_by(L2) %>% summarise(X=median(fullumap_1),Y=median(fullumap_2))
LabelsL1<-MetaData.merged %>% group_by(L1) %>% summarise(X=median(fullumap_1),Y=median(fullumap_2))

# Create color palettes
MetaData.merged$L1<-factor(MetaData.merged$L1,levels=sort(unique(MetaData.merged$L1)))
ColsL1<-paletteer::paletteer_d("ggsci::category10_d3")[1:length(levels(MetaData.merged$L1))]
names(ColsL1)<-levels(MetaData.merged$L1)

MetaData.merged$L2<-factor(MetaData.merged$L2,levels=sort(unique(MetaData.merged$L2)))
ColsL2<-paletteer::paletteer_d("ggsci::default_igv")[1:length(levels(MetaData.merged$L2))]
names(ColsL2)<-levels(MetaData.merged$L2)

GetFilteredSpatialData <- function(PATH, size="008um") {
  if(grepl("/outs$", PATH)) {
    base_path <- PATH
    full_path <- file.path(PATH, paste0("binned_outputs/square_", size))
  } else if(grepl("/outs/binned_outputs/square_", PATH)) {
    full_path <- PATH
    base_path <- dirname(dirname(PATH))
  } else {
    base_path <- file.path(PATH, "outs")
    full_path <- file.path(PATH, paste0("outs/binned_outputs/square_", size))
  }
  cat("Use base path:", base_path, "\n")
  cat("Use full path:", full_path, "\n")
  matrix_path <- file.path(full_path, "filtered_feature_bc_matrix.h5")
  if (!file.exists(matrix_path)) {
    matrix_path <- file.path(full_path, "filtered_feature_bc_matrix/matrix.mtx.gz")
    if (!file.exists(matrix_path)) {
      possible_h5 <- list.files(path = base_path, 
                                pattern = "filtered_feature_bc_matrix.h5$", 
                                recursive = TRUE, 
                                full.names = TRUE)
      if (length(possible_h5) > 0) {
        matrix_path <- possible_h5[1]
      } else {
        stop("The representation matrix file cannot be found. Please check the path.")
      }
    }
  }
  
  cat("Using the expression matrix path:", matrix_path, "\n")
  if (grepl(".h5$", matrix_path)) {
    counts <- Read10X_h5(filename = matrix_path)
  } else {
    mtx_dir <- dirname(matrix_path)
    counts <- Read10X(data.dir = mtx_dir)
  }
  
  expressed_barcodes <- colnames(counts)
  cat("Number of barcodes with expression:", length(expressed_barcodes), "\n")
  spatial_path <- file.path(full_path, "spatial")
  tissue_positions_path <- list.files(spatial_path, 
                                      pattern = "tissue_positions.*", 
                                      full.names = TRUE)
  
  if (length(tissue_positions_path) == 0) {
    stop("Organization location file not found")
  }
  
  cat("Use space file:", tissue_positions_path[1], "\n")
  
  if (grepl(".parquet$", tissue_positions_path[1])) {
    tissue_positions_df <- arrow::read_parquet(tissue_positions_path[1]) %>%
      dplyr::rename_with(~ c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
  } else {
    tissue_positions_df <- read.csv(tissue_positions_path[1], header = FALSE) %>%
      dplyr::rename_with(~ c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
  }
  
  cat("Total number of barcodes:", nrow(tissue_positions_df), "\n")
  
  filtered_positions <- tissue_positions_df %>% 
    dplyr::filter(barcode %in% expressed_barcodes)
  
  cat("Number of barcodes after filtering:", nrow(filtered_positions), "\n")
  
  path_scales <- file.path(spatial_path, "scalefactors_json.json")
  if (file.exists(path_scales)) {
    scales <- rjson::fromJSON(file = path_scales)
    
    bcs <- filtered_positions %>% 
      mutate(
        imagerow_scaled = imagerow * scales$tissue_lowres_scalef,
        imagecol_scaled = imagecol * scales$tissue_lowres_scalef,
        imagerow_scaled_round = round(imagerow * scales$tissue_lowres_scalef),
        imagecol_scaled_round = round(imagecol * scales$tissue_lowres_scalef),
        tissue = as.factor(tissue)
      )
  } else {
    warning("Scale factor file not found, using default value 1")
    bcs <- filtered_positions %>% 
      mutate(
        imagerow_scaled = imagerow * 1,
        imagecol_scaled = imagecol * 1,
        imagerow_scaled_round = round(imagerow * 1),
        imagecol_scaled_round = round(imagecol * 1),
        tissue = as.factor(tissue)
      )
  }
  
  path_clusters <- file.path(full_path, "/clusters.csv")
  if (file.exists(path_clusters)) {
    clusters <- read.csv(path_clusters)
    path_umap <- file.path(full_path, "/projection.csv")
    if (file.exists(path_umap)) {
      umap <- read.csv(path_umap)
      bcs <- bcs %>%
        left_join(clusters, by = c(barcode = "Barcode")) %>%
        left_join(umap, by = c(barcode = "Barcode"))
    }
  }
  
  return(bcs)
}


# Fig1.I ------------------------------------------------------------------
# *Control_Upper -----------------------------------------------------------

table(merged.object$sample)
Sample <- "Control_Upper"
BarcodeData<-split(MetaData.merged,MetaData.merged$sample)
BarcodeData<-BarcodeData[[Sample]]
CU_DF <- GetFilteredSpatialData("./D1-5-3shang/outs", size="008um")
match_indices <- match(paste0(CU_DF$barcode, paste0("_", Sample)), rownames(BarcodeData))
table(is.na(match_indices))
CU_DF$Level5<-BarcodeData$L5[match(paste0(CU_DF$barcode,paste0("_",Sample)),rownames(BarcodeData))]
CU_DF$Level5<-factor(CU_DF$Level5,levels = sort(unique(CU_DF$Level5)))
path_to_scales_CU <- "./D1-5-3shang/outs/binned_outputs/square_008um/spatial/scalefactors_json.json"
scales_CU <- rjson::fromJSON(file = path_to_scales_CU)

plot_data_CU <- CU_DF %>% filter(tissue == 1, !is.na(Level5))
x_range_CU <- range(plot_data_CU$imagecol_scaled, na.rm = TRUE)
y_range_CU <- range(-plot_data_CU$imagerow_scaled, na.rm = TRUE)
bar_length_STANDARD <- (1000 / scales_CU$microns_per_pixel) * scales_CU$tissue_lowres_scalef
margin_x_CU <- (x_range_CU[2] - x_range_CU[1]) * 0.05
margin_y_CU <- (y_range_CU[2] - y_range_CU[1]) * 0.05
bar_x_start_CU <- x_range_CU[1] + margin_x_CU
bar_x_end_CU <- bar_x_start_CU + bar_length_STANDARD 
bar_y_pos_CU <- y_range_CU[1] + margin_y_CU

u1 <- CU_DF %>% 
  filter(tissue == 1, !is.na(Level5)) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  scale_color_manual(values = ColsL5) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_CU,
           xend = bar_x_end_CU, 
           y = bar_y_pos_CU,
           yend = bar_y_pos_CU,
           color = "black",
           linewidth = 1.5)
print(u1)

# *Control_Lower -----------------------------------------------------------

Sample <- "Control_Lower"
BarcodeData<-split(MetaData.merged,MetaData.merged$sample)
BarcodeData<-BarcodeData[[Sample]]
CL_DF <- GetFilteredSpatialData("./A1-5-1xia/outs", size="008um")
match_indices <- match(paste0(CL_DF$barcode, paste0("_", Sample)), rownames(BarcodeData))
table(is.na(match_indices))
CL_DF$Level5<-BarcodeData$L5[match(paste0(CL_DF$barcode,paste0("_",Sample)),rownames(BarcodeData))]
CL_DF$Level5<-factor(CL_DF$Level5,levels = sort(unique(CL_DF$Level5)))
path_to_scales_CL <- "./A1-5-1xia/outs/binned_outputs/square_008um/spatial/scalefactors_json.json"
scales_CL <- rjson::fromJSON(file = path_to_scales_CL)
bar_length_CL_native <- (1000 / scales_CL$microns_per_pixel) * scales_CL$tissue_lowres_scalef
rescale_factor <- bar_length_STANDARD / bar_length_CL_native
CL_DF_rescaled <- CL_DF
CL_DF_rescaled$imagecol_scaled <- CL_DF_rescaled$imagecol_scaled * rescale_factor
CL_DF_rescaled$imagerow_scaled <- CL_DF_rescaled$imagerow_scaled * rescale_factor
plot_data_CL_rescaled <- CL_DF_rescaled %>% filter(tissue == 1, !is.na(Level5))
x_range_CL_rescaled <- range(plot_data_CL_rescaled$imagecol_scaled, na.rm = TRUE)
y_range_CL_rescaled <- range(-plot_data_CL_rescaled$imagerow_scaled, na.rm = TRUE)

margin_x_CL_rescaled <- (x_range_CL_rescaled[2] - x_range_CL_rescaled[1]) * 0.05
margin_y_CL_rescaled <- (y_range_CL_rescaled[2] - y_range_CL_rescaled[1]) * 0.05
bar_x_start_CL_rescaled <- x_range_CL_rescaled[1] + margin_x_CL_rescaled
bar_x_end_CL_rescaled <- bar_x_start_CL_rescaled + bar_length_STANDARD 
bar_y_pos_CL_rescaled <- y_range_CL_rescaled[1] + margin_y_CL_rescaled

d1 <- CL_DF_rescaled %>% 
  filter(tissue == 1, !is.na(Level5)) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  scale_color_manual(values = ColsL5) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_CL_rescaled,
           xend = bar_x_end_CL_rescaled, 
           y = bar_y_pos_CL_rescaled,
           yend = bar_y_pos_CL_rescaled,
           color = "black",
           linewidth = 1.5)

print(d1)

# *HFpEF_Upper -----------------------------------------------------------

table(merged.object$sample)
Sample <- "HFpEF_Upper"
BarcodeData<-split(MetaData.merged,MetaData.merged$sample)
BarcodeData<-BarcodeData[[Sample]]
HU_DF <- GetFilteredSpatialData("./D1_3-shang/outs", size="008um")
match_indices <- match(paste0(HU_DF$barcode, paste0("_", Sample)), rownames(BarcodeData))
table(is.na(match_indices))
HU_DF$Level5<-BarcodeData$L5[match(paste0(HU_DF$barcode,paste0("_",Sample)),rownames(BarcodeData))]
HU_DF$Level5<-factor(HU_DF$Level5,levels = sort(unique(HU_DF$Level5)))
path_to_scales_HU <- "./D1_3-shang/outs/binned_outputs/square_008um/spatial/scalefactors_json.json"
scales_HU <- rjson::fromJSON(file = path_to_scales_HU)

bar_length_HU_native <- (1000 / scales_HU$microns_per_pixel) * scales_HU$tissue_lowres_scalef
rescale_factor_hu <- bar_length_STANDARD / bar_length_HU_native 
HU_DF_rescaled <- HU_DF
HU_DF_rescaled$imagecol_scaled <- HU_DF_rescaled$imagecol_scaled * rescale_factor_hu 
HU_DF_rescaled$imagerow_scaled <- HU_DF_rescaled$imagerow_scaled * rescale_factor_hu 

plot_data_hu <- HU_DF_rescaled %>% filter(tissue == 1, !is.na(Level5))
x_range_hu <- range(plot_data_hu$imagecol_scaled, na.rm = TRUE) 
y_range_hu <- range(plot_data_hu$imagerow_scaled, na.rm = TRUE)
margin_x_hu <- (x_range_hu[2] - x_range_hu[1]) * 0.05 
margin_y_hu <- (y_range_hu[2] - y_range_hu[1]) * 0.05

bar_x_start_hu <- x_range_hu[1] + margin_x_hu 
bar_x_end_hu <- bar_x_start_hu + bar_length_STANDARD 
bar_y_pos_hu <- y_range_hu[1] + margin_y_hu 

u2_hu <- HU_DF_rescaled %>% 
  filter(tissue == 1, !is.na(Level5)) %>% 
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  scale_color_manual(values = ColsL5) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_hu, 
           xend = bar_x_end_hu,
           y = bar_y_pos_hu,
           yend = bar_y_pos_hu,
           color = "black",
           linewidth = 1.5) 
print(u2_hu)

# *HFpEF_Lower -----------------------------------------------------------

Sample <- "HFpEF_Lower"
BarcodeData<-split(MetaData.merged,MetaData.merged$sample)
BarcodeData<-BarcodeData[[Sample]]
HL_DF <- GetFilteredSpatialData("./A1_3-xia/outs", size="008um")
match_indices <- match(paste0(HL_DF$barcode, paste0("_", Sample)), rownames(BarcodeData))
table(is.na(match_indices))  
HL_DF$Level5<-BarcodeData$L5[match(paste0(HL_DF$barcode,paste0("_",Sample)),rownames(BarcodeData))]
HL_DF$Level5<-factor(HL_DF$Level5,levels = sort(unique(HL_DF$Level5)))

path_to_scales_HL <- "./A1_3-xia/outs/binned_outputs/square_008um/spatial/scalefactors_json.json"
scales_HL <- rjson::fromJSON(file = path_to_scales_HL)
bar_length_HL_native <- (1000 / scales_HL$microns_per_pixel) * scales_HL$tissue_lowres_scalef
rescale_factor_hl <- bar_length_STANDARD / bar_length_HL_native 
HL_DF_rescaled <- HL_DF
HL_DF_rescaled$imagecol_scaled <- HL_DF_rescaled$imagecol_scaled * rescale_factor_hl 
HL_DF_rescaled$imagerow_scaled <- HL_DF_rescaled$imagerow_scaled * rescale_factor_hl 

plot_data_hl <- HL_DF_rescaled %>% filter(tissue == 1, !is.na(Level5)) 
x_range_hl <- range(plot_data_hl$imagecol_scaled, na.rm = TRUE) 
y_range_hl <- range(plot_data_hl$imagerow_scaled, na.rm = TRUE) 
margin_x_hl <- (x_range_hl[2] - x_range_hl[1]) * 0.05 
margin_y_hl <- (y_range_hl[2] - y_range_hl[1]) * 0.05 
bar_x_start_hl <- x_range_hl[1] + margin_x_hl 
bar_x_end_hl <- bar_x_start_hl + bar_length_STANDARD 
bar_y_pos_hl <- y_range_hl[1] + margin_y_hl 

d2_hl <- HL_DF_rescaled %>% 
  filter(tissue == 1, !is.na(Level5)) %>% 
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") + ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  scale_color_manual(values = ColsL5) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_hl, 
           xend = bar_x_end_hl,
           y = bar_y_pos_hl, 
           yend = bar_y_pos_hl,
           color = "black",
           linewidth = 1.5) 
print(d2_hl)

# Fig1.J ----------------------------------------------------------------------

data_to_plot <- merged.object@meta.data[, c("L5", "condition")]
l5_order <- levels(merged.object$L5) 
if (is.null(l5_order)) {
  print("merged.object$L5 is not a factor")
}
condition_order <- c("HFpEF", "Control") 
data_to_plot$L5 <- factor(data_to_plot$L5, levels = l5_order)
data_to_plot$condition <- factor(data_to_plot$condition, levels = condition_order)

p <- ggplot(data_to_plot, aes(x = condition, fill = L5)) + 
  geom_bar(position = "fill", width = 0.6) + 
  scale_fill_manual(values = ColsL5) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(
    title = "Cell Type Proportions by Condition",
    x = "Condition", 
    y = "Proportion", 
    fill = "Cell Type (L5)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1), 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)











