
library(viridis)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(scattermore)
library(scales)

# Fig5.A ------------------------------------------------------------------

CL_DF_rescaled <- na.omit(CL_DF_rescaled)
CL_DF_rescaled_bcsHD <- CL_DF_rescaled

# SliceA
SliceA <- GetRectangle(
  Spot = "s_008um_00488_00180-1", 
  WidthMicrons = 1300, 
  HeightMicrons = 1100, 
  BarcodeDF = CL_DF_rescaled_bcsHD
)
sliceA_data <-CL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceA)
DFRectA <- data.frame(
  Xmin = min(sliceA_data$imagecol_scaled),
  Xmax = max(sliceA_data$imagecol_scaled),
  Ymin = min(-sliceA_data$imagerow_scaled),
  Ymax = max(-sliceA_data$imagerow_scaled)
)

# SliceD
SliceD <- GetRectangle(
  Spot = "s_008um_00141_00725-1", 
  WidthMicrons = 1300, 
  HeightMicrons = 1100, 
  BarcodeDF = CL_DF_rescaled_bcsHD
)

sliceD_data <-CL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceD)
DFRectD <- data.frame(
  Xmin = min(sliceD_data$imagecol_scaled),
  Xmax = max(sliceD_data$imagecol_scaled),
  Ymin = min(-sliceD_data$imagerow_scaled),
  Ymax = max(-sliceD_data$imagerow_scaled)
)

HL_DF_rescaled <- na.omit(HL_DF_rescaled)
HL_DF_rescaled_bcsHD <- HL_DF_rescaled

# SliceB
SliceB <- GetRectangle(
  Spot = "s_008um_00507_00662-1", 
  WidthMicrons = 1300, 
  HeightMicrons = 1100, 
  BarcodeDF = HL_DF_rescaled_bcsHD
)
sliceB_data <- HL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceB)
DFRectB <- data.frame(
  Xmin = min(sliceB_data$imagecol_scaled),
  Xmax = max(sliceB_data$imagecol_scaled),
  Ymin = min(-sliceB_data$imagerow_scaled),
  Ymax = max(-sliceB_data$imagerow_scaled)
)


# SliceB

SliceC <- GetRectangle(
  Spot = "s_008um_00169_00128-1", 
  WidthMicrons = 1300, 
  HeightMicrons = 1100, 
  BarcodeDF = HL_DF_rescaled_bcsHD
)
sliceC_data <- HL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceC)
DFRectC <- data.frame(
  Xmin = min(sliceC_data$imagecol_scaled),
  Xmax = max(sliceC_data$imagecol_scaled),
  Ymin = min(-sliceC_data$imagerow_scaled),
  Ymax = max(-sliceC_data$imagerow_scaled)
)

# Scale definition
len_1mm_rescaled   <- bar_length_STANDARD 
len_150um_rescaled <- bar_length_STANDARD * (150 / 1000) 

# Slice A
plot_data_full_A <- CL_DF_rescaled_bcsHD %>% filter(tissue == 1)
plot_data_zoom_A <- sliceA_data

x_rng_A    <- range(plot_data_full_A$imagecol_scaled, na.rm=TRUE)
y_rng_A    <- range(-plot_data_full_A$imagerow_scaled, na.rm=TRUE) 
margin_x_A <- diff(x_rng_A) * 0.05
margin_y_A <- diff(y_rng_A) * 0.05
bar_x_start_A <- x_rng_A[1] + margin_x_A
bar_x_end_A   <- bar_x_start_A + len_1mm_rescaled
bar_y_pos_A   <- y_rng_A[1] + margin_y_A

x_rng_zA    <- range(plot_data_zoom_A$imagecol_scaled, na.rm=TRUE)
y_rng_zA    <- range(-plot_data_zoom_A$imagerow_scaled, na.rm=TRUE)
margin_x_zA <- diff(x_rng_zA) * 0.05
margin_y_zA <- diff(y_rng_zA) * 0.05
bar_x_start_zA <- x_rng_zA[1] + margin_x_zA
bar_x_end_zA   <- bar_x_start_zA + len_150um_rescaled
bar_y_pos_zA   <- y_rng_zA[1] + margin_y_zA

DFRectA <- data.frame(
  Xmin = min(plot_data_zoom_A$imagecol_scaled),
  Xmax = max(plot_data_zoom_A$imagecol_scaled),
  Ymin = min(-plot_data_zoom_A$imagerow_scaled), 
  Ymax = max(-plot_data_zoom_A$imagerow_scaled)
)

# Full
Plot1_A <- ggplot(plot_data_full_A, aes(x = imagecol_scaled, y = -imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = c(2000, 2000)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_A, xend = bar_x_end_A, y = bar_y_pos_A, yend = bar_y_pos_A, 
           color = "black", linewidth = 1.5) +
  scale_color_manual(values = subtype_colors) +
  geom_rect(data = DFRectA, aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.8)

# Zoom
Plot2_A <- ggplot(plot_data_zoom_A, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = Level5)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_zA, xend = bar_x_end_zA, y = bar_y_pos_zA, yend = bar_y_pos_zA, 
           color = "black", linewidth = 1.5) +
  scale_fill_manual(values = subtype_colors)

# Slice B
plot_data_full_B <- HL_DF_rescaled_bcsHD %>% filter(tissue == 1)
plot_data_zoom_B <- sliceB_data

x_rng_B    <- range(plot_data_full_B$imagecol_scaled, na.rm=TRUE)
y_rng_B    <- range(plot_data_full_B$imagerow_scaled, na.rm=TRUE) 
margin_x_B <- diff(x_rng_B) * 0.05
margin_y_B <- diff(y_rng_B) * 0.05

bar_x_start_B <- x_rng_B[1] + margin_x_B
bar_x_end_B   <- bar_x_start_B + len_1mm_rescaled
bar_y_pos_B   <- y_rng_B[1] + margin_y_B 

x_rng_zB    <- range(plot_data_zoom_B$imagecol_scaled, na.rm=TRUE)
y_rng_zB    <- range(plot_data_zoom_B$imagerow_scaled, na.rm=TRUE) 
margin_x_zB <- diff(x_rng_zB) * 0.05
margin_y_zB <- diff(y_rng_zB) * 0.05

bar_x_start_zB <- x_rng_zB[1] + margin_x_zB
bar_x_end_zB   <- bar_x_start_zB + len_150um_rescaled
bar_y_pos_zB   <- y_rng_zB[1] + margin_y_zB

DFRectB <- data.frame(
  Xmin = min(plot_data_zoom_B$imagecol_scaled),
  Xmax = max(plot_data_zoom_B$imagecol_scaled),
  Ymin = min(plot_data_zoom_B$imagerow_scaled),
  Ymax = max(plot_data_zoom_B$imagerow_scaled)  
)

# Full
Plot1_B <- ggplot(plot_data_full_B, aes(x = imagecol_scaled, y = imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = c(2000, 2000)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_B, xend = bar_x_end_B, y = bar_y_pos_B, yend = bar_y_pos_B, 
           color = "black", linewidth = 1.5) +
  scale_color_manual(values = subtype_colors) +
  geom_rect(data = DFRectB, aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.8)

# Zoom
Plot2_B <- ggplot(plot_data_zoom_B, aes(x = imagecol_scaled, y = imagerow_scaled, fill = Level5)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_zB, xend = bar_x_end_zB, y = bar_y_pos_zB, yend = bar_y_pos_zB, 
           color = "black", linewidth = 1.5) +
  scale_fill_manual(values = subtype_colors)


# Slice C
plot_data_zoom_C <- sliceC_data

x_rng_zC    <- range(plot_data_zoom_C$imagecol_scaled, na.rm=TRUE)
y_rng_zC    <- range(plot_data_zoom_C$imagerow_scaled, na.rm=TRUE)
margin_x_zC <- diff(x_rng_zC) * 0.05
margin_y_zC <- diff(y_rng_zC) * 0.05
bar_x_start_zC <- x_rng_zC[1] + margin_x_zC
bar_x_end_zC   <- bar_x_start_zC + len_150um_rescaled
bar_y_pos_zC   <- y_rng_zC[1] + margin_y_zC

DFRectC <- data.frame(
  Xmin = min(plot_data_zoom_C$imagecol_scaled),
  Xmax = max(plot_data_zoom_C$imagecol_scaled),
  Ymin = min(plot_data_zoom_C$imagerow_scaled),
  Ymax = max(plot_data_zoom_C$imagerow_scaled)
)

# Full
Plot1_C <- ggplot(plot_data_full_B, aes(x = imagecol_scaled, y = imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = c(2000, 2000)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_B, xend = bar_x_end_B, y = bar_y_pos_B, yend = bar_y_pos_B, 
           color = "black", linewidth = 1.5) +
  scale_color_manual(values = subtype_colors) +
  geom_rect(data = DFRectC, aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.8)

# Zoom 
Plot2_C <- ggplot(plot_data_zoom_C, aes(x = imagecol_scaled, y = imagerow_scaled, fill = Level5)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_zC, xend = bar_x_end_zC, y = bar_y_pos_zC, yend = bar_y_pos_zC, 
           color = "black", linewidth = 1.5) +
  scale_fill_manual(values = subtype_colors)

# Slice D
plot_data_full_D <- CL_DF_rescaled_bcsHD %>% filter(tissue == 1)
plot_data_zoom_D <- sliceD_data

x_rng_D    <- range(plot_data_full_D$imagecol_scaled, na.rm=TRUE)
y_rng_D    <- range(-plot_data_full_D$imagerow_scaled, na.rm=TRUE) 
margin_x_D <- diff(x_rng_D) * 0.05
margin_y_D <- diff(y_rng_D) * 0.05
bar_x_start_D <- x_rng_D[1] + margin_x_D
bar_x_end_D   <- bar_x_start_D + len_1mm_rescaled
bar_y_pos_D   <- y_rng_D[1] + margin_y_D

x_rng_zD    <- range(plot_data_zoom_D$imagecol_scaled, na.rm=TRUE)
y_rng_zD    <- range(-plot_data_zoom_D$imagerow_scaled, na.rm=TRUE) 
margin_x_zD <- diff(x_rng_zD) * 0.05
margin_y_zD <- diff(y_rng_zD) * 0.05

bar_x_start_zD <- x_rng_zD[1] + margin_x_zD
bar_x_end_zD   <- bar_x_start_zD + len_150um_rescaled
bar_y_pos_zD   <- y_rng_zD[1] + margin_y_zD

# Full
Plot1_D <- ggplot(plot_data_full_D, aes(x = imagecol_scaled, y = -imagerow_scaled, color = Level5)) +
  geom_scattermore(pointsize = 3, pixels = c(2000, 2000)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_D, xend = bar_x_end_D, y = bar_y_pos_D, yend = bar_y_pos_D, 
           color = "black", linewidth = 1.5) +
  scale_color_manual(values = subtype_colors) +
  geom_rect(data = DFRectD, aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.8)

# Zoom 
Plot2_D <- ggplot(plot_data_zoom_D, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = Level5)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_zD, xend = bar_x_end_zD, y = bar_y_pos_zD, yend = bar_y_pos_zD, 
           color = "black", linewidth = 1.5) +
  scale_fill_manual(values = subtype_colors)

print(Plot1_A)
print(Plot2_A)

print(Plot1_B)
print(Plot2_B)

print(Plot1_C)
print(Plot2_C)

print(Plot1_D)
print(Plot2_D)


# Fig5.B ------------------------------------------------------------------

merged.object$celltype_condition <- paste0(
  merged.object$L5, "_", merged.object$condition
)
table(merged.object$celltype_condition)

ligand_target_matrix <- readRDS("./ligand_target_matrix_nsga2r_final_mouse.rds")
lr_network <- readRDS("./lr_network_mouse_21122021.rds")

grouping_column <- "celltype_condition" 
sender_celltype <- "Perivascular Fibroblasts_HFpEF" 
receiver_celltype <- "Stressed Cardiomyocytes_HFpEF"
receiver_reference <- "Stressed Cardiomyocytes_Control"

Idents(merged.object) <- merged.object$celltype_condition

# Sender-Focused Approach
expressed_genes_receiver <- get_expressed_genes(receiver_celltype, merged.object, pct = 0.10)
expressed_receptors <- intersect(unique(lr_network$to), expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
expressed_genes_sender <- get_expressed_genes(sender_celltype, merged.object, pct = 0.10)
expressed_potential_ligands <- intersect(potential_ligands, expressed_genes_sender)

DE_table <- FindMarkers(
  object = merged.object,
  ident.1 = receiver_celltype,    
  ident.2 = receiver_reference,   
  min.pct = 0.10,
  logfc.threshold = 0.25,
  assay = "RNA"
)

geneset_oi <- DE_table %>% filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% rownames()
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = expressed_potential_ligands
)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)

target_ligands <- c("Comp", "Thbs1")

potential_receptors <- lr_network %>% 
  filter(from %in% target_ligands) %>%         
  filter(to %in% background_expressed_genes) %>% 
  distinct(to) %>%                            
  pull(to)

lr_network_df_combined <- weighted_networks$lr_sig %>%
  filter(from %in% target_ligands) %>%       
  filter(to %in% potential_receptors) %>%    
  arrange(from, desc(weight))              

print(head(lr_network_df_combined, 20))

receptor_order <- lr_network_df_combined %>%
  group_by(to) %>%
  summarise(total_weight = sum(weight, na.rm = TRUE)) %>%
  arrange(desc(total_weight)) %>%
  pull(to)

plot_df <- lr_network_df_combined %>%
  mutate(to = factor(to, levels = receptor_order))

p_heatmap_combined <- ggplot(plot_df, aes(x = to, y = from, fill = weight)) +
  geom_tile(color = "white", size = 0.2) + 
  scale_fill_gradientn(
    colors = c("#FFF5F5", "#FFD700", "#FF4500", "#8B0000"), 
    name = "Prior Interaction\nWeight"
  ) + 
  labs(
    title = "Ligand-Receptor Interaction Strength (Comp & Thbs1)",
    subtitle = "Receptors ordered by combined interaction weight",
    x = "Receptor expressed in Receiver Cells",
    y = "Ligand from Sender Cells"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"), 
    axis.text.y = element_text(size = 11, face = "bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  )

print(p_heatmap_combined)






# Fig.C,D -------------------------------------------------------------

stressed_cm_AB <- subset(subset_object, 
                         subset = L5 == "Stressed Cardiomyocytes" & 
                           region %in% c("SliceA", "SliceB"))

degs_up_B <- FindMarkers(stressed_cm_AB, 
                         ident.1 = "SliceB", 
                         ident.2 = "SliceA", 
                         group.by = "region", 
                         only.pos = TRUE,         
                         min.pct = 0.1,           
                         logfc.threshold = 0.25)  

significant_up_genes <- degs_up_B[degs_up_B$p_val_adj < 0.05, ]
genes_to_test <- rownames(significant_up_genes)
if(length(genes_to_test) > 0) {
  gene_map <- bitr(genes_to_test, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)

  ego_result <- enrichGO(gene          = gene_map$ENTREZID,
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",              
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,    
                         qvalueCutoff  = 0.2,     
                         readable      = TRUE)    
  
  if (!is.null(ego_result) && nrow(ego_result@result[ego_result@result$p.adjust < 0.05, ]) > 0) {
    p_dot <- dotplot(ego_result, showCategory = 15) + 
      ggtitle("Up-regulated Genes in Slice B Stress CM") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.minor = element_blank()
      )
    print(p_dot)
    top_pathways <- head(ego_result@result, 50)
    for(i in 1:nrow(top_pathways)) {
      cat(paste0("Pathway [", i, "]: ", top_pathways$Description[i], "\n"))
      cat(paste0("Gene: ", top_pathways$geneID[i], "\n"))
      cat("--------------------------------------------------\n")
    }
  } else {
    print("no significantly enrich")
  }
} else {
  print("No significantly different genes were found.")
}

ego_subset <- ego_result
df_all <- ego_result@result

top15_data <- head(df_all, 15)
lipid_data_all <- df_all[grep("lipid", df_all$Description, ignore.case = TRUE), ]
lipid_data_top <- head(lipid_data_all, 3)

combined_data <- rbind(top15_data, lipid_data_top)
combined_data <- unique(combined_data)

ego_subset@result <- combined_data
total_pathways <- nrow(combined_data)

p_dot <- dotplot(ego_subset, showCategory = total_pathways) + 
  ggtitle("Up-regulated Genes in Slice B Stress CM") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

print(p_dot)



# Fig.E -------------------------------------------------------------

Sample <- "Control_Lower"
# Slice A (Control)
plot_data_zoom_A$LR_Score <- merged.object$LR_Score[match(paste0(plot_data_zoom_A$barcode, "_", Sample), rownames(merged.object@meta.data))]
plot_data_zoom_A$LR_Score2 <- merged.object$LR_Score2[match(paste0(plot_data_zoom_A$barcode, "_", Sample), rownames(merged.object@meta.data))]
# Slice D (Control)
plot_data_zoom_D$LR_Score <- merged.object$LR_Score[match(paste0(plot_data_zoom_D$barcode, "_", Sample), rownames(merged.object@meta.data))]
plot_data_zoom_D$LR_Score2 <- merged.object$LR_Score2[match(paste0(plot_data_zoom_D$barcode, "_", Sample), rownames(merged.object@meta.data))]
Sample <- "HFpEF_Lower"
# Slice B (HFpEF Niche)
plot_data_zoom_B$LR_Score <- merged.object$LR_Score[match(paste0(plot_data_zoom_B$barcode, "_", Sample), rownames(merged.object@meta.data))]
plot_data_zoom_B$LR_Score2 <- merged.object$LR_Score2[match(paste0(plot_data_zoom_B$barcode, "_", Sample), rownames(merged.object@meta.data))]
# Slice C (HFpEF Remote)
plot_data_zoom_C$LR_Score <- merged.object$LR_Score[match(paste0(plot_data_zoom_C$barcode, "_", Sample), rownames(merged.object@meta.data))]
plot_data_zoom_C$LR_Score2 <- merged.object$LR_Score2[match(paste0(plot_data_zoom_C$barcode, "_", Sample), rownames(merged.object@meta.data))]

# *comp -------------------------------------------------------------------

max_score <- max(c(plot_data_zoom_A$LR_Score, plot_data_zoom_B$LR_Score, plot_data_zoom_C$LR_Score), na.rm = TRUE)

# --- Plot A (Control) ---
Plot2_A_score <- ggplot(plot_data_zoom_A, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = LR_Score)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "none") + 
  annotate("segment", x = bar_x_start_zA, xend = bar_x_end_zA, y = bar_y_pos_zA, yend = bar_y_pos_zA, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score), name = "Score") +
  ggtitle("Slice Control")

# --- Plot B (HFpEF Niche)
Plot2_B_score <- ggplot(plot_data_zoom_B, aes(x = imagecol_scaled, y = imagerow_scaled, fill = LR_Score)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "none") + 
  annotate("segment", x = bar_x_start_zB, xend = bar_x_end_zB, y = bar_y_pos_zB, yend = bar_y_pos_zB, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score), name = "Interaction\nPotential") +
  ggtitle("Slice HFpEF")

# --- Plot C (HFpEF Remote) ---
Plot2_C_score <- ggplot(plot_data_zoom_C, aes(x = imagecol_scaled, y = imagerow_scaled, fill = LR_Score)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_zC, xend = bar_x_end_zC, y = bar_y_pos_zC, yend = bar_y_pos_zC, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score), name = "Comp-Cd36 L-R_Score") +
  ggtitle("Slice Distal-HFpEF")

# --- Plot D (HFpEF Remote) ---
Plot2_D_score <- ggplot(plot_data_zoom_D, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = LR_Score)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "right") + 
  annotate("segment", x = bar_x_start_zD, xend = bar_x_end_zD, y = bar_y_pos_zD, yend = bar_y_pos_zD, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score), name = "Comp-Cd36 L-R_Score") +
  ggtitle("Slice Distal-Control")

print(Plot2_A_score)
print(Plot2_B_score)
print(Plot2_C_score)
print(Plot2_D_score)


# *thbs -------------------------------------------------------------------

max_score2 <- max(c(plot_data_zoom_A$LR_Score2, plot_data_zoom_B$LR_Score2, plot_data_zoom_C$LR_Score2), na.rm = TRUE)

# --- Plot A (Control) ---
Plot2_A_score2 <- ggplot(plot_data_zoom_A, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = LR_Score2)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "none") + 
  annotate("segment", x = bar_x_start_zA, xend = bar_x_end_zA, y = bar_y_pos_zA, yend = bar_y_pos_zA, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score2), name = "Score") +
  ggtitle("Slice Control")

# --- Plot B (HFpEF Niche) ---
Plot2_B_score2 <- ggplot(plot_data_zoom_B, aes(x = imagecol_scaled, y = imagerow_scaled, fill = LR_Score2)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "none") + 
  annotate("segment", x = bar_x_start_zB, xend = bar_x_end_zB, y = bar_y_pos_zB, yend = bar_y_pos_zB, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score2), name = "Interaction\nPotential") +
  ggtitle("Slice HFpEF")

# --- Plot C (HFpEF Remote) ---
Plot2_C_score2 <- ggplot(plot_data_zoom_C, aes(x = imagecol_scaled, y = imagerow_scaled, fill = LR_Score2)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "none") +
  annotate("segment", x = bar_x_start_zC, xend = bar_x_end_zC, y = bar_y_pos_zC, yend = bar_y_pos_zC, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score2), name = "Thbs1-Cd36 L-R_Score") +
  ggtitle("Slice Distal-HFpEF")

# --- Plot D (Control Remote) ---
Plot2_D_score2 <- ggplot(plot_data_zoom_D, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = LR_Score)) +
  geom_point(shape = 22, size = 4.5, color = "transparent", stroke = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), 
    axis.title = element_blank(), legend.position = "right") + 
  annotate("segment", x = bar_x_start_zD, xend = bar_x_end_zD, y = bar_y_pos_zD, yend = bar_y_pos_zD, 
           color = "black", linewidth = 1.5) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, max_score), name = "Thbs1-Cd36 L-R_Score") +
  ggtitle("Slice Distal-Control")


print(Plot2_A_score2)
print(Plot2_B_score2)
print(Plot2_C_score2)
print(Plot2_D_score2)


# Fig.F -------------------------------------------------------------

gene1 <- "Pdk4"   
gene2 <- "Nppb" 

limit_g1 <- 4  
limit_g2 <- 4 

expr_dual <- FetchData(merged.object, vars = c(gene1, gene2))
expr_dual$barcode <- rownames(expr_dual)

get_blend_color <- function(df_barcode, sample_suffix, expr_data, g1, g2, lim1, lim2) {
  full_barcodes <- paste0(df_barcode, "_", sample_suffix)
  matched_idx <- match(full_barcodes, expr_data$barcode)
  val1 <- expr_data[[g1]][matched_idx]
  val2 <- expr_data[[g2]][matched_idx]

  val1[is.na(val1)] <- 0
  val2[is.na(val2)] <- 0
  
  norm1 <- scales::squish(val1, c(0, lim1)) / lim1
  norm2 <- scales::squish(val2, c(0, lim2)) / lim2

  rgb(red = norm1, green = norm2, blue = 0)
}

Sample <- "Control_Upper"
CU_DF$blend_color <- get_blend_color(CU_DF$barcode, Sample, expr_dual, 
                                     gene1, gene2, limit_g1, limit_g2)

CU_DF %>%
  filter(tissue == 1) %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = blend_color)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  scale_color_identity() +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  labs(title = paste0(gene1, " (Red) + ", gene2, " (Green)"), 
       subtitle = "Yellow = Co-expression") +
  annotate("segment",
           x = bar_x_start_CU, xend = bar_x_end_CU, 
           y = bar_y_pos_CU, yend = bar_y_pos_CU,
           color = "black", 
           linewidth = 1.5)


Sample <- "Control_Lower"
CL_DF_rescaled$blend_color <- get_blend_color(CL_DF_rescaled$barcode, Sample, expr_dual, 
                                              gene1, gene2, limit_g1, limit_g2)
CL_DF_rescaled %>%
  filter(tissue == 1) %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = blend_color)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  scale_color_identity() + 
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  labs(title = paste0(gene1, " (Red) + ", gene2, " (Green)")) + 
  annotate("segment",
           x = bar_x_start_CL_rescaled, xend = bar_x_end_CL_rescaled, 
           y = bar_y_pos_CL_rescaled, yend = bar_y_pos_CL_rescaled,
           color = "white", 
           linewidth = 1.5)


Sample <- "HFpEF_Upper"

HU_DF_rescaled$blend_color <- get_blend_color(HU_DF_rescaled$barcode, Sample, expr_dual, 
                                              gene1, gene2, limit_g1, limit_g2)

HU_DF_rescaled %>%
  filter(tissue == 1) %>%
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = blend_color)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  scale_color_identity() +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  labs(title = paste0(gene1, " (Red) + ", gene2, " (Green)"),
       subtitle = "HFpEF Upper") +
  annotate("segment",
           x = bar_x_start_hu, 
           xend = bar_x_end_hu, 
           y = bar_y_pos_hu, 
           yend = bar_y_pos_hu,
           color = "white",
           linewidth = 1.5)

Sample <- "HFpEF_Lower"

HL_DF_rescaled$blend_color <- get_blend_color(HL_DF_rescaled$barcode, Sample, expr_dual, 
                                              gene1, gene2, limit_g1, limit_g2)

HL_DF_rescaled %>%
  filter(tissue == 1) %>%
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = blend_color)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  scale_color_identity() +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  labs(title = paste0(gene1, " (Red) + ", gene2, " (Green)"),
       subtitle = "HFpEF Lower") +
  annotate("segment",
           x = bar_x_start_hl, 
           xend = bar_x_end_hl, 
           y = bar_y_pos_hl, 
           yend = bar_y_pos_hl,
           color = "black",
           linewidth = 1.5)
