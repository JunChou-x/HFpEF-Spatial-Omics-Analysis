library(dplyr)
library(Seurat)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scattermore)
library(viridis)

# SFig5.A -----------------------------------------------------------------
# Remove noise and artifacts
slice_data_A <- CL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceA)
slice_data_B <- HL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceB)
slice_data_C <- HL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceC)
slice_data_D <- CL_DF_rescaled_bcsHD %>% filter(barcode %in% SliceD)

slice_data_A  <- slice_data_A  %>% filter(Level5 != "Ventricular Cardiomyocytes")
slice_data_B  <- slice_data_B  %>% filter(Level5 != "Ventricular Cardiomyocytes")
slice_data_C  <- slice_data_C  %>% filter(Level5 != "Ventricular Cardiomyocytes")
slice_data_D  <- slice_data_D  %>% filter(Level5 != "Ventricular Cardiomyocytes")

slice_data_A  <- slice_data_A  %>% filter(Level5 != "Atrial Cardiomyocytes")
slice_data_B  <- slice_data_B  %>% filter(Level5 != "Atrial Cardiomyocytes")
slice_data_C  <- slice_data_C  %>% filter(Level5 != "Atrial Cardiomyocytes")
slice_data_D  <- slice_data_D  %>% filter(Level5 != "Atrial Cardiomyocytes")

slice_data_A  <- slice_data_A  %>% filter(Level5 != "Atrial Fibroblasts")
slice_data_B  <- slice_data_B  %>% filter(Level5 != "Atrial Fibroblasts")
slice_data_C  <- slice_data_C  %>% filter(Level5 != "Atrial Fibroblasts")
slice_data_D  <- slice_data_D  %>% filter(Level5 != "Atrial Fibroblasts")

slice_data_A  <- slice_data_A  %>% filter(Level5 != "Atrial SMCs")
slice_data_B  <- slice_data_B  %>% filter(Level5 != "Atrial SMCs")
slice_data_C  <- slice_data_C  %>% filter(Level5 != "Atrial SMCs")
slice_data_D  <- slice_data_D  %>% filter(Level5 != "Atrial SMCs")

count_A <- slice_data_A %>% 
  group_by(Level5) %>% 
  summarise(Count = n()) %>%
  mutate(Slice = "Slice A")

count_B <- slice_data_B %>% 
  group_by(Level5) %>% 
  summarise(Count = n()) %>%
  mutate(Slice = "Slice B")

count_C <- slice_data_C %>% 
  group_by(Level5) %>% 
  summarise(Count = n()) %>%
  mutate(Slice = "Slice C")

count_D <- slice_data_D %>% 
  group_by(Level5) %>% 
  summarise(Count = n()) %>%
  mutate(Slice = "Slice D")

all_counts <- bind_rows(count_A, count_B, count_C, count_D)

total_A <- sum(count_A$Count)
total_B <- sum(count_B$Count)
total_C <- sum(count_C$Count)
total_D <- sum(count_D$Count)

p1 <- ggplot(all_counts, aes(x = Level5, y = Count, fill = Slice)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  ) +
  labs(
    title = "Cell Type Distribution Across Slices",
    subtitle = paste0("Total cells: A=", total_A, ", B=", total_B, ", C=", total_C, ", D=", total_D),
    x = "Cell Subtype",
    y = "Cell Count"
  )


# SFig5.B -----------------------------------------------------------------

vent_percents <- all_percents %>%
  filter(Level5 == "Ventricular Cardiomyocytes")
p_vent <- ggplot(vent_percents, aes(x = Slice, y = Percentage, fill = Level5)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +  
  scale_fill_manual(values = c("Ventricular Cardiomyocytes" = "#CE3D32FF")) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Ventricular Cardiomyocytes Percentage",
    x = "",
    y = "Percentage (%)"
  )

print(p_vent)


# SFig5.C -----------------------------------------------------------------

percent_A <- count_A %>% 
  mutate(Percentage = Count / sum(Count) * 100)

percent_B <- count_B %>% 
  mutate(Percentage = Count / sum(Count) * 100)

percent_C <- count_C %>% 
  mutate(Percentage = Count / sum(Count) * 100)

percent_D <- count_D %>% 
  mutate(Percentage = Count / sum(Count) * 100)

all_percents <- bind_rows(percent_A, percent_B, percent_C, percent_D)

p2 <- ggplot(all_percents, aes(x = Level5, y = Percentage, fill = Slice)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  ) +
  labs(
    title = "Cell Type Percentage Distribution",
    x = "Cell Subtype",
    y = "Percentage (%)"
  )

# SFig5.D -----------------------------------------------------------------

cellsA_base <- SliceA  
cellsB_base <- SliceB  
cellsC_base <- SliceC  
cellsD_base <- SliceD  

cellsA <- paste0(cellsA_base, "_Control_Lower")  
cellsB <- paste0(cellsB_base, "_HFpEF_Lower")   
cellsC <- paste0(cellsC_base, "_HFpEF_Lower")   
cellsD <- paste0(cellsD_base, "_Control_Lower")

print(sum(cellsA %in% colnames(merged.object)))
print(sum(cellsB %in% colnames(merged.object)))
print(sum(cellsC %in% colnames(merged.object)))
print(sum(cellsD %in% colnames(merged.object)))

cellsA_exist <- cellsA[cellsA %in% colnames(merged.object)]
cellsB_exist <- cellsB[cellsB %in% colnames(merged.object)]
cellsC_exist <- cellsC[cellsC %in% colnames(merged.object)]
cellsD_exist <- cellsD[cellsD %in% colnames(merged.object)]

cells_for_niche <- c(cellsA_exist, cellsB_exist, cellsC_exist, cellsD_exist)
subset_object <- subset(merged.object, cells = cells_for_niche)

subset_object$region <- NA 

subset_object@meta.data[cellsA_exist, "region"] <- "SliceA"
subset_object@meta.data[cellsB_exist, "region"] <- "SliceB"
subset_object@meta.data[cellsC_exist, "region"] <- "SliceC"
subset_object@meta.data[cellsD_exist, "region"] <- "SliceD"

subset_object$region <- as.factor(subset_object$region)
Idents(subset_object) <- "region"

region_markers <- FindAllMarkers(
  object = subset_object,
  assay = "RNA",          
  slot = "data",          
  only.pos = FALSE,      
  min.pct = 0.1,
  logfc.threshold = 0.25
)

blacklist <- c("Hba-a2", "Hbb-bs")
clean_top10 <- region_markers %>%
  filter(!gene %in% blacklist) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() 

genes_SliceA <- clean_top10 %>% filter(cluster == "SliceA") %>% pull(gene)
genes_SliceB <- clean_top10 %>% filter(cluster == "SliceB") %>% pull(gene)
genes_SliceC <- clean_top10 %>% filter(cluster == "SliceC") %>% pull(gene)
genes_SliceD <- clean_top10 %>% filter(cluster == "SliceD") %>% pull(gene)

genes_SliceB <- unique(c(genes_SliceB, "Comp", "Thbs1"))
genes_to_plot <- c(genes_SliceC, genes_SliceB, genes_SliceD, genes_SliceA)
genes_to_plot <- unique(genes_to_plot)
levels_order <- c("SliceC", "SliceB", "SliceD", "SliceA")
Idents(subset_object) <- factor(Idents(subset_object), levels = levels_order)

DotPlot(
  subset_object,
  features = genes_to_plot,
  dot.scale = 8
) +
  scale_color_gradient2(
    low = "#3C7DC4",    
    mid = "white",     
    high = "#E64B35",   
    midpoint = 0
  ) +
  RotatedAxis() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 10, face = "italic"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right"
  ) +
  labs(x = "", y = "", title = " ")


# SFig5.E -----------------------------------------------------------------
# *Upgraded pathway ----------------------------------------------------------------------

region_up <- region_markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster)

genes_A <- region_up %>% filter(cluster == "SliceA") %>% pull(gene)
genes_B <- region_up %>% filter(cluster == "SliceB") %>% pull(gene)
genes_C <- region_up %>% filter(cluster == "SliceC") %>% pull(gene)
genes_D <- region_up %>% filter(cluster == "SliceD") %>% pull(gene)

genes_A_entrez <- bitr(genes_A, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
genes_B_entrez <- bitr(genes_B, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
genes_C_entrez <- bitr(genes_C, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
genes_D_entrez <- bitr(genes_D, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

ego_A <- enrichGO(
  gene          = genes_A_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

ego_B <- enrichGO(
  gene          = genes_B_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

ego_C <- enrichGO(
  gene          = genes_C_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

ego_D <- enrichGO(
  gene          = genes_D_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

dotplot(ego_A, showCategory = 20) + ggtitle("SliceA GO BP")
dotplot(ego_B, showCategory = 20) + ggtitle("SliceB GO BP")
dotplot(ego_C, showCategory = 20) + ggtitle("SliceC GO BP")
dotplot(ego_D, showCategory = 20) + ggtitle("SliceD GO BP")

# SFig5.F -----------------------------------------------------------------
# *Downregulation Pathway --------------------------------------------------------------------

region_down <- region_markers %>%
  filter(avg_log2FC < 0) %>%  
  group_by(cluster)

genes_A_down <- region_down %>% filter(cluster == "SliceA") %>% pull(gene)
genes_B_down <- region_down %>% filter(cluster == "SliceB") %>% pull(gene)
genes_C_down <- region_down %>% filter(cluster == "SliceC") %>% pull(gene)
genes_D_down <- region_down %>% filter(cluster == "SliceD") %>% pull(gene)

to_entrez <- function(genes) {
  if(length(genes) == 0) return(NULL)
  bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)$ENTREZID
}

genes_A_entrez_down <- to_entrez(genes_A_down)
genes_B_entrez_down <- to_entrez(genes_B_down)
genes_C_entrez_down <- to_entrez(genes_C_down)
genes_D_entrez_down <- to_entrez(genes_D_down)

run_go_down <- function(gene_entrez) {
  if(is.null(gene_entrez)) return(NULL)
  enrichGO(
    gene          = gene_entrez,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

ego_A_down <- run_go_down(genes_A_entrez_down)
ego_B_down <- run_go_down(genes_B_entrez_down)
ego_C_down <- run_go_down(genes_C_entrez_down)
ego_D_down <- run_go_down(genes_D_entrez_down)

plot_and_print <- function(ego_obj, title_suffix) {
  if (is.null(ego_obj) || nrow(ego_obj) == 0) {
    message(paste("No significant downregulated pathways for", title_suffix))
    return(NULL)
  }
  p <- dotplot(ego_obj, showCategory = 15) + 
    ggtitle(paste(title_suffix, "Downregulated GO BP")) +
    theme(axis.text.y = element_text(size = 10)) 
  print(p)
}

plot_and_print(ego_A_down, "SliceA")
plot_and_print(ego_B_down, "SliceB")
plot_and_print(ego_C_down, "SliceC")
plot_and_print(ego_D_down, "SliceD")

# SFig5.G -----------------------------------------------------------------

expr_data <- FetchData(merged.object, vars = c("Comp", "Cd36"))
merged.object$LR_Score <- sqrt(expr_data$Comp * expr_data$Cd36)
expr_data <- FetchData(merged.object, vars = c("Thbs1", "Cd36"))
merged.object$LR_Score2 <- sqrt(expr_data$Thbs1 * expr_data$Cd36)

Sample <- "Control_Lower"
CL_DF_rescaled$LR_Score <- merged.object$LR_Score[match(paste0(CL_DF_rescaled$barcode, "_", Sample), rownames(merged.object@meta.data))]
CL_DF_rescaled$LR_Score2 <- merged.object$LR_Score2[match(paste0(CL_DF_rescaled$barcode, "_", Sample), rownames(merged.object@meta.data))]

# *comp -------------------------------------------------------------------

d1_score <- CL_DF_rescaled %>% 
  filter(tissue == 1, !is.na(LR_Score)) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = LR_Score)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right"
  ) +
  scale_color_distiller(palette = "Spectral", name = "COMP-CD36\nInteraction") +
  ggtitle("Control_Lower") +
  annotate("segment",
           x = bar_x_start_CL_rescaled,
           xend = bar_x_end_CL_rescaled, 
           y = bar_y_pos_CL_rescaled,
           yend = bar_y_pos_CL_rescaled,
           color = "white",
           linewidth = 1.5)

print(d1_score)

Sample <- "HFpEF_Lower"
HL_DF_rescaled$LR_Score <- merged.object$LR_Score[match(paste0(HL_DF_rescaled$barcode, "_", Sample), rownames(merged.object@meta.data))]
HL_DF_rescaled$LR_Score2 <- merged.object$LR_Score2[match(paste0(HL_DF_rescaled$barcode, "_", Sample), rownames(merged.object@meta.data))]

d2_score <- HL_DF_rescaled %>% 
  filter(tissue == 1, !is.na(LR_Score)) %>% 
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = LR_Score)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") + ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  scale_color_distiller(palette = "Spectral", name = "COMP-CD36\nInteraction") +
  ggtitle("HFpEF_Lower") +
  annotate("segment",
           x = bar_x_start_hl,
           xend = bar_x_end_hl,
           y = bar_y_pos_hl,
           yend = bar_y_pos_hl,
           color = "white",
           linewidth = 1.5)

print(d2_score)

# *thbs -------------------------------------------------------------------

d1_score2 <- CL_DF_rescaled %>% 
  filter(tissue == 1, !is.na(LR_Score2)) %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = LR_Score2)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") +
  ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right") +
  scale_color_distiller(palette = "Spectral", name = "THBS1-CD36\nInteraction") +
  ggtitle("Control_Lower") +
  annotate("segment",
           x = bar_x_start_CL_rescaled,
           xend = bar_x_end_CL_rescaled, 
           y = bar_y_pos_CL_rescaled,
           yend = bar_y_pos_CL_rescaled,
           color = "white",
           linewidth = 1.5)

print(d1_score2)

d2_score2 <- HL_DF_rescaled %>% 
  filter(tissue == 1, !is.na(LR_Score2)) %>% 
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = LR_Score2)) +
  
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  xlab("") + ylab("") +
  theme_set(theme_bw(base_size = 10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  scale_color_distiller(palette = "Spectral", name = "THBS1-CD36\nInteraction") +
  ggtitle("HFpEF_Lower") +
  annotate("segment",
           x = bar_x_start_hl,
           xend = bar_x_end_hl,
           y = bar_y_pos_hl,
           yend = bar_y_pos_hl,
           color = "white",
           linewidth = 1.5)

print(d2_score2)





