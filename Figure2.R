
library(dplyr)
library(ggplot2)
library(scattermore)
library(Seurat)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(Matrix)
library(reshape2)
library(slingshot)
library(SingleCellExperiment)
library(tradeSeq)
library(RColorBrewer)

# Fig2.A ------------------------------------------------------------------
target_L5_types <- c("Ventricular Cardiomyocytes", "Atrial Cardiomyocytes", 
                     "Stressed Cardiomyocytes", "Conduction System Cardiomyocytes")
plot_colors <- ColsL5[target_L5_types]
plot_colors["Other"] <- "grey90"
prepare_data <- function(df) {
  df %>%
    filter(tissue == 1, !is.na(Level5)) %>%
    mutate(Plotting_Group = if_else(Level5 %in% target_L5_types, as.character(Level5), "Other")) %>%
    mutate(Plotting_Group = factor(Plotting_Group, levels = c(sort(target_L5_types), "Other"))) %>%
    arrange(desc(Plotting_Group == "Other")) 
}
plot_data_CU <- prepare_data(CU_DF)
plot_data_CL <- prepare_data(CL_DF_rescaled)
plot_data_HU <- prepare_data(HU_DF_rescaled)
plot_data_HL <- prepare_data(HL_DF_rescaled)

u1_cardio <- plot_data_CU %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = Plotting_Group)) + 
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
  scale_color_manual(values = plot_colors) + 
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_CU,
           xend = bar_x_end_CU, 
           y = bar_y_pos_CU,
           yend = bar_y_pos_CU,
           color = "black",
           linewidth = 1.5)

print(u1_cardio)

d1_cardio <- plot_data_CL %>% 
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = Plotting_Group)) +
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
  scale_color_manual(values = plot_colors) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_CL_rescaled,
           xend = bar_x_end_CL_rescaled, 
           y = bar_y_pos_CL_rescaled,
           yend = bar_y_pos_CL_rescaled,
           color = "black",
           linewidth = 1.5)

print(d1_cardio)

u2_hu_cardio <- plot_data_HU %>% 
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = Plotting_Group)) + 
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
  scale_color_manual(values = plot_colors) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_hu, 
           xend = bar_x_end_hu, 
           y = bar_y_pos_hu, 
           yend = bar_y_pos_hu,
           color = "black",
           linewidth = 1.5) 

print(u2_hu_cardio)

d2_hl_cardio <- plot_data_HL %>% 
  ggplot(aes(x = imagecol_scaled, y = imagerow_scaled, color = Plotting_Group)) + 
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
  scale_color_manual(values = plot_colors) +
  labs(color = "Subtype") +
  ggtitle("") +
  annotate("segment",
           x = bar_x_start_hl, 
           xend = bar_x_end_hl, 
           y = bar_y_pos_hl, 
           yend = bar_y_pos_hl,
           color = "black",
           linewidth = 1.5) 

print(d2_hl_cardio)

# Fig2.B ------------------------------------------------------------------

target_cms <- c(
  "Ventricular Cardiomyocytes", 
  "Atrial Cardiomyocytes", 
  "Stressed Cardiomyocytes", 
  "Conduction System Cardiomyocytes"
)
Idents(merged.object) <- "L5"
cm_subset <- subset(merged.object, idents = target_cms)
cm_markers <- FindAllMarkers(
  cm_subset,
  only.pos = TRUE,       
  min.pct = 0.25,        
  logfc.threshold = 0.25,
  max.cells.per.ident = 500 
)
cm_markers_sig <- cm_markers %>% 
  filter(p_val_adj < 0.05)
top5_genes <- cm_markers_sig %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top5_genes , n = 40)

genes_to_plot <- unique(top5_genes$gene)
levels_order <- c("Conduction System Cardiomyocytes", 
                  "Stressed Cardiomyocytes", 
                  "Atrial Cardiomyocytes", 
                  "Ventricular Cardiomyocytes")
Idents(cm_subset) <- factor(Idents(cm_subset), levels = levels_order)

DotPlot(
  cm_subset,
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
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  ) +
  labs(x = "", y = "")


# Fig2.C ----------------------------------------------------------------

merged.object$condition <- ifelse(grepl("Control", merged.object$sample), "Control", "HFpEF")
merged.object$condition <- factor(merged.object$condition, levels = c("Control", "HFpEF"))
Idents(merged.object) <- "L5"
target_cms <- c(
  "Ventricular Cardiomyocytes", 
  "Atrial Cardiomyocytes", 
  "Stressed Cardiomyocytes", 
  "Conduction System Cardiomyocytes" 
)
all_markers_list <- list()
for (cell_type in target_cms) {
  message(paste0("calculate: ", cell_type, " (HFpEF vs Control)..."))
  markers <- FindMarkers(
    merged.object,
    ident.1 = "HFpEF",
    ident.2 = "Control",
    group.by = "condition",  
    subset.ident = cell_type, 
    max.cells.per.ident = 1000, 
    min.pct = 0.1,
    logfc.threshold = 0.1
  )
  if (nrow(markers) > 0) {
    markers$gene <- rownames(markers)     
    markers$cell_type <- cell_type        
    all_markers_list[[cell_type]] <- markers 
  }
}
final_DE_results <- do.call(rbind, all_markers_list)
sig_DE_results <- final_DE_results %>% 
  filter(p_val_adj < 0.05)

plot_data <- final_DE_results %>%
  mutate(
    p_val_calc = ifelse(p_val == 0, min(p_val[p_val > 0]) * 0.1, p_val),
    log10P = -log10(p_val_calc),
    change = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.1 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.1 ~ "Down",
      TRUE ~ "NS"
    )
  )

genes_ventricular <- c("Pdk4", "Fabp3", "Cpt1b", "Txnip", "Hadhb", 
                       "Myh6", "Myl2", "Myl3", "Lpl", "mt-Nd4l", "mt-Atp8")
genes_stressed <- c("Fabp3", "Pdk4", "Acta1", "Actc1", "Hadhb",
                    "Myl2", "Myl3", "Lpl", "Tnnc1")
genes_atrial <- c("Pdk4", "Fabp3", "Nppb", "Ankrd1", "Cpt1b",
                  "Myl7", "Nppa", "Tnnt2")
genes_conduction <- c("Nr1d1", "Fabp3", "Txnip", "Acta2",
                      "Gja1", "Nppa", "Cd74")

target_map <- bind_rows(
  data.frame(cell_type = "Ventricular Cardiomyocytes",       gene = genes_ventricular),
  data.frame(cell_type = "Stressed Cardiomyocytes",          gene = genes_stressed),
  data.frame(cell_type = "Atrial Cardiomyocytes",            gene = genes_atrial),
  data.frame(cell_type = "Conduction System Cardiomyocytes", gene = genes_conduction)
)

label_data <- plot_data %>%
  inner_join(target_map, by = c("cell_type", "gene"))

p_volcano <- ggplot(plot_data, aes(x = avg_log2FC, y = log10P)) +
  geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NS" = "grey80")) +
  geom_text_repel(
    data = label_data,         
    aes(label = gene),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf,
    fontface = "italic",
    min.segment.length = 0 
  ) +
  facet_wrap(~cell_type, scales = "free", nrow = 2) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10) 
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  labs(
    x = "log2 Fold Change (HFpEF vs Control)",
    y = "-log10(adj P-value)",
    title = " ",
    color = "Expression Status"
  )

print(p_volcano)

# Fig2.D ------------------------------------------------------------------

target_cms <- c(
  "Ventricular Cardiomyocytes", 
  "Atrial Cardiomyocytes", 
  "Stressed Cardiomyocytes", 
  "Conduction System Cardiomyocytes"
)
gsea_results_list <- list()
ora_results_list <- list()
for (current_type in target_cms) {
  cat(paste0("\n\n======================================================\n"))
  cat(paste0("Analyzing: ", current_type, "\n"))
  cat(paste0("======================================================\n"))
  de_data <- final_DE_results %>% 
    dplyr::filter(cell_type == current_type) %>%  
    dplyr::filter(!gene %in% c("Hba-a2", "Hbb-bs"))
  ids <- tryCatch({
    bitr(de_data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  }, error = function(e) { return(NULL) })
  if (is.null(ids)) { 
    cat("ID conversion failed, skip.\n") 
    next 
  }
  de_data_merged <- merge(de_data, ids, by.x = "gene", by.y = "SYMBOL")
  if(nrow(de_data_merged) == 0) { 
    cat("No valid genes found, skip this step.\n") 
    next 
  }
  de_data_unique <- de_data_merged %>%
    dplyr::arrange(dplyr::desc(abs(avg_log2FC))) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  geneList_gsea <- de_data_unique$avg_log2FC
  names(geneList_gsea) <- de_data_unique$ENTREZID
  geneList_gsea <- sort(geneList_gsea, decreasing = TRUE)
  gsea_res <- tryCatch({
    gseGO(
      geneList     = geneList_gsea,
      OrgDb        = org.Mm.eg.db,
      ont          = "BP",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 1,   
      verbose      = FALSE,
      seed         = 123
    )
  }, error = function(e) { return(NULL) })
  if (!is.null(gsea_res)) {
    gsea_results_list[[current_type]] <- gsea_res
    cat("\n  [GSEA] Top Pathways:\n")
    print(
      head(
        as.data.frame(gsea_res) %>% 
          dplyr::arrange(p.adjust) %>% 
          dplyr::select(ID, Description, NES, p.adjust),
        10
      )
    )
  } else {
    cat("  (No significant GSEA results)\n")
  }
}

plot_list <- list()

for (ctype in names(gsea_results_list)) {
  res <- gsea_results_list[[ctype]]
  if (is.null(res)) next
  df <- as.data.frame(res) %>%
    mutate(cell_type = ctype) %>%
    filter(p.adjust < 0.05) 
  if (nrow(df) == 0) next

  top_up <- df %>% 
    filter(NES > 0) %>% 
    arrange(p.adjust) %>% 
    head(3)
  
  top_down <- df %>% 
    filter(NES < 0) %>% 
    arrange(p.adjust) %>% 
    head(3)
  
  plot_list[[ctype]] <- rbind(top_up, top_down)
}

plot_data_gsea <- do.call(rbind, plot_list)
plot_data_gsea$Description <- str_trunc(plot_data_gsea$Description, width = 50)

p_gsea <- ggplot(plot_data_gsea, aes(x = NES, y = reorder(Description, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), 
               color = "grey70", linewidth = 0.5) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_color_gradientn(colours = c("#3C5488", "#8491B4", "grey90", "#F39B7F", "#E64B35"),
                        values = c(0, 0.4, 0.5, 0.6, 1), 
                        guide = guide_colorbar(reverse = TRUE)) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted"), 
    axis.text.y = element_text(size = 10, color = "black"), 
    strip.background = element_rect(fill = "white", color = "black"), 
    strip.text = element_text(face = "bold", size = 12) 
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = "Adj. P-value",
    size = "Gene Count",
    title = "Functional Divergence of Cardiomyocyte Subtypes in HFpEF"
  )

print(p_gsea)


# Fig2.E ------------------------------------------------------------------

lipid_genes <- list(c("Pdk4", "Fabp3", "Cpt1b", "Acadvl", "Slc27a1"))
oxphos_genes <- list(c("Ndufa1","Ndufa4","Ndufa5","Ndufa10", "Ndufb5", "Cox5a", "Atp5a1","Atp5b", "Atp5j"))
stress_genes <- list(c("Ankrd1", "Xirp2", "Fhl2", "Myh7", "Nppb"))

features_plot <- c(
  "Pdk4", "Slc27a1", "Fabp3",           # lipid intake
  "Cpt1b", "Acadvl",            # Fatty acid oxidation
  "Ndufa1","Ndufa4","Ndufa5","Ndufa10",
  "Ndufb5", "Cox5a", "Atp5a1","Atp5b", "Atp5j",
  "Ankrd1", "Xirp2", "Fhl2", "Myh6","Nppb"
)

features_plot <- intersect(features_plot, rownames(merged.object))
target_subtypes <- c(
  "Ventricular Cardiomyocytes", 
  "Atrial Cardiomyocytes", 
  "Stressed Cardiomyocytes", 
  "Conduction System Cardiomyocytes"
)

de_results_list <- list()
Idents(merged.object) <- "L5"

for (ctype in target_subtypes) {
  message(paste("Calculating:", ctype))
  res <- FindMarkers(
    merged.object,
    ident.1 = "HFpEF",
    ident.2 = "Control",
    group.by = "condition",   
    subset.ident = ctype,     
    features = features_plot, 
    min.pct = 0,              
    logfc.threshold = 0      
  )
  res$Gene <- rownames(res)
  res$CellType <- ctype
  de_results_list[[ctype]] <- res
}

plot_data <- do.call(rbind, de_results_list)

plot_data <- plot_data %>%
  mutate(
    log_p = -log10(p_val_adj + 1e-300),
    log_p_capped = ifelse(log_p > 50, 50, log_p),
    avg_log2FC_capped = case_when(
      avg_log2FC > 2 ~ 2,
      avg_log2FC < -2 ~ -2,
      TRUE ~ avg_log2FC
    )
  )

plot_data$Gene <- factor(plot_data$Gene, levels = rev(features_plot)) 
plot_data$CellType <- factor(plot_data$CellType, levels = rev(target_subtypes))

p <- ggplot(plot_data, aes(x = CellType, y = Gene)) +
  geom_point(aes(size = log_p_capped, color = avg_log2FC_capped)) +
  scale_color_gradient2(
    name = "log2FC\n(HFpEF vs Ctrl)",
    low = "#3C5488",   
    mid = "white",     
    high = "#E64B35",  
    midpoint = 0       
  ) +
  scale_size_continuous(
    name = "-log10(P-adj)", 
    range = c(1, 8)    
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),         
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  ) +
  labs(x = "", y = "", title = " ") +
  coord_flip()

print(p)


# Fig2.F --------------------------------------------------------------------

target_idents <- c("Ventricular Cardiomyocytes", 
                   "Stressed Cardiomyocytes")
cells_to_keep <- merged.object@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  filter(L5 %in% target_idents) %>%
  group_by(L5) %>%
  slice_sample(n = 10000) %>% 
  pull(barcode)

sc_small <- subset(merged.object, cells = cells_to_keep)
sc_small <- NormalizeData(sc_small)
sc_small <- FindVariableFeatures(sc_small, selection.method = "vst", nfeatures = 2000)
sc_small <- ScaleData(sc_small)
sc_small <- RunPCA(sc_small)
sc_small <- RunHarmony(sc_small, group.by = "sample")
sc_small <- RunUMAP(sc_small, dims = 1:15, n.neighbors = 30, min.dist = 0.1)
DimPlot(sc_small, reduction = "umap")
sce <- as.SingleCellExperiment(sc_small)
sce <- slingshot(sce, 
                 clusterLabels = 'L5', 
                 reducedDim = 'UMAP', 
                 start.clus = "Ventricular Cardiomyocytes")
pt_matrix <- slingPseudotime(sce)
pt_vec <- pt_matrix[, 1]
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- rep("grey80", length(pt_vec))
valid_cells <- !is.na(pt_vec) 
plotcol[valid_cells] <- colors[cut(pt_vec[valid_cells], breaks=100)]
plot(reducedDims(sce)$UMAP, 
     col = plotcol, 
     pch = 16, 
     asp = 1, 
     cex = 0.6, 
     main = "Pseudotime")
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
df_cm <- df %>%
  filter(celltype %in% c("Ventricular Cardiomyocytes",
                         "Stressed Cardiomyocytes")) %>%
  filter(!is.na(pt))
ks_res <- ks.test(
  df_cm$pt[df_cm$celltype == "Ventricular Cardiomyocytes"],
  df_cm$pt[df_cm$celltype == "Stressed Cardiomyocytes"]
)
ks_res$p.value
pt_values <- slingPseudotime(sce)[, 1]
cell_types <- colData(sce)$L5 
df <- data.frame(Pseudotime = pt_values, CellType = cell_types)
result_table <- df %>%
  group_by(CellType) %>%
  summarise(
    Mean_Time = mean(Pseudotime, na.rm = TRUE), 
    Median_Time = median(Pseudotime, na.rm = TRUE), 
    Cell_Count = n()
  ) %>%
  arrange(Mean_Time)
print(result_table)

plot_df <- data.frame(
  UMAP_1 = reducedDims(sce)$UMAP[, 1],
  UMAP_2 = reducedDims(sce)$UMAP[, 2],
  Pseudotime = slingPseudotime(sce)[, 1],
  CellType = colData(sce)$L5
)

ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = c("#466EB0", "#57B2AB", "#FDD683", "#F16844", "#9E0142")) +
  theme_classic() +
  labs(title = "Pathological Trajectory")


# Extract the proposed sequence
pt <- slingPseudotime(sce)[, 1]
gene_name <- "Pdk4" 
gene_name <- "Fabp3" 
gene_name <- "Acta1" 
gene_name <- "Txnip" 

expr <- as.numeric(as.matrix(logcounts(sce)[gene_name, ]))
group_info <- colData(sce)$condition
print(length(group_info)) 
plot_df <- data.frame(
  Pseudotime = pt,
  Expression = expr,
  Group = group_info
)
plot_df <- plot_df %>% filter(!is.na(Pseudotime))

my_colors <- c("Control" = "#4E84C4", "HFpEF" = "#E63946") 
ggplot(plot_df, aes(x = Pseudotime, y = Expression, color = Group)) +
  geom_point(alpha = 0.4, size = 0.8, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), size = 1.5, se = TRUE) +
  scale_color_manual(values = my_colors) +
  labs(
    title = " ",
    x = "Pseudotime (Healthy -> Stressed)",
    y = "Normalized Expression"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    axis.line = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )







