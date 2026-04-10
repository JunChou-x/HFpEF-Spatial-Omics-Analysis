
library(dplyr)
library(ggplot2)
library(scattermore)
library(stringr)

# Fig4.A ------------------------------------------------------------------

target_L5_types <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells", 
                     "Lymphatic Endothelial Cells")

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

u1_endo <- plot_data_CU %>% 
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

print(u1_endo)

d1_endo <- plot_data_CL %>% 
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

print(d1_endo)

u2_hu_endo <- plot_data_HU %>% 
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

print(u2_hu_endo)

d2_hl_endo <- plot_data_HL %>% 
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

print(d2_hl_endo)


# Fig4.B ------------------------------------------------------------------

DefaultAssay(merged.object) <- "RNA"
target_ecs <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells",
                "Lymphatic Endothelial Cells")
Idents(merged.object) <- "L5"
ec_subset <- subset(merged.object, idents = target_ecs)

ec_markers <- FindAllMarkers(
  ec_subset,
  only.pos = TRUE,       
  min.pct = 0.25,       
  logfc.threshold = 0.25,
  max.cells.per.ident = 500 
)

ec_markers_sig <- ec_markers %>% 
  filter(p_val_adj < 0.05)
top5_genes <- ec_markers_sig %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
top5_genes$cluster <- factor(top5_genes$cluster, levels = target_ecs)
top5_genes <- top5_genes %>% arrange(cluster)
genes_to_plot <- unique(top5_genes$gene)

levels_order <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells",
                  "Lymphatic Endothelial Cells")
Idents(ec_subset) <- factor(Idents(ec_subset), levels = levels_order)

DotPlot(
  ec_subset,
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


# Fig4.C ------------------------------------------------------------------

merged.object$condition <- ifelse(grepl("Control", merged.object$sample), "Control", "HFpEF")
merged.object$condition <- factor(merged.object$condition, levels = c("Control", "HFpEF"))

Idents(merged.object) <- "L5"

target_ecs <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells",
                "Lymphatic Endothelial Cells")

all_markers_list <- list()
for (cell_type in target_ecs) {
  message(paste0("Calculating: ", cell_type, " (HFpEF vs Control)..."))
  markers <- FindMarkers(
    merged.object,
    ident.1 = "HFpEF",
    ident.2 = "Control",
    group.by = "condition",   
    subset.ident = cell_type, 
    max.cells.per.ident = 1000, 
    min.pct = 0,           
    logfc.threshold = 0,   
    return.thresh = 1      
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
top_changes <- sig_DE_results %>%
  group_by(cell_type) %>%
  top_n(n = 100, wt = abs(avg_log2FC)) 

plot_data <- final_DE_results %>%
  mutate(
    p_adj_calc = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0]) * 0.1, p_val_adj),
    log10P = -log10(p_adj_calc),
    change = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.1 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.1 ~ "Down",
      TRUE ~ "NS"
    )
  )

genes_ves   <- c("Pdk4", "Fabp3", "Cpt1b", "Txnip", "Nr1d1", "Depp1", "Lpl", "Aqp1", "Ptgds")
genes_endo  <- c("Cd36", "Pdk4", "Nr1d1", "Ddit4", "Angpt1", "Clu")
genes_lymph <- c("Col3a1", "Igfbp5", "Txnip", "Nr1d1", "Lpl")

target_map <- bind_rows(
  data.frame(cell_type = "Vascular Endothelial Cells",    gene = genes_ves),
  data.frame(cell_type = "Endocardial Endothelial Cells", gene = genes_endo),
  data.frame(cell_type = "Lymphatic Endothelial Cells",   gene = genes_lymph)
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
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.5,
    min.segment.length = 0, 
    segment.color = "black",
    segment.size = 0.3
  ) +
  facet_wrap(~cell_type, scales = "free", nrow = 2) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj P-value)",
    title = "Endothelial Subtypes Volcano Plot",
    color = "Expression Status"
  )

print(p_volcano)





# Fig4.D ------------------------------------------------------------------

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(stringr)
library(ggplot2)

target_ecs <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells",
                "Lymphatic Endothelial Cells")

gsea_results_list <- list()

for (current_type in target_ecs) {
  
  cat(paste0("\n\n======================================================\n"))
  cat(paste0("Analyzing: ", current_type, "\n"))
  cat(paste0("======================================================\n"))
  
  de_data <- final_DE_results %>% 
    dplyr::filter(cell_type == current_type) %>%  
    dplyr::filter(!gene %in% c("Hba-a2", "Hbb-bs"))
  
  ids <- tryCatch({
    bitr(de_data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  }, error = function(e) { return(NULL) })
  if (is.null(ids)) { cat("Conversion failed\n"); next }
  de_data_merged <- merge(de_data, ids, by.x = "gene", by.y = "SYMBOL")
  if(nrow(de_data_merged) == 0) { cat("No effective genes\n"); next }
  de_data_unique <- de_data_merged %>%
    dplyr::arrange(dplyr::desc(abs(avg_log2FC))) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  geneList_gsea <- de_data_unique$avg_log2FC
  names(geneList_gsea) <- de_data_unique$ENTREZID
  geneList_gsea <- sort(geneList_gsea, decreasing = TRUE)
  
  gsea_res <- tryCatch({
    gseGO(geneList     = geneList_gsea,
          OrgDb        = org.Mm.eg.db,
          ont          = "BP",
          minGSSize    = 10,
          maxGSSize    = 500,
          pvalueCutoff = 1, 
          verbose      = FALSE,
          seed         = 123)
  }, error = function(e) { return(NULL) })
  
  if (!is.null(gsea_res)) {
    gsea_results_list[[current_type]] <- gsea_res
    cat("\n  [GSEA] Calculated successfully.\n")
  } else { cat("No significant GSEA results\n") }
}

target_pathways_dict <- list(
  "Vascular Endothelial Cells" = tolower(c(
    "Natural killer cell proliferation",            
    "Neutrophil-mediated killing of bacterium",     
    "Positive regulation of triglyceride lipase activity", 
    "Platelet morphogenesis",                       
    "Platelet formation",                           
    "Negative regulation of phagocytosis"
  )),
  "Endocardial Endothelial Cells" = tolower(c(
    "Negative regulation of protein processing",     
    "Carnitine metabolic process",                   
    "Lipid digestion",                               
    "Mitochondrial fragmentation involved in apoptotic process", 
    "T cell receptor signaling pathway",             
    "Regulation of gastric acid secretion"
  )),
  "Lymphatic Endothelial Cells" = tolower(c(
    "Retinoic acid metabolic process",               
    "Regulation of lymphocyte chemotaxis",           
    "Cerebrospinal fluid circulation",               
    "Regulation of vesicle fusion",                  
    "Positive regulation of vesicle fusion",         
    "Regulation of synaptic vesicle membrane organization"
  ))
)

extracted_data_list <- list()

for (ctype in names(gsea_results_list)) {
  if (!is.null(gsea_results_list[[ctype]])) {
    res_df <- as.data.frame(gsea_results_list[[ctype]])
    target_paths <- target_pathways_dict[[ctype]]
    extracted_df <- res_df %>%
      dplyr::filter(tolower(Description) %in% target_paths) %>%
      dplyr::select(Description, NES, p.adjust) %>%
      dplyr::mutate(cell_type = ctype)
    extracted_data_list[[ctype]] <- extracted_df
  }
}

plot_data_gsea <- do.call(rbind, extracted_data_list)



plot_data_gsea$Description <- stringr::str_to_sentence(plot_data_gsea$Description)


plot_data_gsea$cell_type <- factor(plot_data_gsea$cell_type, 
                                   levels = c("Vascular Endothelial Cells", 
                                              "Endocardial Endothelial Cells", 
                                              "Lymphatic Endothelial Cells"))

plot_data_gsea$Description <- str_trunc(plot_data_gsea$Description, width = 45)

p_gsea <- ggplot(plot_data_gsea, aes(x = NES, y = reorder(Description, NES))) +
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), 
               color = "grey60", linewidth = 0.8) +
  geom_point(aes(color = NES), size = 5) +
  scale_color_gradientn(colours = c("#4575B4", "#91BFDB", "grey90", "#FC8D59", "#D73027"),
                        values = c(0, 0.4, 0.5, 0.6, 1),
                        name = "NES") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
    axis.text.y = element_text(size = 11, color = "black", face = "bold"), 
    axis.text.x = element_text(size = 10),
    strip.background = element_rect(fill = "white"), 
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    
    legend.position = "right"
  ) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = " "
  )

print(p_gsea)



# Fig4.E ------------------------------------------------------------------


Endo_idents <- c(
  "Vascular Endothelial Cells",
  "Endocardial Endothelial Cells",
  "Lymphatic Endothelial Cells"
)

Endo <- subset(merged.object, subset = L5 %in% Endo_idents)
counts_matrix <- GetAssayData(Endo, assay = "RNA", layer = "counts")
counts_matrix <- as(counts_matrix, "dgCMatrix")
counts <- SeuratObject::GetAssayData(Endo, slot = "counts")
if (is.null(counts) || length(counts) == 0) {
  stop("'counts' data not found")
}

meta_data <- Endo@meta.data
initial_obj <- CreateSeuratObject(
  counts = counts, 
  project = Endo@project.name, 
  meta.data = meta_data
)

cleaned_meta_data <- initial_obj@meta.data
for (col_name in colnames(cleaned_meta_data)) {
  if (is.factor(cleaned_meta_data[[col_name]])) {
    message(paste("Transform factor column", col_name))
    cleaned_meta_data[[col_name]] <- as.character(cleaned_meta_data[[col_name]])
  }
  if (is.character(cleaned_meta_data[[col_name]])) {
    if (any(is.na(cleaned_meta_data[[col_name]]))) {
      message(paste("Processing the character column'", col_name, "'NA values", sep=""))
      cleaned_meta_data[[col_name]][is.na(cleaned_meta_data[[col_name]])] <- "unknown"
    }
  }
}

initial_obj@meta.data <- cleaned_meta_data

GetH5ad(
  initial_obj,
  output_path = "./pysenic_ec.h5ad",
  mode = "sc",
  assay = "RNA"
)

# to Python







# Fig4.F ------------------------------------------------------------------

features_plot <- c(
  "Kdr", "Flt1", "Pecam1", "Tek", "Esm1",
  "Thbs1", "Angptl4","Txnip", "Depp1", "Igfbp5"
)
features_plot <- intersect(features_plot, rownames(merged.object))

target_subtypes <- c(
  "Vascular Endothelial Cells", 
  "Endocardial Endothelial Cells", 
  "Lymphatic Endothelial Cells"
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
  labs(x = "", y = "", title = " ")

print(p)



# Fig4.G,H ------------------------------------------------------------------

cellchat_ctrl <- readRDS("./cellchat_ctrl.rds")
cellchat_hfpef <- readRDS("./cellchat_hfpef.rds")

sources.use <- c("Stressed Cardiomyocytes","Perivascular Fibroblasts", "Resident Fibroblasts","Stress-Response Fibroblasts")
targets.use <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells","Lymphatic Endothelial Cells")

sources.use <- c("Vascular Endothelial Cells", "Endocardial Endothelial Cells","Lymphatic Endothelial Cells")
targets.use <-c("Pericytes","Ventricular-associated SMCs", "Arterial SMCs","Atrial-associated  SMCs")

p <- netVisual_bubble(cellchat, 
                      sources.use = sources.use, 
                      targets.use = targets.use, 
                      comparison = c(1, 2), 
                      angle.x = 45,
                      remove.isolate = FALSE, 
                      title.name = " ")
p + theme(
  plot.margin = unit(c(5, 5, 5, 25), "mm") 
)


# Fig4.I ------------------------------------------------------------------

mural_cells <- subset(merged.object, idents = c("Pericytes", "Ventricular-associated SMCs","Atrial-associated SMCs", "Arterial SMCs"))

features_mural <- c(
  "Pdk4", "Fabp3", "Cpt1b", "Angptl4", "Txnip",
  "Acta2", "Tagln", "Myh11", "Rgs5"
)
features_mural <- intersect(features_mural, rownames(mural_cells))

Idents(mural_cells) <- "L5"
target_subtypes <- c("Pericytes", "Ventricular-associated SMCs","Atrial-associated SMCs", "Arterial SMCs")

de_results_list <- list()
for (ctype in target_subtypes) {
  res <- FindMarkers(
    mural_cells,
    ident.1 = "HFpEF",
    ident.2 = "Control",
    group.by = "condition",
    subset.ident = ctype,
    features = features_mural,
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

plot_data$Gene <- factor(plot_data$Gene, levels = rev(features_mural))
plot_data$CellType <- factor(plot_data$CellType, levels = rev(target_subtypes))

p_mural <- ggplot(plot_data, aes(x = CellType, y = Gene)) +
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
    range = c(1,8)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  ) +
  labs(x = "", y = "", title = "")

p_mural


