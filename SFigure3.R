library(ggplot2)
library(scales)
library(dplyr)

# SFig3.A -----------------------------------------------------------------

gene <- "Dcn, Postn, Col3a1"
expr <- FetchData(merged.object, vars = gene)
expr$barcode <- rownames(expr)
expr$barcode_nosample <- sub(paste0("_", Sample), "", expr$barcode)

Sample <- "Control_Upper"
CU_DF$Gene <- expr[[gene]][match(paste0(CU_DF$barcode, paste0("_", Sample)), expr$barcode)]

CU_DF %>%
  filter(tissue == 1, !is.na(Gene)) %>%
  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled, color = Gene)) +
  geom_scattermore(pointsize = 3, pixels = rep(2000, 2)) +
  coord_fixed(ratio = 1, expand = FALSE) + 
  scale_color_gradientn(
    colors = c("lightgrey", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FDE725"),
    values = seq(0, 1, length.out = 6),  
    limits = c(0, 4.5),                  
    oob = scales::squish               
  ) +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank()
  ) +
  labs(color = gene, title = paste(" "))+
  annotate("segment",
           x = bar_x_start_CU,
           xend = bar_x_end_CU, 
           y = bar_y_pos_CU,
           yend = bar_y_pos_CU,
           color = "black",
           linewidth = 1.5)

# SFig3.B -----------------------------------------------------------------

data_to_plot <- fb_subset@meta.data %>% 
  dplyr::select(L5, condition) %>%
  as.data.frame()
l5_order <- levels(fb_subset$L5)
data_to_plot$L5 <- factor(data_to_plot$L5, levels = l5_order)
data_to_plot$condition <- factor(data_to_plot$condition, levels = c("HFpEF", "Control"))

p <- ggplot(data_to_plot, aes(x = condition, fill = L5)) +
  geom_bar(position = "fill", width = 0.65) +  
  scale_y_continuous(labels = percent) +
  coord_flip() +
  scale_fill_manual(values = ColsL5) +
  labs(
    title = " ",
    x = NULL,
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

p



# SFig3.C ------------------------------------------------------------------


ora_results_directional_list <- list()
for (current_type in target_fbs) {
  
  cat(paste0("\n=== Analyzing: ", current_type, " ===\n"))
  de_data_unique <- final_DE_results %>% 
    filter(cell_type == current_type) %>% 
    filter(!gene %in% c("Hba-a2", "Hbb-bs")) 
  if(nrow(de_data_unique) < 5) { 
    cat(paste0("  [Skip]  (<5)\n"))
    next 
  }
  ids <- tryCatch({
    bitr(de_data_unique$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  }, error = function(e) { return(NULL) })
  if (is.null(ids)) { cat("  [Skip] ID conversion failed\n"); next }
  de_data_unique <- merge(de_data_unique, ids, by.x = "gene", by.y = "SYMBOL") %>%
    distinct(ENTREZID, .keep_all = TRUE)
  sig_genes_up <- de_data_unique %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(ENTREZID)
  sig_genes_down <- de_data_unique %>% filter(p_val_adj < 0.05, avg_log2FC < 0) %>% pull(ENTREZID)
  current_cell_results <- list()
  if(length(sig_genes_up) > 5) { 
    ora_up <- tryCatch({
      enrichGO(gene = sig_genes_up, OrgDb = org.Mm.eg.db, ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    }, error = function(e) NULL)
    
    if (!is.null(ora_up) && nrow(as.data.frame(ora_up)) > 0) {
      current_cell_results[["Up"]] <- as.data.frame(ora_up) %>% 
        mutate(Regulation = "Up-regulated", 
               cell_type = current_type) 
      cat(paste0(" Upregulate pathway: ", nrow(current_cell_results[["Up"]]), " 条\n"))
    } else {
      cat("No significant upregulation pathway\n")
    }
  } else { cat("Upregulation of gene deficiency\n") }
  
  if(length(sig_genes_down) > 5) {
    ora_down <- tryCatch({
      enrichGO(gene = sig_genes_down, OrgDb = org.Mm.eg.db, ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    }, error = function(e) NULL)
    
    if (!is.null(ora_down) && nrow(as.data.frame(ora_down)) > 0) {
      current_cell_results[["Down"]] <- as.data.frame(ora_down) %>% 
        mutate(Regulation = "Down-regulated", 
               cell_type = current_type) 
      cat(paste0(" Downregulation pathway: ", nrow(current_cell_results[["Down"]]), "\n"))
    } else {
      cat("No significant downregulation pathways\n")
    }
  } else { cat("Downregulation of genes is insufficient\n") }
  
  if (length(current_cell_results) > 0) {
    ora_results_directional_list[[current_type]] <- bind_rows(current_cell_results)
  }
}

plot_data_directional <- bind_rows(ora_results_directional_list)

if (nrow(plot_data_directional) == 0) {
  stop("No significant pathways were found in any of the subgroups.")
}

plot_data_filtered <- plot_data_directional %>%
  mutate(Description = str_trunc(Description, width = 50)) %>%
  group_by(cell_type, Regulation) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>% 
  ungroup() 
all_pathway_data <- dplyr::bind_rows(ora_results_directional_list, .id = "cell_type")
top_pathways_summary <- all_pathway_data %>%
  dplyr::group_by(cell_type, Regulation) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 20) %>% 
  dplyr::ungroup() %>%
  dplyr::select(cell_type, Regulation, Description, p.adjust, Count) 
unique_types <- unique(top_pathways_summary$cell_type)
for (ctype in unique_types) {
  cat(paste0("\n>>> Cell subsets: ", ctype, " <<<\n"))
  up_data <- top_pathways_summary %>% 
    dplyr::filter(cell_type == ctype, stringr::str_detect(Regulation, "Up"))
  if(nrow(up_data) > 0) {
    cat(paste0("  [↑ Up-regulated] (Top ", nrow(up_data), ")\n"))
    print_df <- up_data %>% 
      dplyr::mutate(Description = stringr::str_trunc(Description, 40)) %>% 
      dplyr::select(Description, p.adjust, Count) %>% 
      as.data.frame()
    print(print_df, row.names = FALSE)
  } else {
    cat("No significant upregulation pathway\n")
  }
  cat("\n")
  down_data <- top_pathways_summary %>% 
    dplyr::filter(cell_type == ctype, stringr::str_detect(Regulation, "Down"))
  if(nrow(down_data) > 0) {
    cat(paste0("  [↓ Down-regulated] (Top ", nrow(down_data), ")\n"))
    print_df_down <- down_data %>% 
      dplyr::mutate(Description = stringr::str_trunc(Description, 40)) %>% 
      dplyr::select(Description, p.adjust, Count) %>% 
      as.data.frame()
    
    print(print_df_down, row.names = FALSE)
  } else {
    cat("No significant downregulation pathways\n")
  }
  
  cat("--------------------------------------------------------\n")
}


plot_data_directional <- dplyr::bind_rows(ora_results_directional_list, .id = "cell_type")
target_pathways <- c(
  # --- 1. Ventricular-associated Fibroblasts ---
  "fatty acid metabolic process",
  "regulation of cellular ketone metabolic process",
  "cholesterol homeostasis",
  "cardiac muscle tissue morphogenesis",
  "muscle system process",
  "cardiac muscle contraction",
  
  # --- 2. Stress-Response Fibroblasts ---
  "regulation of lipid biosynthetic process",
  "cellular ketone metabolic process",
  "SMAD protein signal transduction",
  "actomyosin structure organization",
  "myofibril assembly",
  "striated muscle cell development",
  
  # --- 3. Fibro-Immune Interface ---
  "fatty acid oxidation",
  "adaptive thermogenesis",
  "lipid modification",
  "heart morphogenesis",
  
  # --- 4. Myofibroblasts ---
  "protein secretion",
  "acute-phase response",
  "response to mercury ion",
  "response to copper ion",
  
  # --- 5. Perivascular Fibroblasts ---
  "regulation of blood circulation",
  "regulation of ATP-dependent activity",
  "skeletal muscle cell differentiation",
  "response to xenobiotic stimulus",
  "miRNA metabolic process",
  
  # --- 6. Atrial-associated Fibroblasts ---
  "regulation of the force of heart contraction",
  "cardiac ventricle morphogenesis",
  
  # --- 7. Resident Fibroblasts ---
  "cellular response to leucine starvation"
)

cell_type_order <- c(
  "Ventricular-associated Fibroblasts",
  "Stress-Response Fibroblasts",
  "Fibro-Immune Interface",
  "Myofibroblasts",
  "Perivascular Fibroblasts",
  "Atrial-associated Fibroblasts",
  "Resident Fibroblasts"
)

plot_data_filtered <- plot_data_directional %>%
  filter(tolower(Description) %in% tolower(target_pathways)) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
  mutate(Description = str_trunc(Description, width = 50)) 


p <- ggplot(plot_data_filtered, aes(x = Count, y = reorder(Description, -p.adjust))) +
  geom_segment(aes(x = 0, xend = Count, y = Description, yend = Description), 
               color = "grey80", linewidth = 0.5) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#E64B35", high = "#3C5488", trans = "log10", name = "Adj. P-value") +
  facet_grid(Regulation ~ cell_type, scales = "free", space = "free_y") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), 
    axis.text.y = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = "#F0F0F0", color = "black"), 
    strip.text = element_text(face = "bold", size = 10), 
    strip.text.x = element_text(angle = 0), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  
  labs(
    title = " ",
    x = "Gene Count",
    y = NULL
  )

print(p)



# SFig3.D -----------------------------------------------------------------

genes_to_plot <- c("Pdk4", "Fabp3", "Cpt1b", "Comp")
fb_subset$L5 <- droplevels(fb_subset$L5)

VlnPlot(fb_subset, 
        features = genes_to_plot, 
        group.by = "L5", 
        split.by = "condition", 
        pt.size = 0, 
        ncol = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# SFig3.E -----------------------------------------------------------------

fibro_idents <- c(
  "Resident Fibroblasts", "Atrial-associated Fibroblasts",
  "Ventricular-associated Fibroblasts", "Perivascular Fibroblasts",
  "Myofibroblasts","Stress-Response Fibroblasts",
  "Fibro-Immune Interface"
)

fibro <- subset(merged.object, subset = L5 %in% fibro_idents)

counts_matrix <- GetAssayData(fibro, assay = "RNA", layer = "counts")
counts_matrix <- as(counts_matrix, "dgCMatrix")

counts <- SeuratObject::GetAssayData(fibro, slot = "counts")
if (is.null(counts) || length(counts) == 0) {
  stop("'counts' data not found")
}

meta_data <- fibro@meta.data

initial_obj <- CreateSeuratObject(
  counts = counts, 
  project = fibro@project.name, 
  meta.data = meta_data
)

cleaned_meta_data <- initial_obj@meta.data

for (col_name in colnames(cleaned_meta_data)) {
  if (is.factor(cleaned_meta_data[[col_name]])) {
    message(paste(" Transform factor column:", col_name))
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
  output_path = "./pysenic_fb.h5ad",
  mode = "sc",
  assay = "RNA"
)

# to Python

# SFig3.F -----------------------------------------------------------------

cellchat_ctrl@idents <- as.factor(cellchat_ctrl@meta$L5)
names(cellchat_ctrl@idents) <- colnames(cellchat_ctrl@data)

cellchat_hfpef@idents <- as.factor(cellchat_hfpef@meta$L5)
names(cellchat_hfpef@idents) <- colnames(cellchat_hfpef@data)

cellchat_ctrl@idents <- as.factor(cellchat_ctrl@meta$L5)
names(cellchat_ctrl@idents) <- colnames(cellchat_ctrl@data)

cellchat_hfpef@idents <- as.factor(cellchat_hfpef@meta$L5)
names(cellchat_hfpef@idents) <- colnames(cellchat_hfpef@data)

group.cellType <- union(levels(cellchat_ctrl@idents), levels(cellchat_hfpef@idents))

cellchat_ctrl <- liftCellChat(cellchat_ctrl, group.new = group.cellType)
cellchat_hfpef <- liftCellChat(cellchat_hfpef, group.new = group.cellType)

object.list <- list(Control = cellchat_ctrl, HFpEF = cellchat_hfpef)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat@idents <- as.factor(cellchat@meta$L5)
names(cellchat@idents) <- rownames(cellchat@meta)

gg1 <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "count")
gg1 <- gg1 + ggtitle("Differential Interactions (Count)")

gg2 <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
gg2 <- gg2 + ggtitle("Differential Interactions (Strength)")

gg1 + gg2



# SFig3.G -----------------------------------------------------------------

sources.use <- c(
  "Ventricular-associated Fibroblasts", "Atrial-associated fibroblasts", "Resident fibroblasts",
  "Myofibroblasts", "Perivascular Fibroblasts", "Stress-Response Fibroblasts",
  "Fibro-Immune Interface"
)

targets.use <- c(
  "Ventricular Cardiomyocytes", "Atrial Cardiomyocytes", 
  "Stressed Cardiomyocytes", "Conduction System Cardiomyocytes",
  
  "Ventricular-associated Fibroblasts", "Atrial-associated Fibroblasts", "Resident fibroblasts",
  "Myofibroblasts", "Perivascular Fibroblasts", "Stress-Response Fibroblasts", "Fibro-Immune Interface"
)

p <- netVisual_bubble(cellchat, 
                      sources.use = sources.use, 
                      targets.use = targets.use, 
                      comparison = c(1, 2), 
                      angle.x = 45)


# 5. 打印最终图像
print(p)


# SFig3.H -----------------------------------------------------------------

pathways.show <- c("THBS") 
group.cell <- levels(cellchat@idents)

# Control 
netVisual_aggregate(cellchat_ctrl, 
                    signaling = pathways.show, 
                    layout = "chord", 
                    title.name = "Control: THBS Signaling",
                    show.legend = FALSE)
# HFpEF 
netVisual_aggregate(cellchat_hfpef, 
                    signaling = pathways.show, 
                    layout = "chord", 
                    title.name = "HFpEF: THBS Signaling",
                    show.legend = FALSE)







                         
