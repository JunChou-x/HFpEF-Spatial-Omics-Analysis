
library(dplyr)
library(ggplot2)
library(scattermore)
library(Seurat)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(CellChat)
library(tibble)

# Fig2.A ----------------------------------------------------------------------

target_L5_types <- c("Resident Fibroblasts", "Ventricular-associated Fibroblasts", 
                     "Atrial-associated Fibroblasts", "Perivascular Fibroblasts",
                     "Myofibroblasts","Stress-Response Fibroblasts",
                     "Fibro-Immune Interface")

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

u1_fibro <- plot_data_CU %>% 
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

print(u1_fibro)

d1_fibro <- plot_data_CL %>% 
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

print(d1_fibro)

u2_hu_fibro <- plot_data_HU %>% 
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

print(u2_hu_fibro)

d2_hl_fibro <- plot_data_HL %>% 
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

print(d2_hl_fibro)


# Fig2.B ------------------------------------------------------------------

DefaultAssay(merged.object) <- "RNA"
target_fbs <- rev(c("Resident Fibroblasts", "Ventricular-associated Fibroblasts", 
                    "Atrial-associated Fibroblasts", "Perivascular Fibroblasts",
                    "Myofibroblasts","Stress-Response Fibroblasts",
                    "Fibro-Immune Interface"))

Idents(merged.object) <- "L5"
fb_subset <- subset(merged.object, idents = target_fbs)
table(fb_subset$L5)
table(fb_subset$condition)

fb_markers <- FindAllMarkers(
  fb_subset,
  only.pos = TRUE,       
  min.pct = 0.25,        
  logfc.threshold = 0.25,
  max.cells.per.ident = 500 
)

fb_markers_sig <- fb_markers %>% 
  filter(p_val_adj < 0.05)

top5_genes <- fb_markers_sig %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top5_genes , n = 100)
top5_genes$cluster <- factor(top5_genes$cluster, levels = target_fbs)
top5_genes <- top5_genes %>% arrange(cluster)

genes_to_plot <- unique(top5_genes$gene)

levels_order <- rev(c("Resident Fibroblasts", "Ventricular-associated Fibroblasts", 
                      "Atrial-associated Fibroblasts", "Perivascular Fibroblasts",
                      "Myofibroblasts","Stress-Response Fibroblasts",
                      "Fibro-Immune Interface"))
Idents(fb_subset) <- factor(Idents(fb_subset), levels = levels_order)


DotPlot(
  fb_subset,
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



# Fig2.C ------------------------------------------------------------------

merged.object$condition <- ifelse(grepl("Control", merged.object$sample), "Control", "HFpEF")
merged.object$condition <- factor(merged.object$condition, levels = c("Control", "HFpEF"))
Idents(merged.object) <- "L5"
target_fbs <- c(
  "Resident Fibroblasts", "Atrial-associated Fibroblasts",
  "Ventricular-associated Fibroblasts", "Perivascular Fibroblasts",
  "Myofibroblasts","Stress-Response Fibroblasts",
  "Fibro-Immune Interface"
)

all_markers_list <- list()
for (cell_type in target_fbs) {
  message(paste0("calculating: ", cell_type, " (HFpEF vs Control)..."))
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
  top_n(n = 50, wt = abs(avg_log2FC)) 


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

genes_resident      <- c("Pdgfd", "Smoc1", "Tspan7", "Bckdhb")
genes_atrial      <- c("Ankrd1", "Bgn", "Myl4", "Ccl21a", "Pdzd4", "Cd209g")
genes_ventricular <- c("Fabp3", "Pdk4", "Nr1d1", "Lpl", "Hadhb", "Clu")
genes_perivascular  <- c("Nr1d1", "Zbtb16", "Dbp", "Isg15", "Ddit4", "Nr4a1", "Fos", "Clu", "Dusp1")
genes_myo         <- c("Cilp", "Comp", "Txnip", "Pdk4", "Upk3b", "Aqp1")
genes_stress      <- c("Bmp10", "Lrtm1", "Nr1d1", "Angptl7", "Prg4", "Mgp", "Chad")
genes_immune      <- c("Ucp2", "Depp1", "Fabp3", "Cxcl13", "Upk3b", "Krt19")

target_map <- bind_rows(
  data.frame(cell_type = "Resident Fibroblasts",          gene = genes_resident),
  data.frame(cell_type = "Atrial-associated Fibroblasts",          gene = genes_atrial),
  data.frame(cell_type = "Ventricular-associated Fibroblasts",     gene = genes_ventricular),
  data.frame(cell_type = "Perivascular Fibroblasts",      gene = genes_perivascular),
  data.frame(cell_type = "Myofibroblasts",              gene = genes_myo),
  data.frame(cell_type = "Stress-Response Fibroblasts", gene = genes_stress),
  data.frame(cell_type = "Fibro-Immune Interface",      gene = genes_immune)
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
    min.segment.length = 0, 
    segment.color = "black",
    segment.size = 0.3
  ) +
  facet_wrap(~cell_type, scales = "free", nrow = 4) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj P-value)",
    title = " ",
    color = "Expression Status"
  )

print(p_volcano)


# Fig2.D ------------------------------------------------------------------


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







# Fig2.E ------------------------------------------------------------------

cellchat <- readRDS("./cellchat.rds")
cellchat_ctrl <- readRDS("./cellchat_ctrl.rds")
cellchat_hfpef <- readRDS("./cellchat_hfpef.rds")

set.seed(1424)
target_samples <- c("Control_Lower", "Control_Upper")
max_cells_per_group <- 3000 

meta_df <- merged.object@meta.data %>% 
  tibble::rownames_to_column("barcode") %>%
  filter(sample %in% target_samples) 

sampled_barcodes <- meta_df %>%
  group_by(sample, L5) %>% 
  slice_sample(n = max_cells_per_group) %>% 
  pull(barcode)

sampled_counts <- table(meta_df$L5[meta_df$barcode %in% sampled_barcodes])
print(table(meta_df$L5[meta_df$barcode %in% sampled_barcodes]))
cat(sum(sampled_counts), "cells\n")

object_ctrl <- subset(merged.object, cells = sampled_barcodes)

meta_ctrl <- object_ctrl@meta.data
meta_ctrl$L5 <- as.factor(as.character(meta_ctrl$L5))

data.input_ctrl <- GetAssayData(object_ctrl, assay = "RNA", layer = "counts")
data.input_ctrl <- as(data.input_ctrl, "dgCMatrix")

cellchat_ctrl <- createCellChat(
  object = data.input_ctrl,
  meta = meta_ctrl,
  group.by = "L5"
)
cellchat_ctrl@DB <- CellChatDB.mouse

cellchat_ctrl <- subsetData(cellchat_ctrl)
cellchat_ctrl <- identifyOverExpressedGenes(cellchat_ctrl)
cellchat_ctrl <- identifyOverExpressedInteractions(cellchat_ctrl)
cellchat_ctrl <- computeCommunProb(cellchat_ctrl, raw.use = TRUE, population.size = TRUE)
cellchat_ctrl <- filterCommunication(cellchat_ctrl, min.cells = 10)
cellchat_ctrl <- computeCommunProbPathway(cellchat_ctrl)
cellchat_ctrl <- aggregateNet(cellchat_ctrl)

set.seed(1424) 

target_samples_hfpef <- c("HFpEF_Lower", "HFpEF_Upper")
max_cells_per_group <- 3000 
meta_df <- merged.object@meta.data %>% 
  tibble::rownames_to_column("barcode") %>%
  filter(sample %in% target_samples_hfpef) # 只筛选 HFpEF 组

sampled_barcodes_hfpef <- meta_df %>%
  group_by(sample, L5) %>% 
  slice_sample(n = max_cells_per_group) %>% 
  pull(barcode)

object_hfpef <- subset(merged.object, cells = sampled_barcodes_hfpef)
meta_hfpef <- object_hfpef@meta.data
meta_hfpef$L5 <- as.factor(as.character(meta_hfpef$L5))
meta_hfpef$L5 <- droplevels(meta_hfpef$L5)

data.input_hfpef <- GetAssayData(object_hfpef, assay = "RNA", layer = "counts")
data.input_hfpef <- as(data.input_hfpef, "dgCMatrix")

cellchat_hfpef <- createCellChat(
  object = data.input_hfpef,
  meta = meta_hfpef,
  group.by = "L5"
)

cellchat_hfpef@DB <- CellChatDB.mouse
cellchat_hfpef <- subsetData(cellchat_hfpef)
cellchat_hfpef <- identifyOverExpressedGenes(cellchat_hfpef)
cellchat_hfpef <- identifyOverExpressedInteractions(cellchat_hfpef)
cellchat_hfpef <- computeCommunProb(cellchat_hfpef, raw.use = TRUE, population.size = TRUE)
cellchat_hfpef <- filterCommunication(cellchat_hfpef, min.cells = 10)
cellchat_hfpef <- computeCommunProbPathway(cellchat_hfpef)
cellchat_hfpef <- aggregateNet(cellchat_hfpef)


# Merge
object.list <- list(Control = cellchat_ctrl, HFpEF = cellchat_hfpef)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

df.net <- data.frame(
  group = names(cellchat@net),
  count = sapply(cellchat@net, function(x) sum(x$weight > 0)),
  weight = sapply(cellchat@net, function(x) sum(x$weight))
)

ggplot(df.net, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  labs(
    x = "",
    y = "Interaction Count"
  ) +
  scale_fill_manual(values = c("#4c9aff", "#ff4c4c")) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14)
  )

ggplot(df.net, aes(x = group, y = weight, fill = group)) +
  geom_bar(stat = "identity") +
  labs(
    x = "",
    y = "Total Weight"
  ) +
  scale_fill_manual(values = c("#4c9aff", "#ff4c4c")) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14)
  )



# Fig2.F -----------------------------------------------------------------

get_signal_strength <- function(object_list) {
  df_all <- data.frame()
  for (name in names(object_list)) {
    cc <- object_list[[name]]
    mat <- cc@net$weight
    outgoing <- rowSums(mat)
    incoming <- colSums(mat)

    df <- data.frame(
      CellType = names(outgoing),
      Outgoing = outgoing,
      Incoming = incoming,
      Group = name
    )
    df_all <- rbind(df_all, df)
  }
  return(df_all)
}

plot_data <- get_signal_strength(object.list)

p1 <- ggplot(plot_data, aes(x = CellType, y = Outgoing, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Outgoing Signal Strength (Sender)", y = "Interaction Strength") +
  scale_fill_manual(values = c("Control" = "#3C5488", "HFpEF" = "#E64B35"))

p2 <- ggplot(plot_data, aes(x = CellType, y = Incoming, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Incoming Signal Strength (Receiver)", y = "Interaction Strength") +
  scale_fill_manual(values = c("Control" = "#3C5488", "HFpEF" = "#E64B35"))

print(p1)
print(p2)

cell_order <- c(
  "Ventricular Cardiomyocytes",
  "Atrial Cardiomyocytes", 
  "Conduction System Cardiomyocytes",
  "Stressed Cardiomyocytes",
  
  "Ventricular-associated Fibroblasts",
  "Atrial-associated Fibroblasts",
  "Resident Fibroblasts",
  "Myofibroblasts",
  "Perivascular Fibroblasts",
  "Stress-Response Fibroblasts",
  "Fibro-Immune Interface",
  
  "Vascular Endothelial Cells",
  "Endocardial Endothelial Cells",
  "Lymphatic Endothelial Cells",
  
  "Ventricular-associated SMCs",
  "Atrial-associated SMCs",
  "Arterial SMCs",
  "Pericytes",
  
  "B Cells",
  
  "Mesothelial Cells",
  "Schwann Cells",
  "Neurons",
  "Proliferating Cells",
  "Brown Adipocytes",
  "White Adipocytes"
)

plot_data$CellType <- factor(plot_data$CellType, levels = cell_order)

p1 <- ggplot(plot_data, aes(x = CellType, y = Outgoing, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt")
  ) +
  labs(title = "Outgoing Signal Strength (Sender)", y = "Interaction Strength") +
  scale_fill_manual(values = c("Control" = "#3C5488", "HFpEF" = "#E64B35"))

p2 <- ggplot(plot_data, aes(x = CellType, y = Incoming, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt")
  ) +
  labs(title = "Incoming Signal Strength (Receiver)", y = "Interaction Strength") +
  scale_fill_manual(values = c("Control" = "#3C5488", "HFpEF" = "#E64B35"))

print(p1)
print(p2)

# Fig2.G -----------------------------------------------------------------

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


# Fig2.H -----------------------------------------------------------------

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

