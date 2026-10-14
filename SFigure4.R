
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(forcats)

# SFig4.A -----------------------------------------------------------------

gene <- "Pecam, Npr3, Ccl21a"

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

# SFig4.B -----------------------------------------------------------------

data_to_plot <- ec_subset@meta.data %>% 
  dplyr::select(L5, condition) %>%
  as.data.frame()

l5_order <- levels(ec_subset$L5)
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


data_to_plot <- mural_cells@meta.data %>% 
  dplyr::select(L5, condition) %>%
  as.data.frame()

l5_order <- levels(mural_cells$L5)
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



# SFig4.C,E ------------------------------------------------------------------

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


# SFig4.D -----------------------------------------------------------------

target_gene <- "Cd36"
gene_data <- FetchData(merged.object, vars = c(target_gene, "condition"))
colnames(gene_data) <- c("expression", "group")

stats <- gene_data %>%
  group_by(group) %>%
  summarise(
    Mean_Exp = mean(expression),            
    Pct_Exp = mean(expression > 0) * 100,   
    Cell_Count = n()                        
  )

data_cd36 <- FetchData(merged.object, vars = c("Cd36", "L5", "condition"))
colnames(data_cd36) <- c("Expression", "CellType", "Condition")

stats_summary <- data_cd36 %>%
  group_by(CellType, Condition) %>%
  summarise(
    Mean = mean(Expression),
    Pct = mean(Expression > 0) * 100,
    .groups = 'drop'
  )

comparison_table <- stats_summary %>%
  dplyr::select(CellType, Condition, Mean) %>%
  pivot_wider(names_from = Condition, values_from = Mean, values_fill = 0) %>%
  mutate(Difference = HFpEF - Control) %>% 
  arrange(desc(Difference))

pct_table <- stats_summary %>%
  dplyr::select(CellType, Condition, Pct) %>%
  pivot_wider(names_from = Condition, values_from = Pct, values_fill = 0) %>%
  mutate(Diff_Pct = HFpEF - Control) %>%
  arrange(desc(Diff_Pct))

print("--- CD36 Percent Expressed (%) ---")
print(pct_table, n = 30)

plot_data <- pct_table %>%
  mutate(CellType = fct_reorder(CellType, Diff_Pct))

p1 <- ggplot(plot_data) +
  geom_segment(aes(y = CellType, yend = CellType, 
                   x = Control, xend = HFpEF), 
               color = "grey", size = 1) +
  geom_point(aes(y = CellType, x = Control, color = "Control"), size = 3) +
  geom_point(aes(y = CellType, x = HFpEF, color = "HFpEF"), size = 3) +
  scale_color_manual(values = c("Control" = "#377EB8", "HFpEF" = "#E41A1C")) +
  labs(title = "Changes in Cd36 Expression Percentage (Control vs HFpEF)",
       x = "Percentage of Cells Expressing Cd36 (%)",
       y = NULL,
       color = "Group") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) 

print(p1)



# SFig4.F ------------------------------------------------------------------

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




