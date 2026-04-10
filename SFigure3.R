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

# SFig3.C -----------------------------------------------------------------

genes_to_plot <- c("Pdk4", "Fabp3", "Cpt1b", "Comp")
fb_subset$L5 <- droplevels(fb_subset$L5)

VlnPlot(fb_subset, 
        features = genes_to_plot, 
        group.by = "L5", 
        split.by = "condition", 
        pt.size = 0, 
        ncol = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# SFig3.D -----------------------------------------------------------------

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

# SFig3.E -----------------------------------------------------------------

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

