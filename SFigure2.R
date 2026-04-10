
library(ggplot2)
library(scales)
library(dplyr)
library(reshape2)
library(SCNT)

# SFig1.A --------------------------------------------------------------------

gene <- "Myl2, Nppa, Acta1, Cntn"
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
  # 添加1mm的bar
  annotate("segment",
           x = bar_x_start_CU,
           xend = bar_x_end_CU, 
           y = bar_y_pos_CU,
           yend = bar_y_pos_CU,
           color = "black",
           linewidth = 1.5)


# SFig1.B --------------------------------------------------------------------

data_to_plot <- cm_subset@meta.data %>% 
  dplyr::select(L5, condition) %>%
  as.data.frame()
l5_order <- levels(cm_subset$L5)
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


# SFig1.C --------------------------------------------------------------------

scObj <-readRDS("./scObj.rds")
expr_data3 <- FetchData(HF, vars = c("Nppb", "Nppa", "Ankrd1",  "Acta1","condition"))
expr_long <- melt(expr_data3, id.vars = "condition", 
                  variable.name = "Gene", value.name = "Expression")

ggplot(expr_long, aes(x = condition, y = Expression, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +  
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.4, color = "black") + 
  facet_wrap(~Gene, scales="free_y", nrow=1)+
  scale_fill_manual(values = c("Control" = "#4E84C4", "HFpEF" = "#E63946")) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill="#F0F0F0", color=NA),
    strip.text = element_text(face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face="bold")
  ) +
  ylab("Expression Level")


# SFig1.D -----------------------------------------------------------------

cells_to_keep <- merged.object@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  filter(L5 %in% c("Ventricular Cardiomyocytes", "Stressed Cardiomyocytes")) %>%
  group_by(sample, L5) %>%
  slice_sample(n = 3750) %>%
  pull(barcode)
cm_subset <- subset(merged.object, cells = cells_to_keep)
counts_matrix <- GetAssayData(cm_subset, assay = "RNA", layer = "counts")
counts_matrix <- as(counts_matrix, "dgCMatrix")

counts <- SeuratObject::GetAssayData(cm_subset, slot = "counts")
if (is.null(counts) || length(counts) == 0) {
  stop("verify whether the subject is correct.")
}

meta_data <- cm_subset@meta.data
initial_obj <- CreateSeuratObject(
  counts = counts, 
  project = cm_subset@project.name, 
  meta.data = meta_data
)


cleaned_meta_data <- initial_obj@meta.data
for (col_name in colnames(cleaned_meta_data)) {
  if (is.factor(cleaned_meta_data[[col_name]])) {
    message(paste("Convert the factor column:", col_name))
    cleaned_meta_data[[col_name]] <- as.character(cleaned_meta_data[[col_name]])
  }
  
  if (is.character(cleaned_meta_data[[col_name]])) {
    if (any(is.na(cleaned_meta_data[[col_name]]))) {
      message(paste("Process the character column'", col_name, "' NA ", sep=""))
      cleaned_meta_data[[col_name]][is.na(cleaned_meta_data[[col_name]])] <- "unknown"
    }
  }
}

initial_obj@meta.data <- cleaned_meta_data

GetH5ad(
  initial_obj,
  output_path = "./pysenic_cm.h5ad",
  mode = "sc",
  assay = "RNA"
)

#to Python









