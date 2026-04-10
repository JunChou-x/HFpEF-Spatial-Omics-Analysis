
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


# SFig4.C -----------------------------------------------------------------

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



