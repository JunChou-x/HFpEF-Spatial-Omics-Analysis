library(pheatmap)
library(grid)
library(SCP)


# SFig1.I------------------------------------------------------------------

Idents(Merge) <- Merge$celltype
allmarkers <- FindAllMarkers(Merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
color=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
        "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
        "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
        "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
DimPlot(Merge, reduction = "umap", cols = color, group.by = "celltype")
CellDimPlot(Merge, 
            group.by = "celltype", 
            reduction = "umap", 
            cols = color,
            show_stat = FALSE)


# SFig1.J ----------------------------------------------------------------

markergenes <- c(
  "Tnnt2", "Myh6", "Myh7",        # Cardiomyocyte
  "Acta2", "Myh11", "Tagln",      # Smooth Muscle
  "Pdgfrb", "Rgs5", "Notch3",     # Pericyte
  "Cd3e", "Cd4","Cd8a",           # T Cell
  "Cd68", "Cd163", "Adgre1"  ,    # Macrophage 
  "Col1a1", "Col3a1","Pdgfra",    # Fibroblast
  "Vwf", "Pecam1", "Plpp1",       # Endothelial
  "Cd79a", "Pax5","Ms4a1",        # B Cell
  "Sox10","Plp1","Mpz"            # Schwann Cell
)
avg_exp <- AverageExpression(
  Merge,
  features = markergenes,
  group.by = "celltype",
  assay = "RNA"
)
plot_matrix <- avg_exp$RNA
plot_matrix <- plot_matrix[markergenes, , drop = FALSE]
plot_matrix <- plot_matrix[, c("Cardiomyocyte", 
                               "Smooth Muscle", 
                               "Pericyte", 
                               "T Cell", 
                               "Macrophage", 
                               "Fibroblast", 
                               "Endothelial", 
                               "B Cell",
                               "Schwann Cell"),
                           drop = FALSE]
gene_groups <- data.frame(
  celltype = rep(
    c("Cardiomyocyte", "Smooth Muscle", "Pericyte", "T Cell", 
      "Macrophage", "Fibroblast", "Endothelial","B Cell", "Schwann Cell"),
    times = c(3, 3, 3, 3, 3, 3, 3,3, 3) 
  )
)
rownames(gene_groups) <- markergenes
celltype_colors <- c(
  Fibroblast    = "#A6CEE3",
  Cardiomyocyte = "#1F78B4",
  Endothelial   = "#B2DF8A",
  Pericyte      = "#33A02C",
  Macrophage    = "#FDBF6F",
  `T Cell`      = "#FF7F00",
  `Smooth Muscle` = "#FB9A99",
  `Schwann Cell`  = "#E31A1C",
  `B Cell`        = "#CAB2D6"
)
anno_colors <- list(celltype = celltype_colors)

ht <- pheatmap(
  t(plot_matrix),            
  scale = "column",           
  cluster_rows = FALSE,        
  cluster_cols = FALSE,         
  annotation_col = gene_groups,       
  annotation_colors = anno_colors,
  gaps_col = cumsum(c(3, 3, 3, 3, 3, 3, 3,3)), 
  color = colorRampPalette(c("#008585", "white", "#e12729"))(100),
  border_color = "#eaeaea",
  fontsize_row = 10,            
  fontsize_col = 8,            
  angle_col = "45",            
  main = "Marker Gene Expression Heatmap (Transposed)",
  annotation_names_col = FALSE, 
  annotation_legend = TRUE,
  show_rownames = TRUE,         
  show_colnames = TRUE         
)


# SFig1.K -----------------------------------------------------------------

conf_mat <- table(merged.object$L5, merged.object$DeconvolutionClass)
conf_mat_prop <- prop.table(conf_mat, margin = 1) 
conf_mat_prop[is.nan(conf_mat_prop)] <- 0
remove_types <- c("Brown Adipocytes", 
                  "White Adipocytes", 
                  "Neurons", 
                  "Mesothelial Cells", 
                  "Fibro-Myocardial Interface",
                  "Proliferating Cells",
                  "T-NK_Cell") 
clean_mat <- conf_mat_prop[!rownames(conf_mat_prop) %in% remove_types, ]
clean_mat <- clean_mat[, colSums(clean_mat) > 0]
plot_data <- t(clean_mat)

print(colnames(plot_data))
print(rownames(plot_data))

ordered_rows <- c(
  "Cardiomyocyte",
  "Fibroblast",
  "EndothelialCell",
  "SmoothMuscleCell",
  "Pericyte",
  "Schwann Cell",
  "B_Cell",
  "Macrophage"
)

ordered_cols <- c(
  # Cardiomyocyte
  "Ventricular Cardiomyocytes", "Atrial Cardiomyocytes", 
  "Stressed Cardiomyocytes", "Conduction System Cardiomyocytes",
  # Fibroblast
  "Ventricular-associated Fibroblasts", "Atrial-associated Fibroblasts", "Resident Fibroblasts",
  "Myofibroblasts", "Perivascular Fibroblasts", "Stress-Response Fibroblasts",
  # EndothelialCell
  "Vascular Endothelial Cells", "Endocardial Endothelial Cells", "Lymphatic Endothelial Cells",
  # SmoothMuscleCell
  "Ventricular-associated SMCs", "Arterial SMCs", "Atrial-associated SMCs",
  # Pericyte
  "Pericytes",
  # Schwann Cell
  "Schwann Cells",
  # Immune
  "B Cells",
  "Fibro-Immune Interface"
)

common_rows <- intersect(ordered_rows, rownames(plot_data))
common_cols <- intersect(ordered_cols, colnames(plot_data))
plot_data_sorted <- plot_data[common_rows, common_cols]

p <- pheatmap::pheatmap(plot_data_sorted, 
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE, 
                        display_numbers = FALSE, 
                        color = colorRampPalette(c("white", "firebrick3"))(100),
                        main = "Cell type verification",
                        angle_col = "45",
                        treeheight_row = 0, 
                        treeheight_col = 0,
                        silent = TRUE          
)

grid.newpage()
pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 1))
grid.draw(p$gtable)






