
library(Seurat)
library(dplyr)
library(ggplot2)

# Fig6.A ------------------------------------------------------------------

subset_object <- readRDS("./subset_object.rds")

Merge1<-readRDS("./scOjb1.rds")
Merge2<-readRDS("./scOjb2.rds")

Merge2_fibro <- subset(Merge2, celltype == "Fibroblasts")
Merge2_fibro$celltype <- "Fibroblast"
Merge2_fibro$condition <- ifelse(Merge2_fibro$group == "ct", "Control", 
                                 ifelse(Merge2_fibro$group == "hfpef", "HFpEF", NA))
Merge <- merge(x = Merge1, y = Merge2_fibro, add.cell.ids = c("M1", "M2"))
Merge <- ifelse(grepl("^M1_", colnames(Merge)), "Dataset_1", "Dataset_2")
Merge <- NormalizeData(Merge, verbose = FALSE)
Merge <- FindVariableFeatures(Merge, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Merge <- ScaleData(Merge, verbose = FALSE)
Merge <- RunPCA(Merge, npcs = 30, verbose = FALSE)
Merge <- RunHarmony(Merge, group.by.vars = "batch", plot_convergence = FALSE)



sobj_sliceB <- subset(subset_object, subset = region == "SliceB")
Idents(sobj_sliceB) <- "L5"

markers_stress_cm <- FindMarkers(
  sobj_sliceB,
  ident.1 = "Stressed Cardiomyocytes",
  only.pos = TRUE,       
  min.pct = 0.25,        
  logfc.threshold = 0.25 
)

top10_stress_cm <- markers_stress_cm %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(10)

markers_perivascular_fib <- FindMarkers(
  sobj_sliceB,
  ident.1 = "Perivascular Fibroblasts",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

sig_genes_stress <- markers_stress_cm %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(50) %>% 
  rownames()

sig_genes_perivascular <- markers_perivascular_fib %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(50) %>% 
  rownames()

sobj_hfpef <- subset(Merge, subset = condition == "HFpEF")

sobj_hfpef_cm <- subset(sobj_hfpef, subset = celltype == "Cardiomyocyte")
sobj_hfpef_cm <- AddModuleScore(
  sobj_hfpef_cm,
  features = list(sig_genes_stress),
  name = "Stress_Score_Robust"
)

sobj_hfpef_fib <- subset(sobj_hfpef, subset = celltype == "Fibroblast")
sobj_hfpef_fib <- AddModuleScore(
  sobj_hfpef_fib,
  features = list(sig_genes_perivascular),
  name = "Perivascular_Score_Robust"
)

print_score_stats <- function(scores, name) {
  cat(paste0("\n================ ", name, " Statistical Overview ================\n"))
  
  cat("--- (Summary) ---\n")
  print(summary(scores))
  
  cat("\n---  (Quantiles) ---\n")
  qs <- quantile(scores, probs = c(0, 0.25, 0.5, 0.75, 0.80, 0.90, 0.95, 0.99, 1))
  print(qs)
  
  cat("\n--- The top 10 cell scores ---\n")
  print(head(sort(scores, decreasing = TRUE), 10))
  
  cat(paste0("\nThreshold reference (Top 10%): > ", round(qs['90%'], 4), "\n"))
  cat(paste0("Threshold reference (Top 5%):  > ", round(qs['95%'], 4), "\n"))
}

top_stress_cells <- colnames(sobj_hfpef_cm)[sobj_hfpef_cm$Stress_Score_Robust1 > quantile(sobj_hfpef_cm$Stress_Score_Robust1, 0.90)]
top_perivascular_cells <- colnames(sobj_hfpef_fib)[sobj_hfpef_fib$Perivascular_Score_Robust1 > quantile(sobj_hfpef_fib$Perivascular_Score_Robust1, 0.90)]

print(paste("Stressed CM:", length(top_stress_cells)))
print(paste("Perivascular Fib:", length(top_perivascular_cells)))

subset_stress <- subset(sobj_hfpef_cm, cells = top_stress_cells)
subset_perivascular <- subset(sobj_hfpef_fib, cells = top_perivascular_cells)

subset_stress$Target_Group <- "Stressed Cardiomyocytes"
subset_perivascular$Target_Group <- "Perivascular Fibroblasts"
target_combined <- merge(x = subset_stress, y = subset_perivascular)
Idents(target_combined) <- "Target_Group"

scObj$celltype_refined <- scObj$celltype
# Stressed cardiomyocytes
scObj$celltype_refined[top_stress_cells] <- "Stressed Cardiomyocyte"

# Perivascular fibroblasts
scObj$celltype_refined[top_Perivascular_cells] <- "Perivascular Fibroblast"

table(scObj$celltype_refined)


DimPlot(
  scObj,
  reduction = "umap",
  group.by = "celltype_refined",
  
  repel = TRUE,
  pt.size = 0.4
)




# Fig6.B ------------------------------------------------------------------

target_combined <- JoinLayers(target_combined)
cellchat <- createCellChat(object = target_combined, group.by = "Target_Group")

CellChatDB <- CellChatDB.mouse 
cellchat@DB <- CellChatDB 

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

df.net <- subsetCommunication(cellchat)
select_pairs <- data.frame(interaction_name = c("THBS2_CD36", "THBS2_CD47"))
source_cells <- "Perivascular Fibroblasts"
target_cells <- "Stressed Cardiomyocytes"

p_bubble <- netVisual_bubble(
  cellchat, 
  sources.use = source_cells, 
  targets.use = target_cells, 
  pairLR.use = select_pairs,      
  remove.isolate = FALSE,
  angle.x = 45
) + ggtitle("Validation: Thbs2 Signaling in scRNA-seq")

print(p_bubble)




# Fig6.C ------------------------------------------------------------------

pathways.show <- c("THBS") 
par(mfrow=c(1,1))
netVisual_aggregate(
  cellchat, 
  signaling = pathways.show,  
  layout = "chord", 
  sources.use = "Perivascular Fibroblasts", 
  targets.use = "Stressed Cardiomyocytes", 
  title.name = "Thbs2 Signaling Flow"
)


