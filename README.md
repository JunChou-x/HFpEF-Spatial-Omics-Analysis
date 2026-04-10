# High-resolution spatial transcriptomics reveals fibroblast-cardiomyocyte metabolic coupling in heart failure with preserved ejection fraction

## 📖 Overview
This repository contains the custom code and analytical pipelines for the study: **"High-resolution spatial transcriptomics reveals fibroblast-cardiomyocyte metabolic coupling in heart failure with preserved ejection fraction"**.

Heart failure with preserved ejection fraction (HFpEF) is driven by complex interactions among metabolic stress, fibrosis, and microvascular dysfunction. In this study, we integrated single-cell-informed deconvolution with Visium HD spatial transcriptomics to map the myocardial landscape in a two-hit murine model of early-stage HFpEF. Our analysis identifies spatially organized fibroblast-cardiomyocyte metabolic coupling as a prominent feature of HFpEF and suggests niche-associated signaling programs (such as COMP/THBS1-CD36) as potential therapeutic targets.

## 🗂️ Repository Structure
The analysis scripts in this repository are named according to the corresponding figures in the manuscript (e.g., `Figure1.R`, `SFigure2.R`, etc.). 

Collectively, these scripts cover the following main analytical workflows:
* **Data Preprocessing & Annotation:** Quality control, clustering, and cell type annotation for both scRNA-seq and Visium HD datasets.
* **Spatial Deconvolution:** Mapping single-cell signatures to spatial microenvironments to identify pathological niches.
* **Transcriptomic Remodeling Analysis:** Differential expression and pathway enrichment analysis evaluating metabolic reprogramming (e.g., fatty acid oxidation) in fibroblasts and endothelial cells.
* **Cell-Cell Communication:** Inferring intercellular signaling networks (specifically focusing on the COMP/THBS-CD36 axis) across matched cell populations.

## 💻 System Requirements

### Hardware Requirements
The spatial transcriptomics and deconvolution analysis require a standard computer/server with at least 64GB RAM.

### Software Dependencies
The analyses were primarily performed using R and Python. Major packages include:

**R (>= 4.2.0)**
* Seurat (>= 4.3.0 or 5.0.0)
* CellChat
* clusterProfiler
* (Add your spatial deconvolution package here, e.g., RCTD)

**Python (>= 3.9)**
* Scanpy
* pySCENIC

## 🚀 Usage
1. Clone this repository:
   ```bash
   git clone [https://github.com/YourUsername/HFpEF-Spatial-Metabolism.git](https://github.com/YourUsername/HFpEF-Spatial-Metabolism.git)
   cd HFpEF-Spatial-Metabolism
