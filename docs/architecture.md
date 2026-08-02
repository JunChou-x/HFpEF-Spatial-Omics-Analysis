# HFpEF Analysis Architecture

## Runtime model

This repository is a stateful research notebook expressed as scripts. The
important boundary is the R session: preparation scripts create objects and
figure scripts consume those objects by name.

```text
Inputs on disk
  |
  +-- scRNA_object.R
  |       `-- scObj.rds
  |
  +-- spatial_object.R
  |       `-- merged.object (in memory)
  |
  `-- external RDS objects
          `-- Figure6.R

merged.object + spatial data frames from Figure1.R
  |
  +-- Figure2.R -> cardiomyocyte analyses
  +-- Figure3.R -> fibroblast analyses + CellChat
  +-- Figure4.R -> endothelial/mural analyses + CellChat
  +-- Figure5.R -> spatial niches + NicheNet
  `-- SFigure1.R ... SFigure5.R

SFigure2.R/SFigure3.R -> H5AD -> external pySCENIC -> Loom
                                      |
                                      `-- python_pyscenic_*.py
```

## Object contracts

| Producer | Object or file | Consumers |
| --- | --- | --- |
| `scRNA_object.R` | `scObj.rds` and object `Merge` | downstream analyses that need the scRNA reference |
| `spatial_object.R` | `merged.object` with sample/cell-type metadata | `Figure1.R` and all main/supplementary figure scripts |
| `Figure1.R` | `CU_DF`, `CL_DF_rescaled`, `HU_DF_rescaled`, `HL_DF_rescaled`, palette objects | Figures 2-5 and supplementary figures |
| `SFigure2.R` | `pysenic_cm.h5ad` | external pySCENIC and `python_pyscenic_cm.py` |
| `SFigure3.R` | `pysenic_fb.h5ad` | external pySCENIC and `python_pyscenic_fb.py` |
| `Figure4.R` | `pysenic_ec.h5ad` | external pySCENIC and `python_pyscenic_ec.py` |
| external preprocessing | `subset_object.rds`, `scOjb1.rds`, `scOjb2.rds` | `Figure6.R` |

## Path conventions

All scripts assume the repository root is the working directory. Relative
paths point to one of three places:

- raw 10x/Visium directories (`Control1`, `D1-5-3shang/outs`, and similar),
- intermediate RDS/H5AD/Loom files at the repository root, or
- `data/spatial_data` used by `spatial_object.R`.

There is no shared configuration layer. Before running a figure, check its
`readRDS`, `read.csv`, `read_parquet`, and `read_h5ad` calls and confirm that
the expected files exist.

## Dependency groups

- Single-cell/spatial core: Seurat, SeuratObject, BPCells, Arrow, Harmony,
  DoubletFinder, and future.
- Figure statistics: tidyverse, clusterProfiler, org.Mm.eg.db, tradeSeq,
  monocle3, slingshot, CellChat, NicheNet, and plotting helpers.
- Regulatory summaries: Scanpy, AnnData, loompy, NumPy, pandas, SciPy,
  Matplotlib, seaborn, and `adjustText`.

The repository does not pin versions. Treat the original analysis environment
as part of the reproducibility record.

## Maintenance guidance

- Keep figure numbering aligned with the manuscript.
- Document any new shared object in this file before adding another figure
  script that depends on it.
- Prefer writing generated objects under a dedicated results directory in new
  code; preserve the legacy root-level paths until all consumers are migrated.
- Do not commit raw sequencing data, private RDS objects, or generated PDFs.
