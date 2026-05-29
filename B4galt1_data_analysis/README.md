# B4galt1 Data Analysis
## B4galt1 single-cell analysis (LSK)

This repo provides **analysis code** for the B4galt1 LSK single‑cell RNA‑seq workflow. **Raw or processed sequencing data are not included**—obtain inputs separately and set paths in each script or Rmd. Main workflow: `B4galt1_publication.Rmd` includes **transfer learning** (ML) for cell type classification. **Multiome** analysis uses `pf4_multiome.R` (10x RNA+ATAC with Seurat/Signac, WNN integration, Pf4/megakaryocyte‑focused markers, and KO vs WT differential analyses; requires your own `filtered_feature_bc_matrix.h5` and ATAC fragment files). Optional **scRNA** add‑ons (main Rmd object, not multiome): `mk_velocyto.R` (RNA velocity) and `trajectory.R` (Monocle3).

## 1) System requirements
- OS: macOS (tested), Linux (expected compatible)
- CPU: Standard desktop/laptop; no GPU required
- RAM: 16+ GB (core); 32+ GB recommended for scRNA analysis
- RAM: 360 GB (12 cores) recommended for scMultiome analysis
- Disk: 200+ GB free

Tested software versions:
- R 4.2–4.3
- Key R packages: Seurat 4.3.0, sctransform 0.4.0, SeuratWrappers, SingleR, celldex, scater, tricycle, cowplot, patchwork, ggplot2, ggrepel, ggpubr, aplot, gridExtra, dplyr, tidyr, data.table, cluster, diptest, **xgboost** (for transfer learning ML), pivottabler, openxlsx, clusterProfiler, org.Mm.eg.db, msigdbr, velocyto.R (optional), monocle3 (optional).

No non‑standard hardware required.

## 2) Installation guide
1. Install R (4.3 recommended) and optionally RStudio.
2. Install packages:
```r
install.packages(c(
  "Seurat","SeuratWrappers","sctransform","SingleR","celldex","scater",
  "tricycle","cowplot","patchwork","ggplot2","ggrepel","ggpubr","aplot",
  "gridExtra","dplyr","tidyr","data.table","cluster","diptest","xgboost",
  "pivottabler","openxlsx","clusterProfiler","msigdbr"
))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
# Optional
install.packages("velocyto.R")
BiocManager::install("monocle3")
```
Typical install time: 10–30 min on a normal desktop.

## 3) Running the analyses
Point the main Rmd at your 10x‑style inputs (e.g. under `filt_copy/`, one folder per sample, plus processed matrices if your setup uses `PROCESSED_DATA/`). Adjust every path in the Rmd to match your layout.

Run main analysis:
```r
rmarkdown::render("B4galt1_publication.Rmd")
```
Outputs are written into an auto‑named folder (starts with the sample list and `mito_20_seurat4.3.0_sct0.4.0_ccGenes/`) containing PDFs and RDS files (e.g., `combined_SCT_finalClusters.Rds`). The analysis includes **transfer learning** using XGBoost for cell type classification (LTHSC, STHSC, MPP2, MPP3, MPP4).

Expected runtime: ~30–60 min for the main Rmd on a normal desktop.

**scRNA** RNA velocity (optional): adjust paths inside `mk_velocyto.R` to your loom files and Seurat object from the main scRNA pipeline, then:
```r
source("mk_velocyto.R")
```
Expected runtime: ~30–90 min.

**scRNA** trajectory (optional): requires objects from the velocity step above, then:
```r
source("trajectory.R")
```
Expected runtime: ~20–60 min.

**Multiome** (10x scRNA+scATAC): edit **hardcoded** paths at the top of `pf4_multiome.R` (`setwd`, sample folders with `filtered_feature_bc_matrix.h5` and ATAC `fragments.tsv.gz`). Install extra packages as needed: **Signac**, **BSgenome.Mmusculus.UCSC.mm10**, **EnsDb.Mmusculus.v79** (`BiocManager::install(...)` and CRAN). Then:
```r
source("pf4_multiome.R")
```
The script merges peaks, runs RNA QC and **SCTransform** (with cell‑cycle regression), **Seurat** RNA integration, **Signac** TF‑IDF/LSI, **WNN** UMAP, marker/feature plots (e.g. Pf4, megakaryocyte/platelet genes), per‑cluster **FindMarkers** (RNA and peaks), and **GSEA** (KEGG/Reactome).

Expected runtime: highly dependent on multiome library size; plan for more wall time and RAM than the scRNA‑only Rmd. ~4-8 h.

## 4) Instructions for use
- Provide your own 10x outputs (barcodes/features/matrix or equivalent) and place them under a directory like `filt_copy/`, one subfolder per sample (or change paths in the Rmd).
- In `B4galt1_publication.Rmd`, update:
  - `samples` and `phenotypes`
  - Path in `input_10x(...)`
  - `mito_cutoff` as desired
- Knit the Rmd to generate integrated objects and figures.

Notes:
- If adding/removing samples (e.g., removed `LSK_KO_2`), ensure any hardcoded lines referencing that sample are removed or replaced by loops.
- Update file paths in `mk_velocyto.R` and `trajectory.R` (scRNA) and in `pf4_multiome.R` (multiome) to your environment before running.

## Reproduction instructions
With **your own data** in place and paths updated:
- Knit `B4galt1_publication.Rmd` to regenerate main figures and RDS outputs under the output folder.
- Run `pf4_multiome.R` after path adjustments to reproduce **multiome** results (paired RNA+ATAC).
- Optionally run `mk_velocyto.R` and `trajectory.R` after path adjustments to reproduce **scRNA** velocity and trajectory outputs (not used for the multiome script).
- Session info files (e.g., `sessionInfo.txt`, `finalsessionInfo.txt`) are produced to capture versions.

## Contact
Questions about running or reproducing results: contact Hoffmeister lab maintainers.
