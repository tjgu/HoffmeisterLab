library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(SeuratDisk)
library(dplyr)
library(SingleR)
library(celldex)
library(gprofiler2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(loomR)


setwd("/scratch/g/lcashdol/Hoffmeister_pf4_multiome/")
### create annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "mm10"
mt_genes <- c("mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Cytb", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6")

### load data
## WT3620 ##
wt3620.counts <- Read10X_h5(filename = "./wt3620/filtered_feature_bc_matrix.h5")
wt3620_fragpath <- "./wt3620/atac_fragments.tsv.gz"
non_mt_genes <- setdiff(rownames(wt3620.counts$`Gene Expression`), mt_genes)
length(rownames(wt3620.counts$`Gene Expression`))
wt3620.counts$`Gene Expression` <- wt3620.counts$`Gene Expression`[non_mt_genes,]
length(rownames(wt3620.counts$`Gene Expression`))

wt3620 <- CreateSeuratObject(counts = wt3620.counts$`Gene Expression`, assay = "RNA", project = "wt3620")
wt3620$condition <- "wt"
# create ATAC assay and add it to the object
wt3620[["ATAC"]] <- CreateChromatinAssay(
  counts = wt3620.counts$Peaks,
  sep = c(":", "-"),
  fragments = wt3620_fragpath,
  annotation = annotation,
  genome = 'mm10'
)

## VAV3702 ##
vav3702.counts <- Read10X_h5(filename = "./vavko_3702/filtered_feature_bc_matrix.h5")
vav3702_fragpath <- "./vavko_3702/atac_fragments.tsv.gz"
length(rownames(vav3702.counts$`Gene Expression`))
non_mt_genes <- setdiff(rownames(vav3702.counts$`Gene Expression`), mt_genes)
vav3702.counts$`Gene Expression` <- vav3702.counts$`Gene Expression`[non_mt_genes,]
length(rownames(vav3702.counts$`Gene Expression`))

vav3702 <- CreateSeuratObject(counts = vav3702.counts$`Gene Expression`, assay = "RNA", project = "vav3702")
vav3702$condition <- "ko"
# create ATAC assay and add it to the object
vav3702[["ATAC"]] <- CreateChromatinAssay(
  counts = vav3702.counts$Peaks,
  sep = c(":", "-"),
  fragments = vav3702_fragpath,
  annotation = annotation
)

## merge 
combined <- merge(wt3620, vav3702)
# rm(wt3620.counts, vav3702.counts)

## find commone peaks 
peaks <- reduce(unlist(as(c(wt3620@assays$ATAC@ranges,
                            vav3702@assays$ATAC@ranges),
                          "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]

counts_atac_merged <- FeatureMatrix(combined@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(combined))
combined[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                                         fragments = combined@assays$ATAC@fragments,
                                         annotation = combined@assays$ATAC@annotation,
                                         sep = c(":","-"),
                                         genome = 'mm10')
rm(counts_atac_merged)

# QC
combined <- PercentageFeatureSet(combined, pattern = "^mt-", col.name = "percent.mt", assay = "RNA")
combined <- NucleosomeSignal(combined, assay = "ATAC")
combined <- TSSEnrichment(combined, assay = "ATAC")

VlnPlot(combined,
        features = c("nFeature_RNA",
                     "percent.mt",
                     "nFeature_ATAC",
                     "TSS.enrichment",
                     "nucleosome_signal"),
        ncol = 5, assay = "RNA")

combined <- subset(combined,
                 subset = nFeature_RNA > 10 &
                   nFeature_RNA < 3000 &
                   nFeature_ATAC > 1000 &
                   nFeature_ATAC < 30000 &
                   TSS.enrichment > 2.5 &
                   nucleosome_signal < 1.5
)
## add cell cycle score 
combined <- SCTransform(combined)
s.genes <- gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
combined$CC.Difference <- combined$S.Score - combined$G2M.Score
## regress cell cycle 
combined <- SCTransform(combined, vars.to.regress = c("CC.Difference"))
combined <- RunPCA(combined)
ElbowPlot(combined)
combined <- RunUMAP(combined, dims = 1:10)
combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.2)
DimPlot(combined, label = TRUE)
DimPlot(combined, label = TRUE, split.by = "orig.ident")
DimPlot(combined, label = TRUE, group.by = "orig.ident")

## integration (anchors)
object_list <- SplitObject(combined, split.by = "orig.ident")
object_list <- lapply(object_list, function(obj) {
  obj <- SCTransform(obj, vars.to.regress = c("CC.Difference"))
  obj <- FindVariableFeatures(obj, selection.method = "sct")
})
int_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = int_features)

# find the anchors for integration
int_anchors <- FindIntegrationAnchors(object.list = object_list, anchor.features = int_features, 
                                      normalization.method = "SCT")

# remove large objects from environment
rm(object_list)

# perform integration
combine_integrated <- IntegrateData(int_anchors, normalization.method = "SCT")
DefaultAssay(combine_integrated) <- "integrated"
combine_integrated <- FindVariableFeatures(combine_integrated)
combine_integrated <- RunPCA(combine_integrated, features = VariableFeatures(combine_integrated))
ElbowPlot(combine_integrated)
combine_integrated <- RunUMAP(combine_integrated, dims = 1:5)
combine_integrated <- FindNeighbors(combine_integrated, dims = 1:5)
combine_integrated <- FindClusters(combine_integrated, resolution = 0.3)
DimPlot(combine_integrated, label = TRUE)
DimPlot(combine_integrated, label = TRUE, split.by = "orig.ident")
DimPlot(combine_integrated, label = T, group.by = "Phase")

p1 <- DimPlot(combined, label = F, split.by = "orig.ident") & NoAxes() & ggtitle("scRNAseq before integration")
p2 <- DimPlot(combine_integrated, label = F, split.by = "orig.ident") & NoAxes() & ggtitle("scRNAseq after integration")
p1 + p2 & theme(plot.title = element_text(hjust = 0.5, size = 20), 
                legend.text = element_text(size = 20),strip.text = element_text(size = 20))

### ATAC 
DefaultAssay(combine_integrated) <- "ATAC"
combine_integrated <- FindTopFeatures(combine_integrated, min.cutoff = 50)
combine_integrated <- RunTFIDF(combine_integrated, method = 1)
combine_integrated <- RunSVD(combine_integrated, n = 50)
p1 <- ElbowPlot(combine_integrated, reduction="lsi")
p2 <- DepthCor(combine_integrated, n = 30)
p1 | p2
combine_integrated <- RunUMAP(combine_integrated,
                              reduction = "lsi",
                              dims = 2:5,
                              reduction.name = "umap_atac",
                              reduction.key = "UMAPATAC_")
combine_integrated <- FindNeighbors(object = combine_integrated, reduction = 'lsi', dims = 2:5)
combine_integrated <- FindClusters(object = combine_integrated, verbose = FALSE, algorithm = 3, resolution = 0.1)
DimPlot(combine_integrated,
        split.by = "orig.ident",
        reduction = "umap_atac") & NoAxes()

combine_integrated <- FindMultiModalNeighbors(combine_integrated, reduction.list = list("pca", "lsi"), dims.list = list(1:5, 2:5))
combine_integrated <- RunUMAP(combine_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
combine_integrated <- FindClusters(combine_integrated, graph.name = "wsnn", algorithm = 3, resolution = 0.2)
DimPlot(combine_integrated, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE)
DimPlot(combine_integrated, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)
DimPlot(combine_integrated, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE)

p1 <- DimPlot(combine_integrated, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE, group.by = "orig.ident") + ggtitle("RNA") & NoAxes()
p2 <- DimPlot(combine_integrated, reduction = "umap_atac", label = TRUE, label.size = 2.5, repel = TRUE, group.by = "orig.ident") + ggtitle("ATAC")& NoAxes()
p3 <- DimPlot(combine_integrated, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE, group.by = "orig.ident") + ggtitle("WNN")& NoAxes()
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))

### all Markers
DefaultAssay(combine_integrated) <- "SCT"
combine_integrated <- PrepSCTFindMarkers(combine_integrated)

features <- c("Itga2b","Pf4", "Gata1", "Vwf", "Fli1", "Mki67", "Gp9", "Gata2")
FeaturePlot(combine_integrated, features = "B4galt1", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Itga2b", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Pf4", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Gata1", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Vwf", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Fli1", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Mki67", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Gp9", split.by = "condition", reduction = "wnn.umap")
FeaturePlot(combine_integrated, features = "Gata2", split.by = "condition", reduction = "wnn.umap")

#DGE across all clusters between wt and ko

all_gene_results <- data.frame()
for(cluster in c(0,1,2,3,4,5,6,7)){
  this_fm <- FindMarkers(object = subset(combine_integrated, idents = cluster), group.by = "condition",
                         ident.1 = "ko", ident.2 = "wt", recorrect_umi = F)
  this_fm$gene <- rownames(this_fm)
  this_fm$cluster <- as.character(cluster)
  all_gene_results <- rbind(all_gene_results, this_fm)
}
sig_gene_results <- all_gene_results[all_gene_results$p_val_adj<0.05,]

###GSEA across all clusters between wt and ko
all_kegg_results <- data.frame()
for(cluster in c(0,1,2,3,4,5,6,7)){
  this_fm <- FindMarkers(object = subset(combine_integrated, idents = cluster), group.by = "condition",
                         ident.1 = "ko", ident.2 = "wt", logfc.threshold = 0, recorrect_umi = F)
  this_lfc <- this_fm$avg_log2FC
  names(this_lfc) <- mapIds(x = org.Mm.eg.db, keys = rownames(this_fm), keytype = "SYMBOL", column = "ENTREZID")
  this_lfc <- this_lfc[order(this_lfc, decreasing = TRUE)]
  kegg_res <- gseKEGG(this_lfc, organism = "mmu", nPerm = 10000, pvalueCutoff = 1)
  kegg_res <- setReadable(kegg_res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kegg_res_df <- kegg_res@result
  kegg_res_df$cluster <- cluster
  all_kegg_results <- rbind(all_kegg_results, kegg_res_df)
}

all_reactome_results <- data.frame()
for(cluster in c(0,1,2,3,4,5,6,7)){
  this_fm <- FindMarkers(object = subset(combine_integrated, idents = cluster), group.by = "condition",
                         ident.1 = "ko", ident.2 = "wt", logfc.threshold = 0, recorrect_umi = F)
  this_lfc <- this_fm$avg_log2FC
  names(this_lfc) <- mapIds(x = org.Mm.eg.db, keys = rownames(this_fm), keytype = "SYMBOL", column = "ENTREZID")
  this_lfc <- this_lfc[order(this_lfc, decreasing = TRUE)]
  react_res <- gsePathway(this_lfc, organism = "mouse", nPerm = 10000, pvalueCutoff = 1)
  react_res <- setReadable(react_res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  react_res_df <- react_res@result
  react_res_df$cluster <- cluster
  all_reactome_results <- rbind(all_reactome_results, react_res_df)
}


DefaultAssay(combine_integrated) <- "ATAC"
all_peak_results <- data.frame()
for(cluster in c(0,1,2,3,4,5,6,7)){
  this_fm <- FindMarkers(object = subset(combine_integrated, idents = cluster), group.by = "condition",
                         ident.1 = "ko", ident.2 = "wt")
  this_fm$gene <- rownames(this_fm)
  this_fm$cluster <- as.character(cluster)
  all_peak_results <- rbind(all_peak_results, this_fm)
}
sig_peak_results <- all_peak_results[all_peak_results$p_val_adj<0.05,]

splitted <- strsplit(sig_peak_results$gene, "-")
sig_peak_results$chr <- sapply(splitted, "[[", 1)
sig_peak_results$start <- sapply(splitted, "[[", 2)
sig_peak_results$end <- sapply(splitted, "[[", 3)
WriteXLS::WriteXLS(x=sig_peak_results, "Sig_peak_vav3702_wt3620_removeMT_results.xlsx", SheetNames = "peak", 
                   AdjWidth = T, row.names = T, col.names = T, FreezeRow = 1)

splitted <- strsplit(all_peak_results$gene, "-")
all_peak_results$chr <- sapply(splitted, "[[", 1)
all_peak_results$start <- sapply(splitted, "[[", 2)
all_peak_results$end <- sapply(splitted, "[[", 3)
chr16_peak <- all_peak_results[all_peak_results$chr == "chr16", ]
WriteXLS::WriteXLS(x=chr16_peak, "chr16_peak_vav3702_wt3620_removeMT_results.xlsx", SheetNames = "peak", 
                   AdjWidth = T, row.names = T, col.names = T, FreezeRow = 1)


WriteXLS::WriteXLS(x = list(sig_gene_results, all_gene_results, all_kegg_results, all_reactome_results), 
                   ExcelFileName = "all_vav3702_wt3620_removeMT_results.xlsx", 
                   SheetNames = c("sig_gene_results", "all_genes_results", "all_kegg_results", "all_reactome_results"), 
                   AdjWidth = T, FreezeRow = 1, row.names = T, col.names = T)

CoveragePlot(
  object = combine_integrated,
  region = "Itga2b",
  features = "Itga2b",
  annotation = T,
  peaks = T,
  group.by = "orig.ident",
  tile = TRUE,
  links = TRUE, expression.assay = "SCT"
)
CoveragePlot(
  object = combine_integrated,
  region = "B4galt1",
  features = "B4galt1",
  annotation = T,
  peaks = T,
  group.by = "condition",
  tile = TRUE,
  links = TRUE, expression.assay = "SCT"
)
CoveragePlot(
  object = combine_integrated,
  region = "Muc13",
  features = "Muc13",
  annotation = T,
  peaks = T,
  group.by = "condition",
  tile = TRUE,
  links = TRUE, expression.assay = "SCT"
)

CoveragePlot(
  object = combine_integrated,
  region = "Pbx1",
  features = "Pbx1",
  annotation = T,
  peaks = T,
  group.by = "condition",
  tile = TRUE,
  links = TRUE, expression.assay = "SCT"
)
## motif analysis
library(TFBSTools)
library(JASPAR2020)
# add motif
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

combine_integrated <- AddMotifs(combine_integrated, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

# find motif
open_peaks <- AccessiblePeaks(combine_integrated)
meta.feature = na.omit(combine_integrated[['ATAC']]@meta.features[open_peaks, ])
query.feature = na.omit(combine_integrated[['ATAC']]@meta.features[sig_peak_results$gene, ])
peaks_matched <- MatchRegionStats(meta.feature = meta.feature,
                                  query.feature = query.feature,
                                  n = 50000)

motif_enrichment_mpp2 <- FindMotifs(combine_integrated,
                                    features = sig_peak_results$gene[sig_peak_results$cluster == "5"],
                                    background = peaks_matched) 
motif_enrichment_lt <- FindMotifs(combine_integrated,
                                    features = sig_peak_results$gene[sig_peak_results$cluster == "0"|sig_peak_results$cluster == "1"],
                                    background = peaks_matched) 
all_motif_enrichment <- data.frame()
for(cluster in c(0:(length(unique(combine_integrated@active.ident)) - 1))){
  motif_res <- FindMotifs(combine_integrated,
                          features = sig_peak_results$gene[sig_peak_results$cluster == cluster],
                          background = peaks_matched)
  motif_res$cluster <- cluster
  all_motif_enrichment <- rbind(all_motif_enrichment, motif_res)
}

sig_motif_enrichment <- all_motif_enrichment[all_motif_enrichment$p.adjust < 0.05, ]

WriteXLS::WriteXLS(x=list(sig_motif_enrichment, all_motif_enrichment), 
                   "all_motif_removeMT_results.xlsx", SheetNames = c("sig_motif", "all_motifg"), 
                   AdjWidth = T, row.names = T, col.names = T, FreezeRow = 1)

MotifPlot(combine_integrated, motifs = c("MA0140.2", "MA0766.2", "MA0037.3", "MA0482.2", "MA0036.3", "MA0035.4", "MA1104.2"), ncol=4)

WriteXLS::WriteXLS(x=motif_enrichment_mpp2, "motif_mpp2_vav_wt_results.xlsx", SheetNames = "motif", AdjWidth = T, row.names = T, col.names = T, FreezeRow = 1)
WriteXLS::WriteXLS(x=motif_enrichment_lt, "motif_lt_vav_wt_results.xlsx", SheetNames = "motif", AdjWidth = T, row.names = T, col.names = T, FreezeRow = 1)

## Trajectory analysis
library(SeuratWrappers)
library(monocle3)
library(Matrix)

DefaultAssay(combine_integrated) <- "SCT"

#### trajectory ####
combine_integrated_cds <- as.cell_data_set(combine_integrated)
combine_integrated_cds

#get cell metadata
colData(combine_integrated_cds)

## get gene metadata 
fData(combine_integrated_cds)
rownames(fData(combine_integrated_cds))[1:10]

##add gene_short_name column 
fData(combine_integrated_cds)$gene_shourt_name <- rownames(fData(combine_integrated_cds))

## get counts 
counts(combine_integrated_cds)

## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(combine_integrated_cds@colData@rownames)))
names(reacreate.partition) <- combine_integrated_cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

combine_integrated_cds@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster <- combine_integrated@active.ident
combine_integrated_cds@clusters$UMAP$clusters <- list_cluster

#assign umap coordinate - cell embeddings 
# combine_integrated_cds@int_colData@listData$reducedDims$UMAP <- combine_integrated@reductions$umap@cell.embeddings
combine_integrated_cds@int_colData@listData$reducedDims$UMAP <- combine_integrated@reductions$wnn.umap@cell.embeddings

#plot
plot_cells(combine_integrated_cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
combine_integrated_cds <- learn_graph(combine_integrated_cds, use_partition = TRUE)

plot_cells(combine_integrated_cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
combine_integrated_cds <- order_cells(combine_integrated_cds, reduction_method = 'UMAP', 
                                      root_cells = colnames(combine_integrated_cds[, clusters(combine_integrated_cds) == 0]))

plot_cells(combine_integrated_cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = T,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(combine_integrated_cds)
combine_integrated_cds$monocle3_pseudotime <- pseudotime(combine_integrated_cds)
data.pseudo <- as.data.frame(colData(combine_integrated_cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(combine_integrated_cds, neighbor_graph = 'principal_graph', cores = 8)

# visualizing pseudotime in seurat 

combine_integrated$pseudotime <- pseudotime(combine_integrated_cds)

FeaturePlot(combine_integrated, features = 'pseudotime', reduction = "wnn.umap", label = T)

##wt
DefaultAssay(combine_integrated) <- "SCT"
wt <- subset(combine_integrated, subset = orig.ident == "wt3620")

#### trajectory ####
wt_cds <- as.cell_data_set(wt)
wt_cds

#get cell metadata
colData(wt_cds)

## get gene metadata 
fData(wt_cds)
rownames(fData(wt_cds))[1:10]


##add gene_short_name column 
fData(wt_cds)$gene_shourt_name <- rownames(fData(wt_cds))

## get counts 
counts(wt_cds)

## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(wt_cds@colData@rownames)))
names(reacreate.partition) <- wt_cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

wt_cds@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster <- wt@active.ident
wt_cds@clusters$UMAP$clusters <- list_cluster

#assign umap coordinate - cell embeddings 
# wt_cds@int_colData@listData$reducedDims$UMAP <- wt@reductions$umap@cell.embeddings
wt_cds@int_colData@listData$reducedDims$UMAP <- wt@reductions$wnn.umap@cell.embeddings


#plot
plot_cells(wt_cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
wt_cds <- learn_graph(wt_cds, use_partition = TRUE)

plot_cells(wt_cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
wt_cds <- order_cells(wt_cds, reduction_method = 'UMAP', 
                                      root_cells = colnames(wt_cds[, clusters(wt_cds) == 0]))

plot_cells(wt_cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = T,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(wt_cds)
wt_cds$monocle3_pseudotime <- pseudotime(wt_cds)
data.pseudo <- as.data.frame(colData(wt_cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(wt_cds, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

wt$pseudotime <- pseudotime(wt_cds)

FeaturePlot(wt, features = 'pseudotime', reduction = "wnn.umap", label = T)

##ko
DefaultAssay(combine_integrated) <- "SCT"
ko <- subset(combine_integrated, subset = orig.ident == "vav3702")

#### trajectory ####
ko_cds <- as.cell_data_set(ko)
ko_cds

#get cell metadata
colData(ko_cds)

## get gene metadata 
fData(ko_cds)
rownames(fData(ko_cds))[1:10]


##add gene_short_name column 
fData(ko_cds)$gene_shourt_name <- rownames(fData(ko_cds))

## get counts 
counts(ko_cds)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(ko_cds@colData@rownames)))
names(reacreate.partition) <- ko_cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

ko_cds@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster <- ko@active.ident
ko_cds@clusters$UMAP$clusters <- list_cluster

#assign umap coordinate - cell embeddings 
# wt_cds@int_colData@listData$reducedDims$UMAP <- wt@reductions$umap@cell.embeddings
ko_cds@int_colData@listData$reducedDims$UMAP <- ko@reductions$wnn.umap@cell.embeddings


#plot
plot_cells(ko_cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
ko_cds <- learn_graph(ko_cds, use_partition = TRUE)

plot_cells(ko_cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
ko_cds <- order_cells(ko_cds, reduction_method = 'UMAP', 
                      root_cells = colnames(ko_cds[, clusters(ko_cds) == 0]))

plot_cells(ko_cds,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = T,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(ko_cds)
ko_cds$monocle3_pseudotime <- pseudotime(ko_cds)
data.pseudo <- as.data.frame(colData(ko_cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(ko_cds, neighbor_graph = 'principal_graph', cores = 8)

# visualizing pseudotime in seurat 

ko$pseudotime <- pseudotime(ko_cds)

FeaturePlot(ko, features = 'pseudotime', reduction = "wnn.umap", label = T)
