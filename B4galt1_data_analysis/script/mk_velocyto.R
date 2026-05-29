library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(stringr)
library(monocle3)

setwd("")
x <- combined_SCT

## b4
b4_seurat <- subset(x, subset = orig.ident == "LSK_KO_8")

b4.loom <- ReadVelocity(file = "possorted_genome_bam_I7379.loom")
b4 <- as.Seurat(x = b4.loom); rm(b4.loom)

## convert barcode 
old_barcode <- rownames(b4@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_3")

b4 <- RenameCells(b4, new.names = new_barode)
b4 <- AddMetaData(b4, b4_seurat@meta.data)
b4 <- b4[, rownames(b4@meta.data[!is.na(b4@meta.data$seurat_clusters),])]


##add umap and pca 
b4@reductions$umap <- b4_seurat@reductions$umap
b4@reductions$pca <- b4_seurat@reductions$pca
##run velocyto
b4 <- RunVelocity(object = b4, deltaT = 1, kCells = 25, fit.quantile = 0.02)
b4 <- SetIdent(b4, value = b4@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = b4)))
names(x = ident.colors) <- levels(x = b4)
cell.colors <- ident.colors[Idents(object = b4)]
names(x = cell.colors) <- colnames(x = b4)
show.velocity.on.embedding.cor(emb = Embeddings(object = b4, reduction = "umap"), 
                               vel = Tool(object = b4, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))


## ko1
x_ko1 <- subset(x, subset = orig.ident == "LSK_KO_1")
ko1.loom <- ReadVelocity(file = "possorted_genome_bam_U8IP3.loom")
ko1 <- as.Seurat(x = ko1.loom)
# ko <- SCTransform(object = ko, assay = "spliced")
rm(ko1.loom)
# ko <- RunPCA(object = ko, verbose = FALSE)

old_barcode <- rownames(ko1@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_1")

ko1 <- RenameCells(ko1, new.names = new_barode)
ko1 <- AddMetaData(ko1, x_ko1@meta.data)
ko1 <- ko1[, rownames(ko1@meta.data[!is.na(ko1@meta.data$seurat_clusters),])]


##add umap and pca 
ko1@reductions$umap <- x_ko1@reductions$umap
ko1@reductions$pca <- x_ko1@reductions$pca
##run velocyto
ko1 <- RunVelocity(object = ko1, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ko1 <- SetIdent(ko1, value = ko1@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = ko1)))
names(x = ident.colors) <- levels(x = ko1)
cell.colors <- ident.colors[Idents(object = ko1)]
names(x = cell.colors) <- colnames(x = ko1)
show.velocity.on.embedding.cor(emb = Embeddings(object = ko1, reduction = "umap"), 
                               vel = Tool(object = ko1, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))

## ko3
x_ko3 <- subset(x, subset = orig.ident == "LSK_KO_3")
ko3.loom <- ReadVelocity(file = "possorted_genome_bam_Q7H8E.loom")
ko3 <- as.Seurat(x = ko3.loom)
# ko <- SCTransform(object = ko, assay = "spliced")
rm(ko3.loom)
# ko <- RunPCA(object = ko, verbose = FALSE)

old_barcode <- rownames(ko3@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_2")

ko3 <- RenameCells(ko3, new.names = new_barode)
ko3 <- AddMetaData(ko3, x_ko3@meta.data)
ko3 <- ko3[, rownames(ko3@meta.data[!is.na(ko3@meta.data$seurat_clusters),])]


##add umap and pca 
ko3@reductions$umap <- x_ko3@reductions$umap
ko3@reductions$pca <- x_ko3@reductions$pca
##run velocyto
ko3 <- RunVelocity(object = ko3, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ko3 <- SetIdent(ko3, value = ko3@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = ko3)))
names(x = ident.colors) <- levels(x = ko3)
cell.colors <- ident.colors[Idents(object = ko3)]
names(x = cell.colors) <- colnames(x = ko3)
show.velocity.on.embedding.cor(emb = Embeddings(object = ko3, reduction = "umap"), 
                               vel = Tool(object = ko3, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))

## wt
wt_seurat <- subset(x, subset = orig.ident == "LSK_WT_8")
wt.loom <- ReadVelocity(file = "possorted_genome_bam_33UY5.loom")
wt <- as.Seurat(x = wt.loom)
# ko <- SCTransform(object = ko, assay = "spliced")
rm(wt.loom)
# ko <- RunPCA(object = ko, verbose = FALSE)

old_barcode <- rownames(wt@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_7")

wt <- RenameCells(wt, new.names = new_barode)
wt <- AddMetaData(wt, wt_seurat@meta.data)
wt <- wt[, rownames(wt@meta.data[!is.na(wt@meta.data$seurat_clusters),])]


##add umap and pca 
wt@reductions$umap <- wt_seurat@reductions$umap
wt@reductions$pca <- wt_seurat@reductions$pca
##run velocyto
wt <- RunVelocity(object = wt, deltaT = 1, kCells = 25, fit.quantile = 0.02)
wt <- SetIdent(wt, value = wt@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = wt)))
names(x = ident.colors) <- levels(x = wt)
cell.colors <- ident.colors[Idents(object = wt)]
names(x = cell.colors) <- colnames(x = wt)
show.velocity.on.embedding.cor(emb = Embeddings(object = wt, reduction = "umap"), 
                               vel = Tool(object = wt, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))

## lsk_wt1
x_wt1 <- subset(x, subset = orig.ident == "LSK_WT_1")
wt1.loom <- ReadVelocity(file = "possorted_genome_bam_F5ILL.loom")
wt1 <- as.Seurat(x = wt1.loom)
# ko <- SCTransform(object = ko, assay = "spliced")
rm(wt1.loom)
# ko <- RunPCA(object = ko, verbose = FALSE)

old_barcode <- rownames(wt1@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_4")

wt1 <- RenameCells(wt1, new.names = new_barode)
wt1 <- AddMetaData(wt1, x_wt1@meta.data)
wt1 <- wt1[, rownames(wt1@meta.data[!is.na(wt1@meta.data$seurat_clusters),])]


##add umap and pca 
wt1@reductions$umap <- x_wt1@reductions$umap
wt1@reductions$pca <- x_wt1@reductions$pca
##run velocyto
wt1 <- RunVelocity(object = wt1, deltaT = 1, kCells = 25, fit.quantile = 0.02)
wt1 <- SetIdent(wt1, value = wt1@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = wt1)))
names(x = ident.colors) <- levels(x = wt1)
cell.colors <- ident.colors[Idents(object = wt1)]
names(x = cell.colors) <- colnames(x = wt1)
show.velocity.on.embedding.cor(emb = Embeddings(object = wt1, reduction = "umap"), 
                               vel = Tool(object = wt1, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))

## lsk_wt2
x_wt2 <- subset(x, subset = orig.ident == "LSK_WT_2")
wt2.loom <- ReadVelocity(file = "possorted_genome_bam_82NDA.loom")
wt2 <- as.Seurat(x = wt2.loom)
# ko <- SCTransform(object = ko, assay = "spliced")
rm(wt2.loom)
# ko <- RunPCA(object = ko, verbose = FALSE)

old_barcode <- rownames(wt2@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_5")

wt2 <- RenameCells(wt2, new.names = new_barode)
wt2 <- AddMetaData(wt2, x_wt2@meta.data)
wt2 <- wt2[, rownames(wt2@meta.data[!is.na(wt2@meta.data$seurat_clusters),])]


##add umap and pca 
wt2@reductions$umap <- x_wt2@reductions$umap
wt2@reductions$pca <- x_wt2@reductions$pca
##run velocyto
wt2 <- RunVelocity(object = wt2, deltaT = 1, kCells = 25, fit.quantile = 0.02)
wt2 <- SetIdent(wt2, value = wt2@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = wt2)))
names(x = ident.colors) <- levels(x = wt2)
cell.colors <- ident.colors[Idents(object = wt2)]
names(x = cell.colors) <- colnames(x = wt2)
show.velocity.on.embedding.cor(emb = Embeddings(object = wt2, reduction = "umap"), 
                               vel = Tool(object = wt2, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))

## lsk_wt4
x_wt4 <- subset(x, subset = orig.ident == "LSK_WT_4")
wt4.loom <- ReadVelocity(file = "possorted_genome_bam_F79DF.loom")
wt4 <- as.Seurat(x = wt4.loom)
# ko <- SCTransform(object = ko, assay = "spliced")
rm(wt4.loom)
# ko <- RunPCA(object = ko, verbose = FALSE)

old_barcode <- rownames(wt4@meta.data)
splited <- strsplit(old_barcode, split = ":")
preffix <- sapply(splited, "[[", 2)
nucleotides_only <- str_sub(sapply(splited, "[[", 2), 1,16) 
new_barode <- paste0(nucleotides_only, "-1", "_6")

wt4 <- RenameCells(wt4, new.names = new_barode)
wt4 <- AddMetaData(wt4, x_wt4@meta.data)
wt4 <- wt4[, rownames(wt4@meta.data[!is.na(wt4@meta.data$seurat_clusters),])]


##add umap and pca 
wt4@reductions$umap <- x_wt4@reductions$umap
wt4@reductions$pca <- x_wt4@reductions$pca
##run velocyto
wt4 <- RunVelocity(object = wt4, deltaT = 1, kCells = 25, fit.quantile = 0.02)
wt4 <- SetIdent(wt4, value = wt4@meta.data$seurat_clusters)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = wt4)))
names(x = ident.colors) <- levels(x = wt4)
cell.colors <- ident.colors[Idents(object = wt4)]
names(x = cell.colors) <- colnames(x = wt4)
show.velocity.on.embedding.cor(emb = Embeddings(object = wt4, reduction = "umap"), 
                               vel = Tool(object = wt4, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, 
                               xlim = c(-15, 15), ylim = c(-15,15))
