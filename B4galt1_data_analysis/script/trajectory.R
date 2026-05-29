library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

setwd("")
x <- combined_SCT

#### b4 ####
cds_b4 <- as.cell_data_set(b4_seurat)
cds_b4

#get cell metadata
colData(cds_b4)

## get gene metadata 
fData(cds_b4)
rownames(fData(cds_b4))[1:10]


##add gene_short_name column 
fData(cds_b4)$gene_shourt_name <- rownames(fData(cds_b4))

## get counts 
counts(cds_b4)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_b4@colData@rownames)))
names(reacreate.partition) <- cds_b4@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_b4@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_b4 <- b4_seurat@active.ident
cds_b4@clusters$UMAP$clusters <- list_cluster_b4

#assign umap coordinate - cell embeddings 
cds_b4@int_colData@listData$reducedDims$UMAP <- b4_seurat@reductions$umap@cell.embeddings

#plot
plot_cells(cds_b4,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_b4 <- learn_graph(cds_b4, use_partition = TRUE)

plot_cells(cds_b4,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_b4 <- order_cells(cds_b4, reduction_method = 'UMAP', root_cells = colnames(cds_b4[, clusters(cds_b4) == 4]))

plot_cells(cds_b4,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_b4)
cds_b4$monocle3_pseudotime <- pseudotime(cds_b4)
data.pseudo <- as.data.frame(colData(cds_b4))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_b4, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

b4_seurat$pseudotime <- pseudotime(cds_b4)

FeaturePlot(b4_seurat, features = 'pseudotime')

b4_umap_coords <- Embeddings(b4_seurat, reduction = "umap")
b4_umap_coords <- as.data.frame(b4_umap_coords)

# Add specific metadata column (e.g., seurat_clusters)
b4_umap_coords$pseudotime <- b4_seurat$pseudotime

write.csv(b4_umap_coords, "b4_pseudotime.csv")


#### ko1 ####

cds_ko1 <- as.cell_data_set(x_ko1)
cds_ko1

#get cell metadata
colData(cds_ko1)

## get gene metadata 
fData(cds_ko1)
rownames(fData(cds_ko1))[1:10]


##add gene_short_name column 
fData(cds_ko1)$gene_shourt_name <- rownames(fData(cds_ko1))

## get counts 
counts(cds_ko1)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_ko1@colData@rownames)))
names(reacreate.partition) <- cds_ko1@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_ko1@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_ko1 <- x_ko1@active.ident
cds_ko1@clusters$UMAP$clusters <- list_cluster_ko1

#assign umap coordinate - cell embeddings 
cds_ko1@int_colData@listData$reducedDims$UMAP <- x_ko1@reductions$umap@cell.embeddings

#plot
plot_cells(cds_ko1,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_ko1 <- learn_graph(cds_ko1, use_partition = TRUE)

plot_cells(cds_ko1,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_ko1 <- order_cells(cds_ko1, reduction_method = 'UMAP', root_cells = colnames(cds_ko1[, clusters(cds_ko1) == 4]))

plot_cells(cds_ko1,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_ko1)
cds_ko1$monocle3_pseudotime <- pseudotime(cds_ko1)
data.pseudo <- as.data.frame(colData(cds_ko1))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_ko1, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_ko1$pseudotime <- pseudotime(cds_ko1)

FeaturePlot(x_ko1, features = 'pseudotime')

ko1_umap_coords <- Embeddings(x_ko1, reduction = "umap")
ko1_umap_coords <- as.data.frame(ko1_umap_coords)
# Add specific metadata column (e.g., seurat_clusters)
ko1_umap_coords$pseudotime <- x_ko1$pseudotime
write.csv(ko1_umap_coords, "ko1_pseudotime.csv")

#### ko3 ####

cds_ko3 <- as.cell_data_set(x_ko3)
cds_ko3

#get cell metadata
colData(cds_ko3)

## get gene metadata 
fData(cds_ko3)
rownames(fData(cds_ko3))[1:10]


##add gene_short_name column 
fData(cds_ko3)$gene_shourt_name <- rownames(fData(cds_ko3))

## get counts 
counts(cds_ko3)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_ko3@colData@rownames)))
names(reacreate.partition) <- cds_ko3@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_ko3@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_ko3 <- x_ko3@active.ident
cds_ko3@clusters$UMAP$clusters <- list_cluster_ko3

#assign umap coordinate - cell embeddings 
cds_ko3@int_colData@listData$reducedDims$UMAP <- x_ko3@reductions$umap@cell.embeddings

#plot
plot_cells(cds_ko3,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_ko3 <- learn_graph(cds_ko3, use_partition = TRUE)

plot_cells(cds_ko3,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_ko3 <- order_cells(cds_ko3, reduction_method = 'UMAP', root_cells = colnames(cds_ko3[, clusters(cds_ko3) == 4]))

plot_cells(cds_ko3,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_ko3)
cds_ko3$monocle3_pseudotime <- pseudotime(cds_ko3)
data.pseudo <- as.data.frame(colData(cds_ko3))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_ko3, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_ko3$pseudotime <- pseudotime(cds_ko3)

FeaturePlot(x_ko3, features = 'pseudotime')

ko3_umap_coords <- Embeddings(x_ko3, reduction = "umap")
ko3_umap_coords <- as.data.frame(ko3_umap_coords)
# Add specific metadata column (e.g., seurat_clusters)
ko3_umap_coords$pseudotime <- x_ko3$pseudotime
write.csv(ko3_umap_coords, "ko3_pseudotime.csv")

#### wt ####

cds_wt <- as.cell_data_set(wt_seurat)
cds_wt

#get cell metadata
colData(cds_wt)

## get gene metadata 
fData(cds_wt)
rownames(fData(cds_wt))[1:10]


##add gene_short_name column 
fData(cds_wt)$gene_shourt_name <- rownames(fData(cds_wt))

## get counts 
counts(cds_wt)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_wt@colData@rownames)))
names(reacreate.partition) <- cds_wt@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_wt@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_wt <- wt_seurat@active.ident
cds_wt@clusters$UMAP$clusters <- list_cluster_wt

#assign umap coordinate - cell embeddings 
cds_wt@int_colData@listData$reducedDims$UMAP <- wt_seurat@reductions$umap@cell.embeddings

#plot
plot_cells(cds_wt,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_wt <- learn_graph(cds_wt, use_partition = TRUE)

plot_cells(cds_wt,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_wt <- order_cells(cds_wt, reduction_method = 'UMAP', root_cells = colnames(cds_wt[, clusters(cds_wt) == 4]))

plot_cells(cds_wt,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_wt)
cds_wt$monocle3_pseudotime <- pseudotime(cds_wt)
data.pseudo <- as.data.frame(colData(cds_wt))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_wt, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

wt_seurat$pseudotime <- pseudotime(cds_wt)

FeaturePlot(wt_seurat, features = 'pseudotime')

wt_umap_coords <- Embeddings(wt_seurat, reduction = "umap")
wt_umap_coords <- as.data.frame(wt_umap_coords)
# Add specific metadata column (e.g., seurat_clusters)
wt_umap_coords$pseudotime <- wt_seurat$pseudotime
write.csv(wt_umap_coords, "wt_pseudotime.csv")

#### wt1 ####

cds_wt1 <- as.cell_data_set(x_wt1)
cds_wt1

#get cell metadata
colData(cds_wt1)

## get gene metadata 
fData(cds_wt1)
rownames(fData(cds_wt1))[1:10]


##add gene_short_name column 
fData(cds_wt1)$gene_shourt_name <- rownames(fData(cds_wt1))

## get counts 
counts(cds_wt1)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_wt1@colData@rownames)))
names(reacreate.partition) <- cds_wt1@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_wt1@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_wt1 <- x_wt1@active.ident
cds_wt1@clusters$UMAP$clusters <- list_cluster_wt1

#assign umap coordinate - cell embeddings 
cds_wt1@int_colData@listData$reducedDims$UMAP <- x_wt1@reductions$umap@cell.embeddings

#plot
plot_cells(cds_wt1,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_wt1 <- learn_graph(cds_wt1, use_partition = TRUE)

plot_cells(cds_wt1,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_wt1 <- order_cells(cds_wt1, reduction_method = 'UMAP', root_cells = colnames(cds_wt1[, clusters(cds_wt1) == 4]))

plot_cells(cds_wt1,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_wt1)
cds_wt1$monocle3_pseudotime <- pseudotime(cds_wt1)
data.pseudo <- as.data.frame(colData(cds_wt1))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_wt1, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_wt1$pseudotime <- pseudotime(cds_wt1)

FeaturePlot(x_wt1, features = 'pseudotime')

wt1_umap_coords <- Embeddings(x_wt1, reduction = "umap")
wt1_umap_coords <- as.data.frame(wt1_umap_coords)
# Add specific metadata column (e.g., seurat_clusters)
wt1_umap_coords$pseudotime <- x_wt1$pseudotime
write.csv(wt1_umap_coords, "wt1_pseudotime.csv")

#### wt2 ####

cds_wt2 <- as.cell_data_set(x_wt2)
cds_wt2

#get cell metadata
colData(cds_wt2)

## get gene metadata 
fData(cds_wt2)
rownames(fData(cds_wt2))[1:10]


##add gene_short_name column 
fData(cds_wt2)$gene_shourt_name <- rownames(fData(cds_wt2))

## get counts 
counts(cds_wt2)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_wt2@colData@rownames)))
names(reacreate.partition) <- cds_wt2@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_wt2@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_wt2 <- x_wt2@active.ident
cds_wt2@clusters$UMAP$clusters <- list_cluster_wt2

#assign umap coordinate - cell embeddings 
cds_wt2@int_colData@listData$reducedDims$UMAP <- x_wt2@reductions$umap@cell.embeddings

#plot
plot_cells(cds_wt2,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_wt2 <- learn_graph(cds_wt2, use_partition = TRUE)

plot_cells(cds_wt2,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_wt2 <- order_cells(cds_wt2, reduction_method = 'UMAP', root_cells = colnames(cds_wt2[, clusters(cds_wt2) == 4]))

plot_cells(cds_wt2,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_wt2)
cds_wt2$monocle3_pseudotime <- pseudotime(cds_wt2)
data.pseudo <- as.data.frame(colData(cds_wt2))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_wt2, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_wt2$pseudotime <- pseudotime(cds_wt2)

FeaturePlot(x_wt2, features = 'pseudotime')

wt2_umap_coords <- Embeddings(x_wt2, reduction = "umap")
wt2_umap_coords <- as.data.frame(wt2_umap_coords)
# Add specific metadata column (e.g., seurat_clusters)
wt2_umap_coords$pseudotime <- x_wt2$pseudotime
write.csv(wt1_umap_coords, "wt2_pseudotime.csv")

#### w4 ####

cds_wt4 <- as.cell_data_set(x_wt4)
cds_wt4

#get cell metadata
colData(cds_wt4)

## get gene metadata 
fData(cds_wt4)
rownames(fData(cds_wt4))[1:10]


##add gene_short_name column 
fData(cds_wt4)$gene_shourt_name <- rownames(fData(cds_wt4))

## get counts 
counts(cds_wt4)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_wt4@colData@rownames)))
names(reacreate.partition) <- cds_wt4@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_wt4@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_wt4 <- x_wt4@active.ident
cds_wt4@clusters$UMAP$clusters <- list_cluster_wt4

#assign umap coordinate - cell embeddings 
cds_wt4@int_colData@listData$reducedDims$UMAP <- x_wt4@reductions$umap@cell.embeddings

#plot
plot_cells(cds_wt4,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_wt4 <- learn_graph(cds_wt4, use_partition = TRUE)

plot_cells(cds_wt4,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_wt4 <- order_cells(cds_wt4, reduction_method = 'UMAP', root_cells = colnames(cds_wt4[, clusters(cds_wt4) == 4]))

plot_cells(cds_wt4,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_wt4)
cds_wt4$monocle3_pseudotime <- pseudotime(cds_wt4)
data.pseudo <- as.data.frame(colData(cds_wt4))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_wt4, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_wt4$pseudotime <- pseudotime(cds_wt4)

FeaturePlot(x_wt4, features = 'pseudotime')

wt4_umap_coords <- Embeddings(x_wt4, reduction = "umap")
wt4_umap_coords <- as.data.frame(wt4_umap_coords)
# Add specific metadata column (e.g., seurat_clusters)
wt4_umap_coords$pseudotime <- x_wt4$pseudotime
write.csv(wt4_umap_coords, "wt4_pseudotime.csv")


# 1. Pseudotime values (vector)
pseudotime <- x_wt4$pseudotime
# 2. Trajectory coordinates (matrix)
trajectory_coords <- reducedDims(cds_wt4)$UMAP
# 3. UMAP coordinates (separate - from Seurat)
umap_coords <- Embeddings(x_wt4, reduction = "umap")
# 4. Metadata (from original object)
metadata <- x_wt4@meta.data
# 5. Combined data frame
trajectory_df <- data.frame(
  cell_id = names(pseudotime),
  pseudotime = pseudotime,
  trajectory_1 = trajectory_coords[,1],
  trajectory_2 = trajectory_coords[,2],
  umap_1 = umap_coords[,1],
  umap_2 = umap_coords[,2],
  cluster = metadata$seurat_clusters
)




#### all_ko ####
cds_all_ko <- as.cell_data_set(x_all_ko)
cds_all_ko

#get cell metadata
colData(cds_all_ko)

## get gene metadata 
fData(cds_all_ko)
rownames(fData(cds_all_ko))[1:10]


##add gene_short_name column 
fData(cds_all_ko)$gene_shourt_name <- rownames(fData(cds_all_ko))

## get counts 
counts(cds_all_ko)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_all_ko@colData@rownames)))
names(reacreate.partition) <- cds_all_ko@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_all_ko@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_all_ko <- x_all_ko@active.ident
cds_all_ko@clusters$UMAP$clusters <- list_cluster_all_ko

#assign umap coordinate - cell embeddings 
cds_all_ko@int_colData@listData$reducedDims$UMAP <- x_all_ko@reductions$umap@cell.embeddings

#plot
plot_cells(cds_all_ko,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_all_ko <- learn_graph(cds_all_ko, use_partition = TRUE)

plot_cells(cds_all_ko,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_all_ko <- order_cells(cds_all_ko, reduction_method = 'UMAP', root_cells = colnames(cds_all_ko[, clusters(cds_all_ko) == 4]))

plot_cells(cds_all_ko,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_all_ko)
cds_all_ko$monocle3_pseudotime <- pseudotime(cds_all_ko)
data.pseudo <- as.data.frame(colData(cds_all_ko))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_all_ko, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_all_ko$pseudotime <- pseudotime(cds_all_ko)

FeaturePlot(x_all_ko, features = 'pseudotime')


#### all_wt ####
cds_all_wt <- as.cell_data_set(x_all_wt)
cds_all_wt

#get cell metadata
colData(cds_all_wt)

## get gene metadata 
fData(cds_all_wt)
rownames(fData(cds_all_wt))[1:10]


##add gene_short_name column 
fData(cds_all_wt)$gene_shourt_name <- rownames(fData(cds_all_wt))

## get counts 
counts(cds_all_wt)


## cluster cells using clustering info from seurat umap 

# assign paritiion 
reacreate.partition <- c(rep(1, length(cds_all_wt@colData@rownames)))
names(reacreate.partition) <- cds_all_wt@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_all_wt@clusters$UMAP$partitions <- reacreate.partition

#assign cluster info 
list_cluster_all_wt <- x_all_wt@active.ident
cds_all_wt@clusters$UMAP$clusters <- list_cluster_all_wt

#assign umap coordinate - cell embeddings 
cds_all_wt@int_colData@listData$reducedDims$UMAP <- x_all_wt@reductions$umap@cell.embeddings

#plot
plot_cells(cds_all_wt,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) + 
  theme(legend.position = 'right')


# learn trajectory graph
cds_all_wt <- learn_graph(cds_all_wt, use_partition = TRUE)

plot_cells(cds_all_wt,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           group_label_size = 5)

# order cells in pseudotime 
cds_all_wt <- order_cells(cds_all_wt, reduction_method = 'UMAP', root_cells = colnames(cds_all_wt[, clusters(cds_all_wt) == 4]))

plot_cells(cds_all_wt,
           color_cells_by = 'pseudotime', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE)

## cells ordered by monocle3 pseudotime 

pseudotime(cds_all_wt)
cds_all_wt$monocle3_pseudotime <- pseudotime(cds_all_wt)
data.pseudo <- as.data.frame(colData(cds_all_wt))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot()

##finding genes that change as a function of pseudotime 
pseudo_test <- graph_test(cds_all_wt, neighbor_graph = 'principal_graph', cores = 8)

# gene_test %>%
#   arrange(q_value) %>%
#   filter(status == 'OK') %>%
#   head()
# FeaturePlot(x, features = c())

# visualizing pseudotime in seurat 

x_all_wt$pseudotime <- pseudotime(cds_all_wt)

FeaturePlot(x_all_wt, features = 'pseudotime')



