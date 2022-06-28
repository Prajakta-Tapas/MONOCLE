library(monocle3)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
library(gdata)
library(readr)

matrix = readMM('matrix.mtx')
cell_metadata <- read.table('barcodes.tsv' , header = F)
gene_metadata = read.table('features.tsv' , header = F)


cds <- new_cell_data_set(matrix,
                         cell_metadata,
                         gene_metadata)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds , reduction_method = "UMAP")

#
plot_cells(cds)

plot_cells(cds, color_cells_by="scMAGIC")

cds <- reduce_dimension(cds, reduction_method="tSNE")

plot_cells(cds, reduction_method="tSNE", color_cells_by="scMAGIC")

plot_cells(cds, color_cells_by="scMAGIC", label_cell_groups=FALSE)

cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

cds_subset <- choose_cells(cds)

pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)

plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

cds_subset = cluster_cells(cds_subset, resolution=1e-2)
plot_cells(cds_subset, color_cells_by="cluster")

colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])
## I have changed 1st 2 numbers only
colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                        "1"="Oligodendrocyte progenitor cell_AdultTemporalLobe",
                                                        "2"="Oligodendrocyte progenitor cell_AdultTemporalLobe",
                                                        "3"="Vulval precursors",
                                                        "4"="Sex myoblasts",
                                                        "5"="Sex myoblasts",
                                                        "6"="Vulval precursors",
                                                        "7"="Failed QC",
                                                        "8"="Vulval precursors",
                                                        "10"="Unclassified neurons",
                                                        "11"="Distal tip cells")
plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")
