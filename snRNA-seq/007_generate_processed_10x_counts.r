#############################################################
## Need to go from counts in a Seurat object to the 10X format?
# Recently found that some tools like Scrublet (doublet detection), require scRNA-seq counts to be in 10X format.
# * barcodes.tsv
# * genes.tsv
# * matrix.mtx
######  Generate features cells mat for geo / to share using rio
library(Seurat)
library(tidyverse)
library(rio)
library(DropletUtils)
library(dplyr)
#save_Dir
save_dir = "./sc-atac/aoa/geo_info/processed_mat/"
## objects
e13_rna_exci <- readRDS("./ao/data_pipeline1/mergee13/exci_lin/e13_exci_lin_celltypes_obj1.rds")
e16_rna_exci <- readRDS("./ao/data_pipeline2/merge/exci_lin/e16_exci_lin_celltypes_obj1.rds")
p0_rna_exci <- readRDS("./ao/data_pipeline1/mergep0/exci_lin/p0_exci_lin_celltypes_obj1.rds")
# Output data #10x formats 
write10xCounts(x = e13_rna_exci@assays$RNA@counts, path =paste0(save_dir, "./e13_exci"), version="3")
write10xCounts(x = p0_rna_exci@assays$RNA@counts, path =paste0(save_dir, "./p0_exci"), version="3")
write10xCounts(x = e16_rna_exci@assays$RNA@counts, path =paste0(save_dir, "./e16_exci"), version="3")
#### other data to share
###################  UMAP EXCI
# Access UMAP coordinates from the reductions slot
umap_coordinates <- e13_rna_exci@reductions$umap@cell.embeddings
# Save UMAP coordinates to a CSV file
export(umap_coordinates, paste0(save_dir, "umap_coord/e13_exci_umap.csv"),row.names=TRUE) # comma-separated values
### Access UMAP coordinates from the reductions slot
umap_coordinates <- e16_rna_exci@reductions$umap@cell.embeddings
# Save UMAP coordinates to a CSV file
export(umap_coordinates, paste0(save_dir, "umap_coord/e16_exci_umap.csv"),row.names=TRUE) # comma-separated values
### Access UMAP coordinates from the reductions slot
umap_coordinates <- p0_rna_exci@reductions$umap@cell.embeddings
# Save UMAP coordinates to a CSV file
export(umap_coordinates, paste0(save_dir, "umap_coord/p0_exci_umap.csv"),row.names=TRUE) # comma-separated values
