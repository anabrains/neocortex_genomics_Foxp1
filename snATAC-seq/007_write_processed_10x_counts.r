###################################################################################
###########################  Get 10x Matrix formats from modified rio write10xCounts
# Load necessary libraries
library(Seurat)
library(Signac)
library(tidyverse)
library(rio)
library(DropletUtils)
library(dplyr)
library(Matrix)
library(data.table)
###
save_dir = "./sc-atac/aoa/geo_info/processed_mat/"
source('./sc-atac/aoa/write10xCountsAtac_1.r')
##
p0_atac_exci_m <- readRDS("./sc-atac/aoa/processed_objs/merged_harmony/p0_merged_exci_lin_m.rds")
e16_atac_exci_m <- readRDS("./sc-atac/aoa/processed_objs/merged_harmony/e16_merged_exci_lin_m.rds")
e13_atac_exci_m <- readRDS("./sc-atac/aoa/processed_objs/merged_harmony/e13_merged_exci_lin_m.rds")

## modified rio function
write10xCountsAtac(x = e13_atac_exci_m@assays$peaks@counts, path =paste0(save_dir, "/atac_mat/e13_atac_exci"), version="3", custom_name="e13_exci_atac")
write10xCountsAtac(x = e16_atac_exci_m@assays$peaks@counts, path =paste0(save_dir, "/atac_mat/e16_atac_exci"), version="3", custom_name="e16_exci_atac")
write10xCountsAtac(x = p0_atac_exci_m@assays$peaks@counts, path =paste0(save_dir, "/atac_mat/p0_atac_exci"), version="3", custom_name="p0_exci_atac")
##############################################################################################
###################  UMAP EXCI ATAC
# Access UMAP coordinates from the reductions slot
umap_coordinates <- e13_atac_exci_m@reductions$umap@cell.embeddings
# Save UMAP coordinates to a CSV file
export(umap_coordinates, paste0(save_dir, "umap_coord/atac_standard/e13_exci_atac_umap.csv"),row.names=TRUE) # comma-separated values
##
# Access UMAP coordinates from the reductions slot
umap_coordinates <- e16_atac_exci_m@reductions$umap@cell.embeddings
# Save UMAP coordinates to a CSV file
export(umap_coordinates, paste0(save_dir, "umap_coord/atac_standard/e16_exci_atac_umap.csv"),row.names=TRUE) # comma-separated values
### Access UMAP coordinates from the reductions slot
umap_coordinates <- p0_atac_exci_m@reductions$umap@cell.embeddings
# Save UMAP coordinates to a CSV file
export(umap_coordinates, paste0(save_dir, "umap_coord/atac_standard/p0_exci_atac_umap.csv"),row.names=TRUE) # comma-separated values
##