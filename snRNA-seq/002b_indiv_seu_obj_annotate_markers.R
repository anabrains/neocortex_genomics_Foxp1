## reads .rdata seurat objects 
## creates a umap plot, then finds all markers ##compares to dev ctx datasets
## hypergro enrichment code

library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(cowplot)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(svglite) #for saving .svg
library(future);
plan(strategy='multicore',workers=22);
#22.67 GiB needed for this; input in bytes
options(future.globals.maxSize = 30000 * 1024^2);
#makes the umap
gen_cluster_plot <- function(seu_obj, curr_mouse_id, results_folder) {
  curr_plot = DimPlot(seu_obj, label = T, label.size = 6, pt.size = 2)
  curr_plot = curr_plot + plot_annotation(paste('Clusters of', curr_mouse_id))
  plot_filename = paste("UMAP_Plot_", curr_mouse_id, ".pdf", sep="")
  ggsave(file.path(results_folder, plot_filename), plot = curr_plot, width=12, height=8, units= "in", dpi=300)
}

#files_i_need
# values aka files we need
zylka_p0=read.table(file="./zylka/zylka-p0-cellBarcodes/zylka_p0_markers_BarcodeAnnotated_forEnrichment.csv",sep=",", header = TRUE)
zylka_markers = read.table('./ao/eao/zylka_e14+p0_markers.csv', sep = ",")
zylka_e= read.table(file="./zylka/zylkaE14_matched/ZYLKA_E14_MARKERS_LIST.csv",sep=",", header = TRUE)
##dibella dev mouse ctx
da_e13= read.table(file="./ao/eao/arlotta_dev_mouse/e13/arlotta_e13_devCtx_makers_for_enrich.csv",sep=",", header = TRUE)
library(readr);
da_dev_ctx <- read_csv("./ao/eao/arlotta_dev_mouse/arlotta_developing_Ctx_makers_for_enrich.csv");
###arlotta p0
da_p1 <- read_csv("./ao/eao/arlotta_dev_mouse/p1/arlotta_p1_devCtx_makers_for_enrich.csv");
##
da_e16 <- read_csv("./ao/eao/arlotta_dev_mouse/e16/arlotta_e16_devCtx_makers_for_enrich.csv")
#
source('./ao/eao/hypergeometric_test.R');

# Allows us to save objects with a new name defined by a string # https://stackoverflow.com/a/50665893

saveit <- function(..., string, file) {
  x <- list(...)
  names(x) <- string
  save(list=names(x), file=file, envir=list2env(x))
}

#data folder is where the .rdata objs live
data_folder=("./ao/data_pipeline1/doublets_removed_seu_objs")
#where to save output .csvs and hypergeo comparisons
results_folder = ("./ao/data_pipeline1/all_markers")
#where to save the umap plot
results_plot = ("./ao/data_pipeline1/plots")

filenames <- list.files(path = data_folder, pattern = "\\.rdata$", full.names = T);
for (i in 1:length(filenames)){
  curr_name <- load(filenames[i])
  ###https://stackoverflow.com/questions/9083907/how-to-call-an-object-with-the-character-variable-of-the-same-name
  ##get lets you get the R object by the sanme character version of its name
  curr <- get(curr_name)
  #create umap plot
  gen_cluster_plot(seu_obj = curr, curr_name, results_plot)
  
  seu_obj.markers <- FindAllMarkers(object = curr, only.pos=TRUE, min.pct = 0.25, logfc.threshold= 0.25);
  seu_obj.markers <- filter(seu_obj.markers, p_val_adj <= 0.05);
  #paste(unlist(strsplit(noquote(as.character(AO1_no_db@meta.data$orig.ident)), '_'))[[1]], "all_mrks", sep="_");
  save_name <- paste(unlist(strsplit(noquote(curr_name), '_'))[[1]], "all_mrks", sep="_");
  print(paste0("[X] Found Markers for ", print(deparse(substitute(curr))))); 
  saveit(seu_obj_sub=seu_obj.markers, string=save_name, file=file.path(results_folder,paste0(save_name,".rdata")));
  setwd(results_folder)
  write.csv(seu_obj.markers,file=paste(
    save_name,"_all_cluster_markers.csv", 
    sep =""), row.names=F);
  
  print(paste0("[X] Subsetting Markers ", print(deparse(substitute(curr))))); 
  
  seu_obj_tbl <- seu_obj.markers[,6:7];
  seu_obj_tbl <- seu_obj_tbl[c(2,1)];
  colnames(seu_obj_tbl) = c("Gene","DEFINITION");
  
  save_name <- paste(unlist(strsplit(noquote(curr_name), '_'))[[1]], "no-dubs", sep="_");
  saveit(object=seu_obj_tbl, string=save_name, file=file.path(results_folder,paste0(save_name,".rdata")));
  
  write.csv(seu_obj_tbl,file=paste(
    save_name,"_cluster_gene_markers.csv", 
    sep =""), row.names=F);
  
  
  print(paste0("[X] Saved both markers ", print(deparse(substitute(curr))))); 
  
  overlap_sets(table1 = zylka_markers, table2 = seu_obj_tbl, background_n=length(rownames(curr)), file_prefix = paste(save_name,"hypergeo_zylkaAll",sep = "_"))
  #overlap_sets(table1 = zylka_e, table2 = seu_obj_tbl, background_n=length(rownames(curr)), file_prefix = paste(save_name,"hypergeo_zylka_e14",sep = "_"))
  overlap_sets(table1 = zylka_p0, table2 = seu_obj_tbl, background_n=length(rownames(curr)), file_prefix = paste(save_name,"hypergeo_zylka_p0",sep = "_"))
  
  overlap_sets(table1 = da_dev_ctx, table2 = seu_obj_tbl, background_n=length(rownames(curr)), file_prefix = paste(save_name,"hypergeo_daAll",sep = "_"))
  #overlap_sets(table1 = da_e13, table2 = seu_obj_tbl, background_n=length(rownames(curr)), file_prefix = paste(save_name,"hypergeo_da_e13",sep = "_"))
  overlap_sets(table1 = da_p1, table2 = seu_obj_tbl, background_n=length(rownames(curr)), file_prefix = paste(save_name,"hypergeo_da_p0",sep = "_"))
  
  print(paste0("[X] Hypergeometric enrichment done for ", print(deparse(substitute(curr))))); 
  
  #end
  #save_name <- paste(unlist(strsplit(noquote(curr_name), '_'))[[1]]);

}

