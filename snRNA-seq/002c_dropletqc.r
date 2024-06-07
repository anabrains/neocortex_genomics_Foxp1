#rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DropletQC)
## Install packages
#devtools::install_github("powellgenomicslab/DropletQC")

# Allows us to save objects with a new name defined by a string
saveit <- function(..., string, file) {
  x <- list(...)
  names(x) <- string
  save(list=names(x), file=file, envir=list2env(x))
}
## LOOP INVARIANTS

## specify gtf path file
gtfPath = "./ao/refdata-gex-mm10-2020-A/genes/genes.gtf"
base_dir= "./ao/"
#the seurat object with doublets removed
seu_path <- paste0(base_dir, "data_pipeline1/doublets_removed_seu_objs")
#where to save the plots of nucFrac plots
results_folder <- paste0(base_dir, "data_pipeline1/nucFrac/plots/")
# obj nuc frac folder 
nucFracData = "./ao/data_pipeline1/nucFrac/"
##what are the folder names? specify from data_pipeline1 where the ones for this timepoint are at #mine are AO(a number)_GENOTYPE_AGE_BATCH
folder_names <- list.files(path = paste0(base_dir, "data_pipeline1"), pattern = "AO.*_[CTRL|CKO]", full.names=F)


for (i in 1:length(folder_names)){
  print(folder_names[i])
  countMatDir <- paste0(base_dir, "data_pipeline1/", folder_names[i], "/cellbender");
  filenames <- list.files(path = countMatDir, pattern = "\\.h5$", full.names = T);
  
  mouse_name <- unlist(strsplit(folder_names[i], "_"))[1]
  
  for (j in 1:length(filenames)){
    print(filenames[j])
    countMat <- Read10X_h5(filenames[j])
    #where are your bams stored? add info to the base dir
    bam_path <- paste0(base_dir, "scrna/", folder_names[i], "/outs/possorted_genome_bam.bam")
    
    seu_obj_list <- list.files(path = seu_path, pattern = paste0(mouse_name, ".*\\.rdata$"), full.names = T);
    
    for (k in 1:length(seu_obj_list)) {
      print(seu_obj_list[k])
      seu_obj = load(seu_obj_list[k]);
      
      ## Calculate nuclear fraction (ratio of reads mapping to introns)
      nucFrac = nuclear_fraction_annotation(
        annotation_path = gtfPath,
        bam = bam_path,
        barcodes = colnames(countMat),
        cores = 22, verbose = T)
      
      ## Load the meta data (contains cell type annotation, UMI count etc.).
      meta = get(seu_obj)@meta.data
      
      # Add nuclear fraction into the meta data
      meta$NuclearFraction = nucFrac[rownames(meta),]
      print("[X] Saving PDF :) ")
      #setwd("./ao/data_pipeline1/nucFrac/plots/")
      ## Plot nuclear fraction per cell type ##
      #pdf(paste0(mouse_name,'NuclearFraction_Per_Cluster.pdf'))
      p1=ggboxplot(meta, x = 'seurat_clusters', y = 'NuclearFraction', color = 'seurat_clusters') + rotate_x_text(90)
      nucFrac_plot_filename = paste("NuclearFraction_Per_Cluster_", mouse_name, ".pdf", sep="")
      ggsave(file.path(results_folder, nucFrac_plot_filename), plot=p1, width=12, height=8, units="in", dpi=300)
      #dev.off()

      print(paste0("[X] PDF Saved for Mouse ", mouse_name))
      
      # solution: Missing value where true/false needed
      #x = NA
      #if(is.na(x)) {x=FALSE} else {if(x) {x}} 
      
      setwd("./ao/data_pipeline1/nucFrac/")
      
      ## save in obj meta data
      seu_obj$nucFract <- c()
      nucFrac1 <- as.vector(meta$NuclearFraction)
      seu_obj$nucFract <- nucFrac1
      #
      saveRDS(seu_obj, file = paste0(nucFracData,mouse_name,'_nucFrac.rds'))
      ###
      
      
    }
  }
}

##rm low quality 

# #ao17
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO17_seu_obj_doublet_removed.rdata")
# ao17= subset(AO17_no_db, idents = c('24'), invert=T)
# save(ao17, file='ao17_seu_obj.rdata')
# #ao15
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO15_seu_obj_doublet_removed.rdata")
# ao15= subset(AO15_no_db, idents = c('16'), invert=T)
# save(ao15, file='ao15_seu_obj.rdata')
# #ao14
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO14_seu_obj_og.rdata")
# ao14= subset(AO14_og, idents = c('13'), invert=T)
# save(ao14, file='ao14_seu_obj.rdata')
# #ao8
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO8_seu_obj_doublet_removed.rdata")
# ao8= subset(AO8_no_db, idents = c('0','3','5','13'), invert=T)
# save(ao8, file='ao8_seu_obj.rdata')
# #
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO7_seu_obj_doublet_removed.rdata")
# ao7= subset(AO7_no_db, idents = c('0','1','19'), invert=T)
# save(ao7, file='ao7_seu_obj.rdata')
# #
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO5_seu_obj_doublet_removed.rdata")
# ao5= subset(AO5_no_db, idents = c('0','1','19'), invert=T)
# save(ao5, file='ao5_seu_obj.rdata')
# #
# load("./ao/data_pipeline1/doublets_removed_seu_objs/AO1_seu_obj_doublet_removed.rdata")
# ao1= subset(AO1_no_db, idents = c('0','3','18'), invert=T)
# save(ao1, file='ao1_seu_obj.rdata')
