 
library(dplyr)

 read_dir = "./sc-atac/aoa/gene_lists/degs_mouse/"
 #below will plot with adj p-value aka FDR in box #overlap_sets
 source('./hypergeometric_test_edit.R')
 ##use loaded exci rna objs 
 ##e13_rna #e16_rna #p0_rna
 ndd_filt_ppsyn <- read.csv("./sc-atac/aoa/gene_lists/ndd_filt_ppsyn.csv")
 ndd_synaptic_gene_list_mouse <- ndd_synaptic_gene_list_mouse %>%   arrange(match(Definition, c("SFARI_ASD_HIGH","Satterstrom_ASD","Fu_ASD", "PreSynaptic","PostSynaptic")))

 e13_degs <- read.csv(paste0(read_dir,"e13_degs_logfc15_renamed.csv"))
 e16_degs <- read.csv(paste0(read_dir,"e16_degs_logfc15_renamed.csv"))
 p0_degs <- read.csv(paste0(read_dir,"p0_degs_logfc15_renamed.csv"))

#hypergeometric enrichment
 overlap_sets(e13_degs, ndd_synaptic_gene_list_mouse, background_n = dim(e13_rna)[1], file_prefix = 'e13_degs_newPlot')
 overlap_sets(e16_degs, ndd_synaptic_gene_list_mouse, background_n = dim(e16_rna)[1], file_prefix = 'e16_degs_newPlot')
 overlap_sets(p0_degs, ndd_synaptic_gene_list_mouse, background_n = dim(p0_rna)[1], file_prefix = 'p0_degs_newPlot')
 
