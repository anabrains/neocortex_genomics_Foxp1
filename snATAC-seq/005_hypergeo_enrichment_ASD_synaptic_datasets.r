 library(dplyr)
 setwd("./sc-atac/aoa/gene_lists/edited_heatmap")

 ##load objects 
 #e13_atac #e16_atac #p0_atac
 #below will plot with adj p-value aka FDR in box #overlap_sets
 source('./sc-atac/aoa/hypergeometric_test_edit.R')
 ###
 ndd_synaptic_gene_list_mouse <- read.csv("./sc-atac/aoa/gene_lists/ndd_synaptic_gene_list_mouse.csv")
 # Reordering values 
 ndd_synaptic_gene_list_mouse <- ndd_synaptic_gene_list_mouse %>%   arrange(match(Definition, c("SFARI_ASD_HIGH","Satterstrom_ASD","Fu_ASD", "PreSynaptic","PostSynaptic")))
 ###
 e13_dars <- read.csv("./sc-atac/aoa/gene_lists/dars_mouse/e13_dars.csv")
 e16_dars <- read.csv("./sc-atac/aoa/gene_lists/dars_mouse/e16_dars.csv")
 p0_dars <- read.csv("./sc-atac/aoa/gene_lists/dars_mouse/p0_dars.csv")
 ##
 overlap_sets(e13_dars, ndd_synaptic_gene_list_mouse, dim(e13_atac@assays$RNA)[1], file_prefix = 'e13_dars_newPlot')
 overlap_sets(e16_dars, ndd_synaptic_gene_list_mouse, dim(e16_atac@assays$RNA)[1], file_prefix = 'e16_dars_newPlot')
 overlap_sets(p0_dars, ndd_synaptic_gene_list_mouse, dim(p0_atac@assays$RNA)[1], file_prefix = 'p0_dars_newPlot')
 