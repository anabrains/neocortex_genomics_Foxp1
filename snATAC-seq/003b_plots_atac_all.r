suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(require(dplyr))
library(scCustomize)
library(gridExtra)
library(magrittr)
library(SingleCellExperiment)
library(cellAlign)
library(ggpubr)
library(destiny)
library(ggthemes)
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(GenomeInfoDb)
library(patchwork)
library(biovizBase) #BiocManager::install("biovizBase")
library(EnsDb.Mmusculus.v79) #bioCManager install if not found
set.seed(1234)
##
p0_all <- readRDS("./sc-atac/aoa/processed_objs/p0_all_merged_celltypes.rds")
e13_all <- readRDS("./sc-atac/aoa/processed_objs/e13_all_merged.rds")
e16_all <- readRDS("./sc-atac/aoa/processed_objs/e16_all_merged.rds")
### lots of colors
#"RG"='#8B93FF', 'IP'='#cb6bee',  'L2-4'='#8792fb', 'L5-6'='#DF67A0'
my_cols <- c("RG"='#D461C2', 'IP'='#8C28E6', "MigNeurons"='#A1EFE7', "MigCells.1"='#A1EFE7',
             'L2-4'='#FFBD33','L2-4.1'='#FFBD33', 'L2-4.2'='#FFBD33', 'L2-4.3'='#FFBD33','L2-4.4'='#FFBD33', 'L2-4.5'='#FFBD33', 'L2-4.6'='#FFBD33','L2-4.7'='#FFBD33', 'L2-4.8'='#FFBD33', 'L2-4.9'='#FFBD33','L2-4.10'='#FFBD33',
             'L4'='#D28BF6','L4.1'='#D28BF6',
             'DL CPN'='#E270CC', 'DL CPN.2'='#E270CC',
             'L5-6'='#74CAE6','L5-6.0'='#74CAE6', 'L5-6.1'='#74CAE6', 'L5-6.2'='#74CAE6','L5-6.3'='#74CAE6',
             'L5-6.4'='#74CAE6','L5-6.5'='#74CAE6','L5-6.6'='#74CAE6','L5-6.7'='#74CAE6','L5-6.8'='#74CAE6', 'L6'='#74CAE6',
             'L1_CR'='#6bbdc1',"L1"='#6bbdc1',
             'Choroid Plexus'='#DB797A','ChoroidPlexus'='#DB797A', 'Endothelial'='#B9BA5F','Endo'='#B9BA5F', 'Pericytes'='#b1f1cc', 'Microglia'='#ce851f',
             'Oligodendrocytes'='#A4DFF2', 'OPCs'='#A4DFF2',
             "G.E." = "#496989",
             'Int'='#6bbd80', 'Int.1'='#6bbd80', 'Int.2'='#6bbd80', 'Int.3'='#6bbd80','Int.4'='#6bbd80', 'Int.5'='#6bbd80', 'Int.6'='#6bbd80','Int.7'='#6bbd80',
             'Stri.Inhib.'='#6cbead','Stri.Inhib'='#6cbead',
             'Int / Thalm'='#559766',"Thalm/Int"='#559766',
             'VLMC'="#ba9623", "cortical hem"="#627254",
             'Nonspec'='#bcb8b1','Nonspec.1'='#bcb8b1','Nonspec.2'='#bcb8b1','Nonspec.3'='#bcb8b1','Nonspec.4'='#bcb8b1',
             'CyclingGliaCells'='#6a994e',
             'Astrocytes'='#1a936f','Astrocytes1'='#1a936f',
             'Oligodendrocytes'='#619b8a')

######################  E13 plots
ao_atac = e13_all
ao_atac$CellType1 = factor(ao_atac$CellType1,levels=c("ChoroidPlexus","Pericytes", "VLMC", "VLMC.1","Microglia","cortical hem","Endothelial",
                                  "G.E.","G.E.1","G.E.2", "Int.3", "Int.5","Int.6","Int.7","Int.8",
                                  "Nonspec", "Stri.Inhib","Stri.Inhib.1",  "Thalm/Int", "Thalm/Int.2","Thalm/Int.3",
                                  "L1", "RG.2","RG.3","RG.4","RG.5","RG.6","RG.8","RG.9", 
                                  "IP.1","IP.2", "MigNeurons","MigNeurons.1", "L5-6.0","L5-6.1","L5-6.2","L5-6.3","L5-6.4"))
DimPlot(ao_atac, group.by = 'CellType1', label = T,label.size = 2)
#now make this new ordered identity the active identity so it saves
ao_atac <- SetIdent(ao_atac, value = "CellType1"); table(ao_atac$CellType1);  table(Idents(ao_atac)) #do they match
###
setwd("./sc-atac/aoa/plots/qc_plots")
Stacked_VlnPlot(ao_atac, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = 0, group.by = 'CellType1_coll', x_lab_rotate = T, ggplot_default_colors = F, colors_use = my_cols, plot_spacing = 0.3)

pdf("aoe13_vlnPlot_QCs.pdf", width=9, height=6.2);
Stacked_VlnPlot(ao_atac, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = 0, group.by = 'CellType1_coll', x_lab_rotate = T, ggplot_default_colors = F, colors_use = my_cols, plot_spacing = 0.3)
dev.off()

plt1 <- ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$CellType1_coll, fill=genotype)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.1));
pdf("aoe13_all_stack_barplot_genotype.pdf", width=8.2, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$CellType1_coll, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.1));
pdf("aoe13_all_stack_barplot_orig-ident.pdf", width=8.2, height=4);
plt1
dev.off();
plt1 <- DimPlot(ao_atac, group.by = 'CellType1_coll', pt.size = 3, cols = my_cols, label = T, label.size = 5, repel = T)
pdf("aoe13_all_dimplot_label.pdf", width=14, height=10);
plt1
dev.off();
###################################
########################### E16
##get atac obj all cells and re-order identities for plotting and collapsed for ease of plotting
ao_atac = e16_all
ao_atac$CellType1 = factor(ao_atac$CellType1,levels=c("Choroid Plexus", "Endo", "Pericytes","Microglia","OPCs", "Nonspec", "Stri.Inhib.", "Int / Thalm","Int", "Int.1", "Int.2","Int.3",
                                                      "RG.1","RG.2", "RG.3", "RG.4","IP.0","IP.1","IP.2",                 
                                                      "Migrating neurons", "Migrating neurons.1","Migrating neurons.2", "Migrating neurons.3",                 
                                                      "L2-4.1", "L2-4.2","L2-4.3" ,"L4-6",
                                                      "L5_6.1","L5_6.2","L5-6","L5-6.0","L5-6.1", "L5-6.2", "L5-6.3","L5-6.4","L5-6.5","L5-6.6","L5-6.7"))

                                                      
pdf("aoe16_vlnPlot_QCs.pdf", width=9, height=6.2);
Stacked_VlnPlot(ao_atac, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = 0, group.by = 'CellType1_coll', x_lab_rotate = T, ggplot_default_colors = F, colors_use = my_cols, plot_spacing = 0.3)
dev.off()

plt1 <- ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$CellType1_coll, fill=genotype)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.1));
pdf("aoe16_all_stack_barplot_genotype.pdf", width=8.2, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$CellType1_coll, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.1));
pdf("aoe16_all_stack_barplot_orig-ident.pdf", width=8.2, height=4);
plt1
dev.off();
plt1 <- DimPlot(ao_atac, group.by = 'CellType1_coll', pt.size = 3, cols = my_cols, label = T, label.size = 5, repel = T)
pdf("aoe16_all_dimplot_label.pdf", width=14, height=10);
plt1
dev.off();
###############################
ao_atac = p0_all
ao_atac <- SetIdent(ao_atac, value = "CellType1")
library(scCustomize)
new_idents_vec = c("VLMC","Endothelial",'Pericytes',
                   "Nonspec", "Nonspec", "Nonspec", "Nonspec","Microglia",
                   "Int", "Int", "Int", "Int", "Int", "Int", "Int", "CyclingGliaCells",
                   "Astrocytes", "Astrocytes","Oligodendrocytes",  "MigCells.1",
                   "L2-4", "L2-4", "L2-4", "L2-4", "L2-4", "L2-4", "L2-4", "L2-4", "L2-4","L2-4", "L2-4", "L2-4",
                   "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6", "L5-6")
ao_atac = Rename_Clusters(ao_atac, new_idents_vec, meta_col_name = 'CellType1')
ao_atac[["CellType1_coll"]] <- Idents(obj = ao_atac)

plt1 <- ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$CellType1_coll, fill=genotype)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =55, vjust = 0.5));
pdf("aop0_all_stack_barplot_genotype.pdf", width=8.2, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$CellType1_coll, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
pdf("aop0_all_stack_barplot_orig-ident.pdf", width=8.2, height=4);
plt1
dev.off();

plt1 <- DimPlot(ao_atac, group.by = 'CellType1_coll', pt.size = 3, cols = my_cols, label = T, label.size = 5, repel = T)
pdf("aop0_all_dimplot_label.pdf", width=14, height=10);
plt1
dev.off();