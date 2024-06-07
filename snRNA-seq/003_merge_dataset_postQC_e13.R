##merge datasets after cellBender, doublets removed, nucFrac cleaning
library(Seurat)
library(harmony)
library(ggplot2)
library(ggpubr)
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
plan(strategy='multicore',workers=26);
#22.67 GiB needed for this; input in bytes
options(future.globals.maxSize = 30000 * 1024^2);

#appropriate e13 objs #automate later  #test now
load("./ao/data_pipeline1/doublets_removed_seu_objs/AO09_seu_obj_og.rdata")
load("./ao/data_pipeline1/doublets_removed_seu_objs/AO10_seu_obj_og.rdata")
load("./ao/data_pipeline1/doublets_removed_seu_objs/AO13_seu_obj_og.rdata")
load("./ao/data_pipeline1/doublets_removed_seu_objs/AO16_seu_obj_doublet_removed.rdata")
load("./ao/data_pipeline1/nucFrac/ao14_seu_obj.rdata")
load("./ao/data_pipeline1/nucFrac/ao15_seu_obj.rdata")
load("./ao/data_pipeline1/nucFrac/ao17_seu_obj.rdata")
load("./ao/data_pipeline1/nucFrac/ao18_seu_obj.rdata")

aoe13 <- merge(AO09_og, y = c(ao14, ao15, AO16_no_db, AO10_og, AO13_og, ao17, ao18 ), add.cell.ids = c("ao09", "ao14", "ao15", 'ao16', 'ao10', 'ao13', 'ao17', 'ao18'), project = "ctx")


# Run the standard workflow for visualization and clustering
ao_combined <- NormalizeData(aoe13);
ao_combined <- FindVariableFeatures(ao_combined);
print("[X] Scale Data Now")
ao_combined <- ScaleData(ao_combined, vars.to.regress = c("percent.mt","nCount_RNA", "batch"), verbose = FALSE);
print("[X] Done with Scale Data")
ao_combined <- RunPCA(ao_combined, npcs = 100, verbose = FALSE);
ao_combined <- RunHarmony(ao_combined, group.by.vars = "batch")
#elbow
ElbowPlot(ao_combined, ndims = 100)
ao_combined <- JackStraw(ao_combined, num.replicate = 100, dims=100);
ao_combined <- ScoreJackStraw(ao_combined, dims=1:100)
JackStrawPlot(ao_combined, dims = 1:100)

# t-SNE and Clustering
ao_combined <- FindNeighbors(ao_combined, reduction = "harmony", dims = 1:79);
ao_combined <- FindClusters(ao_combined, resolution = 2.0);
ao_combined <- RunUMAP(ao_combined, reduction = "harmony", dims = 1:79,spread=0.8);

# Visualization
DimPlot(ao_combined, group.by = c("genotype", "ident", "batch"), ncol = 3)
DimPlot(ao_combined, reduction = "umap", label = TRUE, label.size = 5, repel = T);
#
nbt=BuildClusterTree(ao_combined,reorder = TRUE,reorder.numeric = TRUE,dims = 1:70);
PlotClusterTree(nbt);
DimPlot(nbt, label = T, label.size = 5, repel = T); #
ao_combined <- nbt;
#
#Visualize save
p1 <- DimPlot(ao_combined, reduction = "umap", group.by = "genotype");
p2 <- DimPlot(ao_combined, reduction = "umap", label = TRUE, label.size = 5, repel = T)+NoLegend();
plot_grid(p1, p2);
pdf("aoe13_merged_umap.pdf", width=7.2, height=6);
plot_grid(p1, p2);
dev.off();

pdf("aoe13-merged-umap-split.pdf", width=7.2, height=6);
DimPlot(ao_combined, reduction = "umap", split.by = "genotype", label = T, repel = T,label.size = 5);
dev.off();
#
################################################################## find markers and do hypergeometric enrichment #############################################################3
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
#fnction
source('./ao/eao/hypergeometric_test.R');

# Allows us to save objects with a new name defined by a string
saveit <- function(..., string, file) {
  x <- list(...)
  names(x) <- string
  save(list=names(x), file=file, envir=list2env(x))
}


DefaultAssay(ao_combined) <- "RNA";

results_folder= "./ao/data_pipeline1/mergee13"

seu_obj.markers <- FindAllMarkers(object = ao_combined, only.pos=TRUE, min.pct = 0.25, logfc.threshold= 0.25);
seu_obj.markers <- filter(seu_obj.markers, p_val_adj <= 0.05);
save_name <- paste(unlist(strsplit(noquote(unique(ao_combined$age)), '_'))[[1]], "merge-all_mrks", sep="_");
print(save_name)
print(paste0("[X] Found Markers for ", print(deparse(substitute(ao_combined))))); 
saveit(seu_obj_sub=seu_obj.markers, string=save_name, file=file.path(results_folder,paste0(save_name,".rdata")));
setwd(results_folder)
write.csv(seu_obj.markers,file=paste(
  save_name,"_all_cluster_markers.csv", 
  sep =""), row.names=F);

print(paste0("[X] Subsetting Markers ", print(deparse(substitute(ao_combined))))); 

seu_obj_tbl <- seu_obj.markers[,6:7];
seu_obj_tbl <- seu_obj_tbl[c(2,1)];
colnames(seu_obj_tbl) = c("Gene","DEFINITION");

save_name <- paste(unlist(strsplit(noquote(unique(ao_combined@meta.data$age)), '_'))[[1]], "no-dubs", sep="_");
saveit(object=seu_obj_tbl, string=save_name, file=file.path(results_folder,paste0(save_name,".rdata")));

write.csv(seu_obj_tbl,file=paste(
  save_name,"_cluster_gene_markers.csv", 
  sep =""), row.names=F);


print(paste0("[X] Saved both markers ", print(deparse(substitute(ao_combined))))); 

#overlap_sets(table1 = zylka_markers, table2 = seu_obj_tbl, background_n=length(rownames(ao_combined)), file_prefix = paste(save_name,"hypergeo_zylkaAll",sep = "_"))
overlap_sets(table1 = zylka_e, table2 = seu_obj_tbl, background_n=length(rownames(ao_combined)), file_prefix = paste(save_name,"hypergeo_zylka_e14",sep = "_"))
#overlap_sets(table1 = zylka_p0, table2 = seu_obj_tbl, background_n=length(rownames(ao_combined)), file_prefix = paste(save_name,"hypergeo_zylka_p0",sep = "_"))

#overlap_sets(table1 = da_dev_ctx, table2 = seu_obj_tbl, background_n=length(rownames(ao_combined)), file_prefix = paste(save_name,"hypergeo_daAll",sep = "_"))
overlap_sets(table1 = da_e13, table2 = seu_obj_tbl, background_n=length(rownames(ao_combined)), file_prefix = paste(save_name,"hypergeo_da_e13",sep = "_"))
#overlap_sets(table1 = da_p1, table2 = seu_obj_tbl, background_n=length(rownames(ao_combined)), file_prefix = paste(save_name,"hypergeo_da_p0",sep = "_"))

################################### label transfer
##########################################################################################
###label your dataset
##post p0 merged by integration by genotype
#attempt to label transfer from zylka
#load zylka P0 object where cell.ids have the published cell identity in zylka paper, Loo et al 
load("./zylka/zylka-p0-cellBarcodes/zylka-p0-seuratObject_CellIDsByBarcodes.rdata");
load("./zylka/zylkaE14_matched/zylkaE14_byBarcodes+CellTypesInMetadata.rdata")
load("./ao/eao/arlotta_dev_mouse/p1/arlotta_p1_dev_mouse_ctx.rdata"); #gives #da_p1

load("./ao/eao/arlotta_dev_mouse/arlotta_dev_mouse_ctx.rdata");
#gives da_devCtx

load("./ao/eao/arlotta_dev_mouse/e13/arlotta_e13_dev_mouse_ctx.rdata");
#gives da_e13

###ok since i am trying to label the ao P0 dataset that was integrated
pancreas.query <- ao_combined;
reference_data = da_e13
##zylka p0 reference
pancreas.anchors <- FindTransferAnchors(reference = reference_data, query = pancreas.query, project.query = FALSE, dims = 1:75,npcs = 75);

predictions <- TransferData(anchorset = pancreas.anchors, refdata = reference_data$New_cellType, dims = 1:75);
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions);
DimPlot(pancreas.query, group.by ="predicted.id");
## put the predicted id into a metadata called CellType
pancreas.query[["CellType_da"]] = pancreas.query@meta.data$predicted.id;

hist(pancreas.query$prediction.score.max);
abline(v = 0.5, col = "red");
################### 
#pancreas.query <- ao_combined;
reference_data= zylkae14;
##zylka p0 reference
pancreas.anchors <- FindTransferAnchors(reference = reference_data, query = pancreas.query, project.query = FALSE, dims = 1:75,npcs = 75);

predictions <- TransferData(anchorset = pancreas.anchors, refdata = reference_data$CellType, dims = 1:75);
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions);
DimPlot(pancreas.query, group.by ="predicted.id");
## put the predicted id into a metadata called CellType
pancreas.query[["CellType_zylk"]] = pancreas.query@meta.data$predicted.id;

hist(pancreas.query$prediction.score.max);
abline(v = 0.5, col = "red");
#####
##

saveRDS(pancreas.query, file='aoe13-merge-lb.rds')
###
################################  Plots ######################
#
dot = DotPlot(aoe13, group.by = "tree.ident", features = c("Pax6", "Mki67", "Emx2","Emx1",'Lhx2', "Eomes", "Erbb4","Sst", "Gad2","Foxp1", "Foxp2",'Tle4',"Bcl11b",'Fezf2','Foxg1','Neurod6','Neurod2', "Bcl6","Cux1", "Rorb","Cdh8","Satb2", "Mef2c", "Rbfox3","Reln"));
dotr = dot + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust=1))  + scale_colour_gradient2(low = "#D4E6F1", mid = "#F8F8FF", high = "#CD0045");
dotr;
pdf("aoe13_tree_ident_dot_plot.pdf", width=12, height=6);
dotr
dev.off();
####

plt1 <- ggplot(ao_combined@meta.data, aes(x=ao_combined@meta.data$tree.ident, fill=genotype)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
pdf("aoe13_merged_genotype.pdf", width=8, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_combined@meta.data, aes(x=ao_combined@meta.data$tree.ident, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
pdf("aoe13_merged_orig-ident.pdf", width=8, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_combined@meta.data, aes(x=ao_combined@meta.data$tree.ident, fill=batch)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
pdf("aoe13_merged_batch.pdf", width=8, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_combined@meta.data, aes(x=ao_combined@meta.data$tree.ident, fill=predicted_doublet)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
pdf("aoe13_merged_doublet.pdf", width=8, height=4);
plt1
dev.off();

plt1 <- ggplot(ao_combined@meta.data, aes(x=ao_combined@meta.data$tree.ident, fill=Phase)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
pdf("aoe13_merged_phase.pdf", width=8, height=4);
plt1
dev.off();

pdf("aoe13-merged-by-genotype-umap-split.pdf", width=10, height=6);
DimPlot(ao_combined, reduction = "umap", split.by = "genotype", pt.size = 2, label = T, repel = T,label.size = 5);
dev.off();
#
pdf("aoe13-merged-umap-lab.pdf", width=10, height=6);
DimPlot(ao_combined, reduction = "umap", group.by = "tree.ident", pt.size = 2, label = T, repel = T,label.size = 5);
dev.off();
#
pdf("aoe13-merged-umap-split_by_cluster.pdf", width=10, height=6);
DimPlot(ao_combined, split.by = 'tree.ident', ncol = 6, label = T, repel=T, pt.size = 2)+NoLegend()
dev.off();
#
pdf("aoe13-merged-by-umap-nolab.pdf", width=10, height=6);
DimPlot(ao_combined, reduction = "umap", group.by = "tree.ident", label = F, repel = T,label.size = 5, pt.size=1);
dev.off();
#
#
pdf("aoe13-merged-by-umap-orig.ident.split.pdf", width=10, height=6);
DimPlot(ao_combined, reduction = "umap", group.by = "tree.ident",split.by = 'orig.ident', label = F, repel = T,label.size = 5, pt.size=2, ncol=2
);
dev.off();

pdf("aoe13-merged-by-predID.pdf", width=10, height=6);
DimPlot(pancreas.query, reduction = "umap", group.by = "predicted.id", label = T, repel = T,label.size = 5, pt.size=1);
dev.off();
#
#################################################################################################
#######################################################################################################
#######################################  label clusters now too #######################################
obj <- RenameIdents(object = aoe13, 
                    
                    "1" = "ChoroidPlexus",
                    "2" = "Microglia",
                    "3" = "Pericytes",
                    "4" = "Pericytes1",
                    "5" = "Endothelial",
                    "6" = "Pericytes2",
                    
                    "7" = "L1",
                    "8" = "Int.1",
                
                    
                    "10" = "Int.2",
                    "11" = "Int.3",
                    "12" = "Int.4",
                    "13" = "Int.5",
                    "14" = "Int.6",
                    "15" = "Thalm/Int",
                    "16" = "Thalm/Int.1",
                    "17"= "Thalm/Int.2",
                    "18" = "Thalm/Int.3",
                    "19" = "Int.7", "20" = "Int.8",
                    
                    "21" = "Stri.Inhib",
                    "22" = "Stri.Inhib.1",
                    "25" = "MigNeurons",
                    "26" = "MigNeurons.1",
                    
                    '9'= 'L5-6.0',
                    "23" = "L5-6.1",
                    "24" = "L5-6.2",
                   
                    "27"= "L5-6.3",
                    "28" = "L5-6.4",
                    "29" = "cortical hem", "30" = "RG.1",
                    "31" = "RG.2",
                    "32" = "RG.3",
                    "33" = "RG.4",
                    "34" = "RG.5",
                    "35" = "RG.6",
                    "36" = "RG.7",
                    "37"= "RG.8",
                    "38" = "RG.9",
                    "39" = "IP.1", "40" = "IP.2",
                    "41"="G.E.","42"="G.E.1", "43"="G.E.2"
);
####
obj$CellType1 <- Idents(obj);
# Next, switch the identity class of all cells to reflect replicate ID
obj <- SetIdent(obj, value=obj$CellType1);
####### set idents back to treee.idents to collapsed IDs #########################

#############################################################################################
################################### RE-SET IDENTS ###################################
obj <- SetIdent(obj, value=obj$tree.ident);
#######################################  label clusters now too #######################################
obj <- RenameIdents(object = obj, 
                    
                    "1" = "ChoroidPlexus.E13",
                    "2" = "Microglia.E13",
                    "3" = "Pericytes.E13",
                    "4" = "Pericytes.E13",
                    "5" = "Endothelial.E13",
                    "6" = "Pericytes.E13",
                    
                    "7" = "L1.E13",
                    "8" = "Int.E13",
                    '9'= 'L5-6.E13',
                    "10" = "Int.E13",
                    "11" = "Int.E13",
                    "12" = "Int.E13",
                    "13" = "Int.E13",
                    "14" = "Int.E13",
                    "15" = "Thalm/Int.E13",
                    "16" = "Thalm/Int.E13",
                    "17"= "Thalm/Int.E13",
                    "18" = "Thalm/Int.E13",
                    "19" = "Int.E13", "20" = "Int.E13",
                    
                    "21" = "Stri.Inhib.E13",
                    "22" = "Stri.Inhib.E13",
                    "23" = "L5-6.E13",
                    "24" = "L5-6.E13",
                    "25" = "MigNeurons.E13",
                    "26" = "MigNeurons.E13",
                    "27"= "L5-6.E13",
                    "28" = "L5-6.E13",
                    "29" = "cortical hem.E13", "30" = "RG.E13",
                    "31" = "RG.E13",
                    "32" = "RG.E13",
                    "33" = "RG.E13",
                    "34" = "RG.E13",
                    "35" = "RG.E13",
                    "36" = "RG.E13",
                    "37"= "RG.E13",
                    "38" = "RG.E13",
                    "39" = "IP.E13", "40" = "IP.E13",
                    "41"="G.E.E13","42"="G.E.E13", "43"="G.E.E13"
);
####
obj$CellType1_coll <- Idents(obj);
#################################################################
######################### save the CellTypes for all cells
aoe13 <- obj
saveRDS(aoe13, file='aoe13_merge_celltype.rds')
############################################################################
###########################################################################
#########################    subset  e13 ###################################################
obj$CellType1_coll
obj <- SetIdent(obj, value=obj$tree.ident);
qs <- WhichCells(aoe13m, idents = c('16','17','18','19','20','21'))
DimPlot(aoe13m, cells.highlight = qs, group.by = 'predicted.id', split.by = 'genotype')
int_lin = WhichCells(aoe13, idents =c('41','42','43','8','11','12','13','14','19','20'))
exc_cells = WhichCells(aoe13, idents =c('9','23','24','25','26','27','28','30','31','32','33','34','35','36','37','38','39','40'))
e13_exc = subset(aoe13, cells = exc_cells, invert=F)
e13_int = subset(aoe13, cells = int_lin, invert=F)
saveRDS(exc_cells, file="e13_exci_lineage_barcodes.rds")
saveRDS(int_lin, file="e13_interneuron_lineage_barcodes.rds")
#
saveRDS(e13_exc, file = "e13_exci_lin.rds")
saveRDS(e13_int, file = "e13_interneuron_lin.rds")
#####################################################################################################
#########################  post subset: exci  ##################################################
##########################################################################################################################################
################################################## re-cluster ####################################################
e13_exc <- FindVariableFeatures(e13_exc, selection.method = "vst", nfeatures = 2000);
Sys.time()
e13_exc <- ScaleData(e13_exc, verbose = FALSE,vars.to.regress = c("percent.mt","nCount_RNA", "batch"));
Sys.time()
e13_exc <- RunPCA(e13_exc, npcs = 70, verbose = FALSE); Sys.time(); ElbowPlot(e13_exc, ndims = 70);
#dims
e13_exc <- JackStraw(e13_exc, num.replicate = 100, dims = 70);
e13_exc <- ScoreJackStraw(e13_exc, dims = 1:70); JackStrawPlot(e13_exc, dims = 1:70)
## Harmony
e13_exc <- RunHarmony(e13_exc, group.by.vars = "batch")

# t-SNE and Clustering
e13_exc <- FindNeighbors(e13_exc,reduction = "harmony", dims = 1:48);
e13_exc <- FindClusters(e13_exc, resolution = c(0.8,1.0,1.2,2.0));
Idents(e13_exc) <- "RNA_snn_res.1.2"
e13_exc <- RunUMAP(e13_exc, reduction = "harmony", dims = 1:48,spread = 0.8);
# Visualization
DimPlot(e13_exc, reduction = "umap", label = TRUE, label.size = 5, repel = T);
#
nbt=BuildClusterTree(e13_exc,reorder = TRUE,reorder.numeric = TRUE,dims = 1:48);
PlotClusterTree(nbt);
DimPlot(nbt, label = T, label.size = 5, repel = T); #
DimPlot(nbt, label = T, label.size = 5, repel = T, split.by = 'tree.ident', ncol=5)
#
e13_exc <- nbt;
saveRDS(e13_exc, file = "e13_exci_lin.rds")
## Use FindAllMarkers
## Use hypergeometric enrichment
## annotate 
########################################3
#######################################################################################################
############################################################################################################################
########################################### re-rename clusters tree.ident   #############################################################
obj <- nbt;
obj <- SetIdent(obj, value=obj$tree.ident);
###
obj <- RenameIdents(object = obj, 
                    "10" = "RG",
                    "11"="RG.1",
                    "12"="RG.2",
                    "13" = "RG.3",
                    "14" = "RG.4",
                    "15" = "RG.5",
                    "16" = "RG.6",
                    "7" = "IP",
                    "9" = "IP.1",
                    '8'= 'MigNeurons',
                    
                    "1" = "L5-6.1",
                    "2" = "L5-6.2",
                    "3" = "L5-6.3",
                    "4" = "L5-6.4",
                    "5" = "L5-6.5",
                    "6" = "L5-6.6"
                    );
obj$CellType <- Idents(obj);
# Next, switch the identity class of all cells to reflect replicate ID
obj <- SetIdent(obj, value=obj$CellType);
DimPlot(obj, label = T, label.size = 5, repel = T); #
#############################
obj <- SetIdent(obj, value=obj$tree.ident);
obj <- RenameIdents(object = obj, 
                    "10" = "RG.E13",
                    "11"="RG.1.E13",
                    "12"="RG.2.E13",
                    "13" = "RG.3.E13",
                    "14" = "RG.4.E13",
                    "15" = "RG.5.E13",
                    "16" = "RG.6.E13",
                    "7" = "IP.E13",
                    "9" = "IP.1.E13",
                    '8'= 'MigNeurons',
                    
                    "1" = "L5-6.1.E13",
                    "2" = "L5-6.2.E13",
                    "3" = "L5-6.3.E13",
                    "4" = "L5-6.4.E13",
                    "5" = "L5-6.5.E13",
                    "6" = "L5-6.6.E13"
);
obj$CellType_age <- Idents(obj);
## Next, switch the identity class of all cells to reflect replicate ID
#obj <- SetIdent(obj, value=obj$CellType_age);
DimPlot(obj, label = T, label.size = 5, repel = T, group.by = 'CellType_age');
### next
obj <- SetIdent(obj, value=obj$tree.ident);

obj <- RenameIdents(object = obj, 
                    "10" = "RG.E13",
                    "11"="RG.E13",
                    "12"="RG.E13",
                    "13" = "RG.E13",
                    "14" = "RG.E13",
                    "15" = "RG.E13",
                    "16" = "RG.E13",
                    "7" = "IP.E13",
                    "9" = "IP.E13",
                    '8'= 'MigNeurons.E13',
                    
                    "1" = "L5-6.E13",
                    "2" = "L5-6.E13",
                    "3" = "L5-6.E13",
                    "4" = "L5-6.E13",
                    "5" = "L5-6.E13",
                    "6" = "L5-6.E13"
);
obj$CellType_coll <- Idents(obj);
# Next, switch the identity class of all cells to reflect replicate ID
obj <- SetIdent(obj, value=obj$CellType_coll);
DimPlot(obj, label = T, label.size = 5, repel = T); #
############ save with IDs newest cell type labeled e13 object
obj <- SetIdent(obj, value=obj$CellType);
setwd("./ao/data_pipeline1/mergee13/exci_lin")
objdir = "./ao/data_pipeline1/mergee13/exci_lin"
e13_exc <- obj; saveRDS(e13_exc, file = "e13_exci_lin_celltypes_obj.rds")
e13_rna <- readRDS(paste0(objdir,"/","e13_exci_lin_celltypes_obj.rds"))
##arrange
e13_rna$CellType_coll <- factor(e13_rna$CellType_coll,levels=c("RG.E13","IP.E13","MigNeurons.E13", "L5-6.E13"))
e13_rna <- SetIdent(e13_rna, value = "CellType_coll")
