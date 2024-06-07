suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(require(dplyr))
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
#################
########
########
read_files_dir = "./sc-atac/aoa"
## set file save directory
objdir = "./sc-atac/aoa/processed_objs/"
#set pdf save directory
####
## set merged obj no additional processing save directory
mergedir = "./sc-atac/aoa/processed_objs/merged"
plotdir = "./sc-atac/aoa/plots/removing_batch_plot"
######################################################################################################3
#sep by geno datasets
e13_ko <- readRDS(file=paste0(objdir,"e13_cko_combined.rds") )
e13_ctrl <- readRDS(file=paste0(objdir,"e13_ctrl_combined.rds") )

# make combined obj snice it uses there lsi reductions
combined <- merge(
  x = e13_ctrl,
  y = list(e13_ko),
  add.cell.ids = c("e13ctrl", "e13cko")
)

######################## nornalize
## the ctrl and cko                                                                                     ### 
library(future) 
plan("multicore", workers = 16)
options(future.globals.maxSize = 8500 * 1024^2)

combined <- RunTFIDF(combined) #norm
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined) #dim analysis
#####  with harmony          
combined <- RunHarmony(object = combined, group.by.vars = c('batch', 'orig.ident'), reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
combined <- FindNeighbors(combined, dims = 2:30, reduction = 'harmony')
combined <- FindClusters(combined, algorithm = 3, resolution = 0.6, reduction = 'harmony')#, future.seed=TRUE)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'harmony')
DimPlot(combined, group.by = "genotype", label = F, repel = TRUE, pt.size = 2, label.size = 6)+ ggtitle("E13 snATAC-seq cells - harmonized") 
ggsave(paste0(plotdir,"/","aoe13_atac_UMAP_mergeAll_geno_correction_harmony.pdf"), width = 18, height = 10)
DimPlot(combined, split.by = "batch", label = F, repel = TRUE, pt.size = 2, label.size = 6)+ ggtitle("E13 snATAC-seq cells - harmonized") 
ggsave(paste0(plotdir,"/","aoe13_atac_UMAP_mergeAll_batch_correction_harmony.pdf"), width = 18, height = 10)
nbt=BuildClusterTree(combined,reorder = TRUE,reorder.numeric = TRUE,reduction='harmony')
combined <- nbt
DimPlot(combined, label=T, label.size=6, repel=T) 
ggplot(combined@meta.data, aes(x=combined@meta.data$tree.ident, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
################################################################## start treeIdent plots #######################################################
##############################################################################################################################################
#######################################################  make normalized # bar plots ###########################################################
nbt=BuildClusterTree(combined,reorder = TRUE,reorder.numeric = TRUE,reduction='harmony')
obj <- nbt
##normalize the cell numbers per cell type cluster 
library(sweep)
###get your values from the seurat object
## what to divide each cell count value by, e.g. total cell # per cluster (each column is the total cell count per cluster)
the_div <- table(obj$tree.ident)
#makes a table with rows of tree.ident, divided by Genotype in the column (get count for Genotype or w/e variable)
cellnum <- table(obj$tree.ident, obj$genotype)
#now use the sweep library #MARGIN =1 means do the function by column #divives the cellnum by total cell counts of the_div
normalize_tree.idents <- sweep(cellnum, MARGIN = 1, STATS = the_div, FUN = "/")
############ make ggplot barplot
norm_cell_mat <- as.data.frame.table(normalize_tree.idents);
#give label names
colnames(norm_cell_mat) <- c("tree.ident", "Genotype", "Proportion")  #x-lab, #legend, #y-lab
#make ggplot filled
ggplot(norm_cell_mat, aes(fill=Genotype, y=Proportion, x=tree.ident))+geom_bar(position = "Fill",stat = "identity")+ theme(axis.text.x=element_text(angle =55, vjust = 0.5))+  ggtitle("e13 Proportions of Cells in Cluster by Genotype (normalized) - harmony") + theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12),  axis.title=element_text(size=13));
ggsave(paste0(plotdir,"/","aoe13_atac_barplot-normalized_mergeAll_geno_correction_harmony_treeIdent.pdf"), width = 18, height = 10)
########## orig ident plots normalized
###get your values from the seurat object
## what to divide each cell count value by, e.g. total cell # per cluster (each column is the total cell count per cluster)
the_div <- table(obj$tree.ident)
#makes a table with rows of tree.ident, divided by orig.ident in the column (get count for orig.ident or w/e variable)
cellnum <- table(obj$tree.ident, obj$orig.ident)
#now use the sweep library #MARGIN =1 means do the function by column #divives the cellnum by total cell counts of the_div
normalize_tree.idents <- sweep(cellnum, MARGIN = 1, STATS = the_div, FUN = "/")
############ make ggplot barplot
norm_cell_mat <- as.data.frame.table(normalize_tree.idents);
#give label names
colnames(norm_cell_mat) <- c("tree.ident", "orig.ident", "Proportion")  #x-lab, #legend, #y-lab
#make ggplot filled
ggplot(norm_cell_mat, aes(fill=orig.ident, y=Proportion, x=tree.ident))+geom_bar(position = "Fill",stat = "identity")+ theme(axis.text.x=element_text(angle =55, vjust = 0.5))+  ggtitle("e13 Proportions of Cells in Cluster by orig.ident (normalized) - harmony") + theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12),  axis.title=element_text(size=13));
ggsave(paste0(plotdir,"/","aoe13_atac_barplot-normalized_mergeAll_origIdent_correction_harmony_treeIdent.pdf"), width = 18, height = 10)
################################################################## end normalized section #######################################################
#############################      UMAPS
DimPlot(nbt, group.by = "genotype", label = F, repel = TRUE, pt.size = 2, label.size = 6)+ ggtitle("e13 scATAC-seq cells - harmonized") 
ggsave(paste0(plotdir,"/","aoe13_atac_UMAP_mergeAll_geno_correction_harmony_treeIdent.pdf"), width = 18, height = 10)
DimPlot(nbt, split.by = "batch", label = F, repel = TRUE, pt.size = 2, label.size = 6)+ ggtitle("e13 scATAC-seq cells - harmonized") 
ggsave(paste0(plotdir,"/","aoe13_atac_UMAP_mergeAll_batch_correction_harmony_treeIdent.pdf"), width = 18, height = 10)
##############################################################################################################################################
################################################################## end treeIdent plots #######################################################
##############################################################################################################################################
combined <- nbt
hdir = "./sc-atac/aoa/processed_objs/merged_harmony"
saveRDS(combined, file=paste0(hdir,"/","e13_merged.rds"))
##########################################################################################################################################################################################3
############%%%%%%%%%%%%%%%%%%%%%%% label transfer  %%%%%%%%%%%#######################################################################
###################################################################################################################################################
## set file save directory
plotdir = "./sc-atac/aoa/plots"
objdir = "./sc-atac/aoa/processed_objs/"
hdir = "./sc-atac/aoa/processed_objs/merged_harmony"
combined <- readRDS(file=paste0(hdir,"/","e13_merged.rds"))


############ if combined obj looks good, label transfer from snRNA-seq  ################### 
ao_atac <- combined
gene.activities <- GeneActivity(ao_atac)
## add the gene activity matrix to the Seurat object as a new assay and normalize it
ao_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
detach('package:Signac', unload=TRUE) #to remove future() else normalize gets stuck
library(Signac)

ao_atac <- NormalizeData(
  object = ao_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ao_atac$nCount_RNA)
)
saveRDS(ao_atac, file=paste0(hdir,"/","e13_all_merged.rds"))
################ load your sn-RNA-seq data used for annotation  #CellType1 ################

#load atac data to be labeled
ao_atac <- readRDS(file=paste0(hdir,"/","e13_all_merged.rds"))
DefaultAssay(ao_atac) <- 'RNA'
#check its good
FeaturePlot(
  object = ao_atac,
  features = c('Sst','Pvalb',"Gad2","Neurod6","Rorb","Cux1", "Pax6", "Eomes","Satb2"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
############################################################################
########%%%%%%%%%          LABEL TRANSFER     %%%%%%%%%%%%%%################
#############################################################################
##load my own AO E13 reference RNA
e13_rna <- readRDS('./ao/data_pipeline1/mergee13/aoe13_merge_celltype.rds') #aoe13_all cells

varFeatures <- VariableFeatures(e13_rna);
transfer.anchors <- FindTransferAnchors(reference = e13_rna, query = ao_atac, features = varFeatures, 
                                        reference.assay = "RNA", query.assay = 'RNA', reduction = "cca");

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = e13_rna$CellType1, 
                                     weight.reduction = ao_atac[["lsi"]], dims = 2:30);

ao_atac <- AddMetaData(ao_atac, metadata = celltype.predictions)
hist(ao_atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(ao_atac$prediction.score.max > 0.5)
##Idents(ao_atac) <- 'predicted.id'
ao_atac[['CellType1']] <- ao_atac$predicted.id; Idents(ao_atac) <- ao_atac$CellType1
DimPlot(ao_atac, label=T, repel=T, label.size=8)
saveRDS(ao_atac, file=paste0(hdir,"/","e13_all_merged.rds"))
## RE-LOAD
#e13_atac <- readRDS(file=paste0(hdir,"/","e13_all_merged.rds"))
#ao_atac <- e13_atac
## plots for this
DefaultAssay(ao_atac) <- 'RNA'
dot = DotPlot(ao_atac, features = c("Foxp1","Pax6", "Mki67", "Emx2","Emx1", "Eomes", "Erbb4", "Sst","Gad2","Npy", "Cux1", "Rorb","Satb2", "Bcl11b", "Foxp2","Tle4","Tbr1", "Bcl6","Mef2c", "Rbfox3", "Reln"), group.by='predicted.id');
dotr = dot + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust=1))  + scale_colour_gradient2(low = "#D4E6F1", mid = "#F8F8FF", high = "#CD0045")+ ggtitle("Marker genes for E13 all cells")+ theme(plot.title = element_text(hjust = 0.5));
dotr;
ggsave(paste0(plotdir,"/","ao-e13-merged-dotPlot-celltype-allCells.pdf"), width = 12, height = 10)
VlnPlot(ao_atac, features = c("Emx1", "Gad2", "Eomes"), pt.size = 0, ncol=1);
VlnPlot(ao_atac, features = c("peak_region_fragments","pct_reads_in_peaks", "nucleosome_signal", "TSS.enrichment"), pt.size = 0, ncol=1);
ggsave(paste0(plotdir,"/","ao-e13-merged-VlnPlot-celltype-allCells_qc.pdf"), width = 12, height = 10)
##################################################  sanity checks   ##########################3
#arlotta e13   ## double check
load("./ao/eao/arlotta_dev_mouse/e13/arlotta_e13_dev_mouse_ctx.rdata");
#da_devCtx all cells ages #e13 only
load("./ao/eao/arlotta_dev_mouse/arlotta_dev_mouse_ctx.rdata");
#ao_atac <- pbmc1 #aoe13_ctrl_atac;

varFeatures <- VariableFeatures(da_devCtx);
transfer.anchors <- FindTransferAnchors(reference = da_devCtx, query = ao_atac, features = varFeatures, 
                                        reference.assay = "RNA", query.assay = 'RNA', reduction = "cca");

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = da_devCtx$New_cellType, 
                                     weight.reduction = ao_atac[["lsi"]], dims = 2:30);

ao_atac <- AddMetaData(ao_atac, metadata = celltype.predictions)
hist(ao_atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(ao_atac$prediction.score.max > 0.5)
ao_atac[['CellType_da']] <- ao_atac$predicted.id; Idents(ao_atac) <- ao_atac$CellType_da
DimPlot(ao_atac, group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 1, label.size = 6) + ggtitle("e13 scATAC-seq cells") 
ggsave(paste0(plotdir,"/","ao-e13-merged-UMAP-ao-celltype-allCells_da.pdf"), width = 12, height = 10)
#####################################################################
######### FindAllMarkers as needed; hypergeo
## one group of cells non-spec based on predicitons
nonspec <- WhichCells(ao_atac, idents = c('Low quality cells'))
ao_atac <- SetIdent(ao_atac, value = "CellType1")
Idents(ao_atac, cells = nonspec) <- 'Nonspec'
table(Idents(ao_atac))
ao_atac[["CellType1"]] <- Idents(ao_atac)
##################################################
##########%%%%   CellType_cpll  ordered for plotting
ao_atac = RenameIdents(ao_atac, c("ChoroidPlexus"="ChoroidPlexus","Pericytes"="Pericytes", "VLMC"="VLMC", "VLMC.1"="VLMC","Microglia","cortical hem"="cortical hem","Endothelial"="Endothelial",
                                  "G.E."="G.E.","G.E.1"="G.E.","G.E.2"="G.E.", "Int.3"="Int", "Int.5"="Int","Int.6"="Int","Int.7"="Int","Int.8"="Int",
                                  "Nonspec"="Nonspec", "Stri.Inhib"="Stri.Inhib","Stri.Inhib.1"="Stri.Inhib","Thalm/Int"="Thalm/Int", "Thalm/Int.2"="Thalm/Int","Thalm/Int.3"="Thalm/Int",
                                  "L1"="L1", "RG.2"= "RG","RG.3"="RG","RG.4"="RG","RG.5"="RG","RG.6"="RG","RG.8"="RG","RG.9"="RG", 
                                  "IP.1"="IP","IP.2"="IP", "MigNeurons"="MigNeurons","MigNeurons.1"="MigNeurons", "L5-6.0"="L5-6","L5-6.1"="L5-6","L5-6.2"="L5-6","L5-6.3"="L5-6","L5-6.4"="L5-6"))
ao_atac[["CellType1_coll"]] <- Idents(object = ao_atac)
my_cols <- c("RG"='#D461C2', 'IP'='#8C28E6', "MigNeurons"='#A1EFE7', 'L2-4'='#FFBD33',
             'L5-6'='#74CAE6',    'L1_CR'='#6bbdc1',
            'ChoroidPlexus'='#DB797A', 'Endothelial'='#B9BA5F','Endo'='#B9BA5F', 'Pericytes'='#b1f1cc', 'Microglia'='#ce851f', 'Oligodendrocytes'='#A4DFF2', 'OPCs'='#A4DFF2',
             "G.E." = "#496989", 'Int'='#6bbd80', 'Stri.Inhib.'='#6cbead','Int / Thalm'='#559766',"Thalm/Int"='#559766',
             'VLMC'="#ba9623", "cortical hem"="#627254",'Nonspec'='#bcb8b1',   'CyclingGliaCells'='#6a994e', 'Astrocytes'='#1a936f')#'Oligodendrocytes'='#619b8a')
plt1 <- DimPlot(ao_atac, group.by = 'CellType1_coll', pt.size = 3, cols = my_cols)
pdf("aoe13_all_dimplot.pdf", width=14, height=10);
plt1
dev.off();
#########################################################
obj <- SetIdent(obj, value=obj$CellType1_coll);
saveRDS(e13_all_atac,file = "./sc-atac/aoa/processed_objs/e13_all_merged.rds")
#######################################################################################################################################################
############################################   Get EXCI NEURONL LIN ###################################################################################
#####################################  EXCI LIN ONLY, RECLUSTER     ##################################################################################
######################################################################################################################################################
# get your ao_atac obj and filter only wanted cells, then calculate DARs
## rm non exci cells as list
e13_exci_atac <- subset(e13_exci_atac, cells = c(lq_cells,chpch,non_spec_thalm, inhib_nonspec,ge,ge_cells,potential_interneurons,potential_interneurons1, endo, nonexci_2), invert = T)
e13_exci_atac <- subset(e13_exci_atac, cells = c(nonexci_2), invert = T) #didnt work first time?
#check
#e13_exci_atac <- subset(ao_atac, idents = c('Apical progenitors','Immature neurons', 'Intermediate progenitors','CThPN','Migrating neurons','UL CPN'))
#e13_exci_atac <- subset(e13_exci_atac, cells = c(cmany), invert = T)
##recluster
e13_exci_atac <- FindTopFeatures(e13_exci_atac, min.cutoff = 20)
e13_exci_atac <- RunSVD(e13_exci_atac)
e13_exci_atac <- FindNeighbors(e13_exci_atac, dims = 2:30, reduction = 'lsi')
e13_exci_atac <- FindClusters(e13_exci_atac, algorithm = 3, resolution = 0.6)
#e13_exci_atac <- RunUMAP(e13_exci_atac, dims = 2:50, reduction = 'lsi')
library(harmony)
e13_exci_atac <- RunHarmony(object = e13_exci_atac, group.by.vars = c('batch','orig.ident'), reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
e13_exci_atac <- RunUMAP(e13_exci_atac, dims = 2:30, reduction = 'harmony')
nbt=BuildClusterTree(e13_exci_atac,reorder = TRUE,reorder.numeric = TRUE,reduction='harmony')
e13_exci_atac <- nbt
DimPlot(e13_exci_atac, label=T, label.size=6, repel=T) 
ggplot(e13_exci_atac@meta.data, aes(x=e13_exci_atac@meta.data$tree.ident, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));

saveRDS(e13_exci_atac, file=paste0(hdir,"/","e13_merged_exci_lin.rds")) 

DimPlot(e13_exci_atac, group.by = 'genotype', pt.size = 0.1)
DimPlot(e13_exci_atac, split.by='orig.ident',ncol=3) 
DimPlot(e13_exci_atac, split.by='tree.ident',ncol=5) 
###### re-annotate w/ exci lin
#CellType
e13_exci_atac <- readRDS(file=paste0(hdir,"/","e13_merged_exci_lin.rds"))
e13_exci <- readRDS("./ao/data_pipeline1/mergee13/exci_lin/e13_exci_lin_celltypes_obj.rds")
ao_atac <- e13_exci_atac
varFeatures <- VariableFeatures(e13_exci); #reference dataset
transfer.anchors <- FindTransferAnchors(reference = e13_exci, query = ao_atac, features = varFeatures, 
                                        reference.assay = "RNA", query.assay = 'RNA', reduction = "cca");

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = e13_exci$CellType, 
                                     weight.reduction = ao_atac[["lsi"]], dims = 2:30);
ao_atac <- AddMetaData(ao_atac, metadata = celltype.predictions)
hist(ao_atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(ao_atac$prediction.score.max > 0.5)
ao_atac$CellType <- ao_atac$predicted.id
saveRDS(ao_atac, file=paste0(hdir,"/","e13_merged_exci.rds")) 
################################################  plots  ################################################
DimPlot(ao_atac, group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 1, label.size = 6) + ggtitle("scATAC-seq cells") 
ggsave(paste0(plotdir,"/","aoe13_atac_UMAP_exci_lin_labelled.pdf"), width = 12, height = 10)

DimPlot(ao_atac, group.by = "CellType_coll", label = F, repel = TRUE, pt.size = 3, label.size = 6, cols = c("pink", "#D3B5E5","#D3B5E5", "#BBE7FE", "#8FDDE7")) + ggtitle("E13 scATAC-seq cells") 
ggsave(paste0(plotdir,"/","aoe13_atac_UMAP_exci_lin_unlabelled.pdf"), width = 18, height = 10)

ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$predicted.id, fill=orig.ident)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
ggsave(paste0(plotdir,"/","aoe13_atac_barplot_exci_lin_orig-ident.pdf"), width = 12, height = 10)
ggplot(ao_atac@meta.data, aes(x=ao_atac@meta.data$predicted.id, fill=genotype)) + geom_bar(position = "fill")+  theme(axis.text.x=element_text(angle =90, vjust = 0.5));
ggsave(paste0(plotdir,"/","aoe13_atac_barplot_exci_lin_geno.pdf"), width = 12, height = 10)
VlnPlot(ao_atac, features = c("peak_region_fragments", "pct_reads_in_peaks", "nucleosome_signal", "TSS.enrichment"), pt.size=0, ncol=2, group.by='predicted.id')
ggsave(paste0(plotdir,"/","aoe13_atac_VlnPlot_exci_lin_QC.pdf"), width = 12, height = 10)
#################### Run  ??
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = e13_exci$CellType_coll, 
                                     weight.reduction = ao_atac[["lsi"]], dims = 2:30);
ao_atac <- AddMetaData(ao_atac, metadata = celltype.predictions)
hist(ao_atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(ao_atac$prediction.score.max > 0.5)
ao_atac[['CellType_coll']] <- ao_atac$predicted.id; Idents(ao_atac) <- ao_atac$CellType_coll
#CoveragePlot(object = ao_atac, region = c("chr6-144154829-144209568"),extend.upstream=100,extend.downstream = 100)
#ggsave(paste0(plotdir,"/","ao-e13-merged-CovPlot_Sox5_ctc.pdf"), width = 8, height = 6)
saveRDS(ao_atac, file=paste0(hdir,"/","e13_merged_exci.rds")) 

# ###### for future steps
# # Assuming 'combined' is your Seurat object with a chromatin assay
# chrom_assay <- combined[["peaks"]]
# # Get the peak names
# peak_names <- rownames(chrom_assay)
# # Filter out peaks from unwanted chromosomes chrM
# filtered_peak_names <- peak_names[!grepl("^chrX|^chrY|^chrM|^J", peak_names)]
# # Subset the chromatin assay to include only the filtered peaks
# filtered_chrom_assay <- chrom_assay[filtered_peak_names, ]
# ### Ensure that the cell names match
# ####rownames(filtered_chrom_assay) == rownames(combined)
# # Replace the original chromatin assay with the filtered one
# combined[["peaks"]] <- filtered_chrom_assay
# # Add the filtered chromatin assay to the Seurat object as a new assay
# combined[["filtered_peaks"]] <- filtered_chrom_assay
# ####
# # Assuming 'combined' is your Seurat object with a chromatin assay
# chrom_assay <- combined[["peaks"]]

# # Get the peak names
# peak_names <- rownames(chrom_assay)

# # Filter out peaks from unwanted chromosomes
# filtered_peak_names <- peak_names[!grepl("^chrX|^chrY|^chrM|^J", peak_names)]

# # Subset the chromatin assay to include only the filtered peaks
# filtered_counts <- chrom_assay@counts[filtered_peak_names, ]

# # Create a new ChromatinAssay with the filtered counts
# filtered_chrom_assay <- CreateChromatinAssay(
#   counts = filtered_counts,
#   min.cells = 0,
#   min.features = 0
# )
# ## Add the filtered chromatin assay to the Seurat object as a new assay ##could also put in a new slot #"filt_peaks"
# combined[["peaks"]] <- filtered_chrom_assay
# ####
# combined <- RunTFIDF(combined) #normalize
# combined <- FindTopFeatures(combined, min.cutoff = 20)
# combined <- RunSVD(combined) #dim reduction
# ### for better clustering
# combined <- RunHarmony(object = combined, group.by.vars = c('batch', 'orig.ident'), reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
# combined <- FindNeighbors(combined, dims = 2:30, reduction = 'harmony')
# combined <- FindClusters(combined, algorithm = 3, resolution = 0.6, reduction = 'harmony')#, future.seed=TRUE)
# combined <- RunUMAP(combined, dims = 2:30, reduction = 'harmony')
# combined=BuildClusterTree(combined,reorder = TRUE,reorder.numeric = TRUE,reduction='harmony')
# PlotClusterTree(combined)
# saveRDS(combined, file=paste0(hdir,"/e13_merged_exci.rds")) 
