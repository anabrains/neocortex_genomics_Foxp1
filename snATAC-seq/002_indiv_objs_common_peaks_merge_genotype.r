## intermediate_qc_one_ibj_before_all_merge_join_peaks
library(Seurat)
library(gridExtra)
library(dplyr)
library(magrittr)
library(SingleCellExperiment)
library(cellAlign)
library(ggpubr)
library(destiny)
library(ggthemes)
library(Signac)
require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(GenomeInfoDb)
library(patchwork)
library(biovizBase) #BiocManager::install("biovizBase")
library(EnsDb.Mmusculus.v79) #bioCManager install if not found

########
########
read_files_dir = "./sc-atac/aoa"
## set file save directory
objdir = "./sc-atac/aoa/processed_objs"
#set pdf save directory
plotdir = "./sc-atac/aoa/plots"
#peak directory
peakdir = "./sc-atac/aoa/peaks"
##save pre-subseted objs
savedir_presub = "./sc-atac/aoa/processed_objs/pre-subset"
###############################################################################################################################################################################
####################################################################### CHANGE ME AKA WHICH FILES #############################################################################
## Automate listing dirs, pulling files ##seems like last char is auto wildcar
dir_list = dir(path = read_files_dir, pattern = "AOE13")
## refine dir list ##onnly cKOs at E13
##### CHANGE ME!!
objects <- c(1,5:6)
dir_list = dir_list[objects]

##not in loop bc want parallel in non parallel instances
  animal = dir_list[1]
##### metadata##############################
  orig_ident = unlist(strsplit(animal, split = "_"))[1]
  geno = unlist(strsplit(animal, split = "_"))[2]
  age = unlist(strsplit(animal, split = "_"))[3]
  batch_id = unlist(strsplit(animal, split = "_"))[4]
#########################################################################################################################################
#########################################################################################################################################
## cellranger meta data
##also called id: paste0("meta_",unlist(strsplit(dir_list[1], split = "_"))[1])
meta_e13a1 <- read.table(
  file = paste0(read_files_dir,"/",dir_list[1],"/outs/","singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
#
e13a1_cr_counts <- Read10X_h5(filename = paste0(read_files_dir,"/",dir_list[1],"/outs/","filtered_peak_bc_matrix.h5"));
## generate macs2 peaks
#macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOP03_CKO_P0_1/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2")
## or loacs already generated macs2 peaks ## previously combined per genotype per timepoint
macs2_peaks <- get(load(paste0(peakdir,"/","combined_peaks_e13_ctrl.rdata")))

chrom_assay <- CreateChromatinAssay(
  counts = e13a1_cr_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments =  paste0(read_files_dir,"/",dir_list[1],"/outs/","fragments.tsv.gz"),
#  min.cells = 3,
#  min.features = 200
)

#make seuart obj
cr_seu_obj <- CreateSeuratObject(chrom_assay, assay = "peaks", meta.data = meta_e13a1)
#
library(future)
plan("multicore", workers = 10)
# make counts table with CR cells but MACS2 peaks #takes time
MACS2_COUNTS <- FeatureMatrix(fragments = Fragments(cr_seu_obj), features = macs2_peaks, cells = colnames(cr_seu_obj) )

## Create chromatin assay using MACS2 counts
MACS2_CHRASSAY <- CreateChromatinAssay(counts = MACS2_COUNTS, genome = "mm10", fragments = Fragments(cr_seu_obj))#,min.cells = 3, min.features = 200)

## Create Seurat object using MACS2-based chromatin assay
pbmc <- CreateSeuratObject(MACS2_CHRASSAY, assay = "peaks", meta.data = meta_e13a1)
#check it out
pbmc[['peaks']]

granges(pbmc)
############################################# metadata  #########################################################
# add information to identify dataset of origin
# orig ident
pbmc$orig.ident <- unlist(strsplit(dir_list[1], split = "_"))[1]
# genotype
pbmc$genotype <- unlist(strsplit(dir_list[1], split = "_"))[2]
#age
pbmc$age <- unlist(strsplit(dir_list[1], split = "_"))[3]
# batch
pbmc$batch <- unlist(strsplit(dir_list[1], split = "_"))[4]

##################################################################################################################################################################
###################################################################### add doublet data ##########################################################################
#load the meta data saved from ArchR Doublet enrichmen etc
load("./sc-atac/aoa/run_archr/archr_e13_metadata.rdata")
##keep only specific sample here AOE13A1
#subset dataframe
data_frame_mod <- e13archrmeta[e13archrmeta$Sample=="AOE13A1",]
#remove anything before the # sign in ArchR metadata filtered column so cellIDs match
gsub(".*#", "", rownames(data_frame_mod)) -> rownames(data_frame_mod)
library(dplyr)
df_1 <- dplyr::select(data_frame_mod, DoubletEnrichment)
##########
pbmc <- AddMetaData(
  object = pbmc,
  metadata = df_1,
  col.name = 'doublet_enrichment')
#
#sum(is.na(pbmc$doublet_enrichment))
########################################## save
setwd(savedir_presub)
seu_e13a1 <- pbmc
save(seu_e13a1, file = 'e13a1.rdata')
####################################################################################################################################################################
## extract gene annotations from EnsDb 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79);
## Update annotation
## change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC';
genome(annotations) <- "mm10";

# add the gene information to the object
Annotation(pbmc) <- annotations;

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
#TSS enrichment scores by grouping the cells based on the score and plotting the accessibility signal over all TSS sites.
#takes some time
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

#plotting of the TSS enrichment signal for different groups of cells using the TSSPlot()
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend();
setwd(plotdir);
ggsave("ao-e13a1-tss-enrichment.pdf", width = 12, height = 8)

#fragment length periodicity
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group', region= 'chr1-1-10000000')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group', region= 'chr6-98930090-99163018')

#Vlnplts of QC 
setwd(plotdir)
options(repr.plot.width=12, repr.plot.height=10)
VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.001,
  ncol = 5
)
ggsave("aoe13a1-qc.pdf", width = 12, height = 8)

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.00,
  ncol = 5
)
ggsave("aoe13a1-qc1.pdf", width = 12, height = 8)
#########################
setwd(savedir_presub)
seu_e13a1 <- pbmc
save(seu_e13a1, file = 'e13a1.rdata')
###### do for all CTRL objs / then for all CKO object with the common peak set
###### AOE13A1; AOE13A2; AOE13A3; AOE13A4; AOE13A5; AOE13A6
##################################################################################################################
##############################$$$ process like this for all replicates   $$$########################################
####################################################################################################################
##### now  re-load from pre-subset folder
############################################
### load objects into memory
require(data.table)
#pre-subset in:
filenames <- list.files(path=savedir_presub, pattern = 'AOE13', full.names = T)
for (i in 1:length(filenames)){
  curr_name <- load(filenames[i])
  curr <- get(curr_name) }
#objects now loaded, subset
pbmc <- seu_e13a1  #seu_e13a5 #seu_e13a6
pbmc<- subset(x= pbmc, subset =  doublet_enrichment > 1.5, invert = TRUE)
########################################################################################
low_prf <- quantile(pbmc$peak_region_fragments, 0.02)
upper_prf <- quantile(pbmc$peak_region_fragments, 0.98)
low_prp <- quantile(pbmc$pct_reads_in_peaks, probs = 0.02)
high_blr <- quantile(pbmc$blacklist_ratio, probs = 0.98)
high_ns <- quantile(pbmc$nucleosome_signal, probs = 0.98)
low_ts <- quantile(pbmc$TSS.enrichment, probs = 0.02)
print(low_prf) ; print(upper_prf); print(low_prp); print(high_blr); print(high_ns); print(low_ts)

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > low_prf &
    peak_region_fragments < upper_prf &
    pct_reads_in_peaks > low_prp &
    blacklist_ratio < high_blr & #checkif0
    nucleosome_signal < high_ns &
    TSS.enrichment > low_ts
)
####################################################
### save them, merge within genotype
## then merge with other genotype, that have been processed the same
##### merge all datasets, adding a cell ID to make sure cell names are unique
ct1 <- pbmc
#ct2 <- pbmc
#ct3 <- pbmc

combined <- merge(
  x = ct1,
  y = list(ct2, ct3),
  add.cell.ids = c("E13A1", "E13A5", "E13A6")
)
combined[["ATAC"]]
saveRDS(combined, file=paste0(objdir,"/","e13_ctrl_combined.rds") )
############################################################ do for cko too then move to next file
