######### GENERATE INDIVIDUAL OBJECT WITH JOINT PEAKS
##test reading 10x counts with MACS2 generated peaks
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
read_files_dir = "./sc-atac/aoa"
## set file save directory
#objdir = "./ao/sc-atac/aoa/signac-macs2-cr-counts"
#set pdf save directory
plotdir = "./sc-atac/aoa/plots"
#peak directory
peakdir = "./sc-atac/aoa/peaks"
# macs2 peaks dir
macs2peakdir = "./ao/sc-atac/aoa/signac-macs2/"
###############  Generate MACS2 peaks from wrapped function in Signac
cr_metadata <- read.csv(
  file = paste0(read_files_dir,"/AOE13A1_CTRL_E13_1/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
);
cr_counts <- Read10X_h5(filename = paste0(read_files_dir,"/AOE13A1_CTRL_E13_1/outs/filtered_peak_bc_matrix.h5"));
## generate macs2 peaks
macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOE13A1_CTRL_E13_1/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2", effective.genome.size = 1.87e+9)
save(macs2_peaks, file = paste0(peakdir,"aoe13a1_peaks.rdata"))
### Do for each replicate per time point
## e.g. next
# macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOE13A2_CKO_E13_1/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2", effective.genome.size = 1.87e+9)
# save(macs2_peaks, file = paste0(peakdir,"aoe13a2_peaks.rdata"))
# macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOE13A3_CKO_E13_2/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2", effective.genome.size = 1.87e+9)
# save(macs2_peaks, file = paste0(peakdir,"aoe13a3_peaks.rdata"))
# macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOE13A4_CKO_E13_2/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2", effective.genome.size = 1.87e+9)
# save(macs2_peaks, file = paste0(peakdir,"/","aoe13a4_peaks.rdata"))
# macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOE13A5_CTRL_E13_2/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2", effective.genome.size = 1.87e+9)
# save(macs2_peaks, file = paste0(peakdir,"aoe13a5_peaks.rdata"))
# macs2_peaks <- CallPeaks(object = paste0(read_files_dir,"/AOE13A6_CTRL_E13_2/outs/fragments.tsv.gz"), macs2.path = "./.conda/envs/peaks/bin/macs2", effective.genome.size = 1.87e+9)
# save(macs2_peaks, file = paste0(peakdir,"aoe13a6_peaks.rdata"))
##########################  LOAD INDIV PEAKS TO MAKE COMMON
##save pre-subseted objs
savedir_presub = "./sc-atac/aoa/processed_objs/pre-subset"
###############################################################################################################################################################################
####################################################################### CHANGE ME AKA WHICH FILES #############################################################################
## Automate listing dirs, pulling files ##seems like last char is auto wildcard
dir_list = dir(path = read_files_dir, pattern = "AOE13")
## refine dir list ## ctrl at E13
##### CHANGE ME!!
objects <- c(1,5:6)
dir_list = dir_list[objects]

##not in loop #wants as is sep. instances
  animal = dir_list[1]
##### metadata##############################
  orig_ident = unlist(strsplit(animal, split = "_"))[1]
  geno = unlist(strsplit(animal, split = "_"))[2]
  age = unlist(strsplit(animal, split = "_"))[3]
  batch_id = unlist(strsplit(animal, split = "_"))[4]
#####
###################################
## get peak information
  ####  CHANGE ME
peak_list <- list.files(peakdir, pattern = "aoe13")
  ####  CHANGE ME 
peak_list = peak_list[objects]
###########################################################################################################################################################################
###############################################################################################################################################################################
# read in peak sets #ctrl
peaks_e13a1 <- get(load(paste0(peakdir,peak_list[1])))
peaks_e13a5 <- get(load(paste0(peakdir,peak_list[2])))
peaks_e13a6 <- get(load(paste0(peakdir,peak_list[3])))

# convert to genomic ranges
gr.e13a1 <- makeGRangesFromDataFrame(peaks_e13a1)
gr.e13a5 <- makeGRangesFromDataFrame(peaks_e13a5)
gr.e13a6 <- makeGRangesFromDataFrame(peaks_e13a6)
#  create a common set of peaks
# Create a unified set of peaks to quantify in each dataset ## reduce does merging of intersecting peaks
combined.peaks <- Signac::reduce(x = c(gr.e13a1, gr.e13a5, gr.e13a6))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
save(combined.peaks, file = paste0(peakdir,'combined_peaks_e13_ctrl.rdata'))