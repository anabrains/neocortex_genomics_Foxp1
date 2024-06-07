
### motifs
#BiocManager::install("JASPAR2020")
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
##
motif_dir = "./sc-atac/aoa/processed_objs/motifs/e13_v3/" ;
hdir = "./sc-atac/aoa/processed_objs/merged_harmony/";
##
mouse_brain <- readRDS(file=paste0(hdir,"/","e13_merged_exci_lin_m.rds")) 
DefaultAssay(mouse_brain) <- 'peaks'
###
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
## clean up any peaks in non-standard chromosomes ##otherwise error when run AddMotifs() 
#Error in .getOneSeqFromBSgenomeMultipleSequences(x, name, start, NA, width,  : 
gr <- granges(mouse_brain)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Mmusculus.UCSC.mm10) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
mouse_brain[['peaks']] <- subset(mouse_brain[['peaks']], features = feat.keep)

# add motif information
mouse_brain <- AddMotifs(object = mouse_brain, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
#stored #mouse_brain@assays$peaks@motifs
saveRDS(mouse_brain, file=paste0(hdir,"/","e13_merged_exci_lin_m.rds")) 
#mouse_brain <- readRDS("./sc-atac/aoa/processed_objs/merged_harmony/e13_merged_exci_lin_m.rds")
#mouse_brain <- readRDS(file=paste0(hdir,"/","e13_merged_exci_lin_m.rds")) 
#####################################
############
coembed <- SetIdent(mouse_brain, value = "CellType_coll")
coembed[["CellType_coll"]] <- Idents(obj = coembed)
coembed$CellType_coll <- factor(coembed$CellType_coll,levels=c("RG.E13","IP.E13","L5-6.mig.E13", "L5-6.E13"))
coembed <- SetIdent(coembed, value = "CellType_coll")
table(Idents(coembed))
##
merged <- coembed #query or what to annotate
Idents(merged) <- "CellType_coll"  #make sure it is active ident
merged$celltype.genotype <- paste(Idents(merged), merged$genotype, sep = "_");
Idents(merged) <- "celltype.genotype";

levels(merged) <- c('RG.E13_CTRL','RG.E13_CKO','IP.E13_CTRL','IP.E13_CKO','L5-6.mig.E13_CTRL','L5-6.mig.E13_CKO','L5-6.E13_CTRL', 'L5-6.E13_CKO')

merged[["celltype.genotype"]] <- Idents(merged)
t1=table(merged@meta.data$celltype.genotype);
t2=rownames(table(merged@meta.data$celltype.genotype));

# dap = differentially accessible peak
DefaultAssay(merged) <- 'peaks'
i=1;
while (i < length(table(Idents(merged)))) {
  control_var = rownames(table(Idents(merged)))[i];
  cko_var = rownames(table(Idents(merged)))[i+1];
  
  # want to save each individual element, so we can append it to some list  #E16 all same batch
  dapi <- FindMarkers(merged, ident.1 = cko_var, ident.2 = control_var, only.pos = FALSE,test.use = 'LR', min.pct = 0.05,latent.vars = c('nCount_peaks','batch')); 	 
  names(dapi)[1] <- paste('cluster',unlist(strsplit(control_var,"_"))[1]);
  top.da <- (dapi[dapi$p_val_adj < 0.005, ])
  
  ## names
  varName <- i;    #iteration loop, not really essential, but can be useful for QC/QA based on length of t1 or t2 above #really just counts
  cell <- cko_var;  #the cellType or cluster and the genotype that comes first in comparison group
  age <- "e13";  #age of data
  
  file_name <- paste0(motif_dir,"all_top_raw/",age,"_","dap_",cell,"_",varName,".rds", sep="");
  saveRDS(top.da, file= file_name)
  
  ### get background peaks meta data ##match the overall GC content in the peak set
  meta.feature <- GetAssayData(merged, assay = "peaks", slot = "meta.features")
  # find peaks open in either control or cko
  open.peaks <- AccessiblePeaks(merged, idents = c(control_var,cko_var))
  
  ##
  dapt <- top.da
  dapiuc <- rownames(dapt[dapt$avg_log2FC > 0.25, ])
  dapid <- rownames(dapt[dapt$avg_log2FC < -0.25, ])
  ######## calculate BG
  # match the overall GC content in the peak set, part 2
  ## enriched motics rel open in cko 
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[dapiuc, ],
    n = 10000
  )
  ## motifs enriched relatively up in cKOs
  cko_open <- FindMotifs(object = merged, features = peaks.matched)
  cko_open <- (cko_open[cko_open$pvalue< 0.05, ])
  name_up <- paste0("open_in_cko_", cell);
  
  
  ##enrichmed motifs in relatively closed in cko peaks
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[dapid, ],
    n = 10000)
  cko_closed <- FindMotifs(object = merged, features = peaks.matched)
  cko_closed <- (cko_closed[cko_closed$pvalue< 0.05, ])
  #names(closest_down)[1] <- paste('cluster',unlist(strsplit(control_var,"_"))[1]);
  name_down <- paste0("closed_in_cko_", cell);

  #file_name <- paste(motif_dir,"/",name_up,"_",age,"_","dap_",varName,".rds", sep="");
  file_name <- paste(motif_dir,name_up,"_",varName,".rds", sep="");
  saveRDS(cko_open, file= file_name)
  file_name <- paste(motif_dir,name_down,"_",varName,".rds", sep="");
  saveRDS(cko_closed, file= file_name)
  
  file_name <- paste(motif_dir,name_up,"_",varName,".csv", sep="");
  write.csv(cko_open, file= file_name, row.names=FALSE)
  file_name <- paste(motif_dir,name_down,"_",varName,".csv", sep="");
  write.csv(cko_closed, file= file_name, row.names=FALSE)
  
  
    p1 <- MotifPlot(object = merged, motifs = rownames(cko_open)[1:20]);
    p1 <- p1 + ggtitle(paste0("Top p-val adj Enriched Motifs cKO open"),cell)+ theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12),  axis.title=element_text(size=13));
    file_name <- paste0(motif_dir,"plots/",cell,"_","cko_open",".pdf", sep="");
    pdf(file_name, width=12, height=10);
    p1
    dev.off();
    
    p1 <- MotifPlot(object = merged, motifs = rownames(cko_closed)[1:20]);
    p1 <- p1 + ggtitle(paste0("Top p-val adj Enriched Motifs cKO closed"),cell)+ theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12),  axis.title=element_text(size=13));
    file_name <- paste0(motif_dir,"plots/",cell,"_","cko_closed",".pdf", sep="");
    pdf(file_name, width=14, height=10);
    p1
    dev.off();
  
  i = i+ 2;
}
