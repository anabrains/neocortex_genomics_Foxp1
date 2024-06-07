library(Signac)
library(tidyverse)
library(ggrepel)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(BiocSingular)
library(clusterProfiler)
library(enrichR)
library(dplyr)
library(Seurat)
library(stringr)
###################
#rm(list=setdiff(ls(), c("e16_atac","merged","convert_mouse_to_human","i")))
#######
hdir = "./sc-atac/aoa/processed_objs/merged_harmony/"
dardir = "./sc-atac/aoa/processed_objs/dars_harmony/p0_updates/"
updir = "./sc-atac/aoa/processed_objs/dars_harmony/p0_updates/up_in_cko_rds/"
downdir = "./sc-atac/aoa/processed_objs/dars_harmony/p0_updates/down_in_cko_rds/"
godir = "./sc-atac/aoa/processed_objs/dars_harmony/p0_updates/go/"
##
###
#########################
merged = readRDS(file=paste0(hdir,"e13_merged_exci_lin.rds"))  

#e16_atac = readRDS(file=paste0(hdir,"e16_merged_exci_lin.rds"))   

#merged = readRDS(file=paste0(hdir,"p0_exci_cleaned.rds")) 
########
ao_atac = merged;
DefaultAssay(merged) <- 'peaks'
Idents(merged) <- "CellType_coll"  #make sure it is active ident
#merged <- RenameIdents(merged, "RG.E13" = "aRG.E13","IP.E13"="IP.E13",'L5-6.mig.E13'='MigNeurons.E13',"L5-6.E13"="L5-6.E13");
merged$celltype.genotype <- paste(Idents(merged), merged$genotype, sep = "_");
Idents(merged) <- "celltype.genotype";
t2=rownames(table(merged@meta.data$celltype.genotype)); t2 ##is in order and ckO is first
#merged$celltype.genotype <- factor(merged$celltype.genotype,levels=c("aRG.E13_CTRL","aRG.E13_CKO", "IP.E13_CKO","IP.E13_CTRL","MigNeurons.E13_CKO", "MigNeurons.E13_CTRL","L5-6.E13_CKO","L5-6.E13_CTRL"))
#t2=rownames(table(merged@meta.data$celltype.genotype)); t2
##table((merged$celltype.genotype))
# # only autosomes
t1=table(merged@meta.data$celltype.genotype);
t2=rownames(table(merged@meta.data$celltype.genotype));

# Basic function to convert mouse to human gene names
source('./sc-atac/aoa/convertMouseToHuman.r')
library(future);
plan(strategy='multicore',workers=30);
###22.67 GiB needed for this; input in bytes
options(future.globals.maxSize = 99990 * 1024^2);
#####
table(Idents(merged))
i = 1; #for E13 start at3 to skip IP.1
while (i < length(table(merged$celltype.genotype))) {
  ####bc i ordered them so contrl would be first, otherwise alplabetocal default is cKO first
  cko_var = rownames(table(merged$celltype.genotype))[i];  control_var = rownames(table(merged$celltype.genotype))[i+1];
  #control_var = rownames(table(Idents(merged)))[i+1]; #cko_var = rownames(table(Idents(merged)))[i];
  started = Sys.time();
      print(paste0("We are now on: ", cko_var, " and ", control_var))

  # ident1 = mut ident2 = control, different than before  #now if positive FC = up in CKO bc group 1 is comparison group ## E16 ATAC has no batch diff b/c all same batch
  degi <- FindMarkers(merged, ident.1 = cko_var, ident.2 = control_var, only.pos = FALSE,test.use = "LR", verbose = FALSE,latent.vars = c('peak_region_fragments','batch'), min.pct = 0.05); 
  #names(degi)[1] <- paste0('cluster_',unlist(strsplit(control_var,"_"))[1]);  
  ##add regions detected, name the column
  degi <- cbind(rownames(degi), data.frame(degi, row.names=NULL)); 
  names(degi)[1] <- "dar_region";
  #clean
  ######  names as labels for saving
  count <- i;    #iteration loop, not really essential, but can be useful for QC/QA based on length of t1 or t2 above #really just counts
  cell <- cko_var;  #the cellType or cluster and the genotype that comes first in comparison group
  cell <- gsub(".", "_", cell, fixed=TRUE) #replace period with underscore
  label <- gsub('(.*)_\\w+', '\\1', cell)  # cell name is RG_E13_CTRL; strip everything after last _ to get ride of CTRL or CKO
  age<-gsub("^.*\\_","",label) ## label was IP_E13, remove 1st underscore and everything before
  ###add column label name
  degi$Cell <- label;
  #save raw DARs
  file_name <- paste0(dardir,"raw_dars/",label,"_raw_dars.csv", sep="");
  write.csv(degi, file= file_name, row.names=FALSE);
  # save raw
  file_name <- paste(dardir,label,"_dar_filt.rds", sep="");
  saveRDS(degi, file= file_name)
  ####filter sort 
  filt = dplyr::filter(degi, p_val_adj < 0.05); 
  filt$Cell <- label;
      print(paste0( "the number of down DARs is ",nrow(filt)))

  ##### Get Names of closest genes, open / closed DARs relative
  degt <- filt
  rownames(degt) <- NULL
  rownames(degt) <- degt[,1] #added regions back as row names
  
  degiuc <- rownames(degt[degt$avg_log2FC > 0.25, ])
  degid <- rownames(degt[degt$avg_log2FC < -0.25, ])

    print(paste0( "the number of up DARs is ",nrow(degiuc)))
    print(paste0( "the number of down DARs is ",nrow(degid)))

  #calclate features closest to DA peaks
  name_up <- paste0(cell,"_OPEN");
  closest_up <- ClosestFeature(merged, degiuc)
  closest_up$Cell <- name_up;
  
  #closets features down in first group or up in second
  name_down <- paste0(cell,"_CLOSED");
  closest_down <- ClosestFeature(merged, degid)
  closest_down$Cell <- name_down;
  
  file_name <- paste(updir,name_up,"_dar_",count,".rds", sep="");
  saveRDS(closest_up, file= file_name);
  file_name <- paste(downdir,name_down,"_dar_",count,".rds", sep="");
  saveRDS(closest_down, file= file_name);
  both_dars = rbind(closest_up, closest_down);
  saveRDS(both_dars, file= paste0(dardir,label,"_dars_all_genes_",count,".rds", sep=""));
  
  file_name <- paste0(dardir,"up_csv/",name_up,"_dar_",count,".csv", sep="");
  write.csv(closest_up, file= file_name, row.names=FALSE);
  file_name <- paste0(dardir,"down_csv/",name_down,"_dar_",count,".csv", sep="");
  write.csv(closest_down, file= file_name, row.names=FALSE);
  write.csv(both_dars, file= paste0(dardir,label,"_dars_all_genes_",count,".csv", sep=""), row.names = F);
  ## count of all DARs post filtered
  df = as.data.frame(cbind(nrow(closest_up),name_up));
  df1 = as.data.frame(cbind(nrow(closest_down),name_down));
  nums = cbind(df,df1)
  write.csv(nums, file= paste0(dardir,"raw_dars_count/",label,"_all_dars_filt_count.csv", sep=""), row.names = F); ##count dars collapsed to genes
    print(paste0("We are now on done with DARs Pt1"))

  ####################################################################################################################
  ####################################################################################################################
  ##### make unique gene lists for Up and Down
  #data frame, values to look for unique in by row; if 2 rows in list, then based on both
  #test_all <- both_dars %>% distinct(gene_name,Cell, .keep_all = TRUE)
  up_mod <- closest_up %>% distinct(gene_name, .keep_all = TRUE);
  up_mod <- up_mod %>% dplyr::select(gene_name, Cell); 
  colnames(up_mod) <- c("Gene","Cell")
  saveRDS(up_mod, file= paste0(dardir,"gene_list/raw_all/",name_up,"_dars_filt.rds", sep=""));
  write.csv(up_mod, file= paste0(dardir,"gene_list/csv/raw_all/up/",name_up,"_dars_filt.csv", sep=""), row.names = F);
  
  down_mod <- closest_down %>% distinct(gene_name, .keep_all = TRUE);
  down_mod <- down_mod %>% dplyr::select(gene_name, Cell); 
  colnames(down_mod) <- c("Gene","Cell")
  saveRDS(down_mod, file= paste0(dardir,"gene_list/raw_all/",name_down,"_dars_filt.rds", sep=""));
  write.csv(down_mod, file= paste0(dardir,"gene_list/csv/raw_all/down/",name_down,"_dars_filt.csv", sep=""), row.names = F);
  ## count of all unique DARs
  df = as.data.frame(cbind(nrow(up_mod),name_up));
  df1 = as.data.frame(cbind(nrow(down_mod),name_down));
  nums = cbind(df,df1)
  write.csv(nums, file= paste0(dardir,"gene_list/csv/raw_all/",label,"_dars_all_unique_filt_count.csv", sep=""), row.names = F); ##count dars collapsed to genes
  ## count of joint unique dars (indiv)
  gene_list = rbind(up_mod,down_mod);
  colnames(gene_list) <- c("Gene","Cell");
  saveRDS(gene_list, file= paste0(dardir,"gene_list/raw_all/",label,"_dars_joint_filt.rds", sep=""));
  write.csv(gene_list, file= paste0(dardir,"gene_list/csv/raw_all/joint/",label,"_dars_joint_filt.csv", sep=""), row.names = F);
  ##########################################################################################################################################################################
  # source function hard coded above mouse to human genes
  #%%%%%%%%%%%%%%%%  RUN GO AND MAKE BUBBLE PLOTS    %%%%%#############################################################
  ################################################################################  GO terms ##################################################  
  ######################## create converted gene list ######################################################################
  # clean
  data_frame_mod <- closest_up[closest_up$gene_biotype %in% c('processed_transcript','protein_coding'),]
  data_frame_mod <- data_frame_mod %>% distinct(gene_name, .keep_all = TRUE)
  gene_list <-   data_frame_mod %>% dplyr::select(gene_name, Cell); 
  colnames(gene_list) <- c("Gene","Cell")
  up_gene_list = gene_list;
  saveRDS(gene_list, file= paste0(dardir,"gene_list/",name_up,"_dars_filtP.rds", sep=""));
  write.csv(gene_list, file= paste0(dardir,"gene_list/csv/",name_up,"_dars_filtP.csv", sep=""), row.names = F);
  ######
  # convert
  cl <- convertMouseToHuman(data_frame_mod[,2])
  # make df
  go <- as.data.frame(cl)
  #add column
  #go$Cell <- c(label)
  go$Cell <- c(paste0(cell,"_OPEN"))
  ##########################
  ## name columns as required
  colnames(go) <- c("Gene","Cell")
  #save
  go %>% write.csv(file=paste0(godir,'go_raw/up/',cell,'_up.csv'),row.names = F)
  ## assign variable name for later joining
  up_go_genes <- go;
  
  ##script required objs
  dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
           'GO_Molecular_Function_2018')
  dge <- go;
  collapsed_output <- data.frame()
  for(cur in as.character(unique(dge$Cell))){
    print(cur)
    # select genes
    cur_genes <- dge %>%
      subset(Cell == cur) %>% .$Gene
    
    # run enrichR on different gene sets:
    cur_result <- enrichr(cur_genes, dbs)
    
    # collapse results into one dataframe
    for(db in dbs){
      cur_result[[db]]$cluster <- cur
      cur_result[[db]]$db <- db
      #cur_result[[db]]$Diagnosis <- 'cKO'
      collapsed_output <- rbind(collapsed_output, cur_result[[db]])
      collapsed_output <- filter(collapsed_output, Adjusted.P.value < 0.05 )
    }
  }
  
  collapsed_output %>% write.csv(file=paste0(godir,cur,"_go_up.csv"), row.names = F)
  ###############%%         PLOTS 
  #### bar plots)
  input_bub <- collapsed_output %>% filter(cluster %in% c(cur),db == "GO_Molecular_Function_2018") %>% mutate(log = -log10(P.value)) %>% group_by(cluster);
  if (nrow(input_bub) > 0) {  
    plotEnrich(input_bub, showTerms = 25, numChar = 80, y = "Count", orderBy = "log")+ggtitle(name_up);
    ggsave(paste0(godir,"plots/ENRICHR_GO_MF_",cur,"_up.pdf"), width = 10, height = 5.5)
  }
  rm(input_bub)
  input_bub <- collapsed_output %>% filter(cluster %in% c(cur),db == "GO_Biological_Process_2018") %>% mutate(log = -log10(P.value)) %>% group_by(cluster);
  if (nrow(input_bub) > 0) {  
    plotEnrich(input_bub, showTerms = 25, numChar = 80, y = "Count", orderBy = "log")+ggtitle(name_up);
    ggsave(paste0(godir,"plots/ENRICHR_GO_BP_",cur,"_up.pdf"), width = 10, height = 5.5)
  }
  #######
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Molecular_Function_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(6,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF");
  if (nrow(input_bub) > 0) {
    ggballoonplot(input_bub, x = "cluster", y = "Term2", size = "log", fill = "log") +
      scale_fill_gradientn(colors = colors) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      guides(size = FALSE)
    ggsave(paste0(godir,"plots/ENRICHR_GO_bubblechartMF_",cur,"_up.pdf"), width = 4, height = 5)
  }
  
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(6,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  if (nrow(input_bub) > 0) {
    ggballoonplot(input_bub, x = "cluster", y = "Term2", size = "log", fill = "log") +
      scale_fill_gradientn(colors = colors) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      guides(size = FALSE)
    ggsave(paste0(godir,"plots/ENRICHR_GO_bubblechartBP_",cur,"_up.pdf"), width = 4, height = 5)
  }
  
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Cellular_Component_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(5,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  if (nrow(input_bub) > 0) {
    ggballoonplot(input_bub, x = "cluster", y = "Term2", size = "log", fill = "log") +
      scale_fill_gradientn(colors = colors) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      guides(size = FALSE)
    ggsave(paste0(godir,"plots/ENRICHR_GO_bubblechartCellC_",cur,"_up.pdf"), width = 4, height = 5)
  }
  ############################################################################################################################################
  #rm(input_bub)
  #rm(collapsed_output); rm(go); rm(data_frame_mod);rm(df);rm(df1);rm(gene_list); rm(cl);
  ################################################################################  next ##################################################  
  ######################## create converted gene list #####################################################################################
  ###########################################################################################################################################
  # clean
  data_frame_mod <- closest_down[closest_down$gene_biotype %in% c('processed_transcript','protein_coding'),]
  data_frame_mod <- data_frame_mod %>% distinct(gene_name, .keep_all = TRUE)
  gene_list <-   data_frame_mod %>% dplyr::select(gene_name, Cell); 
  colnames(gene_list) <- c("Gene","Cell")
  down_gene_list = gene_list;
  saveRDS(gene_list, file= paste0(dardir,"gene_list/",name_down,"_dars_filtP.rds", sep=""));
  write.csv(gene_list, file= paste0(dardir,"gene_list/csv/",name_down,"_dars_filtP.csv", sep=""), row.names = F);
  both_dars_clean = rbind(up_gene_list, down_gene_list);
  saveRDS(both_dars_clean, file= paste0(dardir,"gene_list/",label,"_dars_clean_genes_",count,".rds", sep=""));
  write.csv(both_dars_clean, file= paste0(dardir,"gene_list/csv/",label,"_dars_clean_genes_",count,".csv", sep=""), row.names = F);
  ## count of all unique DARs
  df = as.data.frame(cbind(nrow(up_gene_list),name_up));
  df1 = as.data.frame(cbind(nrow(down_gene_list),name_down));
  nums = cbind(df,df1)
  write.csv(nums, file= paste0(dardir,"gene_list/csv/",label,"_dars_count_protein_coding.csv", sep=""), row.names = F); 
  
  # convert
  cl <- convertMouseToHuman(data_frame_mod[,2])
  # make df
  go <- as.data.frame(cl)
  #add column
  #go$Cell <- c(label)
  go$Cell <- c(paste0(cell,"_CLOSED"))
  ##########################
  ## name columns as required
  colnames(go) <- c("Gene","Cell")
  #save
  go %>% write.csv(file=paste0(godir,'go_raw/down/',cell,'_down.csv'),row.names = F)
  ## assign variable name for later joining
  down_go_genes <- go;
  joint_go_genes <- rbind(up_go_genes, down_go_genes)
  saveRDS(joint_go_genes, file= paste0(dardir,"gene_list/human_",label,"_dars.rds", sep=""));
  write.csv(joint_go_genes, file= paste0(dardir,"gene_list/csv/human/human_",label,"_dars.csv", sep=""), row.names = F);
  ## count of all unique DARs
  df = as.data.frame(cbind(nrow(up_go_genes),name_up));
  df1 = as.data.frame(cbind(nrow(down_go_genes),name_down));
  nums = cbind(df,df1)
  write.csv(nums, file= paste0(dardir,"gene_list/csv/human_",label,"_dars_count_go_hu.csv", sep=""), row.names = F);
  
  ######
  ##script required objs
  dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
           'GO_Molecular_Function_2018')
  dge <- go;
  collapsed_output <- data.frame()
  for(cur in as.character(unique(dge$Cell))){
    print(cur)
    # select genes
    cur_genes <- dge %>%
      subset(Cell == cur) %>% .$Gene
    
    # run enrichR on different gene sets:
    cur_result <- enrichr(cur_genes, dbs)
    
    # collapse results into one dataframe
    for(db in dbs){
      cur_result[[db]]$cluster <- cur
      cur_result[[db]]$db <- db
      #cur_result[[db]]$Diagnosis <- 'cKO'
      collapsed_output <- rbind(collapsed_output, cur_result[[db]])
      collapsed_output <- filter(collapsed_output, Adjusted.P.value < 0.05 )
    }
  }
  
  collapsed_output %>%  write.csv(file=paste0(godir,cur,"_go_down",'.csv'), row.names = F)
  
  ###############%%  PLOTS 
  ### bar plot
  input_bub <- collapsed_output %>% filter(cluster %in% c(cur),db == "GO_Molecular_Function_2018") %>% mutate(log = -log10(P.value)) %>% group_by(cluster);
  if (nrow(input_bub) > 0) {  
    plotEnrich(input_bub, showTerms = 25, numChar = 80, y = "Count", orderBy = "log")+ggtitle(name_down);
    ggsave(paste0(godir,"plots/ENRICHR_GO_MF_",cur,"_down.pdf"), width = 10, height = 5.5)
  }
  rm(input_bub);
  input_bub <- collapsed_output %>% filter(cluster %in% c(cur),db == "GO_Biological_Process_2018") %>% mutate(log = -log10(P.value)) %>% group_by(cluster);
  if (nrow(input_bub) > 0) {  
    plotEnrich(input_bub, showTerms = 25, numChar = 80, y = "Count", orderBy = "log")+ggtitle(name_down);
    ggsave(paste0(godir,"plots/ENRICHR_GO_BP_",cur,"_down.pdf"), width = 10, height = 5.5)
  }
  ###
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Molecular_Function_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(5,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF");
  
  if (nrow(input_bub) > 0) {
    ggballoonplot(input_bub, x = "cluster", y = "Term2", size = "log", fill = "log") +
      scale_fill_gradientn(colors = colors) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      guides(size = FALSE)
    ggsave(paste0(godir,"plots/ENRICHR_GO_bubblechartMF_",cur,"_down.pdf"), width = 4, height = 5)
  }
  
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(5,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  if (nrow(input_bub) > 0) {
    ggballoonplot(input_bub, x = "cluster", y = "Term2", size = "log", fill = "log") +
      scale_fill_gradientn(colors = colors) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      guides(size = FALSE)
    ggsave(paste0(godir,"plots/ENRICHR_GO_bubblechartBP_",cur,"_down.pdf"), width = 4, height = 5)
  }
  
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Cellular_Component_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(5,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  if (nrow(input_bub) > 0) {
    ggballoonplot(input_bub, x = "cluster", y = "Term2", size = "log", fill = "log") +
      scale_fill_gradientn(colors = colors) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      guides(size = FALSE)
    ggsave(paste0(godir,"plots/ENRICHR_GO_bubblechartCellC_",cur,"_down.pdf"), width = 4, height = 5)
  }
  
  
  print(Sys.time());
  ended = Sys.time();
  time_passed = ended - started;   
  print(time_passed);
  ############################################################################################################################################
  i = i+ 2;
  print(paste0("The new i is: ", i))
  #invisible(readline(prompt="Press [enter] to continue"))
}




