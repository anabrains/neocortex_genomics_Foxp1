##-------------------------------------------------------
## LOAD MODULES AND LIBRARIES
#-------------------------------------------------------
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.1.1-gccmkl


rm(list = ls())
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Matrix.utils))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MAST))
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations
library(BiocParallel)
library(magrittr)
library(broom)
library(cowplot)
library(BiocSingular)
library(clusterProfiler)
library(enrichR)

basedir = "./sc-atac/aoa/plots_rna/degs_rna/cko_first/e13/"
savedir = "./sc-atac/aoa/plots_rna/degs_rna/cko_first/e13/"
savedir_up = paste0(basedir,"up_in_cko/")
savedir_down = paste0(basedir,"down_in_cko/")
#rawdir = "./ao/data_pipeline1/mergee13/exci_lin/degs/MAST_GLM/raw_table/"
plotrna = "./sc-atac/aoa/plots_rna/volcano/"
dgerna = "./sc-atac/aoa/plots_rna/degs_rna/cko_first/e13/"
godir = paste0(basedir, "go_dir/")
goplot= paste0(basedir, "go_dir/plots/")

##-------------------------------------------------------
### Test Data
###-------------------------------------------------------
e13 <- readRDS("./ao/data_pipeline1/mergee13/exci_lin/e13_exci_lin_celltypes_obj.rds")
#e16 <- readRDS("./ao/data_pipeline2/merge/exci_lin/e16_exci_lin_celltypes_obj.rds")
#p0 <- readRDS("./ao/data_pipeline1/mergep0/exci_lin/p0_exci_lin_cellTypes_obj.rds")
###

##-------------------------------------------------------
### Source Functions
###-------------------------------------------------------
source("./sc-atac/aoa/convertMouseToHuman.r")
##

##-------------------------------------------------------
### DEG | PSEUDOBULK
###-------------------------------------------------------
seuObjFilt <- e13

# seuObjFilt

Idents(seuObjFilt) <- "CellType_coll"
print(table(seuObjFilt@active.ident))


for (cell in unique(sort(seuObjFilt@active.ident)))
{
  ## select a cell-type
  cluSelected <- cell
  print(cluSelected)
  
  ## subset for a cell-type
  celltype.chosen <- subset(x = seuObjFilt, idents = c(cluSelected))
  celltype.chosen <- NormalizeData(celltype.chosen)
  Idents(celltype.chosen) <- "genotype"
  
  ## Prepare data for MAST ##
  mat <- celltype.chosen@assays$RNA@data
  meta <- celltype.chosen@meta.data
  
  ## Keep genes with 10% expression in at least one species
  cells1 <- WhichCells(celltype.chosen, idents = 'CKO') #was ctrl made cko
  cells2 <- WhichCells(celltype.chosen, idents = 'CTRL') #was cko made ctrl
  pass1 = rowSums(mat[,cells1])/length(cells1) > 0.1
  pass2 = rowSums(mat[,cells2])/length(cells2) > 0.1
  mat = mat[pass1 | pass2,]
  
  ## Create MAST object
  sca <- MAST::FromMatrix(exprsArray = as.matrix(x = mat),
                          cData = meta,
                          fData = data.frame(rownames(mat)))
  
  ## Scale number of detected genes (cngeneson. Same as nFeature_RNA only after filtering)
  cdr2 <- colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  ## Calculate fold change
  avg_logfc <- log(rowMeans(expm1(mat[,cells1])) + 1) - log(rowMeans(expm1(mat[,cells2])) + 1)
  
  ## MAST package with covariates and glm
  mastfix = MAST::zlm(~genotype + cngeneson + batch, sca, ebayes = F, method = 'glm')
  
  summaryCond <- summary(object = mastfix, doLRT = 'genotypeCTRL')
  summaryDt <- summaryCond$datatable
  #write.table(summaryDt,file=paste0(rawdir,cell,"_raw_DGE_MAST_GLM.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  #save(summaryDt,file=paste0(rawdir,cell,"_raw_DGE_MAST_GLM.rdata"))
  #get relevant info
  p_val <- summaryDt[summaryDt$component == "H", 4]
  genes.return <- summaryDt[summaryDt$component == "H", 1]
  to.return <- data.frame(p_val, row.names = genes.return$primerid)
  colnames(to.return)[1] = 'p_value'
  fix_res = to.return
  fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')
  fix_res$avg_logfc = avg_logfc
  fix_res <- rownames_to_column(fix_res, var = paste0(cell,"_gene"))  #tibble
  
  #print(sum(abs(fix_res$avg_logfc) >= 0.1375 & fix_res$adj_p_value <= 0.05))
  mast_allcov_glm <- subset(fix_res, adj_p_value <0.05 & abs(avg_logfc) > 0.05)
  
  up_degs <- subset(fix_res, adj_p_value <0.05 & avg_logfc > 0.15)  
  down_degs <- subset(fix_res, adj_p_value <0.05 & avg_logfc < -0.15)
  
  ## write dge table
  write.table(mast_allcov_glm, paste0(savedir,"DGE_MASTGLM_TABLE_",cell,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(up_degs, paste0(savedir_up,"DGE_MASTGLM_TABLE_UP_",cell,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(down_degs, paste0(savedir_down,"DGE_MASTGLM_TABLE_DOWN_",cell,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  ##-------------------------------------------------------
  ### DEG | Volcano Plot
  ###-------------------------------------------------------
  ###### make volcano plot############################33
  df <- mast_allcov_glm;
  #### create label for auto naming
  name <- colnames(df[1])
  #remove anything after last period
  name <- sub("^(.*)[.].*", "\\1", name)
  #age
  age <- gsub("^.*\\.","", colnames(df[1]))
  age <- gsub('(.*)_\\w+', '\\1', age)
  ## Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
  df$diffexpressed <- "NOT_DEG"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  df$diffexpressed[df$avg_logfc < -0.15 & df$adj_p_value < 0.05] <- "DOWN_REG"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  df$diffexpressed[df$avg_logfc > 0.15 & df$adj_p_value < 0.05] <- "UP_REG"
  # Explore a bit
  head(df[order(df$adj_p_value) & df$diffexpressed == 'DOWN_REG', ])
  colnames(df)[1] = "Gene"
  ####ok try to have names ####################################################
  de <- df
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NOT_DEG"] <- de$Gene[de$diffexpressed != "NOT_DEG"]
  ## set legibilitly theme # Biostatsquid theme
  theme_set(theme_classic(base_size = 20) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))
  ##generate
  p1 <- ggplot(data=de, aes(x=avg_logfc, y=-log10(adj_p_value), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    #theme_minimal() +
    geom_text_repel(size=4.5) +
    scale_color_manual(values=c("#6699ff","grey","#ff8c66")) +
    geom_vline(xintercept=c(-0.15, 0.15), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept=-log10(0.05),col = "gray", linetype = 'dashed') +
    labs(color = 'Key', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # geom_text(size = 4) + 
    ggtitle(paste('DEGs in', name, "at",age))# + scale_x_continuous(breaks = seq(-.8, .5, .15))
  ##save it
  p1
  ggsave(file = paste0(plotrna,age,"_",name,"_volcano.pdf"), width = 16, height = 14) # you can change the size of the output file

  ##-------------------------------------------------------
  ### DEG | End Volcano Plot
  ###-------------------------------------------------------
  
  ##-------------------------------------------------------
  ### DEG | Convert to Human for GO
  ###-------------------------------------------------------
  
  ### convert mouse to human, save for GO
  up_degs <- subset(fix_res, adj_p_value <0.05 & avg_logfc > 0.15)  
  down_degs <- subset(fix_res, adj_p_value <0.05 & avg_logfc < -0.15)
  
  gene_list <- up_degs[,1]
  gene_list_mod <- convertMouseToHuman(gene_list) #if there's col name, remove it [-1,]
  gene_list_mod_up <- as.data.frame(gene_list_mod)
  
  gene_list <- down_degs[,1]
  gene_list_mod <- convertMouseToHuman(gene_list) #cuz you dont want column name so had to get rid of it
  gene_list_mod_down <- as.data.frame(gene_list_mod)
  ##create names
  #cell <- gsub('(.*)_\\w+', '\\1', colnames(df[1]))
  #########
  ##########################  prep columes
  #add column
  gene_list_mod_up$Cell <- c(cell)
  ## name columns as required
  colnames(gene_list_mod_up) <- c("Gene","Cell")
  ###down
  #add column
  gene_list_mod_down$Cell <- c(cell)
  ## name columns as required
  colnames(gene_list_mod_down) <- c("Gene","Cell")
  ############# try to just bind up list to down list #######################3
  gene_list_mod_up1 <- gene_list_mod_up;
  gene_list_mod_down1 <- gene_list_mod_down;
  gene_list_mod_up1$Cell <- c(paste0(cell,"_up"))
  gene_list_mod_down1$Cell <- c(paste0(cell,"_down"))
  colnames(gene_list_mod_down1) <- c("Gene","Cell")
  com_list <- rbind(gene_list_mod_up1,gene_list_mod_down1)
  write.csv(com_list,file=paste0(dgerna,name,"_",age,'_degs.csv'), row.names = F)
  write.csv(gene_list_mod_up1,file=paste0(savedir_up,name,"_",age,'_degs_up_cko.csv'), row.names = F)
  write.csv(gene_list_mod_down1,file=paste0(savedir_down,name,"_",age,'_degs_down_cko.csv'), row.names = F)
  ###########################  
  ############################   script required objs
  dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
           'GO_Molecular_Function_2018')
  ########################################  what is your input ???
  dge <- gene_list_mod_up1;
  ################################################### run
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
  ### save
  collapsed_output %>%
    write.csv(file=paste0(godir,cur,'_cko_go.csv'), row.names = F)
  df_up <- collapsed_output;
  ###############%%%%%%%%        plot   %%%%%%%##############################################
  ###########################################################################################
  colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF");
  # input_bub <- collapsed_output %>% 
  #   filter(cluster %in% c(cur),db == "GO_Molecular_Function_2018") %>% 
  #   mutate(log = -log10(P.value)) %>%
  #   group_by(cluster) %>%
  #   top_n(5,log) %>%
  #   mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
  #   mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  # 
  # p2 <- ggballoonplot(input_bub, x = "cluster", y = "Term2",
  #                     size = "log", fill = "log") +
  #   scale_fill_gradientn(colors = colors) +
  #   scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  #   guides(size = FALSE)
  # p2;
  # ggsave(paste0(goplot,"ENRICHR_GO_bubblechartMolF_UP_",cur,".pdf"), width = 4, height = 5)

  ## next
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(5,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  p2 <- ggballoonplot(input_bub, x = "cluster", y = "Term2",
                      size = "log", fill = "log") +
    scale_fill_gradientn(colors = colors) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
    guides(size = FALSE)
  p2;
  ggsave(paste0(goplot,"ENRICHR_GO_bubble_Bio_UP_",cur,".pdf"), width = 4, height = 5);

  ########################################  what is your input ???
  dge <- gene_list_mod_down1;
  ################################################### run
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
  
  ### save
  collapsed_output %>%
    write.csv(file=paste0(godir,cur,'_cko_go.csv'), row.names = F)
  ################################################################
  ###############%%%%%%%%        plot   %%%%%%%##############################################
  ###########################################################################################
  colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF");
  # input_bub <- collapsed_output %>% 
  #   filter(cluster %in% c(cur),db == "GO_Molecular_Function_2018") %>% 
  #   mutate(log = -log10(P.value)) %>%
  #   group_by(cluster) %>%
  #   top_n(5,log) %>%
  #   mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
  #   mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  # 
  # p2 <- ggballoonplot(input_bub, x = "cluster", y = "Term2",
  #                     size = "log", fill = "log") +
  #   scale_fill_gradientn(colors = colors) +
  #   scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  #   guides(size = FALSE)
  # 
  # p2;
  # ggsave(paste0(goplot,"ENRICHR_GO_bubblechartMolF_down_",cur,".pdf"), width = 4, height = 5);
  # 
  ## next
  input_bub <- collapsed_output %>% 
    filter(cluster %in% c(cur),db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(5,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))
  
  p2 <- ggballoonplot(input_bub, x = "cluster", y = "Term2",
                      size = "log", fill = "log") +
    scale_fill_gradientn(colors = colors) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
    guides(size = FALSE)
  p2;
  ggsave(paste0(goplot,"ENRICHR_GO_bubble_Bio_down_",cur,".pdf"), width = 4, height = 5);
 
  
  ##################################################################
  ##########%%% joint DGE for joint list test   %%%#################
  ########################################3 jint list for joint terms
  dge <- com_list;
  ################################################### run
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
  
  ### save
  collapsed_output %>%
    write.csv(file=paste0(godir,cur,'_all_go.csv'), row.names = F)
  #### end
  
}


##-------------------------------------------------------
##-------------------------------------------------------
## END
##-------------------------------------------------------
##-------------------------------------------------------
