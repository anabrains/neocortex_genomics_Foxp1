
## USAGE
#  Generates individual replicate seurat objects
#  does QC 
#  Runs DoubletFinder
#  Removes overlapping doublets between DoubletFinder and Scrublet
# Be sure to run `module load R/3.5.1-gccmkl` before running this script.
#
# Execute it using:
# 
#   Rscript ./seurat_object_creation.r
#
#to trouble shoot you can work in rstudio by: module load rstudio-desktop/1.1.456 , then: rstudio

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(cowplot)
library(DoubletFinder)
library(stringr)
library(tidyverse)

load_seurat_object <- function(target_folder) {
    # Find the H5 file in the target_folder
    f = list.files(path=target_folder, pattern="\\.h5$", full.names=FALSE)
    if (length(f) > 0) {
        curr.counts <- Read10X_h5(file.path(target_folder, f[1]))
        curr <- CreateSeuratObject(counts=curr.counts) 
        curr
    } else {
        stop(paste0("Could not find .h5 file in ", target_folder))
    }
}

gen_qc_plots <- function(seu_obj, curr_mouse_id, results_folder) {
    curr_qcplot = VlnPlot(seu_obj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
    qc_plot_filename = paste("QC_Plot_", curr_mouse_id, ".pdf", sep="")
    ggsave(file.path(results_folder, qc_plot_filename), plot = curr_qcplot, width=26, height=12, units= "in", dpi=300)
}

append_scrublet_results <- function(seu_obj, results_folder) {
    # open scrublet results
    curr_scrublet_results = read.csv(file.path(results_folder, 'scrublet_results.csv'))

    # add the doublet_score
    score = curr_scrublet_results[,1];
    score = t(score);
    score = as.numeric(score);
    seu_obj = AddMetaData(seu_obj, metadata=score, col.name="doublet_score")

    # now the predicted_doublet
    dub = curr_scrublet_results[,2];
    seu_obj = AddMetaData(seu_obj, metadata=dub, col.name="predicted_doublet")

    seu_obj
}

gen_ribo_plot <- function(seu_obj, curr_mouse_id, results_folder) {
    curr_ribo_plot = VlnPlot(seu_obj, c("percent.ribo"), pt.size = 0)
    ribo_plot_filename = paste("Ribo_Plot_", curr_mouse_id, ".pdf", sep='')
    ggsave(file.path(results_folder, ribo_plot_filename), plot=curr_ribo_plot, width=12, height=8, units="in", dpi=300)
}

gen_feature_plot <- function(seu_obj, curr_mouse_id, results_folder) {
    curr_feature_plot = FeatureScatter(seu_obj, feature1="nCount_RNA", feature2="percent.ribo")
    feature_plot_filename = paste("Ribo_Feature_Plots_", curr_mouse_id, ".pdf", sep='')
    ggsave(file.path(results_folder, feature_plot_filename), plot=curr_feature_plot, width=12, height=8, units="in", dpi=300)
}

gen_rqc_plot <- function(seu_obj, curr_mouse_id, results_folder) {
    curr_rqc_plot <- VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size=0)
    rqc_plot_filename = paste("rQC_Plot_", curr_mouse_id, ".pdf", sep='')
    ggsave(file.path(results_folder, rqc_plot_filename), plot=curr_rqc_plot, width=26, height=12, units="in", dpi=300)
}

gen_bin_qc_plot <- function(seu_obj, curr_mouse_id, results_folder) {
    brks <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000)

    umiBins <- cut(seu_obj$nCount_RNA, brks, include.lowest = T, right = FALSE)
    binData <- cbind(summary(umiBins))
    colnames(binData) <- c("nUMI")
    binData <- as.data.frame(binData)
    binData$bins <- row.names(binData)
    binData2 <- melt(binData)
    binData2$bins <- factor(binData2$bins, levels = row.names(binData))
    plotBins <- ggplot(binData2, aes(bins, value, fill=variable)) + geom_bar(stat="identity", position="dodge") +  facet_grid(. ~ variable) + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1), legend.position="none") +  geom_text(aes(label = value), vjust = -1) + ylim(0, 5000) + labs(title=paste("Histogram of", curr_mouse_id), x="Bins", y="Number of Cells") 
    umi_bin_plot_filename = paste("QC", curr_mouse_id, "nUMI_BINS.pdf", sep='_')
    ggsave(file.path(results_folder, umi_bin_plot_filename), plot=plotBins, width=12, height=8, units="in", dpi=300)

    umiBins <- cut(seu_obj$nFeature_RNA, brks, include.lowest = T, right = FALSE)
    binData <- cbind(summary(umiBins))
    colnames(binData) <- c("nGene")
    binData <- as.data.frame(binData)
    binData$bins <- row.names(binData)
    binData2 <- melt(binData)
    binData2$bins <- factor(binData2$bins, levels = row.names(binData))
    plotBins <- ggplot(binData2, aes(bins, value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") +  facet_grid(. ~ variable) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +  geom_text(aes(label = value), vjust = -1) + ylim(0, 5000) + labs(title=paste("Histogram of", curr_mouse_id), x = "Bins", y = "Number of Cells") 
    gene_bin_plot_filename = paste("QC", curr_mouse_id, "Gene_BINS.pdf", sep='_')
    ggsave(file.path(results_folder, gene_bin_plot_filename), plot=plotBins, width=12, height=8, units="in", dpi=300)
}
gen_qc_plot_clust <- function(seu_obj, curr_mouse_id, results_folder) {
  plotlist <- VlnPlot(object = seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3, pt.size = 0, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist, nrow = 3);
  title <- ggdraw() + draw_label(paste("Plot of", curr_mouse_id,"_clusters_",curr_age, curr_genotype), fontface = 'bold');
  curr_qcplot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  
  qc_plot_filename = paste("QC_Plot_", curr_mouse_id, "_clusters.pdf", sep="")
  ggsave(file.path(results_folder, qc_plot_filename), plot = curr_qcplot, width=26, height=12, units= "in", dpi=300)
}
#ribo included
gen_rqc_plot_clust <- function(seu_obj, curr_mouse_id, results_folder) {
  plotlist <- VlnPlot(object = seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 4, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist);
  title <- ggdraw() + draw_label(paste("Plot of ",curr_mouse_id,"_clusters_",curr_age, curr_genotype), fontface = 'bold');
  curr_rqc_plot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  
  rqc_plot_filename = paste("rQC_Plot_", curr_mouse_id, "_clusters.pdf", sep='')
  ggsave(file.path(results_folder, rqc_plot_filename), plot=curr_rqc_plot, width=26, height=12, units="in", dpi=300)
}

gen_qc_plot_sub <- function(seu_obj, curr_mouse_id, results_folder) {
  plotlist <- VlnPlot(object = seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3, pt.size = 0, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist, nrow = 3);
  title <- ggdraw() + draw_label(paste("Plot of", curr_mouse_id,"_subset_doublets_",curr_age, curr_genotype), fontface = 'bold');
  curr_qcplot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  
  qc_plot_filename = paste("QC_Plot_", curr_mouse_id, "_subset_doublets.pdf", sep="")
  ggsave(file.path(results_folder, qc_plot_filename), plot = curr_qcplot, width=26, height=12, units= "in", dpi=300)
}
#ribo included
gen_rqc_plot_sub <- function(seu_obj, curr_mouse_id, results_folder) {
  plotlist <- VlnPlot(object = seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 4, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist);
  title <- ggdraw() + draw_label(paste("Plot of ",curr_mouse_id,"subset_doublets_",curr_age, curr_genotype), fontface = 'bold');
  curr_rqc_plot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  rqc_plot_filename = paste("rQC_Plot_", curr_mouse_id, "_subset_doublets.pdf", sep='')
  ggsave(file.path(results_folder, rqc_plot_filename), plot=curr_rqc_plot, width=26, height=12, units="in", dpi=300)
}
gen_dub_plot <- function(seu_obj, curr_mouse_id, results_folder) {
  pANN_class_list <- colnames(seu_obj@meta.data)[grepl("pANN_", colnames(seu_obj@meta.data))]
  plotlist <- VlnPlot(object = seu_obj, features = c("doublet_score", pANN_class_list[1], "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist);
  title <- ggdraw() + draw_label(paste("Plot of ",curr_mouse_id,"plot_doublets_",curr_age, curr_genotype), fontface = 'bold');
  curr_dub_plot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  
  dub_plot_filename = paste("Doublet_Plot_", curr_mouse_id, "_clusters.pdf", sep='')
  ggsave(file.path(results_folder, dub_plot_filename), plot=curr_dub_plot, width=26, height=12, units="in", dpi=300)
}
gen_dub_plot1 <- function(seu_obj, curr_mouse_id, results_folder) {
  pANN_class_list <- colnames(seu_obj@meta.data)[grepl("pANN_", colnames(seu_obj@meta.data))]
  plotlist <- VlnPlot(object = seu_obj, features = c("doublet_score", pANN_class_list[1], "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist);
  title <- ggdraw() + draw_label(paste("Plot of ",curr_mouse_id,"doublets-post-subset_",curr_age, curr_genotype), fontface = 'bold');
  curr_dub_plot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  
  dub_plot_filename = paste("Doublet_Plot_", curr_mouse_id, "_clusters_post-subset.pdf", sep='')
  ggsave(file.path(results_folder, dub_plot_filename), plot=curr_dub_plot, width=26, height=12, units="in", dpi=300)
}
##you need the list of estimated doublets from doubletFinder in the pANN ##here works if there are two elements
gen_dub_plot2 <- function(seu_obj, curr_mouse_id, results_folder) {
  pANN_class_list <- colnames(seu_obj@meta.data)[grepl("pANN_", colnames(seu_obj@meta.data))]
  plotlist <- VlnPlot(object = seu_obj, features = c("doublet_score", pANN_class_list[2], "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, combine = F);
  p1 <- cowplot::plot_grid(plotlist = plotlist);
  title <- ggdraw() + draw_label(paste("Plot of ",curr_mouse_id,"doublets-post-subset-k2_",curr_age, curr_genotype), fontface = 'bold');
  curr_dub_plot<- cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1));
  
  dub_plot_filename = paste("Doublet_Plot_", curr_mouse_id, "_clusters_post-subset-k2.pdf", sep='')
  ggsave(file.path(results_folder, dub_plot_filename), plot=curr_dub_plot, width=26, height=12, units="in", dpi=300)
}

# Copy from https://stackoverflow.com/a/50665893
# Allows us to save objects with a new name defined by a string
saveit <- function(..., string, file) {
  x <- list(...)
  names(x) <- string
  save(list=names(x), file=file, envir=list2env(x))
}

start_dir <- "./ao/data_pipeline2"
data_path <- "cellbender"
filtered_results_folder <- file.path(start_dir,"doublets_removed_seu_objs");


working_files <- list.files(path=start_dir, recursive=FALSE, full.names=FALSE)

print("WELCOME TO <SEURAT> :D")

# Create a Seurat object from matrix data
for (working_folder_name in working_files){

    # create the full path to check
    data_path_folder = file.path(start_dir, file.path(working_folder_name, data_path))
    if (file.exists(data_path_folder)) {
        print(data_path_folder)
        # Filename format is: <Mouse ID>_<Genotype>_<Age>_<Batch>
        split_name = unlist(strsplit(working_folder_name, '_'))
        curr_mouse_id = split_name[1]
        curr_genotype = split_name[2]
        curr_age = split_name[3]
        curr_batch = split_name[4]

        results_folder <- file.path(start_dir, working_folder_name)

        print("[X] Loading Seurat Object")
        curr <- load_seurat_object(data_path_folder)
        
        # mito and calculate percent ribo
        print("[X] mito and calculating percent ribo")
        curr[["percent.mt"]] <- PercentageFeatureSet(curr, pattern="^mt-")
        curr[["percent.ribo"]] <- PercentageFeatureSet(curr, pattern="^Rp")

        # qcPlots
        print("[X] Generating QC Plots")
        gen_qc_plots(curr, curr_mouse_id, results_folder)

        ## Append metadata
        print("[X] Appending Metadata")
        curr <- append_scrublet_results(curr, results_folder)

        ## Plots
        print("[X] Generating Ribo Plot")
        gen_ribo_plot(curr, curr_mouse_id, results_folder)
        print("[X] Generating Feature Plot")
        gen_feature_plot(curr, curr_mouse_id, results_folder)
        print("[X] Generating RQC Plot")
        gen_rqc_plot(curr, curr_mouse_id, results_folder)

        # Save intermediate response
        print("[X] Saving Intermediate Results")
        save(curr, file=file.path(results_folder, "intermediate_metadatas.rdata"))

        ## QC bins with old variable that doesn't have the ribo data just bc I wanted to save before moving on
        print("[X] Generating Bin QC Plots")
        gen_bin_qc_plot(curr, curr_mouse_id, results_folder)
        ## Bin QC DONE

        # either subset based on info above, exclude above 5% or do Stephs was and exclude about .05 fractional per next section

        curr_subset <- subset(x=curr, subset=nCount_RNA< 10000 & nFeature_RNA >500 & percent.mt<5 & percent.ribo<5); 

        mGenes <- scan("./FILTER_MITO_X_Y/GENCODE_vM17_MM10_ChrM_GENES.txt", what = "", sep = "\n");
        xGenes <- scan("./FILTER_MITO_X_Y/GENCODE_vM17_MM10_ChrX_GENES.txt", "./FILTER_MITO_X_Y/cellranger_mm10_genes_ChrX.txt", what = "", sep = "\n");
        xGenes1 <- scan("./FILTER_MITO_X_Y/cellranger_mm10_genes_ChrX.txt", what = "", sep = "\n");
        yGenes <- scan("./FILTER_MITO_X_Y/GENCODE_vM17_MM10_ChrY_GENES.txt", what = "", sep = "\n");
        mxyGenes <- unique(sort(c(mGenes, xGenes, xGenes1, yGenes))) #4670 genes

        print("[X] Filtering out X,Y, mito genes")
        
        keepGenes <- unique(sort(setdiff(row.names(curr_subset), mxyGenes)));
        curr_filtered <- GetAssayData(object=curr_subset);
        curr_filtered_temp <- as.matrix(curr_filtered);
        keepCells <- colnames(curr_filtered_temp);
        curr_filtered_raw <- GetAssayData(object=curr_subset, slot="counts");
        curr_filtered_raw_temp <- as.matrix(curr_filtered_raw);
        curr.newdata <- curr_filtered_raw_temp[keepGenes,keepCells];
        curr_new <- CreateSeuratObject(counts = curr.newdata, project = paste(curr_mouse_id, "woMT", sep="_"));

        #must add meta data back for % mito and % ribo bc it is lost in the removal of xyz mt genes, so run immediately 
        #after above code b.c. variables are dependent
        metaAll <- as.data.frame(curr@meta.data);
        dim(metaAll);
        AO1pMito <- metaAll[keepCells, "percent.mt"]
        names(AO1pMito) <- row.names(metaAll[keepCells,])
        curr_new$percent.mt <- AO1pMito
        AO1perRibo <- metaAll[keepCells, "percent.ribo"]
        names(AO1perRibo) <- row.names(metaAll[keepCells,])
        curr_new$percent.ribo <- AO1perRibo

        #---- in the filtering step in the pre-process seurat code you need this
        #---- after removing your x,y, etc genes 
        ao1dubScore <- metaAll[keepCells, "doublet_score"]
        names(ao1dubScore) <- row.names(metaAll[keepCells,])
        curr_new$doublet_score <- ao1dubScore

        ao1predDub <- metaAll[keepCells, "predicted_doublet"]
        names(ao1predDub) <- row.names(metaAll[keepCells,])
        curr_new$predicted_doublet <- ao1predDub

        # set custom metadata
        curr_new@meta.data$region <- as.factor("cortex")
        curr_new@meta.data$age <- as.factor(curr_age)
        curr_new@meta.data$genotype <- as.factor(curr_genotype)
        curr_new@meta.data$batch <- as.factor(curr_batch)

        ## Now Run DoubletFinder

        #i pre-process Q/C on seu obj and ran with scrublet, here i use this after adding scrublet in metadata
        ## Pre-process Seurat object -------------------------------------------------------------------------------------------------
        #seu_obj <- CreateSeuratObject(kidney.data)
        print("[X] Running DoubletFinder")
        seu_obj <- curr_new
        seu_obj <- NormalizeData(seu_obj)
        seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
        seu_obj <- ScaleData(seu_obj, vars.to.regress = c("nCount_RNA", "percent_mito"))
        seu_obj <- RunPCA(seu_obj)
        seu_obj <- FindNeighbors(seu_obj, dims = 1:35);
        seu_obj <- FindClusters(seu_obj, resolution = c(1.2));
        seu_obj <- RunUMAP(seu_obj, dims = 1:35)

        ###########################################################################################

        ## pK Identification ---------------------------------------------------------------------------------------------------------
        sweep.res.list_kidney <- paramSweep_v3(seu_obj, PCs = 1:20); #takes a bit
        sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE);
        bcmvn_kidney <- find.pK(sweep.stats_kidney);

        ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
        annotations <- seu_obj@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)## ex: annotations <- seu_obj@meta.data$seurat_clusters
        num_cells = ncol(seu_obj)
        expected_doublet_rate = num_cells/125000 # linear approximation to provided table
        nExp_poi <- round(expected_doublet_rate*nrow(seu_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        #nExp_poi <- round(0.2598*nrow(seu_obj@meta.data)) 

        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

        ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
        bcmvn_kidney #to determine ur pk
        # choose the X with the highest Y
        #pK by animal:  ao1: .001, ao2= .01  ; ao3= 0.05; ao4: 0.05; ao5: 0.23; ao6= .005; ao7 = 0.005; ao8=.005
        ####  auto select you rpK value based on prev function
        pK=as.numeric(as.character(bcmvn_kidney$pK))
        BCmetric=bcmvn_kidney$BCmetric
        pK_choose = pK[which(BCmetric %in% max(BCmetric))]
        print("[X] pK for doubletFinder chosen")
        # plot pK with BCmatric
        pk_name <-paste0(curr_mouse_id,"_pK_choose_plot.pdf")
        pdf(file=pk_name,width=26, height=12)
        par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2);
        plot(x = pK, y = BCmetric, pch = 16,type="b",col = "blue",lty=1)
        abline(v=pK_choose,lwd=2,col='red',lty=2) 
        title(paste("The BCmvn distributions", curr_mouse_id))
        text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
        dev.off()
        
        ######  you need a particular pK, chosen from prev code section
        seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:35, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE);
        
        ##-- now re-use the pANN value that is  added to meta data after running the previous line (optional based on hyp)
        #seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:30, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_obj@meta.data)[grepl("pANN", colnames(seu_obj@meta.data))]);
        #this prev analysis will create a separate pANN slot that can be considered low doublet est

        #optional plotting
        #DF.name = colnames(seu_obj@meta.data)[grepl("DF.classification", colnames(seu_obj@meta.data))];
        #cowplot::plot_grid(ncol = 2, VlnPlot(seu_obj, features = "nFeature_RNA", group.by = scrubDub, pt.size = 0),VlnPlot(seu_obj, features = "nFeature_RNA", group.by = DF.name, pt.size = 0)) 

        ########################################################################################
        ######################---- add information on cell cycle --------#######################
        print("[X] Adding cell cycle information")
        s.genes <- cc.genes$s.genes;
        g2m.genes <- cc.genes$g2m.genes;
        g2m.genes = str_to_title(g2m.genes);
        s.genes = str_to_title(s.genes);
        seu_obj <- CellCycleScoring(seu_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE);
        #the difference between cycling and non cycling #put in le meta data
        seu_obj$CC.Difference <- seu_obj$S.Score - seu_obj$G2M.Score
        print("[X] Cell cycle scores added")
        
        print("[X] Saving object after adding cell cycle info: OG")
        #save this object, move forward
        save_name <- paste(curr_mouse_id, "og", sep="_")
        saveit(seu_obj_sub=seu_obj, string=save_name, file=file.path(results_folder, paste0(curr_mouse_id,"_seu_obj_og.rdata")))
        
        ##############################################################################################
        ##########---------        REMOVE CELLS THAT ARE DETECTED AS DOUBLETS IN BOTH PACAKGES      ---------############
        ########
        print("[X] Removing doublets found in both packages")
        ##cells for scrublet that are predictd to be doublets
        ##scrublet cell witih doublet predicted true
        cells=rownames(seu_obj@meta.data)[which(seu_obj@meta.data$predicted_doublet == "True")]
        
        ##cells that are doublets, in doubletFinder
        #df_cells_remove=rownames(seu_obj@meta.data)[which(seu_obj@meta.data$DF.classifications_0.25_0.29_8442 == "Doublet")]
        ##### function to go through the DF outputs based on pANN, and get the doublets based on distinct expectedDoublet values
        
        df_class_list <- colnames(seu_obj@meta.data)[grepl("DF.classification", colnames(seu_obj@meta.data))];
        
        results_per_df <- list()
        for (cl in df_class_list) {
          temp_result <- rownames(seu_obj@meta.data)[which(seu_obj@meta.data[cl] == "Doublet")]
          results_per_df[[length(results_per_df)+1]] <- list(temp_result)
        }
        
        print(head(results_per_df[2], 1))

        #these should now be saved in results_per_df #one is based of % expected doublet and the other on homotypic type here seurat_clusters
       
        ####Cells To Keep
        #keep <- unique(sort(intersect(cells, cellss)))
        #need to unlist the list you call from results_per_df otherwise error
        keep1 <- unique(sort(intersect(cells, unlist(results_per_df[1]))))
        print(paste0("[X] Found doublets in both packages from keep1 ", print(deparse(substitute(seu_obj))))); 

        keep2 <- unique(sort(intersect(cells, unlist(results_per_df[2]))))       
        print("[X] Keep2 done")
        print(dim(keep2))

        # qcPlots after subsetting
        print("[X] Generating QC Plots before subsetting")
        gen_qc_plot_clust(seu_obj, curr_mouse_id, results_folder)
        #
        gen_rqc_plot_clust(seu_obj, curr_mouse_id, results_folder)
        #
        gen_dub_plot(seu_obj, curr_mouse_id, results_folder)

        if (length(keep1) == 0) {
            print(paste("[X] No overlapping doublets found!",print(curr_mouse_id)))
	}

	if (length(keep1) > 0) {
            print("[X] Subsetting and saving with keep1")
            seu_obj_sub1=subset(seu_obj, cells=keep1, invert=TRUE)
            save_name <- paste(curr_mouse_id, "no_db", sep="_")
            saveit(seu_obj_sub=seu_obj_sub1, string=save_name, file=file.path(filtered_results_folder, paste0(curr_mouse_id,"_seu_obj_doublet_removed.rdata")))
            #save(seu_obj_sub1, file=file.path(results_folder, "seu_obj_sub1.rds"))
            gen_rqc_plot_sub(seu_obj_sub1, curr_mouse_id, results_folder)
            gen_qc_plot_sub(seu_obj_sub1, curr_mouse_id, results_folder)
            gen_dub_plot1(seu_obj_sub1, curr_mouse_id, results_folder)

        }

        if (length(keep2) > 0) {
            print("[X] Subsetting and saving with keep2")
            seu_obj_sub2=subset(seu_obj, cells=keep2, invert=TRUE)
            save_name <- paste(curr_mouse_id, "no_db", sep="_")
            saveit(seu_obj_sub=seu_obj_sub2, string=save_name, file=file.path(filtered_results_folder, paste0(curr_mouse_id,"_seu_obj_doublet_removed.rdata")))
            ##save(seu_obj_sub2, file=file.path(results_folder, "seu_obj_sub2.rds"))
            #gen_rqc_plot_sub(seu_obj_sub2, curr_mouse_id, results_folder)
            #gen_qc_plot_sub(seu_obj_sub2, curr_mouse_id, results_folder)
            gen_dub_plot2(seu_obj_sub2, curr_mouse_id, results_folder)
        }
        

        #break
    }
}
