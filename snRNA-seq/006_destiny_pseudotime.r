rm(list = ls())
library(cellAlign)
library(ggpubr)
library(destiny)
library(ggthemes)
library(Seurat)
library(gridExtra)
library(dplyr)
library(magrittr)
library(SingleCellExperiment)

savedir = "./ao/data_pipeline1/mergee13/exci_lin/pseudotime/destiny_2023/ctrl_collapsed"
e13_exc <- readRDS("./ao/data_pipeline1/mergee13/exci_lin/e13_exci_lin_celltypes_obj.rds")

e13_exc <- SetIdent(e13_exc, value=obj$CellType_coll);
table((e13_exc$CellType_coll))

ctrl_obj <- subset(e13_exc, subset = genotype =="CTRL");
cko_obj <- subset(e13_exc, subset = genotype =="CKO");
#saveRDS(ctrl_obj, file = paste0(savedir,'/e13_control_obj.rds'));
#saveRDS(cko_obj, file = paste0(savedir,'/e13_cko_obj.rds'));
##########################################################
## Feature Selection
hkeepfeat = which((rowSums(as.matrix(ctrl_obj@assays$RNA@counts) > 0) / ncol(ctrl_obj)) > 0.1) %>% names;

## Run Destiny
#hdeng <- ctrl_obj@assays$SCT@data
hdeng <- ctrl_obj@assays$RNA@data;
hdeng2 = hdeng[rownames(hdeng) %in% hkeepfeat,];
sigmas <- find_sigmas(t(as.matrix(hdeng2)));
Sys.time()
hdm <- DiffusionMap(t(as.matrix(hdeng2)), sigma = optimal_sigma(sigmas), k = 1000);
Sys.time()
save(hdm, file = paste0(savedir,"/e13_flox_hdm.rdata"));
##re do for viz
#hdm <- DiffusionMap(t(as.matrix(hdeng2)), n_pcs = 50, k = 1000);

#diffusion map object
hdf <- data.frame(DC1 = eigenvectors(hdm)[,1],
                  DC2 = eigenvectors(hdm)[,2],
                  color = ctrl_obj$CellType_coll[colnames(hdeng2)]);

#pdf('ao_e13_flox_diffusion_map.pdf')
##
ggplot(hdf, aes(x = DC1, y =  DC2, colour = ctrl_obj$CellType_coll[colnames(hdeng2)])) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() + ggtitle("Diffusion Plot E13 flox")+ theme(plot.title = element_text(hjust = 0.5));
ggsave(paste0(savedir,"/ao_e13_flox_destiny_plot_dc13_CellTypeColl.pdf"))

##
ggplot(hdf, aes(x = DC1, y =  eigenvectors(hdm)[,3], colour = ctrl_obj$CellType[colnames(hdeng2)])) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() + ggtitle("Diffusion Plot E13 flox")+ theme(plot.title = element_text(hjust = 0.5));
ggsave(paste0(savedir,"/ao_e13_flox_destiny_plot_dc13_CellType.pdf"))
###
save(hdf, file = paste0(savedir,'/e13_flox_exci_diffusion_hdf_tmp.rdata'));

##
#colors_col <-c("#993399", "#FFCCFF", '#CC0099', "#333333"); colors_col <- colors[as.numeric(hdf$type)]
#scatterplot3d::scatterplot3d(x=eigenvectors(hdm)[,1],y= eigenvectors(hdm)[,2], z= eigenvectors(hdm)[,2], color = colors_col)
# Run pseudotime. #from emre: DPT gives error: stats$g >= gmin is not TRUE#
#used DPT, Diffusion Pseudo Time
hdpt = DPT(hdm);

#pdf('DPT_ao_e13_flox_noInt_20201207.pdf')
plot.DPT(hdpt, root = 1, dcs=c(1,2));
#plot(hdpt, dcs=c(1,2), root=2)
ggsave(paste0(savedir,"/ao_e13_flox_destiny_plot_dpt.pdf"))

#dev.off()
plot(hdpt, root = 1, paths_to = c(2,3), dcs=c(1,2))
ggsave(paste0(savedir,"/ao_e13_flox_destiny_traj_pseudo_plot.pdf"))
#
plot(hdpt, root = 1, paths_to = c(2,3), dcs=c(1,2), col_by = 'branch', pal = viridis::inferno)
ggsave(paste0(savedir,"/ao_e13_flox_destiny_traj_plot.pdf"))
# Find tip cell
tips = find_tips(hdpt);
#tip cell aka start
thetip = which( ctrl_obj[[]][colnames(hdeng2)[tips], 'CellType_coll'] == 'RG.E13') # RG.2#'SVZ2 (migrating) [8-P]' #bc this is the 1st point I see
htc = tips[thetip[1]]
#htc = which.max(hdm@eigenvectors[,3])
##in case you dont know th efull name of your CellTypes, look then go back to the tip
#table(ctrl_obj$CellType)
# if needed, remove outliers #not needed ##remove none
hkeepind = which(hdpt[,htc] >= quantile(hdpt[,htc], 0.00, na.rm = TRUE) & hdpt[,htc] <= quantile(hdpt[,htc], 1.00, na.rm = TRUE))

# Show one dimensional trajectory
hsub = subset(ctrl_obj, cells = colnames(hdeng2)[hkeepind])
pdtime_unnorm = hdpt[hkeepind, htc]
hsub$pdtime = (pdtime_unnorm - min(pdtime_unnorm))/(max(pdtime_unnorm) - min(pdtime_unnorm))

hmeta = hsub[[]]

############
pdf(paste0(savedir,'/DENSITY_ao_e13_flox_FILT.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  theme(legend.position="top")+ theme(legend.position = "right")+ ggtitle("E13 flox pseudotime")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_flox_CellType_coll.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  theme(legend.position="top")+ theme(legend.position = "right")+ ggtitle("E13 flox pseudotime")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_flox_indiv_mice_FILT.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  facet_wrap(~orig.ident)+ theme(legend.position = "right")+ggtitle("E13 flox pseudotime original animal")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_flox_CellType_FILT.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType)) +
  geom_density(alpha = 0.4)+
  theme_classic()  + facet_wrap(~CellType) +
  theme(legend.position = "right")+ggtitle("E13 flox pseudotime by cell-type")+ theme(plot.title = element_text(hjust = 0.5));
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_flox_coll.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  theme(legend.position="top")+ theme(legend.position = "right")+ ggtitle("E13 flox pseudotime")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_flox_indiv_mice_coll.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  facet_wrap(~orig.ident)+ theme(legend.position = "right")+ggtitle("E13 flox pseudotime original animal")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_flox_CellType_coll.pdf'))
ggplot(hmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic()  + facet_wrap(~CellType_coll) +
  theme(legend.position = "right")+ggtitle("E13 flox pseudotime by cell-type")+ theme(plot.title = element_text(hjust = 0.5));
dev.off()

#save variables
#get object you want to save
e13_flox_genes = hkeepfeat
e13_flox_pdtime = hmeta$pdtime
e13_flox_cells = rownames(hmeta)
# flox_list = list(flox_pdtime, flox_genes, flox_cells)
# names(flox_list) = c('pdtime', 'genes', 'cells')
# saveRDS(flox_list, '~/flox_pseudotime_list.rds')

savedir = "./ao/data_pipeline1/mergee13/exci_lin/pseudotime/destiny_2023/ctrl_collapsed/"

e13hmeta = hmeta
save(e13_flox_genes,e13_flox_pdtime, e13_flox_cells, e13hmeta,ctrl_obj, file = paste0(savedir,"/ao_e13_flox_pseudotime_destiny_objects_all_cells.rdata" ))
###################################################################################
###################################################################################
#############  NEXT

## Feature Selection
ckeepfeat = which((rowSums(as.matrix(cko_obj@assays$RNA@counts) > 0) / ncol(cko_obj)) > 0.1) %>% names;


## Run Destiny
#cdeng <- cko_obj@assays$SCT@data
cdeng <- cko_obj@assays$RNA@data;
cdeng2 = cdeng[rownames(cdeng) %in% ckeepfeat,];
sigmas <- find_sigmas(t(as.matrix(cdeng2)));
Sys.time()
cdm <- DiffusionMap(t(as.matrix(cdeng2)), sigma = optimal_sigma(sigmas), k = 1000);
Sys.time()
save(cdm, file = paste0(savedir,"/e13_cko_cdm.rdata"));
##re do for viz
#cdm <- DiffusionMap(t(as.matrix(cdeng2)), n_pcs = 50, k = 1000);

#diffusion map object
cdf <- data.frame(DC1 = eigenvectors(cdm)[,1],
                  DC2 = eigenvectors(cdm)[,3],
                  color = cko_obj$CellType_coll[colnames(cdeng2)]);

#pdf('ao_e13_cko_diffusion_map.pdf')
##
ggplot(cdf, aes(x = DC1, y =  eigenvectors(cdm)[,3], colour = cko_obj$CellType_coll[colnames(cdeng2)])) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() + ggtitle("Diffusion Plot E13 cKO")+ theme(plot.title = element_text(hjust = 0.5));
ggsave(paste0(savedir,"/ao_e13_cko_destiny_plot_dc13_CellTypeColl.pdf"))

##
ggplot(cdf, aes(x = DC1, y =  eigenvectors(cdm)[,3], colour = cko_obj$CellType[colnames(cdeng2)])) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() + ggtitle("Diffusion Plot E13 cKO")+ theme(plot.title = element_text(hjust = 0.5));
ggsave(paste0(savedir,"/ao_e13_cko_destiny_plot_dc13_CellType.pdf"))
###
save(cdf, file = paste0(savedir,'/e13_cko_exci_diffusion_cdf_tmp.rdata'));

##
#colors_col <-c("#993399", "#FFCCFF", '#CC0099', "#333333"); colors_col <- colors[as.numeric(cdf$type)]
#scatterplot3d::scatterplot3d(x=eigenvectors(cdm)[,1],y= eigenvectors(cdm)[,3], z= eigenvectors(cdm)[,2], color = colors_col)
# Run pseudotime. #from emre: DPT gives error: stats$g >= gmin is not TRUE#
#used DPT, Diffusion Pseudo Time
cdpt = DPT(cdm);

#pdf('DPT_ao_e13_cko_noInt_20201207.pdf')
plot.DPT(cdpt, root = 1, dcs=c(1,2));
#plot(cdpt, dcs=c(1,2), root=2)
ggsave(paste0(savedir,"/ao_e13_cko_destiny_plot_dpt.pdf"))

#dev.off()
plot(cdpt, root = 1, paths_to = c(2,3), dcs=c(1,2))
ggsave(paste0(savedir,"/ao_e13_cko_destiny_traj_pseudo_plot.pdf"))
#
plot(cdpt, root = 1, paths_to = c(2,3), dcs=c(1,2), col_by = 'branch', pal = viridis::inferno)
ggsave(paste0(savedir,"/ao_e13_cko_destiny_traj_plot.pdf"))
# Find tip cell
tips = find_tips(cdpt);
#tip cell aka start
thetip = which( cko_obj[[]][colnames(cdeng2)[tips], 'CellType_coll'] == 'RG.E13') 
htc = tips[thetip[1]]
#htc = which.max(cdm@eigenvectors[,3])
##in case you dont know the full name of your CellTypes, look then go back to the tip #table(cko_obj$CellType)
# if needed, remove outliers #not needed
hkeepind = which(cdpt[,htc] >= quantile(cdpt[,htc], 0.00, na.rm = TRUE) & cdpt[,htc] <= quantile(cdpt[,htc], 1.00, na.rm = TRUE))

# Show one dimensional trajectory
csub = subset(cko_obj, cells = colnames(cdeng2)[hkeepind])
pdtime_unnorm = cdpt[hkeepind, htc]
csub$pdtime = (pdtime_unnorm - min(pdtime_unnorm))/(max(pdtime_unnorm) - min(pdtime_unnorm))
cmeta = csub[[]]

############
pdf(paste0(savedir,'/DENSITY_ao_e13_cko_FILT.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  theme(legend.position="top")+ theme(legend.position = "right")+ ggtitle("E13 cko pseudotime")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_cko_CellType_coll.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  theme(legend.position="top")+ theme(legend.position = "right")+ ggtitle("E13 cko pseudotime")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_cko_indiv_mice_FILT.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  facet_wrap(~orig.ident)+ theme(legend.position = "right")+ggtitle("E13 cko pseudotime original animal")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_cko_CellType_FILT.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType)) +
  geom_density(alpha = 0.4)+
  theme_classic()  + facet_wrap(~CellType) +
  theme(legend.position = "right")+ggtitle("E13 cko pseudotime by cell-type")+ theme(plot.title = element_text(hjust = 0.5));
dev.off()

pdf(paste0(savedir,'/DENSITY_ao_e13_cko_coll.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  theme(legend.position="top")+ theme(legend.position = "right")+ ggtitle("E13 cko pseudotime")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf(paste0(savedir,'/DENSITY_ao_e13_cko_indiv_mice_coll.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic() +
  facet_wrap(~orig.ident)+ theme(legend.position = "right")+ggtitle("E13 cko pseudotime original animal")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf(paste0(savedir,'/DENSITY_ao_e13_cko_CellType_coll.pdf'))
ggplot(cmeta, aes(x=pdtime, fill=CellType_coll)) +
  geom_density(alpha = 0.4)+
  theme_classic()  + facet_wrap(~CellType_coll) +
  theme(legend.position = "right")+ggtitle("E13 cko pseudotime by cell-type")+ theme(plot.title = element_text(hjust = 0.5));
dev.off()

#save variables
e13_cko_pdtime = cmeta$pdtime
e13_cko_genes = ckeepfeat
e13_cko_cells = rownames(cmeta)
e13_cko_list = list(cko_pdtime, cko_genes, cko_cells)
names(cko_list) = c('pdtime', 'genes', 'cells')
# 
# saveRDS(cko_list, '~/cko_pseudotime.RDS')

e13cmeta = cmeta

savedir = "./ao/data_pipeline1/mergee13/exci_lin/pseudotime/destiny_2023/cko_redo/plots/"
save(e13_cko_genes,e13_cko_pdtime, e13_cko_cells, e13cmeta,cko_obj, cdm, file = paste0(savedir,"ao_e13_cko_pseudotime_destiny_objects_all_cells.rdata" ))
# ########################
## GET PSEUDOTIME VALUES
e13_ctrl_pt <- e13hmeta$pdtime
e13_cko_pt <- e13cmeta$pdtime
#################### to add to whole object meta data   ############################
##ctrl
# Extract corresponding row names
ctrl_row_names <- rownames(e13hmeta)
# Create a dataframe or list to store both values and corresponding row names
e13_control_pseudo <- data.frame(row_names = ctrl_row_names, column_values = e13_ctrl_pt)
## cko
# Extract corresponding row names
cko_row_names <- rownames(e13cmeta)
# Create a dataframe or list to store both values and corresponding row names
e13_cko_pseudo <- data.frame(row_names = cko_row_names, column_values = e13_cko_pt)
## merge to add to whole exci obj
joint <- rbind(e13_control_pseudo, e13_cko_pseudo)
###
## Set the row names of the new metadata dataframe to match the cell names
cell_barcodes_e13 <- joint$row_names  # row_names has column contains cell names #test match to object
## Remove the cell names column from the metadata dataframe
joint <- as.data.frame(joint[, -1])  
names(joint)[1] <- "pseudotime"  ## name
rownames(joint) <- cell_barcodes_e13 ##make barcodes rownames to facilitate joining
## ##$$$$
e13_exc <- AddMetaData(object = e13_exc, metadata = joint)
#saveRDS(test_e13, "./ao/data_pipeline1/mergee13/exci_lin/e13_exci_lin_celltypes_obj1.rds")
########################
###############  STATISTICAL TEST
# Perform two-sided Mann-Whitney U test
result_e13_pd <- wilcox.test(e13_ctrl_pt, e13_cko_pt, alternative = "two.sided")
#print(result_e13_pd)
#Wilcoxon rank sum test with continuity correction
#data:  e13_ctrl_pt and e13_cko_pt
#W = 55749921, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0
# > summary(e13_cko_pt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.6078  0.6600  0.6872  0.8386  1.0000 
# > summary(e13_ctrl_pt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.2488  0.3251  0.4409  0.6911  1.0000 



#############################################
# plot(hdpt, col_by = 'branch', divide = 3, dcs = c(-1,-3,3), pch = 20)
# plot(cdpt, col_by = 'branch', divide = 3, dcs = c(-1,-3,3), pch = 20)

# plot(hdpt, divide = 3, dcs = c(-1,-3,3), pch = 20)
# ggsave(paste0(savedir,"/ao_e13_ctrl_3d_pseudotime1.pdf"))

# plot(cdpt, divide = 3, dcs = c(-1,-3,3), pch = 20)
# ggsave(paste0(savedir,"/ao_e13_cko_3d_pseudotime1.pdf"))


# plot(hdpt, divide = 3, dcs = c(-1,-4,3), pch = 20)
# ggsave(paste0(savedir,"/ao_e13_ctrl_3d_pseudotime.pdf"))

# plot(cdpt, divide = 3, dcs = c(-1,-4,3), pch = 20)
# ggsave(paste0(savedir,"/ao_e13_cko_3d_pseudotime.pdf"))

