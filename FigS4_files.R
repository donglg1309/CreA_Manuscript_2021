#Figure S4B######################

source("/Users/dongliguo/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")

tiles_2<-make_tiles_for_promoter_ATG("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
GLU_WT_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_HA_1_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")

tiles_2<-make_tiles_for_promoter_ATG("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
GLU_WT_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_HA_1_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")

dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)
colors = colorRampPalette(c('black','black','yellow','yellow','yellow3','yellow3','orange','orange'))(200)
plot_heatmap_B(GLU_HA_1_promoter_single,GLU_HA_1_promoter_single,1500,321)
plot_heatmap_B(GLU_WT_HA_promoter_single,GLU_HA_1_promoter_single,1500,321)
plot_heatmap_B(GLU_S4_1_promoter_single,GLU_S4_1_promoter_single,1500,321)
plot_heatmap_B(GLU_WT_promoter_single,GLU_S4_1_promoter_single,1500,321)

dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)
colors = colorRampPalette(c('black','black','yellow','yellow','yellow3','yellow3','orange','orange'))(200)
plot_heatmap_B(GLU_HA_1_promoter_multiple,GLU_HA_1_promoter_multiple,1500,321)
plot_heatmap_B(GLU_WT_HA_promoter_multiple,GLU_HA_1_promoter_multiple,1500,321)
plot_heatmap_B(GLU_S4_1_promoter_multiple,GLU_S4_1_promoter_multiple,1500,321)
plot_heatmap_B(GLU_WT_promoter_multiple,GLU_S4_1_promoter_multiple,1500,321)

################################

