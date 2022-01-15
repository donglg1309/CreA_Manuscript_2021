

#Figure S3G#####################

#make tiles for summit and get the signal
tiles_2<-make_tiles_for_summit("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/summit_list_files.bed")
GLU_WT<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_S4_2<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_2.bam.bed")
GLU_HA_1<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")
GLU_HA_2<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_2.bam.bed")



dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)

colors = colorRampPalette(c('black','black','yellow','yellow','yellow3','yellow3','orange','orange'))(200)
plot_heatmap_B((GLU_HA_1+GLU_HA_2)/2,((GLU_HA_1+GLU_HA_2)/2)[,130:171],1500,301)
plot_heatmap_B(GLU_WT_HA,GLU_WT_HA,1500,301)
plot_heatmap_B((GLU_S4_1+GLU_S4_2)/2,((GLU_HA_1+GLU_HA_2)/2)[,130:171],1500,301)
plot_heatmap_B(GLU_WT,GLU_WT,1500,301)

################################


#Figure S3H#####################

#The figure is plotted with similiar function of Figure 3B

################################

