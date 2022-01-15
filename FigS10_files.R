#Figure S10A and Figure S10B######################

#make tiles for summit and get the signal
tiles_2<-make_tiles_for_summit("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/summit_list_files.bed")
GLU_WT<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_S4_2<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_2.bam.bed")
GLU_HA_1<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")
GLU_HA_2<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_2.bam.bed")

#plot correlation for the figure
plot(rowMeans(GLU_S4_1[,130:170]),rowMeans(GLU_S4_2[,130:170]),type="p",pch=16,cex=0.4,log="xy",xlim=c(100,5000),ylim=c(100,5000))
plot(rowMeans(GLU_HA_1[,130:170]),rowMeans(GLU_HA_2[,130:170]),type="p",pch=16,cex=0.4,log="xy",xlim=c(100,5000),ylim=c(100,5000))
plot((rowMeans(GLU_HA_1[,130:170])+rowMeans(GLU_HA_2[,130:170]))/2,(rowMeans(GLU_S4_1[,130:170])+rowMeans(GLU_S4_2[,130:170]))/2,type="p",pch=16,cex=0.4,log="xy",xlim=c(100,5000),ylim=c(100,5000))

#calculate correlation values
cor((rowMeans(GLU_S4_1[,130:170])),(rowMeans(GLU_S4_2[,130:170])))
cor((rowMeans(GLU_HA_1[,130:170])),(rowMeans(GLU_HA_2[,130:170])))
cor((rowMeans(GLU_HA_1[,130:170])+rowMeans(GLU_HA_2[,130:170]))/2,(rowMeans(GLU_S4_1[,130:170])+rowMeans(GLU_S4_2[,130:170]))/2)

################################

#Figure S10C#####################

tiles_2<-make_tiles_for_summit("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/combine_from_HA_and_S4_summit_list_files.bed")
GLU_WT_combine<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_combine<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_combine<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_HA_1_combine<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")

tiles_2<-make_tiles_for_summit("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/specific_summits_HA_summit_list_files.bed")
GLU_WT_HA<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_HA<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_HA<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_HA_1_HA<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")

tiles_2<-make_tiles_for_summit("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/specific_summits_S4_summit_list_files.bed")
GLU_WT_S4<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_S4<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_S4<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_HA_1_S4<-normalization_tiles_for_summits("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")

colors = colorRampPalette(c('black','black','yellow','yellow','yellow','yellow','orange','orange'))(200)
plot_heatmap_B(GLU_HA_1_combine,GLU_HA_1_combine[,130:171],1500,301)
plot_heatmap_B(GLU_WT_HA_combine,GLU_HA_1_combine[,130:171],1500,301)
plot_heatmap_B(GLU_S4_1_combine,GLU_HA_1_combine[,130:171],1500,301)
plot_heatmap_B(GLU_WT_combine,GLU_HA_1_combine[,130:171],1700,301)

dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)
colors = colorRampPalette(c('black','black','yellow','yellow','yellow','yellow','orange','orange'))(200)
plot_heatmap_B(GLU_HA_1_HA,GLU_HA_1_HA[,130:171],1200,301)
plot_heatmap_B(GLU_WT_HA_HA,GLU_HA_1_HA[,130:171],1200,301)
plot_heatmap_B(GLU_S4_1_HA,GLU_HA_1_HA[,130:171],1200,301)
plot_heatmap_B(GLU_WT_HA,GLU_HA_1_HA[,130:171],1200,301)

dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)
colors = colorRampPalette(c('black','black','yellow','yellow','yellow3','yellow3','orange','orange'))(200)
plot_heatmap_B(GLU_HA_1_S4,GLU_HA_1_S4[,130:171],1200,301)
plot_heatmap_B(GLU_WT_HA_S4,GLU_HA_1_S4[,130:171],1200,301)
plot_heatmap_B(GLU_S4_1_S4,GLU_HA_1_S4[,130:171],1200,301)
plot_heatmap_B(GLU_WT_S4,GLU_HA_1_S4[,130:171],1200,301)

################################




#Figure S3G#####################
  

dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)

colors = colorRampPalette(c('black','black','yellow','yellow','yellow3','yellow3','orange','orange'))(200)
plot_heatmap_B((GLU_HA_1+GLU_HA_2)/2,((GLU_HA_1+GLU_HA_2)/2)[,130:171],1500,301)
plot_heatmap_B(GLU_WT_HA,GLU_WT_HA,1500,301)
plot_heatmap_B((GLU_S4_1+GLU_S4_2)/2,((GLU_HA_1+GLU_HA_2)/2)[,130:171],1500,301)
plot_heatmap_B(GLU_WT,GLU_WT,1500,301)

################################




