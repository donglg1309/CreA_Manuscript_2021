#Figure 6A######################

#make tiles and calculate binding signals
tiles_2<-make_tiles_for_promoter_ATG("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/combined_HA_CreA_annotated_uniq_genes_v2.xls")
GLU_WT_HA_promoter<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_HA_1_promoter_glucose<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")
GLU_HA_1_promoter_acetate<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/ACE_HA_CreA_1.bam.bed")
GLU_HA_1_promoter_proline<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/PRO_HA_CreA_1.bam.bed")
GLU_HA_1_promoter_cf<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/CF_HA_CreA_1.bam.bed")
GLU_WT_HA_promoter<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_HA_2_promoter_glucose<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_2.bam.bed")
GLU_HA_2_promoter_acetate<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/ACE_HA_CreA_2.bam.bed")
GLU_HA_2_promoter_proline<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/PRO_HA_CreA_2.bam.bed")
GLU_HA_2_promoter_cf<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/CF_HA_CreA_2.bam.bed")

#plot figure

dev.new(width=12, height=6.6)
par(mar=c(1,1,1,1),mfrow=c(1,6),xpd=FALSE)
colors = colorRampPalette(c('black','black','yellow','yellow','orange','orange','red','red'))(200)
plot_heatmap_B(GLU_HA_1_promoter_glucose,GLU_HA_1_promoter_glucose,1500,321)
plot_heatmap_B(GLU_HA_1_promoter_acetate,GLU_HA_1_promoter_glucose,1500,321)
plot_heatmap_B(GLU_HA_1_promoter_proline,GLU_HA_1_promoter_glucose,1500,321)
plot_heatmap_B(GLU_HA_1_promoter_cf,GLU_HA_1_promoter_glucose,1500,321)
plot_heatmap_B(GLU_WT_HA_promoter,GLU_HA_1_promoter_glucose,1500,321)


################################

#Figure 6B######################

boxplot(
(rowMeans(GLU_HA_1_promoter_glucose[72:204,])+rowMeans(GLU_HA_2_promoter_glucose[72:204,]))/2,
(rowMeans(GLU_HA_1_promoter_acetate[72:204,])+rowMeans(GLU_HA_2_promoter_acetate[72:204,]))/2,
(rowMeans(GLU_HA_1_promoter_proline[72:204,])+rowMeans(GLU_HA_2_promoter_proline[72:204,]))/2,
(rowMeans(GLU_HA_1_promoter_cf[72:204,])+rowMeans(GLU_HA_2_promoter_cf[72:204,]))/2,
rowMeans(GLU_WT_HA_promoter[72:204,]),
col=c("tomato","lightblue","orchid","palegreen","saddlebrown"),outline=FALSE,ylim=c(0,800))

#t.test for the figures
t.test((rowMeans(GLU_HA_1_promoter_glucose[72:204,])+rowMeans(GLU_HA_2_promoter_glucose[72:204,]))/2,(rowMeans(GLU_HA_1_promoter_acetate[72:204,])+rowMeans(GLU_HA_2_promoter_acetate[72:204,]))/2)
t.test((rowMeans(GLU_HA_1_promoter_glucose[72:204,])+rowMeans(GLU_HA_2_promoter_glucose[72:204,]))/2,(rowMeans(GLU_HA_1_promoter_proline[72:204,])+rowMeans(GLU_HA_2_promoter_proline[72:204,]))/2)
t.test((rowMeans(GLU_HA_1_promoter_glucose[72:204,])+rowMeans(GLU_HA_2_promoter_glucose[72:204,]))/2,(rowMeans(GLU_HA_1_promoter_cf[72:204,])+rowMeans(GLU_HA_2_promoter_cf[72:204,]))/2)
t.test((rowMeans(GLU_HA_1_promoter_cf[72:204,])+rowMeans(GLU_HA_2_promoter_cf[72:204,]))/2,(rowMeans(GLU_WT_HA_promoter[72:204,])+rowMeans(GLU_WT_HA_promoter[72:204,]))/2)

################################



#Figure 6F######################

#programming is same with Figure 6A

################################

#Figure 6G######################

#programming is same with Figure 6B

################################



