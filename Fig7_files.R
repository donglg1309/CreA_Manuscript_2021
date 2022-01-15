#Figure 7A######################

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

#Figure 7B######################

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



#Figure 7F######################

#get single or multiple CreA binding TF genes
full_list_transcription_factors_v2<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig4/full_list_transcription_factors.xls")
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
signle_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
write.table(mutiple_genes[mutiple_genes%in%full_list_transcription_factors_v2],file="mutiple_binding_TFs.xls",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(signle_genes[signle_genes%in%full_list_transcription_factors_v2],file="signle_binding_TFs.xls",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


#make tiles for the binding and calculate binding density
tiles_2<-make_tiles_for_promoter_ATG("signle_binding_TFs.xls")
GLU_WT_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_S4_2_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_2.bam.bed")
GLU_HA_1_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")
GLU_HA_2_promoter_single<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_2.bam.bed")

#make tiles for the binding and calculate binding density
tiles_2<-make_tiles_for_promoter_ATG("mutiple_binding_TFs.xls")
GLU_WT_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_BIA1_CreA_1.bam.bed")
GLU_WT_HA_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")
GLU_S4_1_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
GLU_S4_2_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_2.bam.bed")
GLU_HA_1_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")
GLU_HA_2_promoter_multiple<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_2.bam.bed")

#plot the figure
plot(colMeans(GLU_S4_1_promoter_single),type="l",ylim=c(150,1000),col="green",lwd=2)
lines(colMeans(GLU_S4_1_promoter_multiple),type="l",ylim=c(150,1000),col="orange",lwd=2)

################################

#Figure 7G######################

#read gene expression and TF gene list
full_list_transcription_factors_v2<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig4/full_list_transcription_factors.xls")
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
signle_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")

#find summits for each binding around TF gene promoters
summit_bed_files<-read.table("CreA_summits_and_annotation_per_annotation_v12.bed",sep="\t")
summit_bed_files_v2<-summit_bed_files[summit_bed_files[,5]%in%"promoter",]
write.table(summit_bed_files_v2[summit_bed_files_v2[,6]%in%full_list_transcription_factors_v2,c(1,2,3,6,2)],file="TF_binding_density_summits.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(summit_bed_files_v2[!summit_bed_files_v2[,6]%in%full_list_transcription_factors_v2,c(1,2,3,6,2)],file="others_binding_density_summits.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

#make tiles for the summits and get the corresponding signals
source("/Users/dongliguo/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")
normalization_tiles_for_summits_counts<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=1
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/3001,ncol=3001,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                       M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]*1000000/(0.05*l)
                                }
                        }
                return(M_0h.p.matrix2)

        }
TF_binding_density_out<-read.table("TF_binding_density_summits.bed")
others_binding_density_out<-read.table("others_binding_density_summits.bed")
tiles_2<-make_tiles_for_summit_percise("TF_binding_density_summits.bed")
CreA_glucose_TFs<-normalization_tiles_for_summits_counts("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
tiles_2<-make_tiles_for_summit_percise("others_binding_density_summits.bed")
CreA_glucose_others<-normalization_tiles_for_summits_counts("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")

#plot the figures
boxplot(rowMaxs(CreA_glucose_TFs[,1450:1550]),rowMaxs(CreA_glucose_others[,1450:1550]),outline=FALSE)
boxplot(rowMeans(CreA_glucose_TFs[,1475:1525]),sample(rowMeans(CreA_glucose_others[,1475:1525]),length(rowMeans(CreA_glucose_TFs[,1475:1525]))),col=c("red","gray"),outline=TRUE,ylim=c(0,800),pch=16)

################################



