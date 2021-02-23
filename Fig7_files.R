
#Figure 7B######################

#initation function for count value
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

#get gene direction info
gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")

#get the relationship for control dataset "AcuK-AcuM"
tiles_2<-make_tiles_for_summit_percise("AcuMHA_Ace_ACTTGA_s_8_sequence.fastq.bam_macs2_summits.bed")
GLU_WT<-normalization_tiles_for_summits_counts("AcuKmyc_Ace_GGCTAC_s_8_sequence.fastq.gz.bam_macs2_summits.bed")
AlcR_distance_plot<-GLU_WT[rowSums(GLU_WT[,500:2500])>0,]
AlcR_distance_plot_names<-as.character(AlcR_summits[rowSums(GLU_WT[,500:2500])>0,11])
AlcR_distance_output<-NULL
for ( i in 1:length(AlcR_distance_plot[,1]))
{       
        AlcR_distance_output[i]<-(which(AlcR_distance_plot[i,500:2500]!=0)-1000)[which.min(abs(which(AlcR_distance_plot[i,500:2500]!=0)-1000))]
}
AlcR_distance_output_v2<-AlcR_distance_output
AlcR_distance_output_v3<-AlcR_distance_output_v2
for ( i in 1:length(AlcR_distance_plot_names_v2))
{
        if (gff0_genome_location[gff0_genome_location[,4]%in%AlcR_distance_plot_names_v2[1],5]==2)
        AlcR_distance_output_v3[i]<-(-AlcR_distance_output_v2[i])
        if (gff0_genome_location[gff0_genome_location[,4]%in%AlcR_distance_plot_names_v2[1],5]==1)
        AlcR_distance_output_v3[i]<-(AlcR_distance_output_v2[i])
}
AcuK_distance_output_v3<-AlcR_distance_output_v3


#get the relationship for "CreA-AlcR"
tiles_2<-make_tiles_for_summit_percise("AlcR_summits_v2.bed")
GLU_WT<-normalization_tiles_for_summits_counts("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/summit_list_files.bed")
AlcR_distance_plot<-GLU_WT[rowSums(GLU_WT[,500:2500])>0,]
AlcR_distance_plot_names<-as.character(AlcR_summits[rowSums(GLU_WT[,500:2500])>0,11])
AlcR_distance_output<-NULL
for ( i in 1:length(AlcR_distance_plot[,1]))
{
        AlcR_distance_output[i]<-(which(AlcR_distance_plot[i,500:2500]!=0)-1000)[which.min(abs(which(AlcR_distance_plot[i,500:2500]!=0)-1000))]
}
AlcR_distance_output_v2<-AlcR_distance_output
AlcR_distance_output_v3<-AlcR_distance_output_v2
for ( i in 1:length(AlcR_distance_plot_names_v2))
{
        if (gff0_genome_location[gff0_genome_location[,4]%in%AlcR_distance_plot_names_v2[1],5]==2)
        AlcR_distance_output_v3[i]<-(-AlcR_distance_output_v2[i])
        if (gff0_genome_location[gff0_genome_location[,4]%in%AlcR_distance_plot_names_v2[1],5]==1)
        AlcR_distance_output_v3[i]<-(AlcR_distance_output_v2[i])
}

install.packages("rgr")
library("rgr")
plot(ecdf(abs(AcuK_distance_output_v3)),cex=0.5,col="black")
lines(ecdf(abs(AlcR_distance_output_v3)),cex=0.5,col="red")

gx.ks.test(abs(AcuK_distance_output_v3),abs(AlcR_distance_output_v3))



################################


#Figure 7C######################

#initation function for count value
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

#get ChIP-seq raw values
tiles_2<-make_tiles_for_summit_percise("AlcR_summits_v2.bed")
AlcR_summits_v2_all<-read.table("AlcR_summits_v2.bed")
GLU_WT<-normalization_tiles_for_summits_counts("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/summit_list_files.bed")
CreA_glucose<-normalization_tiles_for_summits_counts("GLU_CreA_1_CL2016_ChIPmix03_L001_R1.fastq.bam.bed")
CreA_ethanol<-normalization_tiles_for_summits_counts("CreA_ETOH_CCATACACT_ChIPMix51_bt2.bam.bed")
AlcR_glucose<-normalization_tiles_for_summits_counts("AlcR_flag_glu16h_glu_6h_CAGATC_CL_ChIPmix24_R1.fastq.gz.bam.bed")
AlcR_ethanol<-normalization_tiles_for_summits_counts("AlcR_flag_glu16h_ETOH6H_ACTTGA_CL_ChIPmix24_R1.fastq.bam.bed")

#filter the dataset and plot figure
AlcR_CreA_glucose<-CreA_glucose[rowSums(GLU_WT[,1450:1550])>0,]
AlcR_CreA_ethanol<-CreA_ethanol[rowSums(GLU_WT[,1450:1550])>0,]
AlcR_glucose_s<-AlcR_glucose[rowSums(GLU_WT[,1450:1550])>0,]
AlcR_ethanol_s<-AlcR_ethanol[rowSums(GLU_WT[,1450:1550])>0,]
log_2_A<-log2(rowMeans(AlcR_ethanol_s[,1475:1525])/rowMeans(AlcR_glucose_s[,1475:1525]))
log_2_B<-log2(rowMeans(AlcR_CreA_ethanol[,1475:1525])/rowMeans(AlcR_CreA_glucose[,1475:1525]))
plot(log_2_A[log_2_A>0],log_2_B[log_2_A>0],xlim=c(-6,6),ylim=c(-6,6),pch=16)
abline(h=0)
abline(v=0)

################################

#Figure 7H######################

#programming is same with Figure 6A

################################

#Figure 7I######################

#programming is same with Figure 6A

################################



