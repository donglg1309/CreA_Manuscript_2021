
source("/Users/dongliguo/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")
second_mapping_to_summits<-read.table("input_summit_binding.bed_annotation_others_200bp",sep="\t",header=TRUE)
first_mapping_to_summits<-read.table("input_summit_binding.bed_annotation_200bp",sep="\t",header=TRUE)
ATG_utr_regions<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/utr_from_CDS.bed")
UTR_length<-ATG_utr_regions[,3]-ATG_utr_regions[,2]
distance_region_values<-NULL
total_gene_names<-paste(comvert_name(ATG_utr_regions[,4]),"-T",sep="")
end_values_out<-NULL
for ( i in 1:length(second_mapping_to_summits[,11]))
{
        if (sum(total_gene_names==as.character(second_mapping_to_summits[i,11]))>=1)
        distance_region_values[i]<-UTR_length[total_gene_names==as.character(second_mapping_to_summits[i,11])]
        if (sum(total_gene_names==as.character(second_mapping_to_summits[i,11]))==0)
        distance_region_values[i]=0
        if (second_mapping_to_summits[i,10]>200)
        end_values_out[i]=second_mapping_to_summits[i,10]-distance_region_values[i]+1
        if (second_mapping_to_summits[i,10]<=200)
        end_values_out[i]=second_mapping_to_summits[i,10]
}
gene_list_within_regions<-second_mapping_to_summits[end_values_out<200&end_values_out>-1500,11]
gene_list_within_regions_bed<-second_mapping_to_summits[end_values_out<200&end_values_out>-1500,]
gene_list_without_regions<-second_mapping_to_summits[!(end_values_out<200&end_values_out>-1500),]
#gene_list_without_regions_inter<-gene_list_without_regions[as.character(gene_list_without_regions[,11])%in%unique(gene_list_within_regions),]
#gene_list_without_regions_mutiple_second<-gene_list_without_regions_inter[gene_list_without_regions_inter[,10]<c(-1500),]
#bound_genes_A<-rbind(gene_list_within_regions_bed,gene_list_without_regions_mutiple_second)
bound_genes_A<-gene_list_within_regions_bed
second_mapping_genes<-data.frame(bound_genes_A[,2],bound_genes_A[,3],bound_genes_A[,4],bound_genes_A[,1],rep("promoter",length(bound_genes_A[,2])),bound_genes_A[,10],bound_genes_A[,11])
write.table(second_mapping_genes,file="output_files_values_promoter_second.bed",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")
second_mapping_genes<-read.table("output_files_values_promoter_second.bed",sep="\t")
first_mapping_genes<-read.table("output_files_values_promoter_first.bed",sep="\t")
merged_summits_input<-merge(first_mapping_genes,second_mapping_genes,by=c("V1","V2","V3"))
merged_summits_input_combined<-abs(merged_summits_input[,10]+merged_summits_input[,6])/4
merge_part_A<-merged_summits_input[abs(merged_summits_input[,10])>merged_summits_input_combined&abs(merged_summits_input[,6])>merged_summits_input_combined,]
write.table(merge_part_A[merge_part_A[,6]<0&merge_part_A[,10]<0,c(1:3,8:11)],file="merge_part_one.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(merged_summits_input[merged_summits_input[,10]>0,c(1:3,8:11)],file="merge_part_two.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")



