
source("/Users/dongliguo/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")
first_mapping_to_summits<-read.table("input_summit_binding.bed_annotation_200bp",sep="\t",header=TRUE)
ATG_utr_regions<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/utr_from_CDS.bed")
UTR_length<-ATG_utr_regions[,3]-ATG_utr_regions[,2]
distance_region_values<-NULL
total_gene_names<-paste(comvert_name(ATG_utr_regions[,4]),"-T",sep="")
end_values_out<-NULL
for ( i in 1:length(first_mapping_to_summits[,11]))
{
	if (sum(total_gene_names==as.character(first_mapping_to_summits[i,11]))>=1)
	distance_region_values[i]<-UTR_length[total_gene_names==as.character(first_mapping_to_summits[i,11])]
	if (sum(total_gene_names==as.character(first_mapping_to_summits[i,11]))==0)
	distance_region_values[i]=0
	if (first_mapping_to_summits[i,10]>200)
	end_values_out[i]=first_mapping_to_summits[i,10]-distance_region_values[i]+1
	if (first_mapping_to_summits[i,10]<=200)
	end_values_out[i]=first_mapping_to_summits[i,10]
}
gene_list_within_regions<-first_mapping_to_summits[end_values_out<200&end_values_out>-1500,11]
gene_list_within_regions_bed<-first_mapping_to_summits[end_values_out<200&end_values_out>-1500,]
gene_list_without_regions<-first_mapping_to_summits[!(end_values_out<200&end_values_out>-1500),]
gene_list_without_regions_inter<-gene_list_without_regions[as.character(gene_list_without_regions[,11])%in%unique(gene_list_within_regions),]
gene_list_without_regions_mutiple_first<-gene_list_without_regions_inter[gene_list_without_regions_inter[,10]<c(-1500),]

GTF_files_input<-read.table("A_nidulans_FGSC_A4_version_s10-m04-r03_features_with_chromosome_sequences_v2.gtf")
write.table(GTF_files_input[!(as.character(GTF_files_input[,10])%in%unique(gene_list_within_regions)),],file="other_genes_files.gtf",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

bound_genes_A<-rbind(gene_list_within_regions_bed,gene_list_without_regions_mutiple_first)
other_locations<-first_mapping_to_summits[!as.character(first_mapping_to_summits[,1])%in%as.character(bound_genes_A[,1]),]

gene_names_A<-as.character(other_locations[,8])
output_gene_names<-NULL
for ( i in 1:length(gene_names_A))
{
	output_gene_names[i]<-strsplit(gene_names_A," ")[[i]][1]
}

output_gene_names[output_gene_names=="promoter-TSS"]="gene_body"
output_gene_names[output_gene_names=="exon"]="gene_body"
output_gene_names[output_gene_names=="intron"]="gene_body"



write.table(data.frame(bound_genes_A[,2],bound_genes_A[,3],bound_genes_A[,4],bound_genes_A[,1],rep("promoter",length(bound_genes_A[,2])),bound_genes_A[,10],bound_genes_A[,11]),file="output_files_values_promoter_first.bed",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

write.table(data.frame(other_locations[,2],other_locations[,3],other_locations[,4],other_locations[,1],output_gene_names,other_locations[,10],other_locations[,11]),file="output_files_values_others.bed",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")










