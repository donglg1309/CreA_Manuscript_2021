

#Figure 6A######################

library(RColorBrewer)
library("gplots")
my_palette <- colorRampPalette(c("black","black","yellow","yellow"))(n = 90)


#read the conservation info matrix
gene_expression_values_values_v2<-read.table("gene_conservation_values_name.xls",row.names=1)
gene_expression_values_values<-as.matrix(gene_expression_values_values_v2)[,-10]

#read CreA target info
CreA_target<-read.table("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/CreA_summits_and_annotation_v12.xls")
CreA_target2<-as.character(CreA_target$V1)
matrix_one<-matrix(0,nrow=length(gene_expression_values_values[,1]),ncol=length(gene_expression_values_values[1,])+1)
matrix_one[,1:length(gene_expression_values_values[1,])]<-gene_expression_values_values
matrix_one[rownames(gene_expression_values_values)%in%CreA_target2,length(gene_expression_values_values[1,])+1]<-10
colnames(matrix_one)<-c(as.character(files_names)," CreA_targets")

#get S. pombe gene names and the corresponding A. nidulans genes, the similiar programming was run for other species
Spom_Anid_orthologs2<-read.table("/Volumes/LiguoDisk/Methods/Fungi_Motif/Gene_names/FungiDB-3.0_Anidulans_FGSCA4Gene.txt_output_connection/SPBC_Gene_name_output.txt")
Spom_targets<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/conserved_processes/Spom_targets.xls")
Spom_distribute<-Spom_Anid_orthologs2[Spom_Anid_orthologs2[,1]%in%Spom_targets,2]


#function for species relationship order according to CreA binding 
order_files_CreA<-function(matrix_names_out,gene_names_out)
{
        target_part<-matrix_names_out[gene_names_out%in%CreA_target2,]
        non_target_part<-matrix_names_out[!gene_names_out%in%CreA_target2,]
        row.names(target_part)<-gene_names_out[gene_names_out%in%CreA_target2]
        row.names(non_target_part)<-gene_names_out[!gene_names_out%in%CreA_target2]
        output_values<-rbind(
        target_part[order(rowSums(target_part),decreasing=TRUE),],
        non_target_part[order(rowSums(non_target_part),decreasing=TRUE),])
        return(output_values)
}
matrix_one_out<-rbind(
order_files_CreA(matrix_one[rowSums(matrix_one[,1:8])>0,],row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,1:8])>0]),
order_files_CreA(matrix_one[rowSums(matrix_one[,9:15])>=2&rowSums(matrix_one[,1:8])==0,],row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,9:15])>=2&rowSums(matrix_one[,1:8])==0]),
order_files_CreA(matrix_one[rowSums(matrix_one[,9:15])==1&rowSums(matrix_one[,1:8])==0,],row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,9:15])==1&rowSums(matrix_one[,1:8])==0])
)
matrix_one<-matrix_one_out


#plot CreA binding info
output_values_matrix_input_barp_CreA<-matrix(c(
length(row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,1:8])>0][row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,1:8])>0]%in%CreA_target2]),
length(row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,9:15])>=2&rowSums(matrix_one[,1:8])==0][row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,9:15])>=2&rowSums(matrix_one[,1:8])==0]%in%CreA_target2]),
length(row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,9:15])==1&rowSums(matrix_one[,1:8])==0][row.names(gene_expression_values_values_v2)[rowSums(matrix_one[,9:15])==1&rowSums(matrix_one[,1:8])==0]%in%CreA_target2])),nrow=1)
output_values_matrix_input_barp_plot<-rbind(rev(output_values_matrix_input_barp),rev(output_values_matrix_input_barp))
barplot(t(output_values_matrix_input_barp_plot),col=c("whitesmoke","red"))

#plot the species figure
plot_heatmap_E((matrix_one[rowSums(matrix_one[,1:8])>0,]),1,15)
plot_heatmap_E((matrix_one[rowSums(matrix_one[,9:15])>=2&rowSums(matrix_one[,1:8])==0,]),1,15)
plot_heatmap_E(matrix_one[rowSums(matrix_one[,9:15])==1&rowSums(matrix_one[,1:8])==0,],1,15)

#write out the gene name for each cluster
write.table(rownames(gene_expression_values_values)[rowSums(matrix_one[,1:8])>0],file="Aspergillus_nidulans_cluster_one.xls",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(rownames(gene_expression_values_values)[rowSums(matrix_one[,9:15])>=2&rowSums(matrix_one[,1:8])==0],file="Aspergillus_nidulans_cluster_two.xls",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(rownames(gene_expression_values_values)[rowSums(matrix_one[,9:15])==1&rowSums(matrix_one[,1:8])==0],file="Aspergillus_nidulans_cluster_three.xls",row.names=FALSE,col.names=FALSE,quote=FALSE)

#plot binding info for pombe, other species will run similiar script
matrix_pombe_genes<-matrix(0,ncol=2,nrow=length(rownames(matrix_one_out)))
matrix_pombe_genes[rownames(matrix_one_out)%in%Spom_distribute,]<-1
sum(matrix_pombe_genes[1:5333,1])
sum(matrix_pombe_genes[5334:(5333+3978),1])
colors <- colorRampPalette(c("whitesmoke","whitesmoke","red","red"))(n = 90)
plot_heatmap_E(as.matrix(matrix_pombe_genes)+as.matrix(rbind(matrix(0,nrow=1,ncol=2),matrix_pombe_genes[1:10719,]))+
as.matrix(rbind(matrix(0,nrow=2,ncol=2),matrix_pombe_genes[1:10718,]))+
as.matrix(rbind(matrix(0,nrow=3,ncol=2),matrix_pombe_genes[1:10717,]))+
as.matrix(rbind(matrix(0,nrow=4,ncol=2),matrix_pombe_genes[1:10716,]))+
as.matrix(rbind(matrix(0,nrow=5,ncol=2),matrix_pombe_genes[1:10715,])),1,2)


#plot the heatmap for niger data, other species will run similiar script
niger_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig2/conserved_genes/niger_DEGs_nidulans.xls")
matrix_niger_DEGs_genes<-matrix(0,ncol=2,nrow=length(rownames(matrix_one_out)))
matrix_niger_DEGs_genes[rownames(matrix_one_out)%in%niger_DEGs,]<-1
sum(matrix_niger_DEGs_genes[1:5333,1])
sum(matrix_niger_DEGs_genes[5334:(5333+3978),1])
colors <- colorRampPalette(c("whitesmoke","whitesmoke","orange","orange"))(n = 90)
plot_heatmap_E(matrix_niger_DEGs_genes+rbind(matrix(0,nrow=1,ncol=2),matrix_niger_DEGs_genes[1:10719,])+
rbind(matrix(0,nrow=2,ncol=2),matrix_niger_DEGs_genes[1:10718,])+
rbind(matrix(0,nrow=3,ncol=2),matrix_niger_DEGs_genes[1:10717,])+
rbind(matrix(0,nrow=4,ncol=2),matrix_niger_DEGs_genes[1:10716,])+
rbind(matrix(0,nrow=5,ncol=2),matrix_niger_DEGs_genes[1:10715,]),1,2)

################################


#Figure 6B#################################

#similiar plot with Figure 3B

##########################################




