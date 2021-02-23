
#Figure 3A###############

#read up-regulated or downregulated genes, gene expression values, and targets
full_list_files<-read.table("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/Fig1_supplemental_files/glucose_gene_expression/ddsFull_fpkm_values_glu_with_p_values.xls",sep="\t",row.names=1,header=TRUE)
repressed_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")
activated_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
signle_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
all_the_targets<-c(mutiple_genes,signle_genes)

#get the gene expression changes for each class 
input_plot_class_one<-data.frame(full_list_files$log2FoldChange[full_list_files$Row.names%in%(repressed_DEGs[repressed_DEGs%in%all_the_targets])],
full_list_files$pvalue[full_list_files$Row.names%in%(repressed_DEGs[repressed_DEGs%in%all_the_targets])])
input_plot_class_two<-data.frame(full_list_files$log2FoldChange[full_list_files$Row.names%in%(all_the_targets[!all_the_targets%in%c(repressed_DEGs,activated_DEGs)])],
full_list_files$pvalue[full_list_files$Row.names%in%(all_the_targets[!all_the_targets%in%c(repressed_DEGs,activated_DEGs)])])
input_plot_class_three<-data.frame(full_list_files$log2FoldChange[full_list_files$Row.names%in%(c(repressed_DEGs,activated_DEGs)[!c(repressed_DEGs,activated_DEGs)%in%all_the_targets])],
full_list_files$pvalue[full_list_files$Row.names%in%(c(repressed_DEGs,activated_DEGs)[!c(repressed_DEGs,activated_DEGs)%in%all_the_targets])])
input_plot_class_four<-data.frame(full_list_files$log2FoldChange[full_list_files$Row.names%in%(activated_DEGs[activated_DEGs%in%all_the_targets])],
full_list_files$pvalue[full_list_files$Row.names%in%(activated_DEGs[activated_DEGs%in%all_the_targets])])

#plot the figure
plot(input_plot_class_three[,1],-log10(input_plot_class_three[,2]),pch=15,cex=0.5,ylim=c(0,100),xlim=c(-15,15),col="turquoise")
points(input_plot_class_two[,1],-log10(input_plot_class_two[,2]),col='gray',lwd=1,cex=0.5,pch=15)
points(input_plot_class_one[,1],-log10(input_plot_class_one[,2]),col='slateblue1',lwd=1,cex=0.5,pch=15)
points(input_plot_class_four[,1],-log10(input_plot_class_four[,2]),col='pink',lwd=1,cex=0.5,pch=15)
abline(v=log2(2), col="black")
abline(v=log2(1/2), col="black")

########################################################


#Figure 3B#############################################

#read the corresponding genes
repressed_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")
activated_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
signle_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
all_the_targets<-c(mutiple_genes,signle_genes)

#get the genes for each class
class_one_list<-repressed_DEGs[repressed_DEGs%in%all_the_targets]
class_two_list<-all_the_targets[!all_the_targets%in%c(repressed_DEGs,activated_DEGs)]
class_three_list<-(c(repressed_DEGs,activated_DEGs)[!c(repressed_DEGs,activated_DEGs)%in%all_the_targets])
class_four_list<-c(activated_DEGs[activated_DEGs%in%all_the_targets])

#KEGG annotation function for A. nidulans
get_percentage<-function(input_genes)
{       
        pathway_genes<-read.table("/Volumes/LiguoDisk/Methods/iSubpathwayMiner-master/an_xml/pathways_ani_v2.xls")
        pathway_genes_table<-table(as.character(pathway_genes[,1]))
        pathway_genes_table_names<-names(pathway_genes_table)
        output_numbers<-NULL
        output_percentages<-NULL
        for ( i in 1:length(pathway_genes_table_names))
        {       
                gene_list_in_path<-pathway_genes[pathway_genes[,1]%in%pathway_genes_table_names[i],2]
                output_numbers[i]<-length(gene_list_in_path[gene_list_in_path%in%input_genes])
                output_percentages[i]<-output_numbers[i]/pathway_genes_table[i]
        }
        return_data_frame<-data.frame(output_numbers,output_percentages)
        rownames(return_data_frame)<-pathway_genes_table_names
        return(return_data_frame)
}

#KEGG annotation results
class_one_list_A<-get_percentage(class_one_list)
class_two_list_A<-get_percentage(class_two_list)
class_three_list_A<-get_percentage(class_three_list)
class_four_list_A<-get_percentage(class_four_list)

#combine the overall KEGG pathway names
output_names<-unique(c(rownames(class_one_list_A[class_one_list_A[,1]>2,]),
rownames(class_two_list_A[class_two_list_A[,1]>2,]),
rownames(class_three_list_A[class_three_list_A[,1]>2,]),
rownames(class_four_list_A[class_four_list_A[,1]>2,])))
write.table(output_names,file="output_names_processes.xls",row.names=FALSE,col.names=FALSE,quote=FALSE)

#organize the matrix for KEGG heatmap plot
matrix_output_for_plot<-matrix(0,nrow=length(output_names),ncol=4)
for ( i in 1:length(output_names))
{
        matrix_output_for_plot[i,1]<-class_one_list_A[rownames(class_one_list_A)==output_names[i],2]
        matrix_output_for_plot[i,2]<-class_two_list_A[rownames(class_two_list_A)==output_names[i],2]
        matrix_output_for_plot[i,3]<-class_three_list_A[rownames(class_three_list_A)==output_names[i],2]
        matrix_output_for_plot[i,4]<-class_four_list_A[rownames(class_four_list_A)==output_names[i],2]
}

matrix_output_for_plot_v2<-matrix_output_for_plot
for ( i in 1:length(kegg_names_final_out_v2_cluster))
{
        matrix_output_for_plot_v2[i,]<-matrix_output_for_plot[kegg_names_final_out%in%kegg_names_final_out_v2_cluster[i],]
}
ratio_matrix<-matrix_output_for_plot_v2

#plot the figure
my_palette <- colorRampPalette(c("whitesmoke","whitesmoke","red", "red"))(n = 100)
heatmap.2(as.matrix(ratio_matrix),symm=F,scale='none',symkey=F,symbreaks=F,sepcolor="gray",trace=c("none"),cexRow = 0.6,density.info=c("none"),col=my_palette[30:100],Colv=FALSE,Rowv=FALSE,sepwidth=c(0.01,0.0001), colsep=0:ncol(as.matrix(ratio_matrix)),rowsep=0:nrow(as.matrix(ratio_matrix)),breaks=seq(0,0.5,length.out=72))

##########################################


#Figure 3C#################################

#simple barplot figure

##########################################


#Figure 3D#################################

#simple barplot figure

##########################################









