
#Figure 2B#############################

CreA_summit_annotation<-read.table("CreA_summits_and_annotation_v12_summit_files.bed_annotation_200bp",sep="\t",header=TRUE)
hist(CreA_summit_annotation$V10,breaks =200,col="blue",xlim=c(-3000,1000))

########################################

#Figure 2D##############################

#read repressed genes
repressed_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")
activated_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")

#read CreA targets
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
single_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
all_the_targets<-c(mutiple_genes,single_genes)

#direct effects
repressed_targets<-all_the_targets[all_the_targets%in%repressed_DEGs]
activated_targets<-all_the_targets[all_the_targets%in%activated_DEGs]

#plot pie figure
pie(length(repressed_targets)/length(all_the_targets),length(activated_targets)/length(all_the_targets),length(all_the_targets_others)/length(all_the_targets),col=c("tomato","lightblue","whitesmoke"))

########################################


#Figure 2E##############################

source("/Users/dongliguo/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")

#read gene expression files
glucose_gene_expression<-read.table("glucose_gene_expression.xls",row.names=1,header=TRUE)

#input the used gene list
gene_list_for_plot<-c(
"amyR",
"alcR",
"xlnR",
"rhaR")

#get formal A. nidulans gene names
gene_list_for_plot_normal<-comvert_name(gene_list_for_plot)


#find gene expression for the corresponding genes
matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=6)
for ( i in 1:length(gene_list_for_plot_normal))
{
        matrix_values[i,]<-as.numeric(glucose_gene_expression[row.names(glucose_gene_expression)%in%gene_list_for_plot_normal[i],])
}

#get the ratio of gene expression changes, and gene expression mean valus, and plot the figure
ratio_matrix<-log2(cbind(rowMeans(matrix_values[,1:3])/rowMeans(matrix_values[,1:3]),rowMeans(matrix_values[,4:6])/rowMeans(matrix_values[,1:3])))
ratio_matrix[ratio_matrix>=4]=4
options(digits=4)
library("gplots")
ratio_matrix_values<-data.frame(rowMeans(matrix_values[,1:3]),rowMeans(matrix_values[,4:6]))
ratio_matrix_values[ratio_matrix_values>100]=100
ratio_matrix_values_ori<-data.frame(rowMeans(matrix_values[,1:3]),rowMeans(matrix_values[,4:6]))
my_palette <- colorRampPalette(c("blue", "blue","whitesmoke","whitesmoke","orange", "orange"))(n = 100)
heatmap.2(as.matrix(ratio_matrix_values),symm=T,scale='none',symkey=T,symbreaks=T,sepcolor="black",trace=c("none"),cexRow = 1.3,density.info=c("none"),colsep=c(0,1,2),rowsep=0:(nrow(as.matrix(ratio_matrix_values))+1),sepwidth=c(0.001,0.0001),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,cellnote=round(ratio_matrix_values_ori, digits=1),notecol="black",notecex=1.5)

##############################################


