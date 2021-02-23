
#Figure 4A############

#get the input info
CreA_targets<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v8/CreA_summits_and_annotation_v12.xls")
integrate_list<-as.character(read.table("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig3/all_intergrate_genes.xls")$V1)
full_list_files<-read.table("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/Fig1_supplemental_files/glucose_gene_expression/ddsFull_fpkm_values_glu_with_p_values.xls",sep="\t",row.names=1,header=TRUE)
repressed_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")
activated_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
signle_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
all_the_targets<-c(mutiple_genes,signle_genes)

#seperate the classes
class_one_list<-repressed_DEGs[repressed_DEGs%in%all_the_targets]
class_two_list<-all_the_targets[!all_the_targets%in%c(repressed_DEGs,activated_DEGs)]
class_three_list<-(c(repressed_DEGs,activated_DEGs)[!c(repressed_DEGs,activated_DEGs)%in%all_the_targets])
class_four_list<-c(activated_DEGs[activated_DEGs%in%c(all_the_targets)])

#initation of the analyses
SM_gene_A<-NULL
SM_cluster_A<-NULL
for ( i in 1:length(read_files_into_character("secondary_metabolism_cluster.xls")))
{       
        SM_gene_A<-c(SM_gene_A,integrate_list[integrate_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])
        SM_cluster_A<-c(SM_cluster_A,rep(read_files_into_character("secondary_metabolism_cluster.xls")[i],length(integrate_list[integrate_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])))
}
write.table(data.frame(SM_gene_A,convert_gene_list(SM_gene_A),SM_cluster_A,matrix_all_include),file="gene_matrix_all.xls",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)


#put the raw info
matrix_all_include<-matrix(0,nrow=length(SM_gene_A),ncol=7)
matrix_all_include[SM_gene_A%in%all_the_targets,1]<-1
matrix_all_include[SM_gene_A%in%repressed_DEGs,2]<-1
matrix_all_include[SM_gene_A%in%activated_DEGs,3]<-1
matrix_all_include[SM_gene_A%in%class_one_list,4]<-1
matrix_all_include[SM_gene_A%in%class_two_list,5]<-1
matrix_all_include[SM_gene_A%in%class_three_list,6]<-1
matrix_all_include[SM_gene_A%in%class_four_list,7]<-1


#get the percentage info for each class
percentage_number_class_one<-NULL
numbers_all_class_one<-NULL
for ( i in 1:length(read_files_into_character("secondary_metabolism_cluster.xls")))
{
        percentage_number_class_one[i]<-length(class_one_list[class_one_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])/length(unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i])))
        numbers_all_class_one[i]<-length(class_one_list[class_one_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])
}
percentage_number_class_two<-NULL
numbers_all_class_two<-NULL
for ( i in 1:length(read_files_into_character("secondary_metabolism_cluster.xls")))
{
        percentage_number_class_two[i]<-length(class_two_list[class_two_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])/length(unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i])))
        numbers_all_class_two[i]<-length(class_two_list[class_two_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])
}
percentage_number_class_three<-NULL
numbers_all_class_three<-NULL
for ( i in 1:length(read_files_into_character("secondary_metabolism_cluster.xls")))
{
        percentage_number_class_three[i]<-length(class_three_list[class_three_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])/length(unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i])))
        numbers_all_class_three[i]<-length(class_three_list[class_three_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])
}
percentage_number_class_four<-NULL
numbers_all_class_four<-NULL
for ( i in 1:length(read_files_into_character("secondary_metabolism_cluster.xls")))
{
        percentage_number_class_four[i]<-length(class_four_list[class_four_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])/length(unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i])))
        numbers_all_class_four[i]<-length(class_four_list[class_four_list%in%unique(read_files_into_character(read_files_into_character("secondary_metabolism_cluster.xls")[i]))])
}

#summary and plot the figure
percentage_all_matrix<-data.frame(percentage_number_class_one,percentage_number_class_two,percentage_number_class_three,percentage_number_class_four)
barplot(t(as.matrix(percentage_all_matrix))[,order(rowSums(percentage_all_matrix))],col=c("blue","gray","turquoise","pink"))

###############################


#Figure 4B############

#read raw gene expression values
glucose_gene_expression<-read.table("glucose_gene_expression.xls",row.names=1,header=TRUE)

#read DBA gene cluster names
gene_list_for_plot<-c(
"dbaA",
"dbaB",
"dbaC",
"dbaD",
"dbaF",
"dbaG",
"dbaH",
"pkeA")
gene_list_for_plot_normal<-comvert_name(gene_list_for_plot)

#get the dba gene cluster expression
matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=6)
for ( i in 1:length(gene_list_for_plot_normal))
{
        matrix_values[i,]<-as.numeric(glucose_gene_expression[row.names(glucose_gene_expression)%in%gene_list_for_plot_normal[i],])
}

#plot the heatmap figure
ratio_matrix<-log2(cbind(rowMeans(matrix_values[,1:3])/rowMeans(matrix_values[,1:3]),rowMeans(matrix_values[,4:6])/rowMeans(matrix_values[,1:3])))
ratio_matrix[ratio_matrix>=4]=4
options(digits=4)
library("gplots")
ratio_matrix_values<-data.frame(rowMeans(matrix_values[,1:3]),rowMeans(matrix_values[,4:6]))
ratio_matrix_values[ratio_matrix_values>100]=100
ratio_matrix_values_ori<-data.frame(rowMeans(matrix_values[,1:3]),rowMeans(matrix_values[,4:6]))
my_palette <- colorRampPalette(c("blue", "blue","whitesmoke","whitesmoke","orange", "orange"))(n = 100)
heatmap.2(as.matrix(ratio_matrix_values),symm=T,scale='none',symkey=T,symbreaks=T,sepcolor="black",trace=c("none"),cexRow = 1.3,density.info=c("none"),colsep=c(0,1,2),rowsep=0:(nrow(as.matrix(ratio_matrix_values))+1),sepwidth=c(0.001,0.0001),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,cellnote=round(ratio_matrix_values_ori, digits=1),notecol="black",notecex=1.5)

######################################


#Figure 4D############################
sid_genes_all<-comvert_name(read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig4/sid_genes_all.xls"))

#binding info
mutiple_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_mutiple_binding.xls")
signle_genes<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/gene_list_with_single_binding.xls")
all_the_targets<-c(mutiple_genes,signle_genes)

#gene expression info
full_list_files<-read.table("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/Fig1_supplemental_files/glucose_gene_expression/ddsFull_fpkm_values_glu_with_p_values.xls",sep="\t",row.names=1,header=TRUE)
repressed_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")
activated_DEGs<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")

#function for getting gene expression ratio changes
plot_heatmaps_RNA_seq_v2<-function(gene_list_for_plot_normal)
{
        matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=6)
        for ( i in 1:length(gene_list_for_plot_normal))
        {
                matrix_values[i,]<-as.numeric(full_list_files[as.character(full_list_files[,c(1)])%in%gene_list_for_plot_normal[i],c(26,27,28,14,15,16)])
        }
        matrix_values[matrix_values=0]=0.0000000001
        ratio_matrix<-log2(cbind(rowMeans(matrix_values[,4:6])/rowMeans(matrix_values[,1:3])))
        ratio_matrix[ratio_matrix>4]=4
        ratio_matrix[ratio_matrix<c(-4)]=-4
        ratio_matrix_v2<-data.frame(ratio_matrix,ratio_matrix)
        return(as.matrix(ratio_matrix_v2))
}

#function for getting binding info
plot_heatmaps_bindings<-function(gene_list_for_plot_normal)
{
        matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=2)
        matrix_values[gene_list_for_plot_normal%in%all_the_targets,]<-1
        return(as.matrix(matrix_values))
}

##plot for gene expression changes
my_palette <- colorRampPalette(c("blue", "blue","whitesmoke","whitesmoke","red", "red"))(n = 100)
ratio_matrix_plot<-plot_heatmaps_RNA_seq_v2(sid_genes_all)
heatmap.2(ratio_matrix_plot,symm=T,scale='none',symkey=T,symbreaks=T,sepcolor="black",trace=c("none"),cexRow = 1.3,density.info=c("none"),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,breaks=seq(-4,4,length.out=101),colsep=c(0,2),rowsep=0:nrow(as.matrix(ratio_matrix_plot)),sepwidth=c(0.01,0.0001))

#plot for binding info
ratio_matrix_plot<-plot_heatmaps_bindings(sid_genes_all)
my_palette <- colorRampPalette(c("black","black","yellow", "yellow"))(n = 100)
heatmap.2(ratio_matrix_plot,symm=T,scale='none',symkey=T,symbreaks=T,sepcolor="black",trace=c("none"),cexRow = 1.3,density.info=c("none"),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,breaks=seq(0,1,length.out=101),colsep=(0,2),rowsep=0:nrow(as.matrix(ratio_matrix_plot)),sepwidth=c(0.01,0.0001))


######################################


#Figure 4E############################

source("/Users/dongliguo/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")
#plot barplot for plate test data
gene_list_for_plot<-read.table("sexual_plot_data.xls")
dev.new(width=10, height=4,5)
par(mar=c(3,5,3,3),mfrow=c(1,4),xpd=FALSE)
matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=6)
for ( i in 1:length(gene_list_for_plot_normal))
{
        y1<-as.matrix(read.table("sexual_plot_data.xls"))
        y1.means <- apply(y1,2,mean)
        y.means <- max(y1.means)
        y1.sd <- apply(y1,2,sd)
        barx <- barplot(y1.means, col=c("black","gray"), srt=45,axis.lty=1,axes=FALSE,ylab=NULL,ylim=c(0,3000))
        error.bar(barx,y1.means, y1.sd,lwd=0.7)
        axis(2,col.axis="black",at=c(0,1500,3000),labels=c(0,1500,3000),lwd=1,las=2,cex.axis=2)
}
dev.new(width=10, height=4,5)
y1<-as.matrix(read.table("T7_S4_plate_test_sexual.xls"))
y1.means <- apply(y1,2,mean)
y.means <- max(y1.means)
y1.sd <- apply(y1,2,sd)
barx <- barplot(y1.means, col=c("black","gray"), srt=45,axis.lty=1,axes=FALSE,ylab=NULL,ylim=c(0,300))
error.bar(barx,y1.means, y1.sd,lwd=0.7)
axis(2,col.axis="black",at=c(0,150,300),labels=c(0,150,300),lwd=1,las=2,cex.axis=2)

######################################




#Figure 4F############################

#read gene list files
gene_list_with_mutiple_bindings<-as.character(read.table("gene_list_with_mutiple_binding.xls")$V1)
gene_list_with_single_bindings<-as.character(read.table("gene_list_with_single_binding.xls")$V1)

#read gene expression values
RNA_seqs_gene_expression_values<-read.table("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig3/glucose_gene_expression.xls",row.names=1,header=TRUE)
gene_list_all<-row.names(RNA_seqs_gene_expression_values)
all_the_gene_names<-read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/A_nidulan_10_gene_final_annotation2_names.xls")

#get gene binding density values
tiles_2<-make_tiles_for_promoter_ATG("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig2/A_nidulan_10_gene_final_annotation2_names.xls")
CreA_glucose_GFP<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_S4_CreA_1.bam.bed")
CreA_glucose<-normalization_tiles_for_promoter("/Volumes/LiguoDisk/CreA_RNA_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_HA_CreA_1.bam.bed")

#read the file path for the used gene list
all_the_input_files_metablism<-read_files_into_character("all_the_input_files_names_oxygen_plot.xls")

#get binding density values
plot_figures_out<-function(input_list)
{       
        input_matrix<-CreA_glucose_matrix[1:length(input_list),]
        for ( i in 1:length(input_list))
        {       
                input_matrix[i,]=CreA_glucose_matrix[all_the_gene_names%in%input_list[i],]
        }       
return(input_matrix)
}       

#match the colors for the binding values
match_colors2<-function(input_files,color_values,length_input)
{       
        output_colors<-NULL
        the_continous_values<-seq(from=1,to=2000,length.out=length_input)
        for ( i in 1:length(input_files))
        {       
                for ( j in 1:length(the_continous_values) )
                {       
                        if(input_files[i]==2000)
                        {       
                                output_colors[i]=color_values[length_input]
                                break   
                        }       
                        if(input_files[i]<=the_continous_values[j+1]&&input_files[i]>=the_continous_values[j])
                        {       
                                output_colors[i]<-color_values[j]
                                break   
                        }       
                }       
        }       
return(output_colors)
}       

#match color for gene expression values
match_colors<-function(input_files,color_values,length_input)
{
	output_colors<-NULL
	the_continous_values<-seq(from=-3,to=3,length.out=length_input)
        for ( i in 1:length(input_files))
        {       
                for ( j in 1:length(the_continous_values) )
                {       
                        if(input_files[i]==3)
                        {output_colors[i]=color_values[length_input]
                        break}
                        if(input_files[i]<=the_continous_values[j+1]&&input_files[i]>=the_continous_values[j])
                        {
                        output_colors[i]<-color_values[j]
                        break
                        }
                }
        }
return(output_colors)
}

#initiate the plot 
plot (c(0,20000),c(1400,1400),col = "gray80",xlim=c(0,10000),ylim=c(-1,1400),type="l",lwd=1,xaxt="n",axes=FALSE)
#get the binding density values
for ( input_file_names in 1:length(all_the_input_files_metablism))
{
	
	#get the corresponding gene name
	acetate_key_genes2<-read_files_into_character(all_the_input_files_metablism[input_file_names])
	#get the binding values
	CreA_glucose_matrix<-cbind(rowMeans(cbind(CreA_glucose[,72:214],CreA_glucose_GFP[,72:214])),rowMeans(cbind(CreA_glucose[,72:214],CreA_glucose_GFP[,72:214])))
	CreA_glucose_matrix[CreA_glucose_matrix>1000]=2000
	binding_values<-plot_figures_out(acetate_key_genes2)[,1]
        binding_values[!acetate_key_genes2%in%c(gene_list_with_mutiple_bindings,gene_list_with_single_bindings)]=1
        binding_values[acetate_key_genes2%in%c(gene_list_with_mutiple_bindings,gene_list_with_single_bindings)]=2000
	#get the expression values
	output_gene_expression<-matrix(0,nrow=length(acetate_key_genes2),ncol=6)
	output_gene_list<-NULL
	d=1
	for ( i in 1:length(acetate_key_genes2))
	{
        	output_gene_expression[i,]<-as.numeric(RNA_seqs_gene_expression_values[gene_list_all%in%acetate_key_genes2[i],])
	}
	gene_expression_alteration<-rowMeans(output_gene_expression[,4:6])/rowMeans(output_gene_expression[,1:3])
	gene_expression_alteration[is.na(gene_expression_alteration)]=1
	gene_expression_alteration[is.infinite(gene_expression_alteration)]=10
	gene_expression_alteration2<-log2(gene_expression_alteration)
	gene_expression_alteration2[gene_expression_alteration2>3]=3
	gene_expression_alteration2[gene_expression_alteration2<c(-3)]=-3
	#reorder the overall values according to binding and gene expresiion changes	
	binding_values_out_formal<-NULL
	expression_values_out_formal<-NULL
	gene_list_output_formal<-NULL
	for ( i in 1:length(table(binding_values)))
	{       
        	if ( i == 1 )
        	{       
                	binding_values_out<-binding_values[binding_values>=as.numeric(names(rev(table(binding_values)))[1])]
                	expression_values_out<-gene_expression_alteration2[binding_values>=as.numeric(names(rev(table(binding_values)))[1])]
                	gene_list_output_out<-acetate_key_genes2[binding_values>=as.numeric(names(rev(table(binding_values)))[1])]
                	binding_values_out_formal<-binding_values_out[order(expression_values_out,decreasing=TRUE)]
                	expression_values_out_formal<-expression_values_out[order(expression_values_out,decreasing=TRUE)]
                	gene_list_output_formal<-gene_list_output_out[order(expression_values_out,decreasing=TRUE)]
                	binding_values_out<-NULL
                	expression_values_out<-NULL
                	gene_list_output_out<-NULL
        	}
                
        	if ( i != 1 )
        	{       
                	binding_values_out<-binding_values[binding_values>=as.numeric(names(rev(table(binding_values)))[i])&binding_values<as.numeric(names(rev(table(binding_values)))[i-1])]
                	expression_values_out<-gene_expression_alteration2[binding_values>=as.numeric(names(rev(table(binding_values)))[i])&binding_values<as.numeric(names(rev(table(binding_values)))[i-1])]      
                	gene_list_output_out<-acetate_key_genes2[binding_values>=as.numeric(names(rev(table(binding_values)))[i])&binding_values<as.numeric(names(rev(table(binding_values)))[i-1])]
                	binding_values_out_formal<-c(binding_values_out_formal,binding_values_out[order(expression_values_out,decreasing=TRUE)])
                	expression_values_out_formal<-c(expression_values_out_formal,expression_values_out[order(expression_values_out,decreasing=TRUE)])
                	gene_list_output_formal<-c(gene_list_output_formal,gene_list_output_out[order(expression_values_out,decreasing=TRUE)])
                	binding_values_out<-NULL
                	expression_values_out<-NULL
                	gene_list_output_out<-NULL
        	}
	}
	#plot the gene expression changes
	the_color_used_v1<-match_colors(expression_values_out_formal,my_palette[20:80],60)
	for ( i in 1:length(the_color_used_v1))
	{
		points (c(10+50*i,10+50*i),c(1400-100*input_file_names,1400-100*input_file_names),pch=22,col="gray",bg=the_color_used_v1[i],lwd=1.5,cex=1.5)
	}
	#plot binding info
	the_color_used_v2<-match_colors2(binding_values_out_formal,my_palette2[1:100],100)
	for ( i in 1:length(the_color_used_v2))
	{
		points (c(10+50*i,10+50*i),c(1370-100*input_file_names,1370-100*input_file_names),pch=21,col="gray",bg=the_color_used_v2[i],lwd=1.5,cex=1.5)
	}
	#writeout the gene list
	write.table(convert_gene_list(gene_list_output_formal),file=paste(all_the_input_files_metablism[input_file_names],"_ordered.xls",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}

################################

#Figure 4G######################

#programming is same with Figure 4E

################################

#Figure 4H######################

#function for getting gene expression changes
plot_heatmaps_RNA_seq<-function(gene_list_for_plot_normal)
{
	matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=6)
	for ( i in 1:length(gene_list_for_plot_normal))
	{       
        	matrix_values[i,]<-as.numeric(full_list_files[as.character(full_list_files[,c(1)])%in%gene_list_for_plot_normal[i],c(26,27,28,14,15,16)])
	}       
        matrix_values[matrix_values=0]=0.0000000001
        ratio_matrix<-log2(cbind(matrix_values[,1:3]/matrix_values[,1:3],matrix_values[,4:6]/matrix_values[,1:3]))
        combined_DEGs<-c(repressed_DEGs,activated_DEGs)
        ratio_matrix[!gene_list_for_plot_normal%in%combined_DEGs]<-0
        ratio_matrix[ratio_matrix>4]=4
        ratio_matrix[ratio_matrix<c(-4)]=-4
        return(as.matrix(ratio_matrix))
}

#function for getting gene expression changes
plot_heatmaps_RNA_seq_v2<-function(gene_list_for_plot_normal)
{
	matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=6)
	for ( i in 1:length(gene_list_for_plot_normal))
	{
        	matrix_values[i,]<-as.numeric(full_list_files[as.character(full_list_files[,c(1)])%in%gene_list_for_plot_normal[i],c(26,27,28,14,15,16)])
	}
        matrix_values[matrix_values=0]=0.0000000001
        ratio_matrix<-log2(cbind(rowMeans(matrix_values[,4:6])/rowMeans(matrix_values[,1:3])))
        ratio_matrix[ratio_matrix>4]=4
        ratio_matrix[ratio_matrix<c(-4)]=-4
        ratio_matrix_v2<-data.frame(ratio_matrix,ratio_matrix)
        return(as.matrix(ratio_matrix_v2))
}

#function for getting binding info
plot_heatmaps_bindings<-function(gene_list_for_plot_normal)
{
	matrix_values<-matrix(0,nrow=length(gene_list_for_plot_normal),ncol=2)
	matrix_values[gene_list_for_plot_normal%in%all_the_targets,]<-1
	return(as.matrix(matrix_values))
}

#read gene list files
sexual_gene_list_used<-comvert_name(read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/sexual_asexual_genes/sexual_genes_list_specific.xls"))
#plot gene expression changes
ratio_matrix_plot<-plot_heatmaps_RNA_seq_v2(sexual_gene_list_used)[order(plot_heatmaps_RNA_seq_v2(sexual_gene_list_used)[,1]),]
my_palette <- colorRampPalette(c("blue", "blue","whitesmoke","whitesmoke","red", "red"))(n = 100)
heatmap.2(ratio_matrix_plot,symm=T,scale='none',symkey=T,symbreaks=T,trace=c("none"),cexRow = 1.3,density.info=c("none"),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,breaks=seq(-4,4,length.out=101))
ratio_matrix_plot<-original_expression_values_v2
#plot raw gene expression values
my_palette <- colorRampPalette(c("blue","whitesmoke","whitesmoke","whitesmoke", "orange"))(n = 100)
original_expression_values<-RNA_seq_gene_expression(sexual_gene_list_used)[order(plot_heatmaps_RNA_seq_v2(sexual_gene_list_used)[,1]),]
original_expression_values_v2<-data.frame(rowMeans(original_expression_values[,1:3]),rowMeans(original_expression_values[,4:6]))
original_expression_values_v2[original_expression_values_v2>=100]=100
ratio_matrix_plot<-original_expression_values_v2
heatmap.2(as.matrix(ratio_matrix_plot),symm=T,scale='none',symkey=T,symbreaks=T,trace=c("none"),cexRow = 1.3,density.info=c("none"),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,breaks=seq(-100,100,length.out=101),colsep=c(0,1,2),rowsep=0:nrow(as.matrix(ratio_matrix_plot)),sepcolor="gray",sepwidth=c(0.01,0.0001))
write.table(convert_gene_list(sexual_gene_list_used)[order(plot_heatmaps_RNA_seq_v2(sexual_gene_list_used)[,1])],file="sexual_ordered_genes.xls",row.names=FALSE,col.names=FALSE,quote=FALSE)
#plot binding inform
ratio_matrix_plot<-plot_heatmaps_bindings(sexual_gene_list_used)[order(plot_heatmaps_RNA_seq_v2(sexual_gene_list_used)[,1]),]
my_palette <- colorRampPalette(c("blue","blue","blue","blue","blue","whitesmoke","whitesmoke","yellow","yellow", "yellow","yellow", "yellow"))(n = 100)
the_color_used<-match_colors(ratio_matrix_plot[,1],my_palette[1:100],100)
plot (c(0,20000),c(1400,1400),col = "gray80",xlim=c(0,5200),ylim=c(-1,800),type="l",lwd=1,xaxt="n",axes=FALSE)
for ( i in 1:length(the_color_used))
{
points (c(1,1),c(10+20*i,10+20*i),pch=21,col="black",bg=the_color_used[i],lwd=1,cex=1.3)
}

#the asexual program is same with the sexual one, just gene list is different

#################################################################


#Figure 4J#################

#programming is same with Figure 4E

################################



