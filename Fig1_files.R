

###Fig1A volcano plot###############################

source("make_tiles_functions.R")

#read gene expression files
full_list_files<-read.table("ddsFull_fpkm_values_glu_with_p_values.xls",sep="\t",row.names=1,header=TRUE)
glucose_expression<-data.frame(
full_list_files$MH11036.glucose.1_S15_L002_Rhat.bam,
full_list_files$MH11036.glucose.2_S23_L002_Rhat.bam,
full_list_files$MH11036.glucose.3_S31_L002_Rhat.bam,
full_list_files$CreA.deltion.glucose.1_S19_L002_Rhat.bam,
full_list_files$CreA.deltion.glucose.2_S27_L002_Rhat.bam,
full_list_files$CreA.deltion.glucose.3_S34_L002_Rhat.bam)

#plot gene expression as points
plot(log10(rowMeans(glucose_expression[,4:6])/rowMeans(glucose_expression[,1:3])),-log10(full_list_files$padj),pch=15,cex=0.3,ylim=c(0,200),xlim=c(-5,5))

#color the upregulated or downregulated genes
points(
        log10(rowMeans(glucose_expression[,4:6])/rowMeans(glucose_expression[,1:3]))[as.character(full_list_files$Row.names)%in%read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")],
        -log10(full_list_files$padj)[as.character(full_list_files$Row.names)%in%read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/repressed_genes.xls")],
        pch=15,
        col='tomato',
        lwd=1,
        cex=0.4
        )
points(
        log10(rowMeans(glucose_expression[,4:6])/rowMeans(glucose_expression[,1:3]))[as.character(full_list_files$Row.names)%in%read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")],
        -log10(full_list_files$padj)[as.character(full_list_files$Row.names)%in%read_files_into_character("/Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/final_plots/Fig1/activated_genes.xls")],
        pch=15,
        col='skyblue2',
        lwd=1,
        cex=0.4
        )

abline(v=log10(2), col="black")
abline(v=log10(1/2), col="black")

####################################################

###Fig1B barplot for KEGG selected process##########

#up-regulated genes
barplot(value_input_up,space=0.5,col="tomato",ylim=c(0,100))
#down-regulated genes
barplot(value_input_down,space=0.5,col="skyblue2",ylim=c(0,100))

####################################################

