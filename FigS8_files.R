
#Figure S8####################

#The data and functions following other Figure script
full_list_transcription_factors_v2<-integrate_list[integrate_list%in%full_list_transcription_factors]
TF_expression_values<-plot_heatmaps_RNA_seq_v2(full_list_transcription_factors_v2)[,1]
TF_binding_level<-plot_heatmaps_bindings(full_list_transcription_factors_v2)[,2]
my_palette <- colorRampPalette(c("skyblue2", "skyblue2","whitesmoke","whitesmoke","tomato", "tomato"))(n = 100)
ratio_matrix_plot<-t(matrix(TF_expression_values[order(TF_expression_values,decreasing=TRUE)],nrow=4))
input_values<-t(matrix(convert_gene_list(full_list_transcription_factors_v2)[order(TF_expression_values,decreasing=TRUE)],nrow=4))
input_values<-as.expression(lapply(full_list_transcription_factors_v2[order(TF_expression_values)],function(a) bquote(italic(.(a)))))
TF_binding_level_plot<-t(matrix(TF_binding_level[order(TF_expression_values,decreasing=TRUE)],nrow=4))
TF_binding_level_plot[TF_binding_level_plot==1]="_        "
TF_binding_level_plot[TF_binding_level_plot==0]=""

#plot the figures
heatmap.2(ratio_matrix_plot,symm=T,scale='none',symkey=T,symbreaks=T,sepcolor="gray",trace=c("none"),cexRow = 1.3,density.info=c("none"),col=my_palette[1:100],Colv=FALSE,Rowv=FALSE,breaks=seq(-4,4,length.out=101),colsep=0:ncol(as.matrix(ratio_matrix_plot)),rowsep=0:nrow(as.matrix(ratio_matrix_plot)),sepwidth=c(0.001,0.00001),cellnote=TF_binding_level_plot,notecol="black",notecex=1.2)

##############################
