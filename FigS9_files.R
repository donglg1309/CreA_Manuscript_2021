
#Figure S9A######################


#read gene expression values into R
full_list_files<-read.table("ddsFull_fpkm_values_glu_with_p_values.xls",sep="\t",row.names=1,header=TRUE)
glucose_WT_input<-full_list_files[,c(19+7,20+7,21+7,7+7,8+7,9+7)]

#plot the correlation figure
dev.new(width=8, height=6)
par(mar=c(0.5,0.5,0.5,0.5),mfrow=c(6,6),xpd=FALSE)
for ( i in 1:6)
{
        for ( j in 1:6 )
        {
        plot(glucose_WT_input[,i],glucose_WT_input[,j],log="xy",xlim=c(0.01,10000),xaxt="n",yaxt="n",xlab=NA,ylab=NA,ylim=c(0.01,50000),pch=15,cex=0.2)
        text(0.1,4500,labels=round(cor(glucose_WT_input[,i],glucose_WT_input[,j]),digit = 1),cex=1)
        }
}




################################

#Figure S9B######################

full_list_files<-read.table("ddsFull_fpkm_values_glu_with_p_values.xls",sep="\t",row.names=1,header=TRUE)

# do PCA analysis
pcA_RNA_seq<-prcomp(t(full_list_files[,c(19+7,20+7,21+7,7+7,8+7,9+7)]))

# plot the principal component score
plot(pcA_RNA_seq$rotation[,c(1,2)],pch=15,cex=1,xlim=c(-1,1),ylim=c(-1,1),col=c(rep("red",3),rep("blue",3)))


################################

