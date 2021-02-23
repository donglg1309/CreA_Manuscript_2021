
########Colors###########

#colors = colorRampPalette(c('white','blue','blue2','blue4'))(200)
#library(RColorBrewer)
#mypalette<-brewer.pal(9,"Blues")
#colors = colorRampPalette(c('black','black','yellow','orange'))(200)
#

#####get R libraries for plots########
#library(matrixStats)

#source("https://bioconductor.org/biocLite.R")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.11")
#install.packages("matrixStats")
#library("BiocManager")
#biocLite("GenomicRanges")
#biocLite("rtracklayer")
#biocLite("IRanges")
#install.packages("gmp")
#install.packages("gplots")
#install.packages("DescTools")
#install.packages("gplots")
#install.packages("DescTools")

library(matrixStats)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
#library(gmp) 
library(gplots)

#BiocManager::install(c("GenomicRanges", "rtracklayer","IRanges"))


######################################

enrich_pvalue <- function(N, A, B, k)
{
     m <- A + k
     n <- B + k
     i <- k:min(m,n)

     as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}


make_tiles_for_TES<-function(x){
        
	gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)
        
	gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
	chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
		for ( i in 1:length(asd_in_gene_name) )
        	{
                	chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
        	}
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
	A_value<-NULL
        	for( i in 1:length(end(ranges(gene_select_region_CreA))))
        	{
                	A_value[i]=0
        	}
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        	if(gene_select_region_CreA$score[i]=="1")
                	A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]-500,end(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        	else
                	A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+500,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
                        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }

######################################

make_tiles_for_promoter_novel<-function(x)
        {
        gene_select_region_CreA<-import.bed(x)
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }


#####################################

make_tiles_for_promoter<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
	x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }

###############################################

make_tiles_for_promoter_ATG<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/a_nidulans_exon_CDS_files_data_names.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }

##############################################

make_tiles_for_body_with_long_promoters_ATG<-function(x){
       gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/a_nidulans_exon_CDS_files_data_names.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_feature_peaks<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_feature_peaks))))
        {
                A_value[i]=0
        }
        tiles=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i],length.out=201)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i],length.out=201))
        tiles2=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i]-3000,start(ranges(gene_feature_peaks))[i],length.out=300)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i]+3000,end(ranges(gene_feature_peaks))[i],length.out=300))
        tiles3=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i]+500,length.out=60)
                else
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i]-500,length.out=60))
        tiles_1<-NULL
        tiles_1<-matrix(NA, nrow = length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1]), ncol = length(tiles[1,]))
        for ( i in 1:length(tiles2[,1]))
        {
                tiles_1[i,]<-tiles2[i,]
        }
        for ( i in (1+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])))
        {
                tiles_1[i,]<-tiles[i-length(tiles2[,1]),]
        }
        for ( i in (1+length(tiles[,1])+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1])))
        {
                tiles_1[i,]<-tiles3[i-(length(tiles[,1])+length(tiles2[,1])),]
        }
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(gene_feature_peaks$name, each=561), 1:561, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_feature_peaks)),each=561)),
                        ranges=IRanges(start=as.vector(tiles_1),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles_1)))),
                        seqinfo=si)
        return(tiles_2)
}





###############################################

make_tiles_for_promoter_narrow_size<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }


###############################################################################

make_tiles_for_promoter_sc<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/S_cerevisiae_sacCer3/annotation/saccharomyces_cerevisiae_R64-1-1_20110208_without_tRNA_convert.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/S_cerevisiae_sacCer3/genome_size.xls")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }

#####################################################################################

make_tiles_for_promoter_sp<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/S_pombe_ASM294v2.31/annotation/Schizosaccharomyces_pombe.ASM294v2.31_convert.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=321)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=321))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/S_pombe_ASM294v2.31/genome_size.xls")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=321)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }

###################################################################################


make_tiles_for_promoter_precise<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=4501)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=4501))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=4501), 1:4501, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=4501)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=2),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }



################################################################################

make_tiles_for_promoter2<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-3000,start(ranges(gene_select_region_CreA))[i]+1500,length.out=4500)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+3000,end(ranges(gene_select_region_CreA))[i]-1500,length.out=4500))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=4500), 1:4500, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=4500)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=1),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }



#######################################

make_tiles_for_short_promoter<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-500,start(ranges(gene_select_region_CreA))[i]+0,length.out=21)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+500,end(ranges(gene_select_region_CreA))[i]-0,length.out=21))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=21), 1:21, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=21)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }
######################################

make_tiles_for_short_promoter_utr<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-500,start(ranges(gene_select_region_CreA))[i]+300,length.out=21)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+500,end(ranges(gene_select_region_CreA))[i]-300,length.out=21))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=21), 1:21, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=21)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
        }



######################################

make_tiles_for_summit_with_orietation<-function(x,y){
	
	gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        y1<-read.table(y)
        asd_in_gene_name<-as.character(y1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }

	x1<-read.table(x)
	x2<-x1[,1:4]
	output_strand<-data.frame(x2,chr_region[,5])
        write.table(output_strand,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	CreA_glucose_peaks<-import.bed("inter_file.bed")
	gene_select_region_CreA<-CreA_glucose_peaks
	A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-1500,start(ranges(gene_select_region_CreA))[i]+1500,length.out=301)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+1500,end(ranges(gene_select_region_CreA))[i]-1500,length.out=301))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=301), 1:301, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=301)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
}

make_tiles_for_distance_with_orietation<-function(x,y){

        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        y1<-read.table(y)
        asd_in_gene_name<-as.character(y1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }

        x1<-read.table(x)
        x2<-x1[,1:4]
        output_strand<-data.frame(x2,chr_region[,5])
        write.table(output_strand,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        CreA_glucose_peaks<-import.bed("inter_file.bed")
        gene_select_region_CreA<-CreA_glucose_peaks
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
        if(gene_select_region_CreA$score[i]=="1")
                A_value[i]+seq(start(ranges(gene_select_region_CreA))[i]-500,start(ranges(gene_select_region_CreA))[i]+500,length.out=1001)
        else
                A_value[i]+seq(end(ranges(gene_select_region_CreA))[i]+500,end(ranges(gene_select_region_CreA))[i]-500,length.out=1001))
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_3 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=1001), 1:1001, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=1001)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=1),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        return(tiles_3)
}
########################################

make_tiles_for_distance_without_orietation<-function(x,y){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-499,501,length.out=1001))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=1001), 1:1001, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=1001)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)

}

########################################

make_tiles_for_distance_without_orietation_abso<-function(x,y){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-50,50,length.out=101))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=101), 1:101, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=101)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=1),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)

}

########################################

make_tiles_for_distance_without_orietation_abso2<-function(x,y){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-500,500,length.out=1001))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=1001), 1:1001, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=1001)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=1),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)

}




########################################

make_tiles_for_summit<-function(x){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-1500,1500,length.out=301))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=301), 1:301, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=301)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)
}

#########################################


make_tiles_for_summit_percise<-function(x){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-1500,1500,length.out=3001))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=3001), 1:3001, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=3001)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=1),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)
}



##########################################



########################################

make_tiles_for_summit_sc<-function(x){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-1500,1500,length.out=301))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/S_cerevisiae_sacCer3/genome_size.xls")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=301), 1:301, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=301)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)
}

#######################################

make_tiles_for_summit_sp<-function(x){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-1500,1500,length.out=301))
        )
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/S_pombe_ASM294v2.31/genome_size.xls")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=301), 1:301, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=301)),
                        ranges=IRanges(start=as.vector(tiles),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles)))),
                        seqinfo=si)
}

########################################

make_tiles_for_body_with_long_promoters<-function(x){
       gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_feature_peaks<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_feature_peaks))))
        {
                A_value[i]=0
        }
        tiles=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i],length.out=201)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i],length.out=201))
        tiles2=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i]-3000,start(ranges(gene_feature_peaks))[i],length.out=300)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i]+3000,end(ranges(gene_feature_peaks))[i],length.out=300))
        tiles3=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i]+500,length.out=60)
                else
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i]-500,length.out=60))
        tiles_1<-NULL
        tiles_1<-matrix(NA, nrow = length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1]), ncol = length(tiles[1,]))
        for ( i in 1:length(tiles2[,1]))
        {
                tiles_1[i,]<-tiles2[i,]
        }
        for ( i in (1+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])))
        {
                tiles_1[i,]<-tiles[i-length(tiles2[,1]),]
        }
        for ( i in (1+length(tiles[,1])+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1])))
        {
                tiles_1[i,]<-tiles3[i-(length(tiles[,1])+length(tiles2[,1])),]
        }
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(gene_feature_peaks$name, each=561), 1:561, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_feature_peaks)),each=561)),
                        ranges=IRanges(start=as.vector(tiles_1),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles_1)))),
                        seqinfo=si)
        return(tiles_2)
}


########################################

make_tiles_for_body_with_normal_promoters<-function(x){
       gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {       
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_feature_peaks<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_feature_peaks))))
        {       
                A_value[i]=0
        }
        tiles=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i],length.out=201)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i],length.out=201))
        tiles2=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i]-500,start(ranges(gene_feature_peaks))[i],length.out=60)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i]+500,end(ranges(gene_feature_peaks))[i],length.out=60))
        tiles3=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i]+500,length.out=60)
                else
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i]-500,length.out=60))
        tiles_1<-NULL
        tiles_1<-matrix(NA, nrow = length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1]), ncol = length(tiles[1,]))
        for ( i in 1:length(tiles2[,1]))
        {
                tiles_1[i,]<-tiles2[i,]
        }
        for ( i in (1+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])))
        {
                tiles_1[i,]<-tiles[i-length(tiles2[,1]),]
        }
        for ( i in (1+length(tiles[,1])+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1])))
        {
                tiles_1[i,]<-tiles3[i-(length(tiles[,1])+length(tiles2[,1])),]
        }
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(gene_feature_peaks$name, each=321), 1:321, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_feature_peaks)),each=321)),
                        ranges=IRanges(start=as.vector(tiles_1),
                        width=50),
                        strand=Rle(rep("*",length(as.vector(tiles_1)))),
                        seqinfo=si)
        return(tiles_2)
}

make_tiles_for_body_with_normal_promoters_precise<-function(x){
       gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        asd_in_gene_name<-as.character(x1$V1)

        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_feature_peaks<-import.bed("inter_file.bed")
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_feature_peaks))))
        {
                A_value[i]=0
        }
        tiles=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i],length.out=2010)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i],length.out=2010))
        tiles2=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i]-500,start(ranges(gene_feature_peaks))[i],length.out=600)
                else
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i]+500,end(ranges(gene_feature_peaks))[i],length.out=600))
        tiles3=sapply(1:length(end(ranges(gene_feature_peaks))),function(i)
                if(gene_feature_peaks$score[i]=="1")
                        A_value[i]+seq(end(ranges(gene_feature_peaks))[i],end(ranges(gene_feature_peaks))[i]+500,length.out=600)
                else
                        A_value[i]+seq(start(ranges(gene_feature_peaks))[i],start(ranges(gene_feature_peaks))[i]-500,length.out=600))
        tiles_1<-NULL
        tiles_1<-matrix(NA, nrow = length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1]), ncol = length(tiles[1,]))
        for ( i in 1:length(tiles2[,1]))
        {
                tiles_1[i,]<-tiles2[i,]
        }
        for ( i in (1+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])))
        {
                tiles_1[i,]<-tiles[i-length(tiles2[,1]),]
        }
        for ( i in (1+length(tiles[,1])+length(tiles2[,1])):(length(tiles[,1])+length(tiles2[,1])+length(tiles3[,1])))
        {
                tiles_1[i,]<-tiles3[i-(length(tiles[,1])+length(tiles2[,1])),]
        }
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(gene_feature_peaks$name, each=3210), 1:3210, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_feature_peaks)),each=3210)),
                        ranges=IRanges(start=as.vector(tiles_1),
                        width=2),
                        strand=Rle(rep("*",length(as.vector(tiles_1)))),
                        seqinfo=si)
        return(tiles_2)
}


make_tiles_for_narrow_summit<-function(x){
        CreA_glucose_peaks<-import.bed(x)
        tiles=sapply(1:length(end(ranges(CreA_glucose_peaks))),function(i)
        (end(ranges(CreA_glucose_peaks))[i]+seq(-100,100,length.out=1))
        )

Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
tiles_2 = GRanges(tilename = paste( rep(CreA_glucose_peaks$name, each=1), 1:1, sep="_" ),
seqnames = Rle( rep(as.character(seqnames(CreA_glucose_peaks)),each=1)),
ranges=IRanges(start=as.vector(tiles),
width=200),
strand=Rle(rep("*",length(as.vector(tiles)))),
seqinfo=si)
return(tiles_2)
}







#########################################

normalization_tiles_for_narraw_bins<-function(x)
	{
        	M_0h<-import.bed(x)
        	start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
        	width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
		M_0h.p=countOverlaps(tiles_2,M_0h)
        	M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2),ncol=1,byrow=TRUE)
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

#######################################################

normalization_tiles_for_narraw_bins2<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2),ncol=1,byrow=TRUE)
                M_0h.p.matrix2<-(M_0h.p.matrix-25)
		M_0h.p.matrix2[M_0h.p.matrix2<1]=1
		M_0h.p.matrix=M_0h.p.matrix2
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

######################################

normalization_tiles_without_extension<-function(x)
	{
               M_0h<-import.bed(x)
               M_0h.p=countOverlaps(tiles_2,M_0h)
               M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/321,ncol=321,byrow=TRUE)
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


######################################


normalization_tiles_without_extension2<-function(x)
        {
               M_0h<-import.bed(x)
               M_0h.p=countOverlaps(tiles_2,M_0h)
               M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/1001,ncol=1001,byrow=TRUE)
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


#######################################

normalization_tiles_for_promoter<-function(x)
	{
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/321,ncol=321,byrow=TRUE)
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

#########################################


normalization_tiles_for_promoter_CAGE_forward<-function(x)
        {
                M_0h<-import.bed(x)
                #start(ranges(M_0h))[as.character(strand(M_0h))=="*"]=start(ranges(M_0h))[as.character(strand(M_0h))=="*"]-(30-width(ranges(M_0h))[as.character(strand(M_0h))=="*"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="*"]=5
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/321,ncol=321,byrow=TRUE)
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

normalization_tiles_for_promoter_CAGE_reverse<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="*"]=start(ranges(M_0h))[as.character(strand(M_0h))=="*"]-(5-width(ranges(M_0h))[as.character(strand(M_0h))=="*"])
                #width(ranges(M_0h))[as.character(strand(M_0h))=="*"]=5
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/321,ncol=321,byrow=TRUE)
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



#########################################

normalization_tiles_for_promoter2<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/4500,ncol=4500,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                        M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]*1000000/(0.001*l)
                                }
                        }
                return(M_0h.p.matrix2)
        }

#######################################

normalization_tiles_for_promoter_precise<-function(x)
        {
                M_0h<-import.bed(x)
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/4501,ncol=4501,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                        M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]*1000000/(0.001*l)
                                }
                        }
                return(M_0h.p.matrix2)
        }




#################################################

#######################################

normalization_tiles_for_short_promoter<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/21,ncol=21,byrow=TRUE)
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



#######################################

normalization_tiles_for_summits<-function(x)
	{
                M_0h<-import.bed(x)
		start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/301,ncol=301,byrow=TRUE)
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

##############################################

normalization_tiles_for_summits_percise_no_extend<-function(x)
        {
                M_0h<-import.bed(x)
               # start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
               # width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/3001,ncol=3001,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                       M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]*1000000/(l)
                                }
                        }
                return(M_0h.p.matrix2)
        }












#####################################

normalization_tiles_for_summits_distance<-function(x)
        {
                M_0h<-import.bed(x)
         #       start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
          #      width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h,type="any")
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/1001,ncol=1001,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                       M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]/length(M_0h)
                                }
                        }
                return(M_0h.p.matrix2)
        }

########################################

normalization_tiles_for_summits_distance2<-function(x)
        {
                M_0h<-import.bed(x)
         #       start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
          #      width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h,type="any")
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/101,ncol=101,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                       M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]
                                }
                        }
                return(M_0h.p.matrix2)
        }



########################################

normalization_tiles_for_gene_body_normal<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/321,ncol=321,byrow=TRUE)
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
########################################

normalization_tiles_for_gene_body_normal_precise<-function(x)
        {
                M_0h<-import.bed(x)
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/3210,ncol=3210,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                        M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]*1000000/l
                                }
                        }
                return(M_0h.p.matrix2)
        }



########################################
normalization_tiles_for_gene_body_reads<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/321,ncol=321,byrow=TRUE)
                l<-length(M_0h)
                M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                        M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]
                                }
                        }
                return(M_0h.p.matrix2)
        }


########################################

normalization_tiles_for_gene_body_long<-function(x)
        {
                M_0h<-import.bed(x)
                start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
                width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
                M_0h.p=countOverlaps(tiles_2,M_0h)
                M_0h.p.matrix=matrix(M_0h.p,nrow=length(tiles_2)/561,ncol=561,byrow=TRUE)
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


########################################

read_files_to_gene_names<-function(x)
{
	gene_names_signal<-read.table(x)
	gene_names_signal2<-as.character(gene_names_signal$V1)
	return(gene_names_signal2)
}

grep_out_files<-function(x)
{
	output_files<-dir()[grepl(x,dir())]
	return(output_files)
}

##########################################
get_matrix_back<-function(x){
n=1
matrix_for_results<-matrix(0,nrow=200,ncol=2*length(grep_out_files(x)))
for ( i in 1:length(grep_out_files(x)))
{
        for_intersect<-read_files_to_gene_names(grep_out_files(x)[i])
        for( j in c("all_the_genes_without_motif","all_the_genes_with_motifs"))
        {
        output_names_a<-paste(j,grep_out_files(x)[i],sep="_")
        intersect_later<-intersect(for_intersect,get(j))
        assign(output_names_a,intersect_later)
        the_gene_list_used<-intersect_later
        binding_signal_summits2<-binding_numbers4[names(binding_numbers4)%in%the_gene_list_used]
        matrix_for_results[,n]<-binding_signal_summits2[1:200]
        n=n+1
        }
}
return(matrix_for_results)
}
###############################################
get_plots_for_bar<-function(x){
matrix_for_results<-x
cutoff=114
barplot(c(length(matrix_for_results[,8][matrix_for_results[,8]>cutoff]),length(matrix_for_results[,6][matrix_for_results[,6]>cutoff]),length(matrix_for_results[,4][matrix_for_results[,4]>cutoff]),length(matrix_for_results[,2][matrix_for_results[,2]>cutoff]),length(matrix_for_results[,7][matrix_for_results[,7]>cutoff]),length(matrix_for_results[,5][matrix_for_results[,5]>cutoff]),length(matrix_for_results[,3][matrix_for_results[,3]>cutoff]),length(matrix_for_results[,1][matrix_for_results[,1]>cutoff])),col="blue",border=NA,width=0.1,space=0.5,cex.axis=2,axes=FALSE,xaxt="n",ylim=c(0,100))
axis(2, col.axis="black", las=0,lwd=3,cex.axis=2,at=c(0,40,80),labels=c("0","40","80"))
}
#################################################
curve_correlation<-function(x,y){
correlation_values<-NULL
for ( i in 1:length(x[,1]))
{
        correlation_values[i]<-cor(x[i,],y[i,])
#	return(correlation_values)
}
hist(correlation_values,xlim=c(-1,1))
#return(correlation_values)
}

###################################################
calculate_correlation<-function(x,y)
{
correlation_values<-cor(log2(x)[!is.infinite(log2(x))][!is.infinite(log2(y)[!is.infinite(log2(x))])],log2(y)[!is.infinite(log2(x))][!is.infinite(log2(y)[!is.infinite(log2(x))])])
return(correlation_values)
}

quantile_seperate<-function(x)
{
	binding_numbers2<-x
	binding_numbers2[x<=quantile(x,0.25)]=0
	binding_numbers2[x>=quantile(x,0.25)&x<=quantile(x,0.5)]=1
	binding_numbers2[x>=quantile(x,0.5)&x<=quantile(x,0.75)]=2
	binding_numbers2[x>=quantile(x,0.75)&x<=quantile(x,1)]=3
	return(binding_numbers2)
}

quantile_seperate_set<-function(x,a,b,c,d)
{
        binding_numbers2<-x
        binding_numbers2[x<=quantile(x,a)]=0
        binding_numbers2[x>=quantile(x,a)&x<=quantile(x,b)]=10
        binding_numbers2[x>=quantile(x,b)&x<=quantile(x,c)]=20
        binding_numbers2[x>=quantile(x,c)&x<=quantile(x,d)]=30
        return(binding_numbers2)
}

#######################################################

MA_return_the_different<-function(x,y){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
weight_value=10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=10
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-10
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3


names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)])

list_for_return<-list(names_increase,names_decrease)

return(list_for_return)

}

#########################################################

MA_return_the_different75<-function(x,y){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
weight_value=10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=10
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-10
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3


names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values+0.75)][change_mean_values[change_y_axis>(weight_value/change_mean_values+0.75)]>10])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)]>c(10)])

list_for_return<-list(names_increase,names_decrease)

return(list_for_return)

}
########################################################

MA_return_the_different750<-function(x,y){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
weight_value=10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=10
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-10
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3


names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values+0.5)][change_mean_values[change_y_axis>(weight_value/change_mean_values+0.5)]>10])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.5)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.5)]>c(10)])

list_for_return<-list(names_increase,names_decrease)

return(list_for_return)

}


########################################################

MA_return_the_different2<-function(x,y){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
weight_value=10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=10
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-10
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3


names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values)]>10])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1))][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1))]>c(10)])

list_for_return<-list(names_increase,names_decrease)

return(list_for_return)

}


########################################################

MA_plot_figure<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(300,-4,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(300,4,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}
###################################################################
MA_plot_figure_represent<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=1,cex=1,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=1,
        cex=1
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=1,
        cex=1
        )
#try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
#try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
#try_value2<-paste("n=",try_values,sep="")
#text(300,-4,labels=try_value2,cex=2)
#try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
#try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
#try_value2_b<-paste("n=",try_values_b,sep="")
#text(300,4,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=2,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=2,cex.axis=2,at=c(10,z/10,z))
}


MA_plot_figure_low<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(2,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(300,-4,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(300,4,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}

##################################

MA_plot_figure_single<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(x+y)/2
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=5
weight_value2=0.75
weight_value3=3

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(0.1,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)][change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)]>weight_value3],change_y_axis[change_y_axis>(weight_value2+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)]>weight_value3],        
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)]>weight_value3],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)]>weight_value3],
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(75,-4,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(75,4,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(0.1,z/10,z))
}
###################################################################

MA_return_the_different_single<-function(x,y){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/2
#weight_value=10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=10
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-10
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=5
weight_value2=0.75
weight_value3=3


names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)][change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)]>weight_value3])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)]>weight_value3])

list_for_return<-list(names_increase,names_decrease)

return(list_for_return)

}

###################################################################

MA_return_the_different_15<-function(x,y){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/2
#weight_value=10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=10
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-10
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=5
weight_value2=1
weight_value3=20


names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)][change_mean_values[change_y_axis>(weight_value/change_mean_values+weight_value2)]>weight_value3])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-weight_value2)]>weight_value3])

list_for_return<-list(names_increase,names_decrease)

return(list_for_return)

}



####################################################################
MA_plot_figure_histone<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(1,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>0][change_mean_values[change_y_axis>0]>1],change_y_axis[change_y_axis>0][change_mean_values[change_y_axis>0]>1],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<0][change_mean_values[change_y_axis<0]>1],change_y_axis[change_y_axis<0][change_mean_values[change_y_axis<0]>1],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<0][change_mean_values[change_y_axis<0]>1]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(50,-1.5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>0][change_mean_values[change_y_axis>0]>1]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(50,1.5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(1,z/10,z))
}
############################################################

MA_return_different_number<-function(x,y,z,k){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

return_values<-NULL
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
return_values[1]<-try_values_b
return_values[2]<-(length(y)-try_values_b-try_values)
return_values[3]<-try_values

return(return_values)
}



##################################################################

MA_plot_figure_ChIP<-function(x,y,z,k){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(1,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>0],change_y_axis[change_y_axis>0],
      #  ylim=c(-2,2),
      #  xlim=c(0,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<0],change_y_axis[change_y_axis<0],
       # ylim=c(-2,2),
       # xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<0]
try_values<-length(try_values_jj)
try_value2<-paste("n=",try_values,sep="")
text(100,-15,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>0]
try_values_b<-length(try_values_jj)
try_value2_b<-paste("n=",try_values_b,sep="")
text(100,15,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(1,10,z/10,z))
}
###################################################################

MA_plot_figure_ChIP2<-function(x,y,z,k){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(1,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>2],change_y_axis[change_y_axis>2],
      #  ylim=c(-2,2),
      #  xlim=c(0,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<c(-2)],change_y_axis[change_y_axis<c(-2)],
       # ylim=c(-2,2),
       # xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<c(-2)]
try_values<-length(try_values_jj)
try_value2<-paste("n=",try_values,sep="")
text(100,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>2]
try_values_b<-length(try_values_jj)
try_value2_b<-paste("n=",try_values_b,sep="")
text(100,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(1,10,z/10,z))
}
###################################################################
MA_plot_figure2<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(100,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(100,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}
###################################################################

MA_plot_figure15<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+1.5)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1.5)]>10],change_y_axis[change_y_axis>(1.5+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1.5)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1.5)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1.5)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1.5)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1.5)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1.5)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1.5)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(100,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1.5+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1.5)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(100,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}

##############################

MA_plot_figure75<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+0.75)][change_mean_values[change_y_axis>(weight_value/change_mean_values+0.75)]>10],change_y_axis[change_y_axis>(0.75+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+0.75)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-0.75)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(100,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(0.75+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+0.75)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(100,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}


#############################

MA_plot_figure_time<-function(x,y,z,k,gene_list_up,gene_list_down,total_genes){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
total_genes2<-total_genes[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k



change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,
        cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[total_genes2%in%gene_list_up],change_y_axis[total_genes2%in%gene_list_up],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=2,
        col='red',
        lwd=1,
        cex=1
        )
points(change_mean_values[total_genes2%in%gene_list_down],change_y_axis[total_genes2%in%gene_list_down],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=6,
        col='blue',
        lwd=1,
        cex=1
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(1000,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(1000,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}
#######################################################################################

MA_plot_figure_time_up<-function(x,y,z,k,gene_list_up,gene_list_down,total_genes){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
total_genes2<-total_genes[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,
        cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[total_genes2%in%gene_list_up],change_y_axis[total_genes2%in%gene_list_up],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=2,
        col='red',
        lwd=1,
        cex=1
        )
#points(change_mean_values[total_genes2%in%gene_list_down],change_y_axis[total_genes2%in%gene_list_down],
#        ylim=c(-2,2),
#        xlim=c(3,6),
#        pch=6,
#        col='blue',
#        lwd=1,
#        cex=1
#        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(1000,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(1000,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}
######################################################################################

MA_plot_figure_time_down<-function(x,y,z,k,gene_list_up,gene_list_down,total_genes){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
total_genes2<-total_genes[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,
        cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
#points(change_mean_values[total_genes2%in%gene_list_up],change_y_axis[total_genes2%in%gene_list_up],
#        ylim=c(-2,2),
#        xlim=c(3,6),
#        pch=2,
#        col='red',
#        lwd=1,
#        cex=1
#        )
points(change_mean_values[total_genes2%in%gene_list_down],change_y_axis[total_genes2%in%gene_list_down],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=6,
        col='blue',
        lwd=1,
        cex=1
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(1000,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(1000,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}



######################################################################################

MA_plot_figure3<-function(x,y,z,k){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3
weight_value=10

plot(change_mean_values,change_y_axis,ylim=c(-k,k),xlim=c(10,z),pch=16,col='grey',lwd=2,cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values)]>10],change_y_axis[change_y_axis>(weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1))][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1))]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1))][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1))]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1))][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1))]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(1000,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(1000,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}
############################################################################





MA_plot_figure2<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
#change_y_axis=change_y_axis2+0.001
change_y_axis3<-change_y_axis2[!is.na(change_y_axis2)]
change_mean_values3<-change_mean_values2[!is.na(change_y_axis2)]
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]>0]=k
change_y_axis3[is.infinite(change_y_axis3)][change_y_axis3[is.infinite(change_y_axis3)]<0]=-k
change_mean_values<-change_mean_values3
change_y_axis<-change_y_axis3

weight_value=10
plot(change_mean_values,change_y_axis,
       ylim=c(-k,k),
       xlim=c(10,z),
        pch=16,
        col='grey',
        lwd=2,
        cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+2)][change_mean_values[change_y_axis>(weight_value/change_mean_values+2)]>10],change_y_axis[change_y_axis>(2+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+2)]>10],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-2)]>c(10)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-2)]>c(10)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )
try_values_jj<-change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-2)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-2)]>c(10)]
try_values<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2<-paste("n=",try_values,sep="")
text(1000,-5,labels=try_value2,cex=2)
try_values_jj<-change_y_axis[change_y_axis>(2+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+2)]>10]
try_values_b<-length(try_values_jj[!is.na(try_values_jj)][!is.infinite(try_values_jj[!is.na(try_values_jj)])])
try_value2_b<-paste("n=",try_values_b,sep="")
text(1000,5,labels=try_value2_b,cex=2)
axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}


###################################################################
make_scales<-function(x,a,c,d,e){
image(seq(0, max(x), length.out=e), 1,
      matrix(seq(0, max(x), length.out=e),e,1),
      col = colors,
      xlab=c, ylab='',
      main=d, yaxt='n',
      lwd=3, axes=TRUE,
        zlim=c(0,a))
box(col='black', lwd=2)
}

make_scales2<-function(a,b,c,d,e){
image(seq(b, a, length.out=e), 1,
      matrix(seq(b, a, length.out=e),e,1),
      col = colors,
      xlab=c, ylab='',
      main=d, yaxt='n',
      lwd=3, axes=TRUE,
        zlim=c(b,a))
box(col='black', lwd=2)
}


##################################################################

comvert_name<-function(x){

#        name_reference<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/gene_names_convert.tab")
#        gene_name_reference_A<-as.character(name_reference$V2)
#        gene_name_reference_B<-as.character(name_reference$V1)
#        gene_name_used<-x
#        gene_formal<-NULL
#        for ( i in 1:length(gene_name_used))
#        {
#                gene_formal[i]<-gene_name_reference_B[gene_name_used[i]==gene_name_reference_A]
#
#        }
#        return(gene_formal)

        convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/gene_names_convert.tab")
        convertor_genes_first<-as.character(convertor_genes$V2)
        convertor_genes_second<-as.character(convertor_genes$V1)
        output<-NULL
        d=1
        for ( i in 1:length(x))
        {
                if(x[i]%in%convertor_genes_first)
                {
                        output[d]<-convertor_genes_second[x[i]==convertor_genes_first]
                        d=d+1
                }
        }
        return(output)


}

##################################################################

MA_plot_figure_above<-function(x,y,z,k){
change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
change_mean_values<-change_mean_values2[change_mean_values2>200]
change_y_axis<-change_y_axis2[change_mean_values2>200]
weight_value=10
plot(change_mean_values,change_y_axis,
       ylim=c(-k,k),
       xlim=c(10,z),
        pch=16,
        col='grey',
        lwd=2,
        cex=2,
       log="x",
        xaxt="n",
        axes=FALSE,
        xlab=NA,
        ylab=NA
        )
abline(h =0, col = "black", lty = 3, lwd=5)
points(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>1],change_y_axis[change_y_axis>(1+weight_value/change_mean_values)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>1],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='red',
        lwd=2,
        cex=2
        )
points(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(1)],change_y_axis[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(1)],
        ylim=c(-2,2),
        xlim=c(3,6),
        pch=16,
        col='blue',
        lwd=2,
        cex=2
        )

axis(2, col.axis="black", las=2,lwd=6,cex.axis=2)
axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,z/10,z))
}

###############################################################

cut_data_equal<-function(x,a)
{
	split_part<-seq(0, max(x), len = a)
#	
#	for ( i in 2:length(split_part))
#	x<-split_part[i-1]
#
#
#
#
	return(split_part)
}

#################################################################

MA_return_the_different_inter<-function(x,y,a,b){

change_y_axis2=log2(y/x)
change_mean_values2=(y+x)/10
change_mean_values<-change_mean_values2[change_mean_values2>a&change_mean_values2<b]
change_y_axis=change_y_axis2[change_mean_values2>a&change_mean_values2<b]

weight_value=10
names_increase<-names(change_mean_values[change_y_axis>(weight_value/change_mean_values+1)][change_mean_values[change_y_axis>(weight_value/change_mean_values+1)]>1])
names_decrease<-names(change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)][change_mean_values[change_y_axis<(weight_value/change_mean_values*c(-1)-1)]>c(1)])
middle_part1<-names(change_mean_values)[!names(change_mean_values)%in%names_increase]
middle_part2<-middle_part1[!middle_part1%in%names_decrease]
list_for_return<-list(names_increase,middle_part2,names_decrease)

return(list_for_return)

}

####################################################################

select_gene_according_to_list<-function(x,y)
{

        asd_in_gene_name<-y
	gene_location_name<-as.character(rownames(x))
	
        gff0_genome_location2<-as.matrix(x)
       # chr_region<-matrix(0,ncol=ncol(x),nrow=length(asd_in_gene_name))
                #for ( i in 1:length(asd_in_gene_name) )
                #{
			chr_region<-gff0_genome_location2[gene_location_name%in%asd_in_gene_name,]
                #}
	rownames(chr_region)<-gene_location_name[gene_location_name%in%asd_in_gene_name]
	return(chr_region)
}

####################################################################

read_files_into_character<-function(x)
{
	files_input2<-read.table(x)
	files_input<-as.character(files_input2$V1)
	return(files_input)
}

####################################################################

quantile_seperate_advance<-function(x,y)
{
	binding_numbers2<-x
	for ( i in 2:length(seq(0,1,length.out=y)))
        {
	binding_numbers2[x>=quantile(x,seq(0,1,length.out=y)[i-1])&x<=quantile(x,seq(0,1,length.out=y)[i])]=i-1
        }
	return(binding_numbers2)
}

make_scales<-function(x,a,c,d){
image(seq(0, max(x), length.out=100), 1,
      matrix(seq(0, max(x), length.out=100),100,1),
      col = colors,
      xlab=c, ylab='',
      main=d, yaxt='n',
      lwd=3, axes=TRUE,
        zlim=c(0,a))
}

###########################################################

#correlation_plot<-function(values_set_one,value_set_two)
#{
#	plot(values_set

plot_heatmap<-function(x,k,l){
        M_0h.p.matrix2<-x[,107:267]
        M_0h.p.matrix2[M_0h.p.matrix2>k]<-k
        image(  x=seq(1, l, length.out=l),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2[order(rowSums(M_0h.p.matrix2)),]),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,k),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
#        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
#        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
#        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
#        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
#        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)

        #axis(1, at=c(1:60), col.axis="black", labels=FALSE, lwd=3, las=0, tck=0)
        #axis(1, at=c(1:1), col.axis="black",labels=c("TSS-500bp"),las=0,tck=0,lwd=0,cex.axis=1)
        #axis(1, at=c(60:60), col.axis="black",labels=c("TSS"),las=0,tck=0,lwd=0,cex.axis=1)
        #axis(1, at=c(261:261), col.axis="black",labels=c("TTS"),las=0,tck=0,lwd=0,cex.axis=1)
        #axis(1, at=c(321:321), col.axis="black",labels=c("TTS+500bp"),las=0,tck=0,lwd=0,cex.axis=1)
        #axis(1, at=c(60:261), col.axis="black", labels=FALSE, lwd=8, las=0, tck=0)
        #axis(1, at=c(261:321), col.axis="black", labels=FALSE, lwd=3, las=0, tck=0)
        #axis(2, col.axis="black", las=2,lwd=8,cex.axis=2)
}

search_gene_lists_values<-function(x,y){
        gene_expression_values_inter<-y[rownames(y)%in%x,]
        rownames(gene_expression_values_inter)<-rownames(y)[rownames(y)%in%x]
        return(gene_expression_values_inter)
}

plot_line<-function(x,y){
	               x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

	plot(colMeans(y)/colMeans(x),t="l",lwd=3,col="red",ylim=c(0.5,1.5),xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
	axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
	abline(h=1,col="blue",lwd=2)
}

plot_line_C<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMeans(y[gene_names2%in%FacB_target_genenames,])/colMeans(x[gene_names2%in%FacB_target_genenames,]),t="l",lwd=3,col="red",ylim=c(0.5,1.5),xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        abline(h=1,col="blue",lwd=2)
}


plot_line2<-function(x,y){
	               x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMedians(y)/colMedians(x),t="l",lwd=3,col="red",ylim=c(0.5,1.5),xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        abline(h=1,col="blue",lwd=2)
}

plot_line_2<-function(x,y){
	x[is.na(x)]<-0.01
	y[is.na(y)]<-0.01
	x[is.infinite(x)]<-10
	y[is.infinite(y)]<-10       
	plot(colMeans(y)/colMeans(x),t="l",lwd=3,col="red",ylim=c(0,10),xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        abline(h=1,col="blue",lwd=2)
}

plot_line2_2<-function(x,y){
	       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10
 
       plot(colMedians(y)/colMedians(x),t="l",lwd=3,col="red",ylim=c(0,10),xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        abline(h=1,col="blue",lwd=2)
}

plot_line_3<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMeans(x),t="l",lwd=3,col=y,xlab=NA,ylim=c(50,200),
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        abline(h=150,col="blue",lwd=2)
}

plot_line_32<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMedians(x),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        abline(h=1,col="blue",lwd=2)
}

plot_line_first<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMeans(x),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=2)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=2)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=2)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}

plot_line_first_compare<-function(x,y,cols,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMeans(y)/colMeans(x),t="l",lwd=3,col=cols,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}



plot_line_first_gene_body<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10
	
        plot(colMeans(x[order(rowMeans(x),decreasing=TRUE),][1:1000,]),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:60), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(61:261), col.axis="black", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(261:321), col.axis="black",labels=FALSE,las=0,tck=0,lwd=4)
	axis(1, at=c(261:261), col.axis="black",labels=c("1"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("0.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(10:10), col.axis="black",labels=c("-0.5k"),las=0,lwd=0,cex.axis=1.5)
	axis(1, at=c(61:61), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}

plot_line_first_gene_body3<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10
        raw_values<-colMeans(x[order(rowMeans(x),decreasing=TRUE),][1:1000,])
        plot_values<-(raw_values-min(raw_values))/(max(raw_values)-min(raw_values))
	k=0
	m=1.4
        plot(plot_values,t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:60), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(61:261), col.axis="black", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(261:321), col.axis="black",labels=FALSE,las=0,tck=0,lwd=4)
        axis(1, at=c(261:261), col.axis="black",labels=c("1"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("0.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(10:10), col.axis="black",labels=c("-0.5k"),las=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(61:61), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}

plot_line_first_gene_body_long<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10
        raw_values<-colMeans(x[order(rowMeans(x),decreasing=TRUE),][1:100,])
        plot_values<-(raw_values-min(raw_values))/(max(raw_values)-min(raw_values))
        k=0
        m=1.4
        plot(plot_values,t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:300), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(301:501), col.axis="black", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(501:561), col.axis="black",labels=FALSE,las=0,tck=0,lwd=4)
        axis(1, at=c(501:501), col.axis="black",labels=c("1"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(551:551), col.axis="black",labels=c("0.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(301:301), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}

plot_line_first_gene_body_long2<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10
        raw_values<-colMeans(x[order(rowMeans(x),decreasing=TRUE),][1:100,])
        plot_values<-(raw_values-min(raw_values))/(max(raw_values)-min(raw_values))
        #k=0
        #m=1.4
        plot(raw_values,t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:300), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(301:501), col.axis="black", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(501:561), col.axis="black",labels=FALSE,las=0,tck=0,lwd=4)
        axis(1, at=c(501:501), col.axis="black",labels=c("1"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(551:551), col.axis="black",labels=c("0.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(301:301), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}


plot_line_first_gene_body2<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMedians(x[order(rowMeans(x),decreasing=TRUE),][1:1000,]),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:60), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(61:261), col.axis="black", labels=FALSE, lwd=8, las=0, tck=0)
        axis(1, at=c(261:321), col.axis="black",labels=FALSE,las=0,tck=0,lwd=4)
        axis(1, at=c(261:261), col.axis="black",labels=c("1"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("0.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(10:10), col.axis="black",labels=c("-0.5k"),las=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(61:61), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
       # abline(h=1,col="blue",lwd=2)
}
plot_line_add_gene_body3<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10
        raw_values<-colMeans(x[order(rowMeans(x),decreasing=TRUE),][1:100,])
        plot_values<-(raw_values-min(raw_values))/(max(raw_values)-min(raw_values))
        k=0
        m=1


        lines(plot_values,t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n")
   #     axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
   #     axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
   #     axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
   #     axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
   #     abline(h=1,col="blue",lwd=2)
}
plot_line_add_gene_body2<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        lines(colMedians(x[order(rowMeans(x),decreasing=TRUE),][1:1000,]),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n")
   #     axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
   #     axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
   #     axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
   #     axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
   #     abline(h=1,col="blue",lwd=2)
}


plot_line_add_gene_body<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        lines(colMeans(x[order(rowMeans(x),decreasing=TRUE),][1:100,]),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n")
   #     axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
   #     axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
   #     axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
   #     axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
   #     abline(h=1,col="blue",lwd=2)
}



plot_line_first_summits<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMeans(x),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:151), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(151:301), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-1.5k"),las=0,tck=0,lwd=0,cex.axis=2)
        axis(1, at=c(291:291), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=2)
        axis(1, at=c(151:151), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=2)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
#        abline(h=1,col="blue",lwd=2)
}

plot_line_first_summits_median<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(colMedians(x),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:151), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(151:301), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-1.5k"),las=0,tck=0,lwd=0,cex.axis=2)
        axis(1, at=c(291:291), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=2)
        axis(1, at=c(151:151), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=2)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
#        abline(h=1,col="blue",lwd=2)
}


plot_line_first_summits_distance<-function(x,y,k,m){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        plot(x,t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, ylim=c(k,m),xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:length(x)), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
#        abline(h=1,col="blue",lwd=2)
}


plot_line_add_summits<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        lines(colMeans(x),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
   #     axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
   #     axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
   #     axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
   #     axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
   #     abline(h=1,col="blue",lwd=2)
}

plot_line_add_summits_distance<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        lines(x,t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
   #     axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
   #     axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
   #     axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
   #     axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
   #     abline(h=1,col="blue",lwd=2)
}


plot_line_add<-function(x,y){
                       x[is.na(x)]<-0.01
        y[is.na(y)]<-0.01
        x[is.infinite(x)]<-10
        y[is.infinite(y)]<-10

        lines(colMeans(x),t="l",lwd=3,col=y,xlab=NA,
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
   #     axis(1, at=c(1:214), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
   #     axis(1, at=c(214:321), col.axis="blue4", labels=FALSE, lwd=8, las=0, tck=0)
   #     axis(1, at=c(10:10), col.axis="black",labels=c("-3k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
   #     axis(1, at=c(214:214), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
   #     axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
   #     abline(h=1,col="blue",lwd=2)
}


plot_heatmap_B<-function(x,y,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
	M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        y_input<-y
        y_input[y_input>z]<-z
	y_input[is.na(y_input)]<-0
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2[order(rowMeans(y_input),decreasing=FALSE),]),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}

plot_heatmap_B2<-function(x,y,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        y_input<-y
        #y_input[y_input>z]<-z
        #y_input[is.na(y_input)]<-0
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2[order(rowMeans(y_input),decreasing=FALSE),]),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}



plot_heatmap_C<-function(x,y,z,k){
        M_0h.p.matrix2<-x[gene_names2%in%FacB_target_genenames,]
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        y_input<-y[gene_names2%in%FacB_target_genenames,]
        y_input[y_input>z]<-z
        y_input[is.na(y_input)]<-0
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2[order(rowMeans(y_input),decreasing=FALSE),]),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}

#install.packages("DescTools")
library(DescTools)
plot_heatmap_E<-function(x,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(Rev(M_0h.p.matrix2,1)),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}

plot_heatmap_E2<-function(x,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(Rev(M_0h.p.matrix2,1)),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}


plot_heatmap_E3<-function(x,y,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(Rev(M_0h.p.matrix2,1)),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(y,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}



plot_correlation_plot<-function(x,y,start_values,end_values,colors_values){
	plot(
		x[order(x,decreasing=TRUE)][1:500],
		y[order(x,decreasing=TRUE)][1:500],
		col=colors_values,
		cex=0.5,
		pch=15,
		type="p",
		xlim=c(start_values,end_values),
		ylim=c(start_values,end_values),	
#		log="xy",
		xaxt="n",
		axes=FALSE,
		xlab=NA,
		ylab=NA
		)
	abline(a=0,b=1,col=4,lwd=3)
	axis(2, col.axis="black", las=2,lwd=6,cex.axis=2,at=c(10,100,1000))
	axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(10,100,1000))
		}

plot_correlation_plot2<-function(x,y,start_values,end_values,colors_values){
        plot(
                x[order(x+y,decreasing=TRUE)][1:500],
                y[order(x+y,decreasing=TRUE)][1:500],
                col=colors_values,
                cex=0.5,
                pch=15,
                type="p",
                xlim=c(start_values,end_values),
                ylim=c(start_values,end_values),
               log="xy",
                xaxt="n",
                axes=FALSE,
                xlab=NA,
                ylab=NA
                )
        abline(a=0,b=1,col=4,lwd=3)
        axis(2, col.axis="black", las=2,lwd=6,cex.axis=2,at=c(50,500,5000))
        axis(1, col.axis="black", las=0,lwd=6,cex.axis=2,at=c(50,500,5000))
                }




######################################################

plot_line_summits<-function(x,y){
        #               x[is.na(x)]<-0.01
        #y[is.na(y)]<-0.01
        #x[is.infinite(x)]<-10
        #y[is.infinite(y)]<-10

        plot(colMeans(x),t="l",lwd=3,col=y,xlab=NA,ylim=c(50,200),
                ylab=NA, xaxt="n",yaxt="n",axes=FALSE)
        axis(1, at=c(1:151), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(151:301), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
        axis(1, at=c(10:10), col.axis="black",labels=c("-1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
        axis(1, at=c(151:151), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
        axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        #abline(h=150,col="blue",lwd=2)
}

plot_line_summits_add<-function(x,y){
       #                x[is.na(x)]<-0.01
       # y[is.na(y)]<-0.01
       # x[is.infinite(x)]<-10
       # y[is.infinite(y)]<-10

        lines(colMeans(x),t="l",lwd=3,col=y,xlab=NA,ylim=c(50,200),
                ylab=NA, xaxt="n",yaxt="n")
 #       axis(1, at=c(1:151), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
 #       axis(1, at=c(151:301), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
 #       axis(1, at=c(10:10), col.axis="black",labels=c("-1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
 #       axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
 #       axis(1, at=c(151:151), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
 #       axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        #abline(h=150,col="blue",lwd=2)
}

plot_line_summits_add_median<-function(x,y){
       #                x[is.na(x)]<-0.01
       # y[is.na(y)]<-0.01
       # x[is.infinite(x)]<-10
       # y[is.infinite(y)]<-10

        lines(colMedians(x),t="l",lwd=3,col=y,xlab=NA,ylim=c(50,200),
                ylab=NA, xaxt="n",yaxt="n")
 #       axis(1, at=c(1:151), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
 #       axis(1, at=c(151:301), col.axis="black", labels=FALSE, lwd=4, las=0, tck=0)
 #       axis(1, at=c(10:10), col.axis="black",labels=c("-1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
 #       axis(1, at=c(311:311), col.axis="black",labels=c("1.5k"),las=0,tck=0,lwd=0,cex.axis=1.5)
 #       axis(1, at=c(151:151), col.axis="black",labels=c("0"),las=0,lwd=0,cex.axis=1.5)
 #       axis(2, col.axis="black", las=2,lwd=4,cex.axis=2)
        #abline(h=150,col="blue",lwd=2)
}


###################################################

get_bed_regionsi_below<-function(x,y,value_input){
value_input<-0
change_y_axis<-y
bed_files<-x
for ( i in 1:length(bed_files[,1]))
{
        chromsnames<-as.character(bed_files[,1][change_y_axis<value_input])
        start_point<-bed_files[,2][change_y_axis<value_input]
        end_point<-bed_files[,3][change_y_axis<value_input]
        names_file<-as.character(bed_files[,4][change_y_axis<value_input])
        bedvalues<-bed_files[,5][change_y_axis<value_input]
}
output_bed_files<-data.frame(chromsnames,start_point,end_point,names_file,bedvalues)
write.table(output_bed_files,file="outboundary_bed_files.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

###################################################


plot_heatmap_cluster<-function(x,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}

plot_heatmap_cluster2<-function(x,z,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2),
                col=my_palette,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}

plot_heatmap_cluster3<-function(x,z,a,k){
        M_0h.p.matrix2<-x
        M_0h.p.matrix2[M_0h.p.matrix2>z]<-z
	M_0h.p.matrix2[M_0h.p.matrix2<a]<-a
        image(  x=seq(0, k, length.out=k),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(a,z),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}


get_read_according<-function(gene_files,strands_genes,input_beds){
import_seq<-import.bed(input_beds)
gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(gene_files)
        asd_in_gene_name<-as.character(x1$V1)
        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA2<-import.bed("inter_file.bed")
        gene_select_region_CreA<-gene_select_region_CreA2[gene_select_region_CreA2$score==strands_genes]
	
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
          {      A_value[i]+seq(start(ranges(gene_select_region_CreA))[i],end(ranges(gene_select_region_CreA))[i],length.out=200)})
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=200), 1:200, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=200)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=20),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        M_0h<-import_seq
        start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
        width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
        M_0h.p=countOverlaps(tiles_2,M_0h)
        M_0h.p.matrix=matrix(M_0h.p,nrow=length(gene_select_region_CreA),byrow=TRUE)
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
######################################################################

get_read_according2<-function(gene_files,strands_genes,input_beds,orignal_bed_files){
import_seq<-import.bed(input_beds)
gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(gene_files)
        asd_in_gene_name<-as.character(x1$V1)
        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA2<-import.bed("inter_file.bed")
        gene_select_region_CreA<-gene_select_region_CreA2[gene_select_region_CreA2$score==strands_genes]
	gene_width_values<-width(ranges(gene_select_region_CreA))
        A_value<-NULL
        for( i in 1:length(end(ranges(gene_select_region_CreA))))
        {
                A_value[i]=0
        }
        tiles2=sapply(1:length(end(ranges(gene_select_region_CreA))),function(i)
          {      A_value[i]+seq(((start(ranges(gene_select_region_CreA))[i]+end(ranges(gene_select_region_CreA))[i])/2)-3000,((start(ranges(gene_select_region_CreA))[i]+end(ranges(gene_select_region_CreA))[i])/2)+3000,length.out=300)})
        Chromatin_length<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/genome.size")
        si<-Seqinfo(seqnames=as.character(Chromatin_length$V1),seqlengths=Chromatin_length$V2,genome="A_nidulans")
        tiles_2 = GRanges(tilename = paste( rep(gene_select_region_CreA$name, each=300), 1:300, sep="_" ),
                        seqnames = Rle( rep(as.character(seqnames(gene_select_region_CreA)),each=300)),
                        ranges=IRanges(start=as.vector(tiles2),
                        width=20),
                        strand=Rle(rep("*",length(as.vector(tiles2)))),
                        seqinfo=si)
        M_0h<-import_seq
        start(ranges(M_0h))[as.character(strand(M_0h))=="-"]=start(ranges(M_0h))[as.character(strand(M_0h))=="-"]-(200-width(ranges(M_0h))[as.character(strand(M_0h))=="-"])
        width(ranges(M_0h))[as.character(strand(M_0h))=="+"]=200
        M_0h.p=countOverlaps(tiles_2,M_0h)
        M_0h.p.matrix=matrix(M_0h.p,nrow=length(gene_select_region_CreA),byrow=TRUE)
        l<-length(import.bed(orignal_bed_files))
        M_0h.p.matrix2<-M_0h.p.matrix
                for ( i in 1:length(M_0h.p.matrix[,1]) )
                        {
                                for ( j in 1:length(M_0h.p.matrix[i,]))
                                {
                                        M_0h.p.matrix2[i,j]<-M_0h.p.matrix[i,j]*1000000/l
                                }
                        }
	mask_signal_values<-M_0h.p.matrix2
	for( i in 1:length(gene_width_values))
	{
	if(ceiling(gene_width_values[i]/(2*20))<150)
	{
	mask_signal_values[i,1:(150-ceiling(gene_width_values[i]/(2*20)))]<-0
	mask_signal_values[i,(150+ceiling(gene_width_values[i]/(2*20))):300]<-0
	}
	}
        return(mask_signal_values)

}

plot_gene_length_heatmap<-function(gene_files,strands_genes,matrix_input_values,max_values,Bin_size)
{
	gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(gene_files)
        asd_in_gene_name<-as.character(x1$V1)
        gene_location_name<-as.character(gff0_genome_location$V4)
        gff0_genome_location2<-as.matrix(gff0_genome_location)
        chr_region<-matrix(0,ncol=5,nrow=length(asd_in_gene_name))
                for ( i in 1:length(asd_in_gene_name) )
                {
                        chr_region[i,]<-gff0_genome_location2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4),]
                }
        write.table(chr_region,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA2<-import.bed("inter_file.bed")
        gene_select_region_CreA<-gene_select_region_CreA2[gene_select_region_CreA2$score==strands_genes]
        gene_width_values<-width(ranges(gene_select_region_CreA))
        M_0h.p.matrix2<-matrix_input_values
        M_0h.p.matrix2[M_0h.p.matrix2>max_values]<-max_values
        M_0h.p.matrix2[is.na(M_0h.p.matrix2)]<-0
        image(  x=seq(0, Bin_size, length.out=Bin_size),
                y=1:length(M_0h.p.matrix2[,1]),
                z=t(M_0h.p.matrix2[order(gene_width_values,decreasing=TRUE),]),
                col=colors,
                xlab=NA,
                ylab=NA,
                lwd=3,
                zlim=c(0,max_values),
                xaxt="n",
                yaxt="n",
                axes=FALSE)
}
		

######################################################################


str_seperate<-function(input_files,delite)
{
output_files<-input_files
for ( i in 1:length(input_files))
{
	output_files[i]<-strsplit(input_files[i],"-")[[1]]
}
return(output_files)
}


######################################################################

pvalues_get<-function(x,y)
{
	x[is.na(x)]=0
	x[is.infinite(x)]=10
	y[is.na(y)]=0
        y[is.infinite(y)]=10
	return(t.test(x,y))



}


######################################################################


order_ranks_files<-function(x)
{

output_character<-x[order(x,decreasing=TRUE)]
return(output_character)

}

#####################################################################

plot_point<-function(x,y,z){
plot(log(x),
     log(y),
     pch=15,
     col='blue',
     lwd=3,
     type="p",
     xlim=c(4,9),
     ylim=c(4,9),
     xlab="repeat two",
     ylab="repeat one"
)
axis(1, col.axis="black", labels=FALSE, lwd=3, las=0)
axis(2, col.axis="black", labels=FALSE, lwd=3,las=2)
legend("topleft", legend=z, col=c("blue"), lwd=5,bty="n",cex=1.5)
}


#####################################################################

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


#####################################################################


totalzscore<-function(x){
CreA_polII2_z<-x
CreA_polII3<-x
pop_sd<-sd(CreA_polII2_z)*sqrt((length(CreA_polII2_z)-1)/(length(CreA_polII2_z)))
pop_mean<-mean(CreA_polII2_z)
for ( i in 1:length(CreA_polII2_z[,1]))
{
#        pop_sd<-sd(CreA_polII2_z[i,])*sqrt((length(CreA_polII2_z[i,])-1)/(length(CreA_polII2_z[i,])))
#        pop_mean<-mean(CreA_polII2_z[i,])
        for ( j in 1:length(CreA_polII2_z[i,]))
        {
                CreA_polII3[i,j]<-(CreA_polII2_z[i,j]-pop_mean)/pop_sd
        }
}
return(CreA_polII3)
}

#####################################################################


rowzscore<-function(x){
CreA_polII2_z<-x
CreA_polII3<-x
for ( i in 1:length(CreA_polII2_z[,1]))
{
        pop_sd<-sd(CreA_polII2_z[i,])*sqrt((length(CreA_polII2_z[i,])-1)/(length(CreA_polII2_z[i,])))
        pop_mean<-mean(CreA_polII2_z[i,])
        for ( j in 1:length(CreA_polII2_z[i,]))
        {
                CreA_polII3[i,j]<-(CreA_polII2_z[i,j]-pop_mean)/pop_sd
        }
}
return(CreA_polII3)
}

###################################################################

rowzscore_revise<-function(x,revised){
CreA_polII2_z<-x
CreA_polII3<-x
for ( i in 1:length(CreA_polII2_z[,1]))
{       
	if (mean(CreA_polII2_z[i,])>revised)
	{
        	pop_sd<-sd(CreA_polII2_z[i,])*sqrt((length(CreA_polII2_z[i,])-1)/(length(CreA_polII2_z[i,])))       
        	pop_mean<-mean(CreA_polII2_z[i,])
        	for ( j in 1:length(CreA_polII2_z[i,]))
        	{       
                	CreA_polII3[i,j]<-(CreA_polII2_z[i,j]-pop_mean)/pop_sd
        	}
	}
	if (mean(CreA_polII2_z[i,])<=revised)
	{
		CreA_polII3[i,]<-seq(0,0,length.out = length(CreA_polII2_z[i,]))
	}
}
return(CreA_polII3)
}




#################################################################

colzscore<-function(x){
CreA_polII2_z<-x
CreA_polII3<-x
for ( i in 1:length(CreA_polII2_z[1,]))
{
        pop_sd<-sd(CreA_polII2_z[,i])*sqrt((length(CreA_polII2_z[,i])-1)/(length(CreA_polII2_z[,i])))
        pop_mean<-mean(CreA_polII2_z[,i])
        for ( j in 1:length(CreA_polII2_z[,i]))
        {
                CreA_polII3[j,i]<-(CreA_polII2_z[j,i]-pop_mean)/pop_sd
        }
}
return(CreA_polII3)
}

##################################################################

library("gplots")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 40)
heatmap_function<-function(x,col_values){
heatmap.2(x,margins=c(5,20),notecex=0,notecol="black",symm=FALSE,trace=c("none"),cexRow = 1.3,density.info=c("none"),col=col_values,Colv=FALSE,Rowv=FALSE)
}

################################################################

get_genes_for_promoter<-function(x){
        gff0_genome_location<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulan_10_gene_final_annotation2.bed")
        x1<-read.table(x)
        x2<-x1$V1
        asd_in_gene_name<-as.character(x2)
        chr_region<-NULL
        chr_start<-NULL
        chr_end<-NULL
        chr_gene_name<-NULL
        chr_strands<-NULL
        gene_location_name<-as.character(gff0_genome_location$V4)
        for ( i in 1:length(asd_in_gene_name) )
        {
                chr_region[i]<-as.character(gff0_genome_location$V1[asd_in_gene_name[i]==as.character(gff0_genome_location$V4)])
                chr_start[i]<-gff0_genome_location$V2[asd_in_gene_name[i]==as.character(gff0_genome_location$V4)]
                chr_end[i]<-gff0_genome_location$V3[asd_in_gene_name[i]==as.character(gff0_genome_location$V4)]
                chr_strands[i]<-gff0_genome_location$V5[asd_in_gene_name[i]==as.character(gff0_genome_location$V4)]
                chr_gene_name[i]<-asd_in_gene_name[i]
        }
        output_bed<-data.frame(region=chr_region,start=chr_start,end=chr_end,gene_name=chr_gene_name,strand=chr_strands)
        write.table(output_bed,file="inter_file.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
        gene_select_region_CreA<-read.table("inter_file.bed")
        gene_select_region_CreA2<-gene_select_region_CreA
	A_value<-NULL
        for( i in 1:length(gene_select_region_CreA[,1]))
        {
                A_value[i]=0
        }
	gene_select_region_CreA2[gene_select_region_CreA[,5]==1,3]=gene_select_region_CreA[gene_select_region_CreA[,5]==1,2]        
	gene_select_region_CreA2[gene_select_region_CreA[,5]==2,2]=gene_select_region_CreA[gene_select_region_CreA[,5]==2,3]	
	gene_select_region_CreA2[gene_select_region_CreA[,5]==1,2]=gene_select_region_CreA[gene_select_region_CreA[,5]==1,2]-500
	gene_select_region_CreA2[gene_select_region_CreA[,5]==2,3]=gene_select_region_CreA[gene_select_region_CreA[,5]==2,3]+500
	return(gene_select_region_CreA2)
        }

############################################

order_matrix<-function(matrix_input,gene_list1,gene_list2)
{
matrix_input2<-matrix(0,nrow=length(gene_list2),ncol=length(matrix_input[1,]))
for ( i in 1:length(gene_list2))
{
	matrix_input2[i,]<-matrix_input[gene_list1==gene_list2[i],]
}
return(matrix_input2)
}


order_matrix_names<-function(matrix_input,gene_list1,gene_list2)
{
matrix_input2<-matrix(0,nrow=length(gene_list2),ncol=length(matrix_input[1,]))
for ( i in 1:length(gene_list2))
{
	matrix_input2[i,]<-matrix_input[gene_list1==gene_list2[i],]
}
colnames(matrix_input2)<-colnames(matrix_input)
rownames(matrix_input2)<-gene_list2

return(matrix_input2)
}

###############################################

normalize_percentage<-function(x)
{
output_matrix<-x
output_matrix<-as.matrix(x)/rowSums(as.matrix(x))
return(output_matrix)
}

plot_image<-function(input_matrix,color_values,numbers_values)
{
        image(  x=seq(1, numbers_values, length.out=numbers_values),
                y=1:length(input_matrix[,1]),
                z=t(input_matrix),
                col=color_values,
                xlab=NA,
                ylab=NA,
                lwd=3,
               # zlim=c(0,k),
                xaxt="n",
                yaxt="n",
                axes=FALSE)

}

col_sums_step<-function(input_matrix)
{
output_matrix<-matrix(0,nrow=(length(input_matrix[,1])/10),ncol=length(input_matrix[1,]))
for(i in 1:(length(input_matrix[,1])/10))
{

	output_matrix[i,]<-colSums(AcuK_binding_level_order[(10*(i-1)+1):(10*(i-1)+10),])

}

return(output_matrix)
}



#######################################

plot_motif<-function(x,y,z,w,a,names){
no_col <- max(count.fields(x, sep = "\t"))
motif_position_list<-read.table(x,sep="\t",fill=TRUE,col.names=1:no_col)
motif_position_list2<-as.matrix(motif_position_list)
motif_position_list_value<-matrix(3000,nrow=length(motif_position_list[,1]),ncol=(no_col-5)/3)

for( i in 1:length(motif_position_list_value[,1]))
{       
        for ( j in 1:((no_col-5)/3))
        {       
                if(!is.na(motif_position_list2[i,3*j-2+5]))
                {       
                        motif_position_list_value[i,j]<-motif_position_list2[i,3*j-2+5]         
                }
        }
}

mytex<-read.table(names)
mytext<-as.character(mytex$V1)
#length_values<-NULL
#for ( k in 1:length(motif_position_list_value[,1]))
#{
#       length_values[k]<-length(motif_position_list_value[k,]) 
#}

#pdf("landscape_YGGRG.pdf", width = 4, height = 60)
#plot (c(0,1700),c(length(motif_position_list_value[,1])+1,length(motif_position_list_value[,1])+1),col = "gray80",xlim=c(0,1700),ylim=c(-1,length(motif_position_list_value[,1])),type="l",lwd=0.01)
plot (c(0,1700),c(length(motif_position_list_value[,1]),length(motif_position_list_value[,1])),col = "gray80",xlim=c(0,2200),ylim=c(-1,length(motif_position_list_value[,1])),type="l",lwd=1,xaxt="n",axes=FALSE)
n=1

for ( i in 1:(length(motif_position_list_value[,1])))
{       
        lines (c(0,1700),c(1*(length(motif_position_list_value[,1])-n),1*(length(motif_position_list_value[,1])-n)),col = "gray80",xlim=c(0,1700),ylim=c(1,length(motif_position_list_value[,1])),type="l",lwd=w)
        text(2100, length(motif_position_list_value[,1])-i+0.5, mytext[i])
        
        for ( j in 1:length(motif_position_list_value[1,]))
        {
        if(as.numeric(motif_position_list_value[i,j])<1700)
        {
        if(motif_position_list[i,5]==1)
        {
        points (c(as.numeric(motif_position_list_value[n,j]),as.numeric(motif_position_list_value[n,j])),c(1*length(motif_position_list_value[,1])-n+a,1*(length(motif_position_list_value[,1])-n+a)),pch=15,col=y,cex=z)
        }
        else if(motif_position_list[i,5]==2)
        {
        points (c(1700-as.numeric(motif_position_list_value[n,j]),1700-as.numeric(motif_position_list_value[n,j])),c(1*length(motif_position_list_value[,1])-n+a,1*(length(motif_position_list_value[,1])-n+a)),pch=15,col=y,cex=z)
        }
        }
        }
        n=n+1

}
abline(v=1500,col="red")
abline(v=1000,col="blue")
#dev.off()
}

plot_motif_add<-function(x,y,z,w,a){
no_col <- max(count.fields(x, sep = "\t"))
motif_position_list<-read.table(x,sep="\t",fill=TRUE,col.names=1:no_col)
motif_position_list2<-as.matrix(motif_position_list)
motif_position_list_value<-matrix(3000,nrow=length(motif_position_list[,1]),ncol=(no_col-5)/3)

for( i in 1:length(motif_position_list_value[,1]))
{
        for ( j in 1:((no_col-5)/3))
        {
                if(!is.na(motif_position_list2[i,3*j-2+5]))
                {
                        motif_position_list_value[i,j]<-motif_position_list2[i,3*j-2+5]
                }
        }
}


#length_values<-NULL
#for ( k in 1:length(motif_position_list_value[,1]))
#{
#        length_values[k]<-length(motif_position_list_value[k,])
#}

#pdf("landscape_YGGRG.pdf", width = 4, height = 60)
n=1
for ( i in 1:(length(motif_position_list_value[,1])))
{

        for ( j in 1:length(motif_position_list_value[1,]))
        {
        if(as.numeric(motif_position_list_value[i,j])<1700)
        {
        if(motif_position_list[i,5]==1)
        {
        points (c(as.numeric(motif_position_list_value[n,j]),as.numeric(motif_position_list_value[n,j])),c(1*length(motif_position_list_value[,1])-n+a,1*(length(motif_position_list_value[,1])-n+a)),pch=15,col=y,cex=z)
        }

        else if(motif_position_list[i,5]==2)
        {
        points (c(1700-as.numeric(motif_position_list_value[n,j]),1700-as.numeric(motif_position_list_value[n,j])),c(1*length(motif_position_list_value[,1])-n+a,1*(length(motif_position_list_value[,1])-n+a)),pch=15,col=y,cex=z)
        }

        }
        }
        n=n+1
}
abline(v=1500,col="red")
abline(v=1000,col="blue")
#dev.off()
}

plot_motif_add_significant<-function(x,y,z,w,a){
no_col <- max(count.fields(x, sep = "\t"))
motif_position_list<-read.table(x,sep="\t",fill=TRUE,col.names=1:no_col)
motif_position_list2<-as.matrix(motif_position_list)
motif_position_list_value<-matrix(3000,nrow=length(motif_position_list[,1]),ncol=(no_col-5)/3)

for( i in 1:length(motif_position_list_value[,1]))
{
        for ( j in 1:((no_col-5)/3))
        {
                if(!is.na(motif_position_list2[i,3*j-2+5]))
                {
                        motif_position_list_value[i,j]<-motif_position_list2[i,3*j-2+5]
                }
        }
}



#length_values<-NULL
#for ( k in 1:length(motif_position_list_value[,1]))
#{
#        length_values[k]<-length(motif_position_list_value[k,])
#}

#pdf("landscape_YGGRG.pdf", width = 4, height = 60)
n=1
for ( i in 1:(length(motif_position_list_value[,1])))
{

        for ( j in 1:length(motif_position_list_value[1,]))
        {
        if(as.numeric(motif_position_list_value[i,j])<1700)
        {
                if(length(peak_regions[(length(motif_position_list_value[,1])-i)==peak_regions[,3],2])>=1)
                {
                        for ( number in 1:length(peak_regions[(length(motif_position_list_value[,1])-i)==peak_regions[,3],2]))
                        {
        #if(motif_position_list_value[i,j]<=peak_regions[(28-i)==peak_regions[,3],2][number]&&motif_position_list_value[i,j]>=peak_regions[(28-i)==peak_regions[,3],1][number])
                        if(motif_position_list[i,5]==1)
                        {
                        if(as.numeric(motif_position_list_value[i,j])<=peak_regions[(length(motif_position_list_value[,1])-i)==peak_regions[,3],2][number]&&as.numeric(motif_position_list_value[i,j])>=peak_regions[(length(motif_position_list_value[,1])-i)==peak_regions[,3],1][number])
                        {
                        points (c(as.numeric(motif_position_list_value[n,j]),as.numeric(motif_position_list_value[n,j])),c(1*length(motif_position_list_value[,1])-n+a,1*(length(motif_position_list_value[,1])-n+a)),pch=15,col=y,cex=z)
                        }
                        }
                if(motif_position_list[i,5]==2)
                {
                        if((1700-as.numeric(motif_position_list_value[i,j]))<=peak_regions[(length(motif_position_list_value[,1])-i)==peak_regions[,3],2][number]&&(1700-as.numeric(motif_position_list_value[i,j]))>=peak_regions[(length(motif_position_list_value[,1])-i)==peak_regions[,3],1][number])
                        points (c(1700-as.numeric(motif_position_list_value[i,j]),1700-as.numeric(motif_position_list_value[i,j])),c(1*length(motif_position_list_value[,1])-i+a,1*(length(motif_position_list_value[,1])-i+a)),pch=15,col=y,cex=z)

                }

        }
        }
        }
        }
        n=n+1
}
abline(v=1500,col="red")
abline(v=1000,col="blue")
#dev.off()
}

#############################

read_motif<-function(x)
{
no_col <- max(count.fields(x, sep = "\t"))
motif_position_list<-read.table(x,sep="\t",fill=TRUE,col.names=1:no_col)
return(motif_position_list)
}

#############################

convert_gene_list<-function(x)
{
#convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/gene_list_for_convert.xls")
#convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/gene_names_convert.tab")
#convertor_genes_first<-as.character(convertor_genes$V2)
#convertor_genes_second<-as.character(convertor_genes$V1)
#return(convertor_genes_first[convertor_genes_second%in%x])

        convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/gene_names_convert.tab")
        convertor_genes_first<-as.character(convertor_genes$V1)
        convertor_genes_second<-as.character(convertor_genes$V2)
        output<-NULL
        d=1
        for ( i in 1:length(x))
        {
                if(x[i]%in%convertor_genes_first)
                {
                        output[d]<-convertor_genes_second[x[i]==convertor_genes_first]
                        d=d+1
                }
        }
        return(output)


}

convert_gene_list_uniq<-function(x)
{
convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/gene_list_for_convert.xls")
convertor_genes_first<-as.character(convertor_genes$V1)
convertor_genes_second<-as.character(convertor_genes$V2)
output<-NULL
for ( i in 1:length(x))
{
output[i]<-convertor_genes_first[convertor_genes_second==x[i]]
}
return(output)
}



#############################

convert_gene_list_sc<-function(x)
{
        convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/S_cerevisiae_sacCer3/annotation/the_convert_gene_lists2.xls")
        convertor_genes_first<-as.character(convertor_genes$V2)
        convertor_genes_second<-as.character(convertor_genes$V1)
        return(convertor_genes_first[convertor_genes_second%in%x])
}

convert_gene_list_sc_order<-function(x)
{
        convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/S_cerevisiae_sacCer3/annotation/the_convert_gene_lists2.xls")
        convertor_genes_first<-as.character(convertor_genes$V2)
        convertor_genes_second<-as.character(convertor_genes$V1)
	output<-NULL
	d=1	
	for ( i in 1:length(x))
	{
		if(x[i]%in%convertor_genes_second)
		{
			output[d]<-convertor_genes_first[x[i]==convertor_genes_second]	
			d=d+1
		}
	}
        return(output)
}

convert_gene_list_sc_order_rev<-function(x)
{
        convertor_genes<-read.table("/Users/dongliguo/Documents/reference_sequences/S_cerevisiae_sacCer3/annotation/the_convert_gene_lists2.xls")
        convertor_genes_first<-as.character(convertor_genes$V2)
        convertor_genes_second<-as.character(convertor_genes$V1)
        output<-NULL
        d=1
        for ( i in 1:length(x))
        {
                if(x[i]%in%convertor_genes_first)
                {
                        output[d]<-convertor_genes_second[x[i]==convertor_genes_first]
                        d=d+1
                }
        }
        return(output)
}




splite_chars<-function(x,y)
{
output<-NULL
for ( i in 1:length(x))
	{
		output[i]<-strsplit(as.character(outputfiles$TF_binding_nums.V2.order.ratio_numbers..)[i],y)[[1]][1]
	}
return(output)
}


#############################

#library("TFBSTools")
#library(seqPattern)
#library(seqinr)
#library(BSgenome.Drerio.UCSC.danRer7)

