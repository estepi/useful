#extract genomic features from annotation file
library(GenomicRanges)
library(GenomicFeatures)
genomeTxDb <- makeTxDbFromGFF( "Hsa38.gtf")
genomeTxDb
columns(genomeTxDb)
###########################################
exons <- exonsBy(genomeTxDb, by="gen") #extract exons by gen 
exonsByGene <- as.data.frame(exons@partitioning)[,4:3]
head(exonsByGene)
colnames(exonsByGene)<-c("ENS","NumOfEx")
head(exonsByGene[order(exonsByGene$NumOfEx, decreasing = T),])
head(exonsByGene[order(exonsByGene$NumOfEx),])
monoexonic<-exonsByGene[exonsByGene$NumOfEx==1,]
dim(monoexonic)#22488
dim(exonsByGene)
########################################
write.table(exonsByGene,
            file="NumberOfExonsByGene.tab", 
            sep="\t", col.names = NA)