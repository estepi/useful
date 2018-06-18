#
library(biomaRt)
listMarts()
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")  
filters = listFilters(ensembl); filters #  ensembl_gene_id
attributes = listAttributes(ensembl); attributes
write.table(attributes, file="ensembl_attribute.txt", sep="\t", quote=F, col.names = NA)
ensIds<-read.csv("SampleClassificationFull.csv", stringsAsFactors = F)
#query:
getBM(attributes = c("start_position","end_position"),
      filters    = "ensembl_gene_id",
      values     = "ENSG00000003756", 
      mart       = ensembl)
#######################################################
df1<-getBM(attributes = c("ensembl_gene_id","start_position","end_position","transcript_count"),
           filters    = "ensembl_gene_id",
           values     = ensIds$Ensemble.Gene.ID, 
           mart       = ensembl)
head(df1)
#######################################################
df2<-getBM(attributes = 
             c("ensembl_gene_id",
               "ensembl_exon_id"),
      filters    = "ensembl_gene_id",
      values     = ensIds$Ensemble.Gene.ID, 
      mart       = ensembl)
#number of exons:
numOfEx<-table(df2$ensembl_gene_id)
#######################################################
ii<-match(df1$ensembl_gene_id, names(numOfEx))
Glength<-df1$end_position - df1$start_position
finalDf<-data.frame(df1, numOfEx[ii], Glength)
#######################################################
vtnames<-read.csv("completeConversionTable_ENSEMBL_VAST.txt", sep="\t")
head(vtnames)
vv<-match(finalDf$ensembl_gene_id, vtnames$FGeneID)
length(which(is.na(vv)))#0
finalDf$vtools<-as.character(vtnames$VT[vv]); head(finalDf)
write.table(finalDf,"NumOfExons.tab", sep="\t", quote=FALSE, col.names = NA)


