library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
filters #  ensembl_gene_id
attributes = listAttributes(ensembl)
attributes
write.table(
  attributes,
  file = "ensembl_attribute.txt",
  sep = "\t",
  quote = F,
  col.names = NA
)
ensIds <-
  read.csv("SampleClassificationFull.csv", stringsAsFactors = F)
#query:
getBM(
  attributes = "ensembl_gene_id",
  filters    = "ensembl_transcript_id",
  values     = "ENST00000348261",
  mart       = ensembl
)
#######################################################
#benjamin predicion
ben3 <-
  read.table("ben/output-ensg-3-sigma-input-estefi-ensgs.txt", header = T)
head(ben3)
ben4 <-
  read.table("ben/output-ensg-4-sigma-input-estefi-ensgs.txt", header = T)
head(ben4)
ben5 <-
  read.table("ben/output-ensg-5-sigma-input-estefi-ensgs.txt", header = T)
head(ben5)
ben6 <-
  read.table("ben/output-ensg-6-sigma-input-estefi-ensgs.txt", header = T)
head(ben6)
#########################################################
ENSGID3 <-
  getBM(
    attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
    filters    = "ensembl_transcript_id",
    values     = unique(ben3$transcript_enst),
    mart       = ensembl
  )
#conversion table
dfBenInt <-
  data.frame(source = ben3$rbp_ensg,
             target = ben3$transcript_enst)
dim(dfBenInt)
head(dfBenInt)
ii <- match(dfBenInt$target, ENSGID3$ensembl_transcript_id)
dfBenInt$targetENS <- ENSGID3$ensembl_gene_id[ii]
head(dfBenInt)
############################################################
conversion <-
  read.table(
    "~ensembl/completeConversionTable_ENSEMBL_VAST.txt",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )

head(conversion)
table(dfBenInt$targetENS)
source2vast <- match(as.character(dfBenInt$source),
                     as.character(conversion$FGeneID))
length(which(is.na(source2vast)))#0
target2vast <- match(as.character(dfBenInt$targetENS),
                     as.character(conversion$FGeneID))
length(which(is.na(target2vast)))#7873

dfBenInt$sourceVT <- conversion$VT[source2vast]
dfBenInt$targetVT <- conversion$VT[target2vast]
dfNet <- data.frame(dfBenInt$sourceVT, dfBenInt$targetVT)
dim(dfNet)
dfNet <- na.omit(dfNet)
dim(dfNet)
rbpN <- graph_from_data_frame(dfNet, directed = F)
