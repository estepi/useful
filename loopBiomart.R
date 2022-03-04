library(biomaRt)
setwd("~/Documents/biomart/")
beddata <- read.table("0720431_Covered_copy.txt", header = F)
colnames(beddata) <- c("chr", "start", "end")
head(beddata)
beddata$chr <- gsub("chr", "", beddata$chr)

filter_crit <- gsub("chr", "", beddata$V1)
filters <- list(filter_crit, beddata$V2, beddata$V3)

filter.start <- list(beddata$V2)
listMarts(
  "http://www.ensembl.org:80/biomart/martservice?type=registry&requestid=biomaRt"
)
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Funciona
bed.query.chr <- getBM(
  attributes = c("external_gene_name"),
  filters =  c("chromosome_name"),
  values =  filter_crit ,
  mart = ensembl
)

# no funciona
bed.query.chr.range <- getBM(
  attributes = c("external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = filters ,
  mart = ensembl
)



  
  