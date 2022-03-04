library(GOstats)
library(biomaRt)
library(data.table)
library("org.Hs.eg.db")
source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
dPSI <- read.csv("df1_dPSI15t.txt", sep = "\t", header = T)
head(dPSI)
##############################################################
vtNames<-as.character(dPSI[,1])
ensembleIds <-
  read.csv("fullConversionTable_vtools_ens.csv",
           sep = ",",
           header = T)
small <- ensembleIds[, 1:2]
notNa <- small[!is.na(small$GENEID), ]
duplicated(notNa)
unique <- notNa[!duplicated(notNa),]
ii <- match(vtNames, as.character(unique$GENE))
length(which(is.na(ii)))
head(notNa)
write.table(vtNames[is.na(ii)], file = "GeneNames_without_ENSIds.txt", quote =
              F)
#queryes for GO
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
df6 <- getBM(
  attributes =
    c(
      "ensembl_gene_id",
      "go_id",
      "name_1006",
      "namespace_1003",
      "hgnc_id",
      "hgnc_symbol"
    ),
  filters    = "ensembl_gene_id",
  values     = as.character(ensembleIds$GENEID),
  mart       = ensembl
)

#query:
bb <- match(df6$ensembl_gene_id, unique$GENEID)
df7 <-
  cbind(df6,
        vtName = unique$GENE[bb],
        equal = df6$hgnc_symbol == unique$GENE[bb])

write.table(df7, file = "Ensembl_GO_hgnc_vtools.txt", quote = F)
##############################################################
vv <- match(as.character(dPSI$Group.1), as.character(df7$vtName))
length(which(is.na(vv)))#159
length(unique(dPSI$Group.1))

hh <- match(as.character(dPSI$Group.1), as.character(df7$hgnc_symbol))
length(which(is.na(hh)))#709
length(which(is.na(hh)))#
dPSI[1, 1:10]
dPSI_Ens <- dPSI
dPSI_Ens$Group.1 <- NULL
dPSI_Ens[1:2, 1:10]
hh <-
  match(as.character(rownames(dPSI_Ens)), as.character(ensembleIds$GENE))
dPSI_Ens$geneID <- as.character(ensembleIds$GENEID[hh])
head(dPSI_Ens$geneID)

gg <- dPSI_Ens[, 8] == 1
colnames(dPSI_Ens[8])
TestGeneIds <- as.character(dPSI_Ens$geneID)[gg]
length(TestGeneIds)#564
universeGeneIds <- as.character(dPSI_Ens$geneID)
length(universeGeneIds)#5873
xx <- as.list(org.Hs.egENSEMBL2EG)
length(xx) #28589
# Gets the entrez gene IDs for the first five Ensembl IDs
entrez <- unlist(xx)
head(entrez)
names(entrez)
length(entrez)
universeEntrezI <-
  match(as.character(universeGeneIds), names(entrez))

length(universeEntrezI)
length(which(is.na(universeEntrezI)))#240
universe <- entrez[universeEntrezI]
length(universe)
uc <- universe[!is.na(universe)]
length(uc)#5633
###########################################################################
selectedEntrezI <- match(as.character(TestGeneIds), names(entrez))

length(selectedEntrezI)#564
length(which(is.na(selectedEntrezI)))#25
selected <- entrez[selectedEntrezI]
sc <- selected[!is.na(selected)]
length(sc)#539
params <- new(
  "GOHyperGParams",
  geneIds = sc,
  universeGeneIds = uc,
  annotation = org.Hs.eg.db,
  ontology = "BP",
  pvalueCutoff = 0.05,
  conditional = FALSE,
  testDirection = "over"
)
hgOver <- hyperGTest(params)
df <- summary(hgOver)
df$padjust <- p.adjust(df$Pvalue)
write.table(df, "AC008073.5_GO.txt", sep="\t")
head(df)
######################################################
test2 <- dPSI_Ens[, colnames(dPSI_Ens) == "SF3B1"]
sum(test2)
test2 == 1
sf3b1Selected <- as.character(dPSI_Ens$geneID[test2 == 1])
length(sf3b1Selected)
selectedEntrezI <-
  match(as.character(sf3b1Selected), names(entrez))

entrezSF3 <- entrez[selectedEntrezI]
sc <- entrezSF3[!is.na(entrezSF3)]
length(sc)
length(entrezSF3)#2427
#######################################################
params <- new(
  "GOHyperGParams",
  geneIds = sc,
  universeGeneIds = uc,
  annotation = org.Hs.eg.db,
  ontology = "BP",
  pvalueCutoff = 1,
  conditional = FALSE,
  testDirection = "over"
)
hgOver <- hyperGTest(params)
df2 <- summary(hgOver)
df2$padjust <- p.adjust(df2$Pvalue)
write.table(df2, "SF3B1_GO_over_classical.txt", sep = "\t")
#######################################################
params <- new(
  "GOHyperGParams",
  geneIds = sc,
  universeGeneIds = uc,
  annotation = org.Hs.eg.db,
  ontology = "BP",
  pvalueCutoff = 1,
  conditional = FALSE,
  testDirection = "under"
)
hgOver <- hyperGTest(params)
df3 <- summary(hgOver)
df3$padjust <- p.adjust(df3$Pvalue)
write.table(df3, "SF3B1_GO_under_classical.txt", sep = "\t")
#######################################################
params <- new(
  "GOHyperGParams",
  geneIds = sc,
  universeGeneIds = uc,
  annotation = org.Hs.eg.db,
  ontology = "BP",
  pvalueCutoff = 1,
  conditional = TRUE,
  testDirection = "over"
)
hgOver <- hyperGTest(params)
df4 <- summary(hgOver)
df4$padjust <- p.adjust(df4$Pvalue)
write.table(df4, "SF3B1_GO_over_conditional.txt", sep = "\t")
dim(df4)
