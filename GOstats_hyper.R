library(GOstats)
library("org.Hs.eg.db")
library(dplyr)
library(ggplot2)
dPSI <- read.csv("df1_dPSI15t.txt", sep = "\t", header = T)
head(dPSI)
##############################################################
vtNames <- as.character(dPSI[, 1])
ensembleIds <-
  read.csv("fullConversionTable_vtools_ens.csv",
           sep = ",",
           header = T)
small <- ensembleIds[, 1:2]
notNa <- small[!is.na(small$GENEID), ]
duplicated(notNa)
unique <- notNa[!duplicated(notNa),]
ii <- match(vtNames, as.character(unique$GENE))
length(which(is.na(ii)))#41
#queryes for GO
dPSI_Ens <- dPSI
dPSI_Ens$Group.1 = NULL
rownames(dPSI_Ens)[1] <- "MISSING"
hh <-
  match(as.character(rownames(dPSI_Ens)), as.character(ensembleIds$GENE))
dPSI_Ens$geneID <- as.character(ensembleIds$GENEID[hh])
head(dPSI_Ens$geneID)
#aca convertir a ENTREZ tb
xx <- as.list(org.Hs.egENSEMBL2EG)
length(xx) #28589
# Gets the entrez gene IDs for the first five Ensembl IDs
#edfine the UNIVERSE, the same for everyKD
entrez <- unlist(xx)
##############################
CatSize <-
  data.frame(sort(table(gos$go_id), decreasing = T))
head(CatSize)
dim(CatSize)
write.table(CatSize, "CatSize.txt", sep = "\t")
write.table(gos, "allGos.txt", sep = "\t")
gos[gos$go_id == "GO:0034470", ]
CatSize[CatSize$Var1 == "GO:0034470", ]
#############################################################
universeGeneIds <- as.character(dPSI_Ens$geneID)
universeEntrezI <-
  match(as.character(universeGeneIds), names(entrez))
head(entrez)
names(entrez)
length(entrez)
dPSI_Ens$Entrez <- entrez[universeEntrezI]
uc <- dPSI_Ens$Entrez[!is.na(dPSI_Ens$Entrez)]
##########################################################
testDF <- dPSI_Ens[, 1:312]

ll <- apply(testDF, 2, function(x) {
  se <- dPSI_Ens$Entrez[x == 1]
  sse <- se[!is.na(se)]
  return(sse)
})

res1 <- lapply(ll, function(x) {
  if (length(x) > 0)
  {
    params <- new(
      "GOHyperGParams",
      geneIds = x,
      universeGeneIds = uc,
      annotation = org.Hs.eg.db,
      ontology = "BP",
      pvalueCutoff = 1,
      conditional = FALSE,
      testDirection = "over"
    )
    hgOver <- hyperGTest(params)
    df <- summary(hgOver)
    df$padjust <- p.adjust(df$Pvalue, method = "fdr")
    rownames(df) <- df$GOBPID
    return(df)
  }
})
#rescue all the categories that  satisfies at least in one KD pval <0.05
######################################################################

