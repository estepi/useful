library(GOstats)
library("org.Hs.eg.db")
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(limma)
genesLRT <-
  read.table("lrt_table_hg19_HF.tab",
             header = T,
             row.names = 1)
dim(genesLRT)#6947
genesLRT[1:5, 306:308]
genesLRT$padjust <- p.adjust(genesLRT$PValue)
onlyValues <- genesLRT[, 1:305]
dim(onlyValues)#6947 filas
affected <- onlyValues
kdsNames <- sub("logFC.group", "", colnames(affected))
length(kdsNames)
colnames(affected) <- kdsNames


affected[abs(affected) < 1.5] <- 0
affected[abs(affected) > 1.5] <- 1
affected[1:5, 1:5]
dim(affected)#6947

####################################
xx <- as.list(org.Hs.egENSEMBL2EG)
length(xx) #28589
# Gets the entrez gene IDs for the first five Ensembl IDs
#edfine the UNIVERSE, the same for everyKD
entrez <- unlist(xx)
universeGeneIds <- as.character(rownames(affected))
universeEntrezI <-
  match(as.character(universeGeneIds), names(entrez))
head(entrez)
names(entrez)
length(entrez)
affected$Entrez <- entrez[universeEntrezI]
uc <- affected$Entrez[!is.na(affected$Entrez)]
length(uc)#6700 com ENTREZ ID
affected[1:10, 300:306]
dim(affected)
head(affected)
testDF <- affected[, 1:305]
colnames(testDF) <- kdsNames
df <- as.data.frame(colSums(testDF))
colnames(df) <- "number"
df$KD <- rownames(df)
####################################
setwd("GO/TestGO/")
pdf("GO_KDs_GEHF.pdf", width = 11.69,  height = 8.27)
ggplot(df, aes(x = reorder(KD, number), y = number)) + geom_bar(stat = "identity", width =
                                                                  0.5) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 3
  ))
dev.off()
####################################
ll <- apply(testDF, 2, function(x) {
  se <- affected$Entrez[x == 1]
  sse <- se[!is.na(se)]
  return(sse)
})
length(ll)
length(uc)

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
head(res1[[1]])
save(res1, file = "res15GenesHF.RData")
rn <- lapply(res1, rownames)
all <- unique(unlist(rn))
length(all)#8
######################################################
pvals <-
  lapply(res1, function(x)
  {
    if (!is.null(dim(x)))
    {
      x %>% select(padjust)
    }
  })

oddsRatios <-
  lapply(res1, function(x)
  {
    if (!is.null(dim(x)))
    {
      x %>% select(OddsRatio)
    }
  })


kds  <- lapply(oddsRatios, function(x) {
  kd <- dim(x)[1]
  return(kd)
})

length(oddsRatios)
vv <- unlist(kds)
nn <- rep(names(vv), vv)
length(nn)
values <- unlist(oddsRatios)
length(values)
pvalsul <- unlist(pvals)
length(pvalsul)#1031845
values[values > 20] <- 20
summary(values)
values[pvalsul < 0.1] <- 1
length(values)#1031845
rrn <- lapply(oddsRatios, rownames)
length(values)
head(values)
dfPlot <- data.frame(unlist(rn), values, KD = nn)
head(dfPlot)
colnames(dfPlot) <- c("GO", "EF", "KD")
str(dfPlot)

setwd("GO/TestGO/")
categFull <-
  read.csv("../GO_all_hg19.cleaned.txt",
           sep = ",",
           header = F)

head(categFull)
colnames(categFull) <- c("GO", "desc")
categ <- read.table("../categories.txt", sep = "\t")
colnames(categ) <- "cat"
head(categ)
selected <- match(categ$cat, categFull[, 2])
goes <- categFull[selected, ]
head(goes)
dim(goes)
selectedCat <- match(dfPlot$GO, goes$GO)
finaldfPlot <- dfPlot[which(!is.na(selectedCat)), ]
head(finaldfPlot)

pdf("GO_GENES_selected_cat_log15_pval.pdf",
    width = 11.69,
    height = 8.27)
ggplot(finaldfPlot, aes(x = KD, y = GO, fill = EF)) +  geom_tile() +  scale_fill_gradient(low = "white", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 3
  ),
  axis.text.y = element_text(size = 5))
dev.off()
getwd()
########################################################
dfPlotSub <- filter(dfPlot, GO %in% goes$GO)
dim(dfPlotSub)
dfSpread <- spread(dfPlotSub, KD, value = EF, fill = 1)
head(dfSpread)
dim(dfSpread)
dfSpread[1:5, 1:5]
rownames(dfSpread) <- dfSpread$GO
dfSpread$GO <- NULL
pheatmap(as.matrix(t(dfSpread)))
pheatmap(
  as.matrix(t(dfSpread)),
  cellwidth = 10,
  cellheight = 1,
  fontsize = 5
)
#chequear list of genes analyzed by AS and by GE:
genesAff <- as.data.frame(row.names(genesLRT))
dim(genesAff)
head(genesAff)

write.table(data.frame(row.names(affected),
                       affected$Entrez),
            file = "Affected_GENES.tab",
            sep = "\t")
#######################vennPlot:
TotalAffected <-
  union(genesAff$`row.names(genesLRT)`, ASAff$dPSI_Ens.geneID)
auxdf <- data.frame(
  TotalAffected %in% genesAff$`row.names(genesLRT)`,
  TotalAffected %in% ASAff$dPSI_Ens.geneID
)
