library(GOstats)
library("org.Hs.eg.db")
library(dplyr)
library(ggplot2)
dPSI<-read.csv("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/ensembl/forGO/df1_dPSI15t.txt", sep="\t", header=T)
head(dPSI)
##############################################################
vtNames<-as.character(dPSI[,1])
ensembleIds<-read.csv("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/ensembl/inputFiles/fullConversionTable_vtools_ens.csv", sep=",", header=T)
small<-ensembleIds[,1:2]
notNa<-small[!is.na(small$GENEID),]
duplicated(notNa)
unique<-notNa[!duplicated(notNa), ]
ii<-match(vtNames, as.character(unique$GENE))
length(which(is.na(ii)))#41
#queryes for GO
dPSI_Ens<-dPSI
dPSI_Ens$Group.1=NULL
rownames(dPSI_Ens)[1]<-"MISSING"
hh<-match(as.character(rownames(dPSI_Ens)), as.character(ensembleIds$GENE) )
dPSI_Ens$geneID<-as.character(ensembleIds$GENEID[hh])
head(dPSI_Ens$geneID)
#aca convertir a ENTREZ tb
xx <- as.list(org.Hs.egENSEMBL2EG); length(xx) #28589
# Gets the entrez gene IDs for the first five Ensembl IDs
#edfine the UNIVERSE, the same for everyKD
entrez<-unlist(xx)
universeGeneIds<-as.character(dPSI_Ens$geneID)
universeEntrezI<-match(as.character(universeGeneIds), names(entrez)); head(entrez);names(entrez);length(entrez)
dPSI_Ens$Entrez<-entrez[universeEntrezI]
uc<-dPSI_Ens$Entrez[!is.na(dPSI_Ens$Entrez)]
##########################################################
dPSI_Ens[1:10,300:314]
dim(dPSI_Ens)
testDF<-dPSI_Ens[,1:312]
colSums(testDF)
ll<-apply(testDF, 2, function(x) {
    se<-dPSI_Ens$Entrez[x==1]
    sse<-se[!is.na(se)]
    return(sse)} )
res1<-lapply(ll, function(x){
    if( length(x) > 0 ) 
  {
    params<-new("GOHyperGParams",
                geneIds=x,
                universeGeneIds=uc,
                annotation=org.Hs.eg.db,
                ontology="BP",
                pvalueCutoff=1,
                conditional=FALSE,
                testDirection="over")
    hgOver<-hyperGTest(params)
    df<-summary(hgOver)
    df$padjust<-p.adjust(df$Pvalue,method = "fdr" )
    rownames(df)<-df$GOBPID
    return(df)  
  }
})
save(res1, file="res1AS.RData")
#rescue all the categories that  satisfies at least in one KD pval <0.05
rn<-lapply(res, rownames)
all<-unique(unlist(rownames));length(all)
rownames1<-lapply(res1, rownames);length(res1)
length(unique(unlist(rownames1)))
######################################################
head(res1[[9]])
oddsRatios <- 
  lapply(res1, function(x)  
      {
      if(!is.null(dim(x)))
         { 
        x%>% select(OddsRatio)  
          }
    }
    )
######################################################
kds  <-lapply(oddsRatios, function(x) {
  kd<-dim(x)[1]
  return(kd)
}
)
length(oddsRatios)
######################################################
padjust <- 
  lapply(res1, function(x)  
  {
    if(!is.null(dim(x)))
    { 
      x%>% select(padjust)  
    }
  }
  )
####################################
kds2  <-lapply(padjust, function(x) {
  kd<-dim(x)[1]
  return(kd)
}
)
######################################################
vv<-unlist(kds)
vv2<-unlist(kds2)
nn<-rep(names(vv), vv); length(nn)
nn2<-rep(names(vv2), vv2); length(nn2)

values<-unlist(oddsRatios)
valuesC<-unlist(padjust)
length(values)
length(valuesC)
values[values>20]<-20
values[values>20]<-20
values[valuesC>0.1]<-1
rrn<-lapply(oddsRatios, rownames)
length(values)
head(values)
dfPlot<-data.frame(unlist(rrn), values,KD=nn)
head(dfPlot)
colnames(dfPlot)<-c("GO","EF","KD")
#############################################################
######################################################################
#select the preferred ones:
setwd("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/GO/TestGO/")
categFull<-read.csv("../GO_all_hg19.cleaned.txt", sep=",", header=F);
colnames(categFull)<-c("GO","desc")
categ<-read.table("../categories.txt", sep="\t")
colnames(categ)<-"cat"
selected<-match(categ$cat, categFull[,2])
goes<-categFull[selected,]
head(goes)
selectedCat<-match(dfPlot$GO,goes$GO)
finaldfPlot<-dfPlot[which(!is.na(selectedCat)),]
head(finaldfPlot)
#####################################################################
dat.mA <- mutate(finaldfPlot, variable = reorder(GO, EF))
pdf("GO_selected_cat_pvalue.pdf", width = 11.69,  height = 8.27)
dat.mA, 
ggplot(data=filter(dat.mA, EF = >2), aes(x = KD, y = GO, fill = EF)) +  geom_tile() + 
  scale_fill_gradient(low = "pink", high = "red", breaks = seq(1:20))+
  theme_classic()+
  theme(  axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
            axis.text.y = element_text(size = 3)
            )
  dev.off()
##################################
#old plot with df only affected
affected<-sort(unlist(lapply(res, function(x) { dim(x)[1] })))
  barplot(affected) #OK
  df<-as.data.frame(affected);head(df)
  colnames(df)<-"number"
  df$KD<-rownames(df)
  pdf("GO_KDs.pdf", width = 11.69,  height = 8.27)
  ggplot(df, aes(x = reorder(KD, number), y = number)) + geom_bar(stat = "identity", width=0.5) +
    theme(  axis.text.x = element_text(angle = 90, hjust = 1, size = 3))
  dev.off()
  getwd()
  