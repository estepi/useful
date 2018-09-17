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
uc<-dPSI_Ens$Entrez[!is.na(dPSI_Ens$Entrez)]; length(uc)
##########################################################
dPSI_Ens[1:10,300:314]
dim(dPSI_Ens)
testDF<-dPSI_Ens[,1:312]
sort(colSums(testDF))
##############################
testDF<-data.frame(CWC22=testDF[,colnames(testDF)=="CWC22"],SF3B1=testDF[,colnames(testDF)=="SF3B1"])

colnames(testDF)

ll<-apply(testDF, 2, function(x) {
    se<-dPSI_Ens$Entrez[x==1]
    sse<-se[!is.na(se)]
    return(sse)} )
length(ll)

res1Cond<-lapply(ll, function(x){
    if( length(x) > 0 ) 
  {
    params<-new("GOHyperGParams",
                geneIds=x,
                universeGeneIds=uc,
                annotation=org.Hs.eg.db,
                ontology="BP",
                pvalueCutoff=1,
                conditional=TRUE,
                testDirection="over")
    hgOver<-hyperGTest(params)
    df<-summary(hgOver)
    df$padjust<-p.adjust(df$Pvalue,method = "fdr" )
    rownames(df)<-df$GOBPID
    return(df)  
  }
})
save(res1Cond, file="Sf3b1_cwc22_cond.RData")
save(res1, file="ASall.RData")
res1<-res1Cond
#rescue all the categories that  satisfies at least in one KD pval <0.05
######################################################################
rn<-lapply(res1, rownames)
all<-unique(unlist(rn));length(all)#9200
######################################################
term <- 
  lapply(res1, function(x)  
  {
    if(!is.null(dim(x)))
    { 
      x%>% select(Term)  
    }
  }
  )

pvals <- 
  lapply(res1, function(x)  
  {
    if(!is.null(dim(x)))
    { 
      x%>% select(padjust)  
    }
  }
  )

oddsRatios <- 
  lapply(res1, function(x)  
  {
    if(!is.null(dim(x)))
    { 
      x%>% select(OddsRatio)  
    }
  }
  )

Count <- 
  lapply(res1, function(x){
    if(!is.null(dim(x)))    { 
      x%>% select(Count)  
    }
  })

Size <- 
  lapply(res1, function(x){
    if(!is.null(dim(x)))    { 
      x%>% select(Size)  
    }
  })




kds  <-lapply(oddsRatios, function(x) {
  kd<-dim(x)[1]
  return(kd)
}
)
length(oddsRatios)
vv<-unlist(kds)
nn<-rep(names(vv), vv); length(nn)
values<-unlist(oddsRatios)
length(values)
pvalsul<-unlist(pvals)

rrn<-lapply(oddsRatios, rownames)
length(values);head(values)

sizetot<-unlist(Size)
countsTot<-unlist(Count)
terms<-unlist(term)

dfPlot<-data.frame(unlist(rn), values,KD=nn); head(dfPlot)
colnames(dfPlot)<-c("GO","EF","KD");dfPlot[1:5,]
dfSpreadEF<-spread(dfPlot, KD, value=EF, fill = 1);head(dfSpreadEF)

dfPlot2<-data.frame(unlist(rn), pvalsul,KD=nn);head(dfPlot2)
colnames(dfPlot2)<-c("GO","Pval","KD")
dfSpreadPv<-spread(dfPlot2, KD, value=Pval, fill = 1);head(dfSpreadPv)

dfPlot3<-data.frame(unlist(rn), sizetot,KD=nn);head(dfPlot3)
colnames(dfPlot3)<-c("GO","Size","KD")
dfSpreadSize<-spread(dfPlot3, KD, value=Size, fill = 1)

dfPlot4<-data.frame(unlist(rn), countsTot,KD=nn);head(dfPlot4)
colnames(dfPlot4)<-c("GO","Count","KD")
dfSpreadCount<-spread(dfPlot4, KD, value=Count, fill = 1)

dfPlot5<-data.frame(unlist(rn), terms,KD=nn);head(dfPlot5)
colnames(dfPlot5)<-c("GO","Term","KD")
dfSpreadTerms<-spread(dfPlot5, KD, value=Term, fill = 1)

length(which(dfSpreadEF$GO == dfSpreadCount$GO))

dfFinal<-data.frame(dfSpreadEF, dfSpreadPv, dfSpreadCount, dfSpreadSize)
ii<-match(dfFinal$GO, dfPlot5$GO)
dfFinal$description<-dfPlot5$Term[ii]
dim(dfFinal)
head(dfFinal)
getwd()
dfFinal$universe<-length(uc)
dfFinal$totCWC22<-length(ll[[1]])
dfFinal$totSF3B1<-length(ll[[2]])
write.table(dfFinal,"GO_SF3B1_CWC22_EF_description.txt", sep="\t")
############################################
