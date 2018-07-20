library(GOstats)
library("org.Hs.eg.db")
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
dPSI<-read.csv("", sep="\t", header=T)
##############################################################
vtNames<-as.character(dPSI[,1])
ensembleIds<-read.csv("fullConversionTable_vtools_ens.csv", sep=",", header=T)
small<-ensembleIds[,1:2]
summary(!is.na(small$GENEID))
notNa<-small[!is.na(small$GENEID),]
unique<-notNa[!duplicated(notNa), ]
ii<-match(vtNames, as.character(unique$GENE))
dPSI_Ens<-dPSI
dPSI_Ens$Group.1= NULL
rownames(dPSI_Ens)[1]<-"MISSING"
hh<-match(as.character(rownames(dPSI_Ens)), as.character(ensembleIds$GENE) )
dPSI_Ens$geneID<-as.character(ensembleIds$GENEID[hh])
#aca convertir a ENTREZ tb
xx <- as.list(org.Hs.egENSEMBL2EG); length(xx) #28589
# Gets the entrez gene IDs for the first five Ensembl IDs
#edfine the UNIVERSE, the same for everyKD
entrez<-unlist(xx)
universeGeneIds<-as.character(dPSI_Ens$geneID)
universeEntrezI<-match(as.character(universeGeneIds), names(entrez)); head(entrez);names(entrez);length(entrez)
dPSI_Ens$Entrez<-entrez[universeEntrezI]
uc<-dPSI_Ens$Entrez[!is.na(dPSI_Ens$Entrez)]

head(dPSI_Ens)
ASAff<-data.frame(dPSI_Ens$geneID ,
                  dPSI_Ens$Entrez)
dim(ASAff)#5873
write.table(ASAff, file="Affected_AS.tab", sep="\t")
##########################################################
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
}
)

#rescue all the categories that  satisfies at least in one KD pval <0.05
rn<-lapply(res, rownames)
all<-unique(unlist(rownames));length(all)
rownames1<-lapply(res1, rownames);length(res1)
length(unique(unlist(rownames1)))
######################################################
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
########################################################################################
oddsRatios <- 
  lapply(res1, function(x)  
      {
      if(!is.null(dim(x)))
         { 
        x%>% select(OddsRatio)  
          }
    }
    )

kds  <-lapply(oddsRatios, function(x) {
  kd<-dim(x)[1]
  return(kd)
}
)

vv<-unlist(kds)
nn<-rep(names(vv), vv); length(nn)
values<-unlist(oddsRatios)
length(values)
values[values>20]<-20
summary(values)
rrn<-lapply(oddsRatios, rownames)
dfPlot<-data.frame(unlist(rrn), values,KD=nn)

colnames(dfPlot)<-c("GO","EF","KD")
#############################################################
pvalues<-
    lapply(res1, function(x)  
    {
      if(!is.null(dim(x)))
      { 
        x%>% select(padjust)  
      }
    }
    )

kds2<-lapply(pvalues, function(x) {
  kd<-dim(x)[1]
  return(kd)})
vv2<-unlist(kds2)
nn2<-rep(names(vv2), vv2); length(nn2)
pvaluestotal<-unlist(pvalues)
length(pvaluestotal)
pvaluestotal[pvaluestotal>0.1]<-1
table(pvaluestotal)
rrn2<-lapply(pvalues, rownames)
dfPlotPVAL<-data.frame(unlist(rrn2), pvaluestotal, KD=nn2)

colnames(dfPlotPVAL)<-c("GO","PV","KD")

######################################################################
#select the preferred ones:
categFull<-read.csv("../GO_all_hg19.cleaned.txt", sep=",", header=F);
colnames(categFull)<-c("GO","desc")
categ<-read.table("../categories.txt", sep="\t")
colnames(categ)<-"cat"
selected<-match(categ$cat, categFull[,2])
goes<-categFull[selected,]

selectedCat<-match(dfPlot$GO,goes$GO)
finaldfPlot<-dfPlot[which(!is.na(selectedCat)),]
#####################################################################
pdf("GO_selected_cat.pdf", width = 11.69,  height = 8.27)
ggplot(dat.mA, aes(x = KD, y = GO, fill = EF)) +  geom_tile() +  scale_fill_gradient(low = "white", high = "red")+
  theme_classic()+
  theme(  axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
            axis.text.y = element_text(size = 5)
            )
  dev.off()
##################################
dfPlotSub<-filter(dfPlot, GO%in%goes$GO)
dfSpread<-spread(dfPlotSub, KD, value=EF, fill = 1)
rownames(dfSpread)<-dfSpread$GO
dfSpread$GO<-NULL
pheatmap(as.matrix(t(dfSpread)), cellwidth = 8, cellheight = 2, fontsize = 2)
###################PVALUES##################
dfPlotSub<-filter(dfPlotPVAL, GO%in%goes$GO)
dfSpread<-spread(dfPlotSub, KD, value=PV, fill = 1)
rownames(dfSpread)<-dfSpread$GO
dfSpread$GO<-NULL
pheatmap(log(t(as.matrix(dfSpread))))
colors<-colorRampPalette(rev(brewer.pal(n=5,name="RdYlBu")))(255)
pheatmap(my_matrix,cluster_cols = FALSE,cellwidth = 30,fontsize = 7,height = 40,show_rownames = FALSE, col=colors)
pheatmap(log(t(as.matrix(dfSpread))), fontsize = 4,height = 50, col=colors)