library(GOstats)
library("org.Hs.eg.db")
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(limma)
#load the objects and produce the matrixes
setwd("TestGO/")
#gO categories:
categFull<-read.csv("../GO_all_hg19.cleaned.txt", sep=",", header=F);
colnames(categFull)<-c("GO","desc")
categ<-read.table("../categories.txt", sep="\t")
colnames(categ)<-"cat"
selected<-match(categ$cat, categFull[,2])
goes<-categFull[selected,]
dir()
######################################
prepareMatrix<-function (rlist, goes, name) 
{
  rn<-lapply(rlist, rownames)
  all<-unique(unlist(rn));
  message(length(all))
  
  pvals <- lapply(rlist, function(x)   {
    if(!is.null(dim(x))){ 
        x%>% select(padjust)  }})
  oddsRatios <- 
  lapply(rlist, function(x)      {
      if(!is.null(dim(x)))       { 
        x%>% select(OddsRatio)  }})
  kds  <-lapply(oddsRatios, function(x) {
        kd<-dim(x)[1] 
        return(kd)})
        
    vv<-unlist(kds)
    nn<-rep(names(vv), vv)
    values<-unlist(oddsRatios)
    pvalsul<-unlist(pvals)
    values[values>20]<-20
    values[pvalsul<0.1]<-1
    rrn<-lapply(oddsRatios, rownames)
    dfPlot<-data.frame(unlist(rn), values,KD=nn)
    colnames(dfPlot)<-c("GO","EF","KD")
    selectedCat<-match(dfPlot$GO,goes$GO)
    labelsGO<-paste(goes$GO, goes$desc, sep=":")
    finaldfPlot<-dfPlot[which(!is.na(selectedCat)),]
#plot1
    file<-paste(name, ".png", sep="")
    message(file)
    png(file)
    gp<-ggplot(finaldfPlot, aes(x = KD, y = GO, fill = EF)) +  geom_tile() +  scale_fill_gradient(low = "white", high = "red")+
    theme_classic()+
    theme(  axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
            axis.text.y = element_text(size = 5)) +
    scale_y_discrete(labels=goes$desc) 
    plot(gp)
    dev.off()
    file<-paste(name, ".pdf", sep="")
    pdf(file, width = 11.69,  height = 8.27)
    plot(gp)
    dev.off()
    ###################################################
    dfPlotSub<-filter(dfPlot, GO%in%goes$GO)
    dfSpread<-spread(dfPlotSub, KD, value=EF, fill = 1)
    #faltan rownames
    file2<-paste(name, "matrix.tab", sep="")
    write.table(dfSpread, file2, sep="\t", col.names = NA)}
#Genes
load("res15GenesHF.RData")
rlistAll<-res1;rm(res1)
name<-"genesAll"
prepareMatrix(rlistAll, goes, name) 
#Genes Up
load("resGenesHFUp.RData")
rlistUp<-res1; rm(res1)
name<-"genesUp"
prepareMatrix(rlistUp, goes, name) 
#Genes Down
load("resGenesHFDown.RData")
rlistDown<-res1; rm(res1)
name<-"genesDown"
prepareMatrix(rlistDown, goes, name) 
#AS
load("res1AS.RData")
rlistAS<-res1; rm(res1)
name<-"AS"
prepareMatrix(rlistAS, goes, name) 

load("resOverlappingAS.RData")
rlistOverlap<-res1; rm(res1)
name<-"GenesAS"
prepareMatrix(rlistOverlap, goes, name) 

load("resOverlappingASUp.RData")
rlistOverlapUp<-res1; rm(res1)
name<-"GenesASUp"
prepareMatrix(rlistOverlapUp, goes, name) 

load("resOverlappingAS.RData")
rlistOverlapDown<-res1; rm(res1)
name<-"GenesASDown"
prepareMatrix(rlistOverlapDown, goes, name) 
