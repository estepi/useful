library(limma)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape) 
##################################################################
#SORT by complex
classKD<-read.table("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/classificationFiltered.txt", sep="\t", header = T)
head(classKD); dim(classKD)#312
head(classKD)
allNodes<-as.character(classKD$VTb)
table(classKD$CLASS)
#################################################################
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/NumOfM100/")
edgeListLC<-read.table("edgelist_048_labchip.tab");head(edgeListLC); dim(edgeListLC)
head(edgeListLC)
edgesC<-paste(edgeListLC$source, edgeListLC$target, sep="->"); length(edgesC)#547
#################################################################
edgeListLR<-read.table("edgelist_055_rnaseq.tab");head(edgeListLR); dim(edgeListLR)
head(edgeListLR)
edgesR<-paste(edgeListLR$source, edgeListLR$target, sep="->"); length(edgesR)#389
#################################################################
PPi<-read.table("ppi_all_names.txt", sep="\t", header=T); head(PPi)
PPi_DF<-PPi[,7:8]; head(PPi_DF)
dim(PPi)
io<-match(classKD$VTb, PPi_DF$geneVastO)
id<-match(classKD$VTb, PPi_DF$geneVastD)
io<-io[!is.na(io)]; length(io)#237
id<-id[!is.na(id)]; length(id)#237
ii<-union(io, id)
kdPPI<-PPi_DF[ii,]; head(kdPPI); dim(kdPPI)#466
edgesPPI<-paste(kdPPI$geneVastO, kdPPI$geneVastD, sep="->"); length(edgesPPI)#466 
#################################################################
#solo seleccionar las interacciones de nuestros KDs
int<-read.table("../interactions/edgelist_intN", sep="\t", header=T)
colnames(int)<-c("number","source","target"); head(int)
head(int)
edgesINT<-paste(int$source, int$target, sep="->"); length(edgesINT)#1930
#################################################################
iclip<-read.table("../interactions/edgelist_iclip", sep="\t", header=T)
colnames(iclip)<-c("number","source","target")
head(iclip); dim(iclip)
edgesICLIP<-paste(iclip$source, iclip$target, sep="->"); length(edgesICLIP)#82
#all available nodes:
##################################################################
#edges todos:
n1<-edgesC
#rnaseq
n2<-edgesR
#ppi
n3<-edgesPPI
#iclip
n4<-edgesICLIP
#############UNION#############
edgesALL<-unique(c(n1,n2,n3,n4))
length(edgesALL)#1351
edgesALL
auxN<-data.frame(labchip=edgesALL%in%n1,
                 rnaseq=edgesALL%in%n2,
                 ppi=edgesALL%in%n3,
                 iclip=edgesALL%in%n4)
vennDiagram(auxN)
#finalMatrixCOmparison
forPlot<-auxN;
head(forPlot)
rownames(forPlot)<- edgesALL
forPlot$edges<-edgesALL
forPlot<-forPlot[rowSums(forPlot[,1:4])>=2,]; dim(forPlot); #119
forPlot[forPlot==TRUE]<-1
forPlot[forPlot== FALSE]<-0
head(forPlot)
forPlotMelt = melt(forPlot, id.vars = "edges")
forPlotMelt$NumOfEdge<-rep(1:119, 4)
head(forPlotMelt)
#################################################################
pdf("presence_ausence_edges.pdf", width = 8.27,  height = 11.69)
ggplot(forPlotMelt, aes(x = variable, y = edges, fill = factor(value)))  +
  theme_minimal() +
  geom_tile(position = "identity",  color = "grey") +
  coord_fixed(ratio = 0.5)+
  #coord_equal()+
  scale_fill_manual(values=c("white","grey"))+
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0,-0.5))+
  theme(
    axis.text.x = element_text(angle = 90, hjust=0, vjust=-5,size=3),
    axis.text.y = element_text(size=3),
    axis.ticks.length=unit(0,"cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border=element_blank(),
    legend.position="none",
    axis.line=element_blank()  )  
dev.off()
getwd()

