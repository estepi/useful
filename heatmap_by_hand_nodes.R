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
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/NumOfM100/topology/")
#################################################################
edgeListLC<-read.table("labchip048_edgelist.tab");head(edgeListLC); dim(edgeListLC)#547
#################################################################
edgeListLR<-read.table("RNA-Seq_compatible_edgelist.tab");head(edgeListLR); dim(edgeListLR)#389
#################################################################
PPi<-read.table("ppi_all_names.txt", sep="\t", header=T); head(PPi)
PPi_DF<-PPi[,7:8]; head(PPi_DF); dim(PPi)
io<-match(classKD$VTb, PPi_DF$geneVastO)
id<-match(classKD$VTb, PPi_DF$geneVastD)
io<-io[!is.na(io)]; length(io)#237
id<-id[!is.na(id)]; length(id)#237
ii<-union(io, id)
kdPPI<-PPi_DF[ii,]; head(kdPPI); dim(kdPPI)#466
#################################################################
int<-read.table("../interactions/edgelist_intN", sep="\t", header=T)
colnames(int)<-c("number","source","target"); head(int)
#################################################################
iclip<-read.table("../interactions/edgelist_iclip", sep="\t", header=T)
colnames(iclip)<-c("number","source","target")
head(iclip); dim(iclip)
#################################################################
exones30all<-read.table("30best_exons_allKDs_edgelist.tab", sep="\t", header=T)
colnames(exones30all)<-c("number","source","target")
head(exones30all)
#all available nodes:
exones30compatible<-read.table("30Best_exons_compatible_edgelist.tab", sep="\t", header=T)
colnames(exones30compatible)<-c("number","source","target")
head(exones30compatible)

##################################################################
#nodes:
#labchip
n1<-unique( c(as.character(edgeListLC$source),
              as.character(edgeListLC$target)));length(n1)
#rnaseq
n2<-unique(c(as.character(edgeListLR$source),
             as.character(edgeListLR$target))); length(n2)
#ppi
n3<-unique(c(as.character(kdPPI$geneVastO),
             as.character(kdPPI$geneVastD))); length(n3)#367
#iclip
n4<-unique(c(as.character(iclip$source),
             as.character(iclip$target))); length(n4)#57
#exones all
n5<-unique(c(as.character(exones30all$source),
             as.character(exones30all$target))); length(n5)#244
#exones compatibles
n6<-unique(c(as.character(exones30compatible$source),
             as.character(exones30compatible$target))); length(n6)#211

auxN<-data.frame(labchip=as.character(classKD$VTb)%in%n1,
                 rnaseq=as.character(classKD$VTb)%in%n2,
                 ppi=as.character(classKD$VTb)%in%n3,
                 iclip=as.character(classKD$VTb)%in%n4,
                 exonesAll=as.character(classKD$VTb)%in%n5,
                 exonesComp=as.character(classKD$VTb)%in%n6,
                 order=as.numeric(classKD$ORDER))

#finalMatrixCOmparison
forPlot<-auxN
head(forPlot)
#colnames(forPlot)<-c("labchip","rnaseq")
rownames(forPlot)<- as.character(classKD$VTb)
forPlot$kds<-as.character(classKD$VTb)
forPlot<-forPlot[rowSums(forPlot[,1:4])>=2,]; dim(forPlot); #227
forPlot[forPlot==TRUE]<-1
forPlot[forPlot== FALSE]<-0
head(forPlot)
forPlotMelt = melt(forPlot, id.vars = c("kds","order"))
head(forPlotMelt)
#################################################################
getwd()
pdf("presence_ausence_nodes.pdf", width = 8.27,  height = 11.69)
ggplot(forPlotMelt, aes(x = variable, y = reorder(kds, -order), fill = factor(value)))  +
  theme_minimal() +
  geom_tile(position = "identity",  color = "grey") +
  coord_fixed(ratio = 0.5)+
  #coord_equal()+
  scale_fill_manual(values=c("white","grey"))+
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0,-0.5))+
  theme(
    axis.text.x = element_text(angle = 90, hjust=0, vjust=-5,size=5),
    axis.text.y = element_text(size=5),
    axis.ticks.length=unit(0,"cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border=element_blank(),
    legend.position="none",
    axis.line=element_blank()  )  
dev.off()
getwd()

##################################################################################
edgesC<-paste(edgeListLC$source, edgeListLC$target, sep="->"); length(edgesC)#547
edgesR<-paste(edgeListLR$source, edgeListLR$target, sep="->"); length(edgesR)#389
edgesPPI<-paste(kdPPI$geneVastO, kdPPI$geneVastD, sep="->"); length(edgesPPI)#466 
edgesICLIP<-paste(iclip$source, iclip$target, sep="->"); length(edgesICLIP)#82

edgesAllExons<-paste(exones30all$source, exones30all$target, sep="->")
edgesCompExons<-paste(exones30compatible$source, exones30compatible$target, sep="->")

n1<-edgesC
#rnaseq
n2<-edgesR
#ppi
n3<-edgesPPI
#iclip
n4<-edgesICLIP
#allExons
n5<-edgesAllExons
#comp
n6<-edgesCompExons
######################################
edgesALL<-unique(c(n1,n2,n3,n4,n5,n6))
length(edgesALL)#2249
length(edgesALL)
auxN<-data.frame(labchip=edgesALL%in%n1,
                 rnaseq=edgesALL%in%n2,
                 ppi=edgesALL%in%n3,
                 iclip=edgesALL%in%n4,
                 allExons=edgesALL%in%n5,
                 allCompExons=edgesALL%in%n6)
forPlot<-auxN;
head(forPlot)
rownames(forPlot)<- edgesALL
forPlot$edges<-edgesALL
length(edgesALL)
head(forPlot)
table(forPlot$edges)
forPlot<-forPlot[rowSums(forPlot[,1:6])>=2,]; dim(forPlot); #119
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
    axis.text.x = element_text(angle = 90, hjust=0, vjust=-5,size=5),
    axis.text.y = element_text(size=5),
    axis.ticks.length=unit(0,"cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border=element_blank(),
    legend.position="none",
    axis.line=element_blank()  )  
dev.off()
setwd("topology/")

