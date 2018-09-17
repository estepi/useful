library(limma)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape) 
##################################################################
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/NumOfM100/")
edgeListLC<-read.table("edgelist_048_labchip.tab");head(edgeListLC); dim(edgeListLC)
labchipG<-graph_from_data_frame(edgeListLC, directed = F)
length(E(labchipG))#547
length(V(labchipG))#196
##################################################################
edgeListLR<-read.table("edgelist_055_rnaseq.tab");head(edgeListLR); dim(edgeListLR)
networkR<-graph_from_data_frame(edgeListLR, directed = F)
length(E(networkR))#389
length(V(networkR))#197
##################################################################
resuting<-intersection(networkR, labchipG,keep.all.vertices = F)
length(E(resuting))#59
length(V(resuting))#141
#links:
l11<-paste(edgeListLC$source, edgeListLC$target, sep="-")
l12<-paste(edgeListLC$target, edgeListLC$source, sep="-")
l1<-c(l11, l12); length(l1); length(unique(l1))

l21<-paste(edgeListLR$source, edgeListLR$target, sep="-")
l22<-paste(edgeListLR$target, edgeListLR$source, sep="-")
l2<-c(l21,l22); length(unique(l21))

allLinks<-union(l11, l21)
aux<-data.frame(allLinks%in%l11, allLinks%in%l21)
vennDiagram(aux)
##################################################################
#nodes:
n1<-unique( c(as.character(edgeListLC$source),
              as.character(edgeListLC$target)));length(n1)
n2<-unique(c(as.character(edgeListLR$source),
             as.character(edgeListLR$target))); length(n2)
allNodes<-union( n1, n2);length(allNodes)#138
auxN<-data.frame(allNodes%in%n1,
                 allNodes%in%n2)
vennDiagram(auxN)
#finalMatrixCOmparison
forPlot<-auxN
head(forPlot)
colnames(forPlot)<-c("labchip","rnaseq")
rownames(forPlot)<- allNodes
forPlot$dataset<-allNodes
forPlot[forPlot==TRUE]<-1
forPlot[forPlot== FALSE]<-0
head(forPlot)
forPlot$value<-NULL
head(forPlot)
forPlotMelt = melt(forPlot, id.vars = "dataset") 
head(forPlotMelt)
#################################################################
#SORT by complex
classKD<-read.table("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/classificationFiltered.txt", sep="\t", header = T)
head(classKD)
table(classKD$CLASS)
#################################################################
pdf("presence_ausence.pdf", width = 8.27,  height = 11.69)
ggplot(forPlotMelt, aes(x = variable, y = dataset, fill = factor(value)))  +
  theme_minimal() +
  geom_tile(position = "identity",  color = "grey") +
  coord_fixed(ratio = 0.8)+
  #coord_equal()+
  scale_fill_manual(values=c("white","grey"))+
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels=as.character(allNodes),  expand = c(0,-0.5))+
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

