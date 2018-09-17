#finalMatrixCOmparison
library(data.table)
library(edgeR)
library(ggplot2)
library(stringr)
install.packages("WGCNA")
library(ape)
install.packages("ape")

setwd("/home/emancini/Dropbox (CRG ADV)/Personal_Gosia/Shared/GosiaAndEstefi/GE_331/")
#using default parameter
GE <- read.csv("geneCounts_331.tab", 
                      header=T, 
                      sep="\t", 
                      stringsAsFactors = F);dim(GE)
GEn<-as.matrix(GE[,-1]); head(GEn); dim(GEn)
##############################################
dim(GEn)#52465
#onlyControls:
finalTable<-read.table("all_lrtTable_331_union_names.txt", header=T, row.names = 1); head(finalTable)
finalTable[1:5,1:5]
dim(ff)
####################################################
ff<-finalTable
colnames(ff)[1]<-"symbol"
table1<-ff; head(ff)
hist(table1$logCPM)
##################################
kdsNames<-sub("logFC.group","",colnames(GE)); length(kdsNames)
kdsNamesI<-match(kdsNames,ff$symbol);kdsNamesI
kdsNamesII<-kdsNamesI[!is.na(kdsNamesI)];kdsNamesII
length(kdsNamesII)
color<-rep("grey", nrow(ff))
color[kdsNamesII]<-"tomato"
table(color)
shape<-rep(20, nrow(ff))
shape[kdsNamesII]<-19
length(GE[kdsNamesII,1])
length(kdsNamesII)
lbls<-as.character(ff$symbol[kdsNamesII])
length(lbls)
rm(kdsLabesl)
####################################################################
getwd()
dir.create("onlyKds")
dir.create("onlyKds2")
setwd("..")
dir.create("all")
setwd("all")

colnames(ff)[2:323]
for (i in 2:323)  
for (i in 2:323)  
{
  file<-paste(sub("logFC.group","",colnames(ff)[i]) ,"all.png", sep=".")
  #file<-paste(sub("logFC.group","",colnames(ff)[i]) , "plot.pdf", sep=".")
  print (file)
  png(file,width=1500,height=1000,res=150)
  #pdf(file, width =11.69, height =  8.27)
  mainG<-sub("logFC.group","",colnames(ff[i]))
  plot(ff$logCPM, ff[,i],
       main=mainG, 
       pch=shape, 
       xlim=c(0,15), 
       ylim = c(-15,15), 
       col=color)
  abline(h=c(-1.5,0,1.5), 
         v=c(0,5))
  text(ff$logCPM[kdsNamesII], 
       ff[kdsNamesII,i],
       labels=lbls,
       cex= 0.4, 
       offset =2,
       family="sans"
  ) 
  dev.off()  
}

##########################################################
head(GEn)
load("lrt_GE.RData")
tt1000<-topTags(lrt, n=1000)
tt5000<-topTags(lrt,n=5000)
head(ff)
lrt$dispersion
ff$dispersion<-lrt$dispersion
png("dispersionALL.png")
plot(ff$dispersion, ff[,4])

abline(v=0, h=0)
dev.off()
fforder <- ff[order(ff$logCPM),]
########################################################
getwd()
setwd("dispersion/")
dir.create("dispersion")
for (i in 2:323)  
{
  file<-paste(sub("logFC.group","",colnames(ff)[i]) ,"dispersion.png", sep=".")
  #file<-paste(sub("logFC.group","",colnames(ff)[i]) , "plot.pdf", sep=".")
  print (file)
  png(file,width=1500,height=1000,res=150)
  #pdf(file, width =11.69, height =  8.27)
  mainG<-sub("logFC.group","",colnames(ff[i]))
  plot(ff$dispersion,ff[,i],
       main=mainG, 
       pch=shape, 
       xlim=c(0,1), 
       ylim = c(-15,15), 
       col=color)
  abline(h=c(-1.5,0,1.5), 
         v=c(0,0.5))
  text(ff$dispersion[kdsNamesII], 
       ff[kdsNamesII,i],
       labels=lbls,
       cex= 0.4, 
       offset =2,
       family="sans"
  ) 
  dev.off()  
}
