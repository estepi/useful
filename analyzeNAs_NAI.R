#finalMatrixCOmparison
setwd("~/useful/")
library(data.table)
library(stringr)
INCFile<-"/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc/HCa_NUMB9.tab"
OUTFile<-"/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc/HCa_NUMB9.tab_notNAs.tab"
source("replaceNA.R")
replaceNA(INCFile, OUTFile)
#########################################
INCFile<-"/home/emancini/Dropbox (CRG)/Personal_Estefania/Arun/vtools/INCLUSION_LEVELS_FULL-Merge.tab"
OUTFile1<-"/home/emancini/Dropbox (CRG)/Personal_Estefania/Arun/vtools/INCLUSION_LEVELS_FULL-Merge_NA.tab"
OUTFile2<-"/home/emancini/Dropbox (CRG)/Personal_Estefania/Arun/vtools/INCLUSION_LEVELS_FULL-Merge_NA_NAI.tab"
source("replaceNAI_corrected.R")
replaceNAI(INCFile, OUTFile1,OUTFile2)
#########################################
infile<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc/HCa_NUMB9.tab_notNAs.tab", 
                     header=T, sep="\t",na.strings = "NA" , row.names = 1)
head(infile)
dim(infile)
table(infile$COMPLEX)
nasC<-which(is.na(infile$COMPLEX))
write.table(infile[nasC,1:6], "~/Dropbox (CRG)/Personal_Estefania/rferrari/COMPLEX_NAs_HS2.tab", sep="\t", col.names = NA, quote = F)
rownames(infile)  <-paste(infile$EVENT, 1:nrow(infile), sep="_")
rownames(infile)[1:10]

