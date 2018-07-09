#finalMatrixCOmparison
library(data.table)
library(stringr)

INCFile<-"/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/default/INCLUSION_LEVELS_FULL-Hsa331-hg19.tab"
OUTFile<-"/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/default/test_NA_NewN3.tab"
source("replaceNA.R")
replaceNA(INCFile, OUTFile)
#########################################
INCFile<-"test.tab"
OUTFile1<-"test_NA_NewN3.tab"
OUTFile2<-"test_NA_NewN3_NAI.tab"
source("replaceNAI.R")
replaceNAI(INCFile, OUTFile1,OUTFile2)
#########################################
