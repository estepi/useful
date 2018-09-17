dPSI<-read.csv("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/ensembl/forGO/df1_dPSI15t.txt", sep="\t", header=T)
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
rownames(dPSI_Ens)
dPSI_Ens[1:5,1:5]

samples<-colnames(dPSI_Ens); length(samples)
samples<-samples[-c(1:7)]; length(samples)
samples<-samples[-306];length(samples)
ii<-match(samples, small$GENE)
length(samples[is.na(ii)])#9 #not expressed

#########################################################################
genesLRT<-read.table("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/GE_hg19_noBadSamples/lrt_table_hg19_HF.tab", header=T, row.names = 1)
genesLRT<-read.table("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/GE_hg19_noBadSamples/lrt_table_hg19.tab", header=T, row.names = 1)
kdsNames<-sub("logFC.group","",colnames(genesLRT));
kdsNames<-kdsNames[-c(306:308)]; length(kdsNames)
vectorLFC<-vector(length = length(kdsNames))
names(vectorLFC)<-kdsNames
dd<-match(names(vectorLFC), genesCounts$X); length(which(is.na(dd)))
rm(df)
df<-data.frame(FC=vectorLFC, symbol=genesCounts$X[dd], ens=row.names(genesCounts)[dd]);head(df)
length(which(is.na(match(df$ens, rownames(genesLRT)))))#59 #19
##############################################################################
dPSI[1:5,1:8]
colnames(dPSI)
NumOfEvents<-colSums(dPSI[,9:ncol(dPSI)])
kdsAS<-match(names(NumOfEvents), df$symbol)
df$numOfE<-NumOfEvents[kdsAS]
head(df)
#chequeamos en GE si son genes poco expresados:
genesCounts<-read.table("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/GE_hg19_noBadSamples/geneCounts_319.tab", header=T, row.names = 1)
cc<-match(kdsNames[is.na(kdsNamesI)], genesCounts$X)
#los genes estan en la tabla genecounts originales --> poca expresion
#########################################################################
#double match:
colnames(genesLRT)<-sub("logFC.group","",colnames(genesLRT));
sI<-match(df$symbol, colnames(genesLRT)); length(sI)
gI<-match(df$ens, rownames(genesLRT)); length(gI)
names(gI)<-df$symbol
names(sI)<-df$symbol
df[1:4,1:4]

for (i in 1:length(gI))
{
  
  print(paste("gene(row)",names(gI[i]) , sep="."   ))   
  print(paste("sample",names(sI[i]), sep=":" ))
  df$FC[i]<-genesLRT[gI[i],sI[i]]
  
}

df[1:4,1:4]
setwd("/home/emancini/Dropbox (CRG ADV)/Personal_Gosia/Shared/GosiaAndEstefi/GE_kds_plots/")
write.table(df,"logFC_KDs_2.tab")
df$FC[is.na(df$FC)]<-0
plot(df$numOfE, df$FC)
abline(h=c(-1,0,1))
####################################
#plot the most affecting ones:
df50<-df[df$numOfE>50,]
plot(df50$numOfE, df50$FC)
abline(h=c(-1,0,1))
df100<-df[df$numOfE>100,]
plot(df100$numOfE, df100$FC)

png("Allevents.png")
plot(df$numOfE, df$FC, 
     main= "Num of Events vs KD efficiency",
     xlab= "Number of events",
     ylab= "KD efficiency (logFC)",
     col= "blue", pch = 19, cex = 0.5, lty = "solid", lwd = 2)
text(df$numOfE, df$FC, labels=df$symbol, cex= 0.7)
abline(h=c(-1,0,1))
dev.off()

png("50events.png")
plot(df50$numOfE, df50$FC, 
     main= "Num of Events vs KD efficiency",
     xlab= "Number of events",
     ylab= "KD efficiency (logFC)",
     col= "blue", pch = 19, cex = 0.5, lty = "solid", lwd = 2)
text(df50$numOfE, df50$FC, labels=df50$symbol, cex= 0.7)
abline(h=c(-1,0,1))
dev.off()

png("100events.png")
plot(df100$numOfE, df100$FC, 
     main= "Num of Events vs KD efficiency",
     xlab= "Number of events",
     ylab= "KD efficiency (logFC)",
     col= "blue", pch = 19, cex = 0.5, lty = "solid", lwd = 2)
text(df100$numOfE, df100$FC, labels=df100$symbol, cex= 0.7)
abline(h=c(-1,0,1))
dev.off()

