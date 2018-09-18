#finalMatrixCOmparison
library(edgeR)
setwd("/home/emancini/Dropbox (CRG ADV)/Personal_Gosia/Shared/GosiaAndEstefi/luisaData/")
#using default parameter
GE <- read.csv("geneCounts_drugs.tab", 
                      header=T, 
                      sep="\t", 
                      stringsAsFactors = F);dim(GE)
head(GE); dim(GE)
##############################################
group<-rep(c("DMSO","SSA", "SUDC", "SUDKB"), each=2)
y <- DGEList(counts=GE, group = group)
design=model.matrix(~group, data=y$samples)
write.table(design, file="designGE.txt")
#filter
cpms<-cpm(y, prior.count = 0.01, log = F); head(cpms)
keep<-which(rowSums(cpms >1) >= 6); dim(y); #13481
y<-y[keep,]
dim(y)
y<-calcNormFactors(y)
y<-estimateGLMCommonDisp(y, design, verbose=TRUE)
fit<-glmFit(y,design)
lrt<-glmLRT(fit, coef=2:4)
fdr<-p.adjust(lrt$table$PValue, method="BH")
tt<-cbind(lrt$table, fdr)
plotMDS(y)
plotSmear(y)
#change NAMES:
names<-read.table("Hsa.ID.names.txt", header = F, sep="\t")
head(names)
ii<-match(rownames(tt), names$V1)
sum(is.na(ii))
final<-cbind(names$V2[ii], tt)
write.table(final, "Drugs_vs_DMSO.tab", sep="\t", col.names = NA)
#################### SSA
lrt<-glmLRT(fit, coef=2)
fdr<-p.adjust(lrt$table$PValue, method="BH")
tt<-cbind(lrt$table, fdr)
#change NAMES:
ii<-match(rownames(tt), names$V1)
sum(is.na(ii))
final<-cbind(names$V2[ii], tt)
write.table(final, "SSA_vs_DMSO.tab", sep="\t", col.names = NA)
#################### SUDC
lrt<-glmLRT(fit, coef=3)
fdr<-p.adjust(lrt$table$PValue, method="BH")
tt<-cbind(lrt$table, fdr)
#change NAMES:
ii<-match(rownames(tt), names$V1)
sum(is.na(ii))
final<-cbind(names$V2[ii], tt)
write.table(final, "SUDC_vs_DMSO.tab", sep="\t", col.names = NA)
#################### SUDKB
lrt<-glmLRT(fit, coef=4)
fdr<-p.adjust(lrt$table$PValue, method="BH")
tt<-cbind(lrt$table, fdr)
#change NAMES:
ii<-match(rownames(tt), names$V1)
sum(is.na(ii))
final<-cbind(names$V2[ii], tt)
write.table(final, "SUDKB_vs_DMSO.tab", sep="\t", col.names = NA)
#################
t1<-read.table("Drugs_vs_DMSO.tab", header=TRUE)
t2<-final

head(t1)
head(t2)

ii<-match(t1$X, rownames(t2))
forPlot<-cbind(t1$logFC.groupSUDKB, t1$PValue, t2$logFC[ii], t2$PValue[ii])
head(forPlot)
cor(forPlot[,1], forPlot[,3])
plot(log(forPlot[,2]),
     log(forPlot[,4]),
     xlab="All vs DMSO", ylab="SUDK vs DMSO", xlim=c(-400,0), ylim=c(-400,0)      )
abline(0,1)

plot(forPlot[,2],
     forPlot[,4],
     xlab="All vs DMSO", ylab="SUDK vs DMSO", 
     xlim=c(0,0.01), ylim=c(0,0.01)      
     )
abline(0,1)
cor(log(forPlot[,2]),
     log(forPlot[,4]))



####analizar con vast-tools para rna all:
getwd()
all<-read.table("cRPKM_AND_COUNTS-Hsa37.tab", header = T, row.names = 1)
cn<-colnames(all)
onlyCounts<-all[,grep("Counts", cn)]
head(onlyCounts)
write(colnames(onlyCounts), file="OnlyCountsSamples.txt")
#check controls:
controls<-grep("DMSO", colnames(onlyCounts))
ssa<-grep("SSA",  colnames(onlyCounts))
length(ssa)
sud<-grep("Sud",  colnames(onlyCounts))
length(sud)
subset<-onlyCounts[,c(controls, ssa, sud)]
head(subset)
rownames(subset)
subset[is.na(subset)]<-0
colnames(subset)<-c("DMSO_1","DMSO_2","DMSO_IP_1","DMSO_IP_2",
                    "SSA_1","SSA_2","SSA_IP_1","SSA_IP_2",
                    "SUDK_1","SUDK_2",
                    "SUDC1_1","SUDC1_2", 
                    "SUDC1_IP_1","SUDC1_IP_2",
                     "SUDK_IP_1","SUDK_IP_2")
#########################################################
group<- c("DMSO","DMSO","DMSO_IP","DMSO_IP",
    "SSA","SSA","SSA_IP","SSA_IP",
    "SUDK","SUDK",
    "SUDC1","SUDC1", 
    "SUDC1_IP","SUDC1_IP",
    "SUDK_IP","SUDK_IP")
dim(subset)#19847
nonzero<-as.matrix(subset[rowSums(subset)>0,]); dim(nonzero)#16834
y <- DGEList(counts=nonzero, group = group)
plotMDS(y)
####divide en IP vs ALL
#group 1: IP
ip<-subset[,grep("IP", colnames(subset))]; dim(ip)
groupip<-c("DMSO_IP", "DMSO_IP",
           "SSA_IP",  "SSA_IP",
           "SUDC1_IP","SUDC1_IP",
           "SUDK_IP", "SUDK_IP")
nonzero<-ip[rowSums(ip)>0,]; dim(nonzero)#15240
class(nonzero)
rownames(nonzero)
y <- DGEList(counts=nonzero, group = groupip)
plotMDS(y)
design=model.matrix(~groupip, data=y$samples)
#filter
cpms<-cpm(y, prior.count = 0.01, log = F); head(cpms)
keep<-which(rowSums(cpms >1) >= 6)
y<-y[keep,]
dim(y)#12138
y<-calcNormFactors(y)
y<-estimateGLMCommonDisp(y, design, verbose=TRUE)
fit<-glmFit(y,design)
lrt<-glmLRT(fit, coef=2:4)
fdr<-p.adjust(lrt$table$PValue, method="BH")
tt<-cbind(lrt$table, fdr)
order_IP<-tt[order(-tt$logFC.groupipSSA_IP),]
head(order_IP)
write.table(tt, "all_vs_dmso_IP", sep="\t", col.names = NA)
top500IP_SSA<-order_IP[1:500,]
#########################################################
##group 2: all
total<-subset[,-grep("IP", colnames(subset))]; dim(ip)
head(total)
grouptotal<-c("DMSO", "DMSO",
           "SSA",  "SSA",
           "SUDK", "SUDK",
           "SUDC1","SUDC1")
nonzero<-as.matrix(total[rowSums(total)>0,]); dim(nonzero)#16686
class(nonzero)
y <- DGEList(counts=nonzero, group = grouptotal); plotMDS(y)
design=model.matrix(~grouptotal, data=y$samples)
cpms<-cpm(y, prior.count = 0.01, log = F); head(cpms)
keep<-which(rowSums(cpms >1) >= 6)
y<-y[keep,]; dim(y)#11125
y<-calcNormFactors(y)
y<-estimateGLMCommonDisp(y, design, verbose=TRUE)
fit<-glmFit(y,design)
lrt<-glmLRT(fit, coef=2:4)
fdr<-p.adjust(lrt$table$PValue, method="BH")
ttotal<-cbind(lrt$table, fdr)
write.table(ttotal, "all_vs_dmso_total", sep="\t", col.names = NA)
head(ttotal)
order_total<-ttotal[order(-ttotal$logFC.grouptotalSSA),] 
head(order_total)
top500total_SSA<-order_total[1:500,]

ii<-match(rownames(top500IP_SSA),
          rownames(top500total_SSA))
#ssa
plot(top500IP_SSA$logFC.groupipSSA_IP[ which(!is.na(ii))],
       top500total_SSA$logFC.grouptotalSSA[ii[!is.na(ii)]] )
abline(0,1)
cor(top500IP_SSA$logFC.groupipSSA_IP[ which(!is.na(ii))],
     top500total_SSA$logFC.grouptotalSSA[ii[!is.na(ii)]] )#081
#sud
plot(top500IP_SSA$logFC.groupipSUDC1_IP[ which(!is.na(ii))],
     top500total_SSA$logFC.grouptotalSUDC1[ii[!is.na(ii)]] )
abline(0,1)
cor(top500IP_SSA$logFC.groupipSUDC1_IP[ which(!is.na(ii))],
     top500total_SSA$logFC.grouptotalSUDC1[ii[!is.na(ii)]] )#074

#sudk
plot(top500IP_SSA$logFC.groupipSUDK_IP[ which(!is.na(ii))],
     top500total_SSA$logFC.grouptotalSUDK[ii[!is.na(ii)]] )
abline(0,1)

cor(top500IP_SSA$logFC.groupipSUDK_IP[ which(!is.na(ii))],
     top500total_SSA$logFC.grouptotalSUDK[ii[!is.na(ii)]] )#069

###add names:
tablaIP<-read.table("all_vs_dmso_IP.txt", header=T, row.names = 1)
head(tablaIP)
head(all)
ii<-match(rownames(tablaIP), rownames(all))
tablaIPfinal<-cbind(all$NAME[ii], tablaIP)
colnames(tablaIPfinal)[1]<-"NAME"
write.table(tablaIPfinal, file="all_vs_dmso_IP.tab", quote=F, sep="\t")
#################################
tablaTotal<-read.table("all_vs_dmso_total.txt", header=T, row.names = 1)
head(tablaTotal)
head(all)
ii<-match(rownames(tablaTotal), rownames(all))
tablaTotalfinal<-cbind(all$NAME[ii], tablaTotal)
colnames(tablaTotalfinal)[1]<-"NAME"
head(tablaTotalfinal)
write.table(tablaTotalfinal, file="all_vs_dmso_total.tab", quote=F, sep="\t")
