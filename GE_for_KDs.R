#gene expression using edgeR
library(edgeR)
setwd("")
table <-
  read.table("geneCounts_319.tab",
             header = T,
             row.names = 1)
dim(table)#52465
##############
colnames(table) <- sub("_geneCounts", "", colnames(table))
condition <- colnames(table)
#rename conditions:
combined <- grep("con$", condition)
condition[combined]
combined
condition[143] <- "IK"
condition[204] <- "PRPF8"
condition[210] <- "RBM17"
condition[232] <- "SF3B1"
condition[250] <- "SMU1"
#############################################################
brep <- grep("_b", condition)
condition[32] <- "CCDC12"
condition[36] <- "CDC5L"
condition[53] <- "CWC22"
condition[114] <- "HFM1"
condition[155] <- "LENG1"
condition[274] <- "SRPK2"
condition[310] <- "XAB2"
#############################################################
controls <- grep("^AA", condition)
condition[controls] <- "AA"
length(controls)#8
table(condition)[table(condition) > 1]
#############################################################
#nonzero:
nonzero <- table[rowSums(table[, -1]) > 0, -1]
dim(nonzero)#44275
length(condition)
y <- DGEList(counts = nonzero, group = condition[-1])
#filter: 100% smaples more than 1cpm
keep <- rowSums(cpm(y) > 1) == 319
length(which(keep))#6947
write.table(rownames(y$counts),
            file = "filteredGenes_1cpm_100.txt",
            quote = F,
            sep = "\t")
y <- y[keep, , keep.lib.sizes = FALSE]
dim(y)
design <- model.matrix( ~ group, y$samples)

write.table(design, file = "matrix_design_319.txt", sep = "\t")
y <- calcNormFactors(y)
y <- estimateDisp(y, design, verbose = TRUE)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
getwd()
save(y, file = "yHF.RData")
save(fit, file = "fitHF.RData")
lrt <- glmLRT(fit, coef = 2:ncol(fit$coefficients))
save(lrt, file = "lrtHF.RData")
###################################################
write.table(lrt$table,
            file = "lrt_table_hg19_HF.tab",
            sep = "\t",
            col.names = NA)
dim(y)
