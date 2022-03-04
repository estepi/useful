library(limma)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape)
##################################################################
#SORT by complex
classKD <-
  read.table("classificationFiltered.txt",
             sep = "\t",
             header = T)
head(classKD)
dim(classKD)#312
classKD[classKD$VTb == "B2", ]

classKD[270:290, ]
allNodes <- as.character(classKD$VTb)
table(classKD$CLASS)
#################################################################
setwd("NumOfM100/topology/")
#################################################################
edgeListLC <-
  read.table("labchip048_edgelist.tab")
head(edgeListLC)
dim(edgeListLC)#547
centLC <-
  read.table(
    "labchip048_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
dim(centLC)#196 nodos
#################################################################
edgeListLR <-
  read.table("rnaseq_edgelist.tab")
head(edgeListLR)
dim(edgeListLR)#583
centR <-
  read.table(
    "rnaseq_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
head(centR)
dim(centR)#197 nodos
#################################################################
exones30all <-
  read.table("30best_exons_allKDs_edgelist.tab",
             sep = "\t",
             header = T)
colnames(exones30all) <- c("number", "source", "target")
centExAll <-
  read.table(
    "30best_exons_allKDs_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#Alt5
alt5all <- read.table("Alt5_op_edgelist.tab", sep = "\t", header = T)
colnames(alt5all) <- c("number", "source", "target")
head(alt5all)
cent5AltAll <-
  read.table(
    "Alt5_op_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#IR
IRall <- read.table("IR_op_edgelist.tab", sep = "\t", header = T)
colnames(IRall) <- c("number", "source", "target")
centIRAll <-
  read.table(
    "IR_op_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#Alt3
alt3all <- read.table("Alt3_op_edgelist.tab", sep = "\t", header = T)
colnames(alt3all) <- c("number", "source", "target")
cent3AltAll <-
  read.table(
    "Alt3_op_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#random30
random30 <- read.table("random_edgelist.tab", sep = "\t", header = T)
colnames(random30) <- c("number", "source", "target")
centRandomCR <-
  read.table(
    "random_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#nodes:
#labchip
n1 <- unique(c(
  as.character(edgeListLC$source),
  as.character(edgeListLC$target)
))
length(n1)
#rnaseq
n2 <- unique(c(
  as.character(edgeListLR$source),
  as.character(edgeListLR$target)
))
length(n2)
#exones all
n5 <- unique(c(
  as.character(exones30all$source),
  as.character(exones30all$target)
))
length(n5)#244
#alt5 all
n7 <- unique(c(as.character(alt5all$source),
               as.character(alt5all$target)))
length(n7)#75
#IR all
n9 <- unique(c(as.character(IRall$source),
               as.character(IRall$target)))
length(n9)#46
#3alt all
n11 <- unique(c(as.character(alt3all$source),
                as.character(alt3all$target)))
length(n11)#126
#random30
n13 <- unique(c(
  as.character(random30$source),
  as.character(random30$target)
))
length(n13)#158

#########################################################################
auxN <- data.frame(
  labchip = as.character(classKD$VTb) %in% n1,
  rnaseq = as.character(classKD$VTb) %in% n2,
  exonesAll = as.character(classKD$VTb) %in% n5,
  alt5All = as.character(classKD$VTb) %in% n7,
  IRAll = as.character(classKD$VTb) %in% n9,
  alt3all = as.character(classKD$VTb) %in% n11,
  random30 = as.character(classKD$VTb) %in% n13,
  order = as.numeric(classKD$ORDER)
)
#finalMatrixCOmparison
forPlot <- auxN
#colnames(forPlot)<-c("labchip","rnaseq")
rownames(forPlot) <- as.character(classKD$VTb)
forPlot$kds <- as.character(classKD$VTb)
dim(forPlot)#312*10
forPlot[forPlot == TRUE] <- 1
forPlot[forPlot == FALSE] <- 0
forPlotMelt = melt(forPlot, id.vars = c("kds", "order"))
#################################################################
pdf(
  "presence_ausence_nodes_AS_ex_alt5_ir_alt3_random_MD.pdf",
  width = 8.27,
  height = 11.69
)
ggplot(forPlotMelt, aes(
  x = variable,
  y = reorder(kds,-order),
  fill = factor(value)
))  +
  theme_minimal() +
  geom_tile(position = "identity",  color = "grey") +
  coord_fixed(ratio = 0.7) +
  #coord_equal()+
  scale_fill_manual(values = c("white", "grey")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0, -0.5)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = -5,
      size = 4
    ),
    axis.text.y = element_text(size = 3),
    axis.ticks.length = unit(0, "cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    axis.line = element_blank()
  )
dev.off()
#stats about centrality rankking:
auxCC <- auxN
rownames(auxCC) <- as.character(classKD$VTb)
head(auxCC)

ii1 <- match(rownames(auxCC), rownames(centLC))#labchip
ii2 <- match(rownames(auxCC), rownames(centR))#rnaseq
ii4 <- match(rownames(auxCC), rownames(centExAll))#AllExons
#####
ii6 <- match(rownames(auxCC), rownames(cent5AltAll))#5altAll
ii8 <- match(rownames(auxCC), rownames(centIRAll)) #IR All
#####
ii10 <- match(rownames(auxCC), rownames(cent3AltAll)) #3AltAll
ii11 <- match(rownames(auxCC), rownames(centRandomCR)) #30random

auxCC$labchip <- centLC$DG[ii1]
auxCC$rnaseq <- centR$DG[ii2]
auxCC$exonesAll <- centExAll$DG[ii4]
auxCC$alt5All <- cent5AltAll$DG[ii6]
auxCC$IRAll <- centIRAll$DG[ii8]
auxCC$alt3all <- cent3AltAll$DG[ii10]
auxCC$random30 <- centRandomCR$DG[ii11]
write.table(auxCC, file = "ComparisonDegree_MD.tab", sep = "\t")
##################################################################
auxCC <- auxN
rownames(auxCC) <- as.character(classKD$VTb)
auxCC$labchip <- centLC$RANK[ii1]
auxCC$rnaseq <- centR$RANK[ii2]
auxCC$exonesAll <- centExAll$RANK[ii4]
auxCC$alt5All <- cent5AltAll$RANK[ii6]
auxCC$IRAll <- centIRAll$RANK[ii8]
auxCC$alt3all <- cent3AltAll$RANK[ii10]
auxCC$random30 <- centRandomCR$RANK[ii11]
write.table(auxCC, file = "ComparisonRANK.tab", sep = "\t")

##################################################################
#lets try to plot in colors absense / presence according degree:
forPlot <- auxN
#colnames(forPlot)<-c("labchip","rnaseq")
rownames(forPlot) <- as.character(classKD$VTb)
forPlot$kds <- as.character(classKD$VTb)
forPlot <- forPlot[rowSums(forPlot[, 1:4]) >= 2, ]
#227
forPlot[forPlot == TRUE] <- 1
forPlot[forPlot == FALSE] <- 0
head(forPlot)
forPlotMelt = melt(forPlot, id.vars = c("kds", "order"))
#################################################################
pdf("presence_ausence_nodes_AS.pdf",
    width = 8.27,
    height = 11.69)
ggplot(forPlotMelt, aes(
  x = variable,
  y = reorder(kds,-order),
  fill = factor(value)
))  +
  theme_minimal() +
  geom_tile(position = "identity",  color = "grey") +
  coord_fixed(ratio = 0.5) +
  #coord_equal()+
  scale_fill_manual(values = c("white", "grey")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0, -0.5)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = -5,
      size = 5
    ),
    axis.text.y = element_text(size = 4),
    axis.ticks.length = unit(0, "cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    axis.line = element_blank()
  )
dev.off()
##################################################################################
edgesC <-
  paste(edgeListLC$source, edgeListLC$target, sep = "->")
length(edgesC)#547
edgesCR <-
  paste(edgeListLC$target, edgeListLC$source, sep = "->")
length(edgesCR)#547
##################################################################################
edgesR <-
  paste(edgeListLR$source, edgeListLR$target, sep = "->")
length(edgesR)#389
edgesRR <-
  paste(edgeListLR$target, edgeListLR$source, sep = "->")
length(edgesR)#389
##################################################################################
edgesAllExons <-
  paste(exones30all$source, exones30all$target, sep = "->")
length(edgesAllExons)#652
edgesAllExonsR <-
  paste(exones30all$target, exones30all$source, sep = "->")
length(edgesAllExonsR)#652
##################################################################################
edgesCompExons <-
  paste(exones30compatible$source, exones30compatible$target, sep = "->")
length(edgesCompExons)#552
edgesCompExonsR <-
  paste(exones30compatible$target, exones30compatible$source, sep = "->")
length(edgesCompExonsR)#552

n1 <- c(edgesC, edgesCR)
length(n1) / 2
#rnaseq
n2 <- c(edgesR, edgesRR)
length(n2) / 2
#allExons
n5 <- c(edgesAllExons, edgesAllExonsR)
length(n5) / 2
#comp
n6 <- c(edgesCompExons, edgesCompExonsR)
length(n6) / 2
######################################
edgesALL <- unique(c(n1, n2, n5, n6))
length(edgesALL)#3078
auxN <- data.frame(
  labchip = edgesALL %in% n1,
  rnaseq = edgesALL %in% n2,
  allExons = edgesALL %in% n5,
  allCompExons = edgesALL %in% n6
)
forPlot <- auxN
rownames(forPlot) <- edgesALL
forPlot$edges <- edgesALL
length(edgesALL)
head(forPlot)
forPlot <- forPlot[rowSums(forPlot[, 1:4]) >= 3,]
forPlot[forPlot == TRUE] <- 1
forPlot[forPlot == FALSE] <- 0
forPlotMelt = melt(forPlot, id.vars = "edges")
forPlotMelt$NumOfEdge <- rep(1:nrow(forPlot), 4)
#################################################################
pdf("presence_ausence_edges_AS.pdf",
    width = 8.27,
    height = 11.69)
ggplot(forPlotMelt, aes(x = variable, y = edges, fill = factor(value)))  +
  theme_minimal() +
  geom_tile(position = "identity",  color = "grey") +
  coord_fixed(ratio = 0.5) +
  #coord_equal()+
  scale_fill_manual(values = c("white", "grey")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0,-0.5)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = -5,
      size = 5
    ),
    axis.text.y = element_text(size = 5),
    axis.ticks.length = unit(0, "cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    axis.line = element_blank()
  )
dev.off()
