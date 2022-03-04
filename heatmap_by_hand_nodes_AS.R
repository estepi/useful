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
head(classKD)
allNodes <- as.character(classKD$VTb)
table(classKD$CLASS)
#################################################################
setwd("Network/NumOfM100/topology/")
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
head(centLC)
#################################################################
edgeListLR <-
  read.table("RNA-Seq_compatible_edgelist.tab")
head(edgeListLR)
dim(edgeListLR)#389
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
head(exones30all)
centExAll <-
  read.table(
    "30best_exons_allKDs_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
head(centExAll)
#all available nodes:
exones30compatible <-
  read.table("30Best_exons_compatible_edgelist.tab",
             sep = "\t",
             header = T)
colnames(exones30compatible) <- c("number", "source", "target")
head(exones30compatible)
centExComp <-
  read.table(
    "30Best_exons_compatible_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
head(centExComp)
##################################################################
#Alt5
alt5all <-
  read.table("best30Alt5_all_edgelist.tab",
             sep = "\t",
             header = T)
colnames(alt5all) <- c("number", "source", "target")
head(alt5all)
cent5AltAll <-
  read.table(
    "best30Alt5_all_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

#all available nodes:
alt5compatible <-
  read.table("best30Alt5_comp_edgelist.tab",
             sep = "\t",
             header = T)
colnames(alt5compatible) <- c("number", "source", "target")
head(alt5compatible)
cent5AlComp <-
  read.table(
    "best30Alt5_comp_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#IR
IRall <- read.table("best30IR_all_edgelist.tab",
                    sep = "\t",
                    header = T)
colnames(IRall) <- c("number", "source", "target")
head(IRall)
centIRAll <-
  read.table(
    "best30IR_all_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

#all available nodes:
IRcompatible <-
  read.table("best30IR_comp_edgelist.tab",
             sep = "\t",
             header = T)
colnames(IRcompatible) <- c("number", "source", "target")
head(IRcompatible)
centIRComp <-
  read.table(
    "best30IR_comp_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

##################################################################
#Alt3
alt3all <-
  read.table("best30Alt3_all_edgelist.tab",
             sep = "\t",
             header = T)
colnames(alt3all) <- c("number", "source", "target")
head(alt3all)
cent3AltAll <-
  read.table(
    "best30Alt3_all_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

#all available nodes:
alt3compatible <-
  read.table("best30Alt3_comp_edgelist.tab",
             sep = "\t",
             header = T)
colnames(alt3compatible) <- c("number", "source", "target")
head(alt3compatible)
cent3AlComp <-
  read.table(
    "best30Alt3_comp_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#random30
random30 <-
  read.table("best30random_edgelist.tab",
             sep = "\t",
             header = T)
colnames(random30) <- c("number", "source", "target")
head(random30)
centRandomCR <-
  read.table(
    "best30random_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

#all available nodes:
random30to <-
  read.table("best30random_tog_edgelist.tab",
             sep = "\t",
             header = T)
colnames(random30to) <- c("number", "source", "target")
head(random30to)
random30togCR <-
  read.table(
    "best30random_tog_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )
##################################################################
#120
best120 <- read.table("best120_edgelist.tab", sep = "\t", header = T)
colnames(best120) <- c("number", "source", "target")
head(best120)
best120CR <-
  read.table(
    "best120_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

#all available nodes:
best120to <-
  read.table("best120together_edgelist.tab",
             sep = "\t",
             header = T)
colnames(best120to) <- c("number", "source", "target")
head(best120to)
best120togCR <-
  -read.table(
    "best120together_CentralityRanking.tab",
    sep = "\t",
    header = T,
    row.names = 1
  )

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
#exones compatibles
n6 <- unique(c(
  as.character(exones30compatible$source),
  as.character(exones30compatible$target)
))
length(n6)#211

#alt5 all
n7 <- unique(c(as.character(alt5all$source),
               as.character(alt5all$target)))
length(n7)#301
#alt5 comp
n8 <- unique(c(
  as.character(alt5compatible$source),
  as.character(alt5compatible$target)
))
length(n8)#301
#IR all
n9 <- unique(c(as.character(IRall$source),
               as.character(IRall$target)))
length(n9)#291
#Ir comp
n10 <- unique(c(
  as.character(IRcompatible$source),
  as.character(IRcompatible$target)
))
length(n10)#291

#3alt all
n11 <- unique(c(as.character(alt3all$source),
                as.character(alt3all$target)))
length(n11)#293

n12 <- unique(c(
  as.character(alt3compatible$source),
  as.character(alt3compatible$target)
))
length(n12)#293
#random30
n13 <- unique(c(
  as.character(random30$source),
  as.character(random30$target)
))
length(n13)#214

n14 <- unique(c(
  as.character(random30to$source),
  as.character(random30to$target)
))
length(n14)#214
#120
n15 <- unique(c(as.character(best120$source),
                as.character(best120$target)))
length(n15)#279

n16 <- unique(c(
  as.character(best120to$source),
  as.character(best120to$target)
))
length(n16)#233
#########################################################################
auxN <- data.frame(
  labchip = as.character(classKD$VTb) %in% n1,
  rnaseq = as.character(classKD$VTb) %in% n2,
  exonesAll = as.character(classKD$VTb) %in% n5,
  exonesComp = as.character(classKD$VTb) %in% n6,
  alt5All = as.character(classKD$VTb) %in% n7,
  alt5Comp = as.character(classKD$VTb) %in% n8,
  IRAll = as.character(classKD$VTb) %in% n9,
  IRComp = as.character(classKD$VTb) %in% n10,
  alt3all = as.character(classKD$VTb) %in% n11,
  alt3comp = as.character(classKD$VTb) %in% n12,
  random30 = as.character(classKD$VTb) %in% n13,
  random30to = as.character(classKD$VTb) %in% n14,
  best120 = as.character(classKD$VTb) %in% n15,
  best120to = as.character(classKD$VTb) %in% n16,
  order = as.numeric(classKD$ORDER)
)

#finalMatrixCOmparison
forPlot <- auxN
head(forPlot)
#colnames(forPlot)<-c("labchip","rnaseq")
rownames(forPlot) <- as.character(classKD$VTb)
forPlot$kds <- as.character(classKD$VTb)
dim(forPlot)#312*10
forPlot <- forPlot[rowSums(forPlot[, 1:4]) >= 2, ]
dim(forPlot)
#227
forPlot[forPlot == TRUE] <- 1
forPlot[forPlot == FALSE] <- 0
head(forPlot)
forPlotMelt = melt(forPlot, id.vars = c("kds", "order"))
head(forPlotMelt)
#################################################################
pdf(
  "presence_ausence_nodes_AS_ex_alt5_ir_alt3_random.pdf",
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
#stats about centrality rankking:
dim(auxN)#312 * 15
auxCC <- auxN
rownames(auxCC) <- as.character(classKD$VTb)
head(auxCC)

ii1 <- match(rownames(auxCC), rownames(centLC))#labchip
ii2 <- match(rownames(auxCC), rownames(centR))#rnaseq
ii3 <- match(rownames(auxCC), rownames(centExComp))#compExons
ii4 <- match(rownames(auxCC), rownames(centExAll))#AllExons
#####
ii5 <- match(rownames(auxCC), rownames(cent5AlComp))#5altComp
ii6 <- match(rownames(auxCC), rownames(cent5AltAll))#5altAll
ii7 <- match(rownames(auxCC), rownames(centIRComp))#IRcom
ii8 <- match(rownames(auxCC), rownames(centIRAll)) #IR All
#####
ii9 <- match(rownames(auxCC), rownames(cent3AlComp)) #3AltComp
ii10 <- match(rownames(auxCC), rownames(cent3AltAll)) #3AltAll

ii11 <- match(rownames(auxCC), rownames(centRandomCR)) #30random
ii12 <- match(rownames(auxCC), rownames(random30togCR)) #

ii13 <- match(rownames(auxCC), rownames(best120CR)) #30random
ii14 <- match(rownames(auxCC), rownames(best120togCR)) #

auxCC$labchip <- centLC$DG[ii1]
auxCC$rnaseq <- centR$DG[ii2]
auxCC$exonesAll <- centExAll$DG[ii4]
auxCC$exonesComp <- centExComp$DG[ii3]
##################################
auxCC$alt5Comp <- cent5AlComp$DG[ii5]
auxCC$alt5All <- cent5AltAll$DG[ii6]
auxCC$IRComp <- centIRComp$DG[ii7]
auxCC$IRAll <- centIRAll$DG[ii8]
##################################
auxCC$alt3comp <- cent3AlComp$DG[ii9]
auxCC$alt3all <- cent3AltAll$DG[ii10]
auxCC$random30 <- centRandomCR$DG[ii11]
auxCC$random30to <- random30togCR$DG[ii12]
auxCC$best120 <- best120CR$DG[ii13]
auxCC$best120to <- best120togCR$DG[ii14]

write.table(auxCC, file = "ComparisonDegree.tab", sep = "\t")
###################################################################
#lets try to plot in colors absense / presence according degree:
forPlot <- auxN
#colnames(forPlot)<-c("labchip","rnaseq")
rownames(forPlot) <- as.character(classKD$VTb)
forPlot$kds <- as.character(classKD$VTb)
forPlot <- forPlot[rowSums(forPlot[, 1:4]) >= 2, ]
#227
forPlot[forPlot == TRUE] <- 1
forPlot[forPlot == FALSE] <- 0
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

head(forPlot)
rownames(forPlot) <- edgesALL
forPlot$edges <- edgesALL
length(edgesALL)
head(forPlot)
forPlot <- forPlot[rowSums(forPlot[, 1:4]) >= 3, ]
dim(forPlot)
#142
forPlot[forPlot == TRUE] <- 1
forPlot[forPlot == FALSE] <- 0
head(forPlot)
forPlotMelt = melt(forPlot, id.vars = "edges")
forPlotMelt$NumOfEdge <- rep(1:nrow(forPlot), 4)
head(forPlotMelt)
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
  scale_y_discrete(expand = c(0, -0.5)) +
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
###############################################################
