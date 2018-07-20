library(GO.db)
zz = Ontology(GOTERM)
table(unlist(zz))

BP<-names(as.list(GOBPOFFSPRING)); length(BP)
CC<-names(as.list(GOCCOFFSPRING)); length(CC)
MF<-names(as.list(GOMFOFFSPRING)); length(MF)
