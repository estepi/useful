setwd("/home/emancini/Dropbox (CRG ADV)/Personal_Gosia/Shared/GosiaAndEstefi/phyper")
#phyper(q,m,n,k,lower.tail=F)
#q = overlap
#m = affected in group 1
#n = total analyzed - total affected group 1
#k = affected group 2
#input matrix:
#a matrix with columns q, m, n k
toy<-read.table("toy.txt")
apply(toy, 1, function(x) {pvavSig<-phyper(x[1],x[2],x[3],x[4],lower.tail=F)} )


  
