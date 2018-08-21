args=(commandArgs(TRUE))
input_file<-as.character(unlist(strsplit(args[1],"="))[2])
output_file<-as.character(unlist(strsplit(args[2],"="))[2])
clust.size=2
library(gplots)
library(ape)
library(RColorBrewer)

expdata = read.table(input_file)
expdata = as.matrix(expdata)

cormx1 = 1-cor(expdata,method="spearman")
cormx.dist1 = as.dist(cormx1)
cluster1 = hclust(cormx.dist1, method = "complete")

noRep<-c("berzamino_rep1", "cabernet.franc_rep1", "carignano_rep1", "chaouch.blanc_rep1",	"chasselas_rep1", "garganega_rep1" , "glera_rep1", "plechistik_rep1", "raboso.piave_rep1", "sahibi.safid_rep1", "traminer_rep1", "V278_rep1", "verduzzo_rep1", "vernaccia_rep1")
x1<-(cluster1$labels[cluster1$order] %in% noRep)
rep_col<- rep(1, length(colnames(expdata)))
rep_col[x1]<-2
rep_col<-as.integer(rep_col)
names(rep_col)<- cluster1$labels[cluster1$order]

labelColors=c("black", "red")                         #funzione che applica i colori alle labels (leafs) per i dendrogrammi semplici
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[rep_col[which(names(rep_col) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
h1<-(cluster1$labels %in% noRep)
rep_col_h1 = rep(1, length(colnames(expdata)))
rep_col_h1[h1] = 2

pdf(output_file)
  #clustering tutti i geni, heatmap+dendrogram+as.phylo
    heatmap.2(cormx1, main = "Spearman+Complete all genes", trace="none", dendrogram="row", Rowv=as.dendrogram(cluster1), Colv=as.dendrogram(cluster1), colRow=rep_col_h1, colCol=rep_col_h1, cexCol=0.3, cexRow=0.3, col=rev(brewer.pal(9, "RdYlGn")))
  par(cex=0.25)
  dendrogram = dendrapply(as.dendrogram(cluster1), colLab)
  plot(dendrogram, horiz=TRUE)
  par(cex=1)
  plot(as.phylo(cluster1), type = "fan", cex=0.3, tip.color=rep_col_h1)
dev.off()

q()
