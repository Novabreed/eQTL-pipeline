args = commandArgs(TRUE)
conn_data =as.character(unlist(strsplit(args[1],"="))[2])
h2_data = as.character(unlist(strsplit(args[2],"="))[2])
results_cis =as.character(unlist(strsplit(args[3],"="))[2])
results_trans = as.character(unlist(strsplit(args[4],"="))[2])
all_genes = as.character(unlist(strsplit(args[5],"="))[2])
outputFile = as.character(unlist(strsplit(args[6],"="))[2])

library(data.table)
library(ggplot2)

########################### 1. Load data ###########################
c = fread(conn_data, fill=TRUE, data.table=FALSE)
rownames(c) = c[,1]
c[,1] = NULL
colnames(c) = c("conn")

h = fread(h2_data, fill=TRUE, data.table=FALSE)
rownames(h) = h[,1]
h[,1] = NULL
colnames(h) = c("h^2")


cisR = fread(results_cis, fill=TRUE, data.table=FALSE)
cisR = cisR[,2]
transR = fread(results_trans, fill=TRUE, data.table=FALSE)
transR = transR[,2]
eGenes = as.data.frame( unique( c(cisR, transR) ) )
eGenes[,2] = "eGene"
colnames(eGenes) = c("gene_ID", "group")

non_eGenes = fread(all_genes, fill=TRUE, data.table=FALSE)
non_eGenes = as.data.frame(non_eGenes[,1])
non_eGenes[,2] = "not_eGene"
colnames(non_eGenes) = c("gene_ID", "group")
non_eGenes = non_eGenes[!(non_eGenes[,1] %in% eGenes[,1]) , , drop=FALSE]

allG = rbind(eGenes, non_eGenes)
rownames(allG) = allG[,1]
allG[,1] = NULL

allG = merge(allG, c, by="row.names")
rownames(allG) = allG[,1]
allG[,1] = NULL
allG = merge(allG, h, by="row.names")
rownames(allG) = allG[,1]
allG[,1] = NULL

allG[,2] = log10(allG[,2])

boxplotC<-ggplot( data=allG, aes(x=allG[,1], y=allG[,2], fill=allG[,1]) ) +
   geom_boxplot(alpha=0.8, outlier.alpha = 0 ) +
   geom_jitter(position=position_jitter(0.1), alpha=0.03) +
   ggtitle(  "Connectivity"  ) +
   xlab("Group") + ylab("Connectivity") +
   scale_fill_brewer(palette="Reds")


boxplotH<-ggplot( data=allG, aes(x=allG[,1], y=allG[,3], fill=allG[,1]) ) +
  geom_boxplot(alpha=0.8, outlier.alpha = 0 ) +
  geom_jitter(position=position_jitter(0.3), alpha=0.01) +
  ggtitle(  "Heritability"  ) +
  xlab("Group") + ylab("Heritability") +
  scale_fill_brewer(palette="Greens")

boxplotX = ggplot( data=allG, aes(x=allG[,2], y=allG[,3]) ) +
  geom_point(alpha=0.05) +
  ggtitle(  "Connectivity"  ) +
  xlab("Connectivity") + ylab("Heritability") +
  scale_fill_brewer(palette="Greens")

pdf(outputFile)
boxplotH
boxplotC
boxplotX
dev.off()

q()
