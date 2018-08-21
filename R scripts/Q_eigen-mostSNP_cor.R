args=(commandArgs(TRUE))
transFile = as.character(unlist(strsplit(args[1],"="))[2])
cisFile = as.character(unlist(strsplit(args[2],"="))[2])
networkEigengenes = as.character(unlist(strsplit(args[3],"="))[2])
snpFile = as.character(unlist(strsplit(args[4],"="))[2])
outputFolder = as.character(unlist(strsplit(args[5],"="))[2])
nSNP = as.numeric(unlist(strsplit(args[6],"="))[2])
library(data.table)
library(reshape2)
library(ggplot2)


eigen = fread(networkEigengenes, data.table=FALSE)
snp = fread(snpFile, data.table=FALSE, fill=TRUE)
rownames(snp) = snp[,1]; colnames(snp) = c("id", colnames(snp)[1: (length(colnames(snp))-1) ]); snp[,1] = NULL
trans = fread(transFile, data.table=FALSE)
cis = fread(cisFile, data.table=FALSE)
allRes = rbind(cis,trans)

mostSNP = sort(table(allRes$SNP), decreasing=TRUE)[1:nSNP]
mostSNP = names(mostSNP)

eigenCor = data.frame(matrix(NA, nrow = length(colnames(eigen)), ncol = length(mostSNP)))
rownames(eigenCor) = colnames(eigen)
for (i in 1:length(mostSNP)) {
  nameSNP = mostSNP[i]
  thisSNP = unlist(snp[nameSNP, , drop=TRUE])
  eigenCor[,i] = as.vector(cor(thisSNP,eigen, use = "complete.obs", method="pearson"))
  colnames(eigenCor)[i] = mostSNP[i]
  }
eigenCor = round(eigenCor, 2)
write.table(eigenCor, paste(outputFolder, "/eigengenes_topSNP_corr.txt", sep=""), quote=FALSE)

pdf(paste(outputFolder, "/eigengenes_topSNP_corr.pdf", sep=""))
  heatMap = melt(as.matrix(eigenCor))
  ggplot(heatMap, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed()
dev.off()

q()
