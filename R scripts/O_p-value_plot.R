args = commandArgs(TRUE)
gene = as.character(unlist(strsplit(args[1],"="))[2])
chr = as.character(unlist(strsplit(args[2],"="))[2])
center = as.numeric(unlist(strsplit(args[3],"="))[2])
range = as.numeric(unlist(strsplit(args[4],"="))[2])

resCis = as.character(unlist(strsplit(args[5],"="))[2])
resTrans = as.character(unlist(strsplit(args[6],"="))[2])
outputFolder = as.character(unlist(strsplit(args[7],"="))[2])
genePos = as.character(unlist(strsplit(args[8],"="))[2])
snpPos = as.character(unlist(strsplit(args[9],"="))[2])
library(data.table)


resC = fread(resCis, data.table=FALSE, fill=TRUE, header=TRUE)
resT = fread(resTrans, data.table=FALSE, fill=TRUE, header=TRUE)
allRes<-rbind(resC,resT)
allRes = allRes[allRes$gene == gene,]
gPos = fread(genePos, data.table=FALSE)
sPos = fread(snpPos, data.table=FALSE)

posres<-merge(allRes,gPos,by.x="gene",by.y="geneid",all=F)
posres<-merge(posres,sPos,by.x="SNP",by.y="snp",all=F)

start = center - range
end = center + range
posResFilt = posres[
    (posres$chr.y == paste("chr", chr, sep="")) &
    (posres$pos >= start) &
    (posres$pos <= end)
    ,]
posResFilt$"p-value" = -log10(posResFilt$"p-value")
posResFilt$pos = (posResFilt$pos / 1000000)
posResFilt$s1 = (posResFilt$s1 / 1000000)
posResFilt$s2 = (posResFilt$s2 / 1000000)
threshold= max(posResFilt$"p-value")/2

pdf(paste(outputFolder, "/" , gene, "_","chr", chr , "_", start/1000000, "-", end/1000000, ".pdf",  sep=""))
plot(posResFilt$pos,posResFilt$"p-value",  xlab="SNP position [Mbp]", ylab="-log10(p-value)", cex=1, col="red", main=paste(gene, " in Chr", chr , " between ", start/1000000, "-", end/1000000, "Mbp", sep="" ))
if (unique(posResFilt$chr.x) == paste("chr", chr, sep="")) {
  if (unique(posResFilt$s2) > unique(posResFilt$s1)) {code = 2} else {code = 1}
  arrows(x0 = unique(posResFilt$s1) , y0 = min(posResFilt$"p-value") + 0.1, x1 = unique(posResFilt$s2) , y1 = min(posResFilt$"p-value") + 0.1, length=0.10, code = code)
}
abline(a=threshold, b=0, col=rgb(20, 20, 20, max = 255, alpha = 100), lwd=0.5)
dev.off()
q()
