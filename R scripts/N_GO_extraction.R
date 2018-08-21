args = commandArgs(TRUE)
GOterm = as.character(unlist(strsplit(args[1],"="))[2])
GOannot = as.character(unlist(strsplit(args[2],"="))[2])
resCis = as.character(unlist(strsplit(args[3],"="))[2])
resTrans = as.character(unlist(strsplit(args[4],"="))[2])
output_folder = as.character(unlist(strsplit(args[5],"="))[2])
genePos = as.character(unlist(strsplit(args[6],"="))[2])
snpPos = as.character(unlist(strsplit(args[7],"="))[2])
GOtermOutput = as.character(unlist(strsplit(args[8],"="))[2])
library(data.table)

########### estrazione geni in un file
GOlist = fread(GOannot, data.table=FALSE, fill=TRUE, header=TRUE)

resC = fread(resCis, data.table=FALSE, fill=TRUE, header=TRUE)
resT = fread(resTrans, data.table=FALSE, fill=TRUE, header=TRUE)

extG = c()
for (i in 1:length(rownames(GOlist)) ) {
  print(i)
  if (grepl(GOterm, GOlist[i,2])) {
    extG = c(extG, GOlist[i,1])
    }
  }
extCis = resC[ resC[,2] %in% extG, ]
extTrans = resT[ resT[,2] %in% extG, ]
write.table(extCis, paste(output_folder, "/", GOtermOutput, "_cis.txt" ,sep="") ,quote=FALSE, sep='\t', row.names=FALSE)
write.table(extTrans, paste(output_folder, "/", GOtermOutput, "_trans.txt" ,sep="") ,quote=FALSE, sep='\t', row.names=FALSE)


########### grafici per ogni cromosoma

gPos = fread(genePos, data.table=FALSE)
sPos = fread(snpPos, data.table=FALSE)
allRes<-rbind(extCis,extTrans)

posres<-merge(allRes,gPos,by.x="gene",by.y="geneid",all=F)
posres<-merge(posres,sPos,by.x="SNP",by.y="snp",all=F)
intGenes = unique(posres[,2])

mychr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
"chr16","chr17","chr18","chr19","chrUn")
posResFilt = posres[posres$chr.x %in% mychr, ]

for (gene in 1:length(intGenes)) {
  thisGene = intGenes[gene]
  print(thisGene)

  pdf(paste(output_folder, "/", GOtermOutput, "_", thisGene ,".pdf" ,sep=""),width=20,height=15)
  for (chr in 1:length(mychr)) {
    print(chr)
    thisResults = posResFilt[(posResFilt$chr.y == mychr[chr]) & (posResFilt$gene == thisGene), ]
    if ( dim(thisResults)[1] == 0 ) { next }
    thisResults$"p-value" = -log10(thisResults$"p-value")
    thisResults$pos = (thisResults$pos)/1000000
    threshold= max(thisResults$"p-value")/2
    plot(thisResults$pos, thisResults$"p-value",  xlab="SNP position [Mbp]", ylab="-log10(p-value)", cex=1, col="red", main=paste(thisGene, " in Chr ", chr , sep="" ))

    if (unique(thisResults$chr.x) == mychr[chr]) {
      if (unique(thisResults$s2) > unique(thisResults$s1)) {code = 2} else {code = 1}
      arrows(x0 = unique(thisResults$s1)/1000000 , y0 = min(thisResults$"p-value") + 0.1, x1 = unique(thisResults$s2)/1000000 , y1 = min(thisResults$"p-value") + 0.1, length=0.10, code = code)
    }
    abline(a=threshold, b=0, col="black")
  }
  dev.off()
}

q()
