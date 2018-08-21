args=(commandArgs(TRUE))
transFile<-as.character(unlist(strsplit(args[1],"="))[2])
cisFile<-as.character(unlist(strsplit(args[2],"="))[2])
genePos<-as.character(unlist(strsplit(args[3],"="))[2])
snpPos = as.character(unlist(strsplit(args[4],"="))[2])
outputFile = as.character(unlist(strsplit(args[5],"="))[2])
networkMod = as.character(unlist(strsplit(args[6],"="))[2])
action = as.character(unlist(strsplit(args[7],"="))[2])
library(data.table)

trans = fread(transFile, data.table=FALSE)
cis = fread(cisFile, data.table=FALSE)
gPos = fread(genePos, data.table=FALSE)
sPos = fread(snpPos, data.table=FALSE)
allRes<-rbind(cis,trans)
posres<-merge(allRes,gPos,by.x="gene",by.y="geneid",all=F)
posres<-merge(posres,sPos,by.x="SNP",by.y="snp",all=F)
if (action == "moduleColors") {
  netMod = fread(networkMod, data.table=FALSE)
  posres = merge(posres,netMod,by.x="gene",by.y="geneNames",all=F)
  }

mychr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
"chr16","chr17","chr18","chr19","chrUn")
posResFilt = posres[posres$chr.x %in% mychr, ]
mylength<-rep(0,length(mychr))
for(aaa in 1:length(mychr))   {
    mylength[aaa]<-max(  c(posResFilt$s1[posResFilt$chr.x==mychr[aaa]],    posResFilt$s2[posResFilt$chr.x==mychr[aaa]],    posResFilt$pos[posResFilt$chr.y==mychr[aaa]])  )
    if (mylength[aaa] == "-Inf") {mylength[aaa] = 10000000}
}
genpos<-c(0,cumsum(mylength))
genpos<-genpos[1:(length(genpos)-1)]
gpos<-data.frame(mychr,genpos)

posResFilt$gs0<-posResFilt$gs1<-posResFilt$gs2<-0
for(aaa in 1:nrow(gpos))
  {
    #SNP pos
    posResFilt$gs0[posResFilt$chr.y==gpos$mychr[aaa]]   <-   posResFilt$pos[posResFilt$chr.y == gpos$mychr[aaa]]   +   gpos$genpos[aaa]
    #start gene
    posResFilt$gs1[posResFilt$chr.x==gpos$mychr[aaa]]   <-   posResFilt$s1[posResFilt$chr.x==gpos$mychr[aaa]]   +   gpos$genpos[aaa]
    #end gene
    posResFilt$gs2[posResFilt$chr.x==gpos$mychr[aaa]]   <-   posResFilt$s2[posResFilt$chr.x==gpos$mychr[aaa]]   +   gpos$genpos[aaa]
  }

if (action == "normal") {
  mychreven<-c("chr2","chr4","chr6","chr8","chr10","chr12","chr14","chr16","chr18","chrUn")
  posResFilt$col<-"black"
  posResFilt$col[posResFilt$chr.y %in% mychreven]<-"red"
  }

if (action == "moduleColors") {
  posResFilt$col<- posResFilt$mergedDynamicColors
  }

posResFilt$gs0<-posResFilt$gs0/1000000
posResFilt$gs1<-posResFilt$gs1/1000000
posResFilt$gs2<-posResFilt$gs2/1000000
genomelength<-c(0,max(c(posResFilt$gs0,posResFilt$gs1,posResFilt$gs2)))
posResFilt<-posResFilt[order(posResFilt$gs0),]

png(outputFile,width=23,height=17,units="cm",res=1000,type="cairo")
plot(posResFilt$gs0, posResFilt$gs1, pch=19,  xlab="SNP position [Mbp]", ylab="Gene start [Mbp]", cex=0.07, col=posResFilt$col)
for (i in 1:length(rownames(gpos)) ) {
  print(gpos$genpos[i])
  abline(v = gpos$genpos[i]/1000000, col=rgb(237, 41, 57, max = 255, alpha = 160),  lty= 3, lwd=0.5)
  abline(h = gpos$genpos[i]/1000000, col=rgb(237, 41, 57, max = 255, alpha = 160),  lty= 3, lwd=0.5)
  }
abline(0,1, col=rgb(20, 20, 20, max = 255, alpha = 160), lwd=0.5)
dev.off()

q()
