args<-commandArgs(TRUE)
input_file = as.character(unlist(strsplit(args[1],"="))[2])
filtered_file = as.character(unlist(strsplit(args[2],"="))[2])
removed_varieties = as.character(unlist(strsplit(args[3],"="))[2])
output_folder = as.character(unlist(strsplit(args[4],"="))[2])

filteredGeneNames = rownames(read.table(filtered_file))
removeList = unlist(strsplit(removed_varieties,","))
expData0 = read.table(input_file)
expData = expData0[,!(colnames(expData0) %in% removeList)]
expData = expData[filteredGeneNames,]


filteredSampleNames = gsub( "_rep1", "", colnames(expData) )
filteredSampleNames = gsub( "_rep2", "", filteredSampleNames)
filteredSampleNames = unique(filteredSampleNames)


######### calcolo Vg #########
gTable <- matrix(0, nrow = length(rownames(expData)), ncol = length(filteredSampleNames))
gTable <-data.frame(gTable)
rownames(gTable)<-rownames(expData)
colnames(gTable)<-filteredSampleNames

#calcolo medie fra stessi genotipi
for (ns in 1:length(colnames(gTable))) {
  print( (colnames(gTable))[ns] )
  pos = grep(filteredSampleNames[ns], colnames(expData))

  if ( length(pos) == 1 ) {
      gTable[, ns] = expData[, pos]
  }

  if ( length(pos) == 2 ) {
      a<-data.frame(expData[,pos[1]], expData[,pos[2]])
      mean<-rowMeans(a)
      gTable[,ns]<-mean
  }

}

#calcolo Vg
Vg = c(rep(0, length(rownames(expData))))
Vg = ( rowSums( (gTable - rowMeans(gTable))^2 ) ) / (dim(gTable)[2])


######### calcolo Ve #########
eTable <- matrix(0, nrow = length(rownames(expData)), ncol = length(filteredSampleNames))
eTable <-data.frame(eTable)
rownames(eTable)<-rownames(expData)
colnames(eTable)<-filteredSampleNames

v1 = c()
for (ns in 1:length(colnames(eTable)) ) {
  print(ns)
  pos = grep(filteredSampleNames[ns], colnames(expData))

  if ( length(pos) == 1 ) {
    v1 = c(v1, ns)
  }

  if ( length(pos) == 2 ) {
    a = data.frame(expData[,pos[1]], expData[,pos[2]])
    varE = rowSums((a - rowMeans(a))^2) / 1 #punto critico
    eTable[,ns] = varE
  }
}

eTable[,v1] = NULL


Ve = c(rep("NA", length(colnames(eTable))))
names(Ve) = colnames(eTable)
Ve = rowMeans(eTable)

######### calcolo H2 #########

Vtot = Vg+Ve
H2 = data.frame(Vg / Vtot)
write.table(H2, paste(output_folder, "H2.txt", sep=""), quote=FALSE, sep='\t', col.names=FALSE)
pdf(paste(output_folder, "H2_hist.pdf", sep=""))
  hist(H2[,1])
dev.off()



#calcolo scarto quadratico medio corretto per la media
media = rowMeans(gTable)
sqmCorr = sqrt(Vg)/media


finalTable = data.frame(
  media, Vtot, Vg, Ve, sqrt(Vg), sqmCorr, H2
)
colnames(finalTable) = c("media", "Vtot", "Vg", "Ve", "sqmG", "sqmG/media", "H^2")
write.table(finalTable, paste(output_folder, "var+H2.txt", sep=""), quote=FALSE, sep='\t', col.names=TRUE)

q()
