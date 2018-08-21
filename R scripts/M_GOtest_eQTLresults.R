args = commandArgs(TRUE)
geneToGO = as.character(unlist(strsplit(args[1],"="))[2])
inputIntGenes = as.character(unlist(strsplit(args[2],"="))[2])
inputAllGenes = as.character(unlist(strsplit(args[3],"="))[2])
outputPrefix = as.character(unlist(strsplit(args[4],"="))[2])
firstGO = as.numeric(unlist(strsplit(args[5],"="))[2])

library(topGO)
library(data.table)

############################ Setting steps ############################

#list of all genes use, NB: it's good practice to use the same initial set of genes used to obtain your interesting genes (format = one column with gene names)
universeGenes = fread(inputAllGenes, data.table=FALSE, fill=TRUE, header=TRUE)
universeGenes = universeGenes[, 1, drop=FALSE]
#list of interesting genes (format = one column with gene names)
interestingGenes = fread(inputIntGenes, data.table=FALSE, fill=TRUE, header=TRUE)

#table with all genes, marked as interesting or not
geneList = (factor(as.integer(universeGenes[,1] %in% interestingGenes[,2])))
names(geneList) = (universeGenes[,1])

#mapping file between genes and GO terms (NB: not every gene has annotated GO terms)
geneID2GO <- readMappings(file = geneToGO)

#construct the GOdata object required for the analysis with topGO
GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)


algorithm="classic" #weight, elim, classic
######################## Run the test ########################
resultFisherMF <- runTest(GOdataMF, algorithm = algorithm, statistic = "fisher")
resultFisherBP <- runTest(GOdataBP, algorithm = algorithm, statistic = "fisher")
resultFisherCC <- runTest(GOdataCC, algorithm = algorithm, statistic = "fisher")
#or run the test with different algorithm/statistic test e.g.:
  #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

#run the analysis only for one GO term


######################## Visualize results ########################
resMF <- GenTable(GOdataMF, classicFisher = resultFisherMF, orderBy = "classicFisher", topNodes = firstGO, ranksOf = "classic")
resBP <- GenTable(GOdataBP, classicFisher = resultFisherBP, orderBy = "classicFisher", topNodes = firstGO, ranksOf = "classic")
resCC <- GenTable(GOdataCC, classicFisher = resultFisherCC, orderBy = "classicFisher", topNodes = firstGO, ranksOf = "classic")

#visualize used GO terms and relative count of associated genes
sel.termsMF = usedGO(GOdataMF)
num.ann.genesMF <- countGenesInTerm(GOdataMF, sel.termsMF)
num.ann.genesMF

sel.termsBP = usedGO(GOdataBP)
num.ann.genesBP <- countGenesInTerm(GOdataBP, sel.termsBP)
num.ann.genesBP

sel.termsCC = usedGO(GOdataCC)
num.ann.genesCC <- countGenesInTerm(GOdataCC, sel.termsCC)
num.ann.genesCC

#print GO graph built using the first n significant nodes and print results tables

pdf(paste(outputPrefix, "_graphs.pdf", sep= ""))

showSigOfNodes(GOdataMF, score(resultFisherMF), firstSigNodes = firstGO, useInfo = 'all')
title(main=paste( "Molecular Function of the most", firstGO, "enriched GO terms", sep=" "))
showSigOfNodes(GOdataBP, score(resultFisherBP), firstSigNodes = firstGO, useInfo = 'all')
title(main=paste( "Biological Process of the most", firstGO, "enriched GO terms", sep=" "))
showSigOfNodes(GOdataCC, score(resultFisherCC), firstSigNodes = firstGO, useInfo = 'all')
title(main=paste( "Cellular Component of the most", firstGO, "enriched GO terms", sep=" "))
dev.off()

write.table(resMF, paste(outputPrefix, "_MF.txt", sep=""), quote=FALSE, sep='\t')
write.table(resBP, paste(outputPrefix, "_BP.txt", sep=""), quote=FALSE, sep='\t')
write.table(resCC, paste(outputPrefix, "_CC.txt", sep=""), quote=FALSE, sep='\t')
q()
