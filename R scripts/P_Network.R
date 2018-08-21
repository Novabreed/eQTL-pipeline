args<-commandArgs(TRUE)
exp_data<-as.character(unlist(strsplit(args[1],"="))[2])
output_folder<-as.character(unlist(strsplit(args[2],"="))[2])

library(data.table)
library(WGCNA)
library(ggplot2)
library(gplots)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

########################### 1. Load data ###########################
dataExpr = fread(exp_data, fill=TRUE, data.table=FALSE)
rownames(dataExpr) = dataExpr[,1]
colnames(dataExpr) = c("id", colnames(dataExpr[1:91]))
dataExpr[,1] = NULL
dataExpr = t(dataExpr)

##################### 2. Analyze network connectivity (k) to choose the right soft-thresholding power (beta) ##### for type "unsigned": adjacency = |cor|^power #####################
powers = c(c(1:10), seq(from = 12, to=20, by=2))                # generate different powers
# test the connectivity with different powers
# nb: the network type is set to "unsigned" i.e. we want to keep both strong negative and positive correlations
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5, networkType = "unsigned")

pdf(paste(output_folder, "/connAnalysis.pdf", sep=""))
    # plot scale-free topology fit index as a function of the soft-thresholding power (we choose the first power that fits scale-free topology with R^2>=0,85)
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=0.5, col="red")
      # this line corresponds to using an R^2 cut-off of h
  abline(h=0.9,col="red")
  # plot mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.5, col="red")
dev.off()

########################### 3. Construct the network & find modules ###########################
#dataExpr <- dataExpr[,sample(1:ncol(dataExpr), 1000, replace=FALSE)]      # select a subset of 500 genes to test the script
softPower = sft$powerEstimate                                             # its better to check with the plot if sft$powerEstimate is the first value that surpass the R^2 cut-off of 0,85.

#adjacency matrix based on correlation between genes, correlation function=biweight midcorrelation
weighted_adjacency_matrix = adjacency(dataExpr, power = softPower, type = "unsigned", corFnc = "bicor")
    geneConnectivity = as.data.frame(rowSums(weighted_adjacency_matrix)-1)
    write.table(geneConnectivity, paste(output_folder, "/gene_conn.txt", sep=""), col.names=FALSE, quote=FALSE)
#use network topology and adjacency matrix to calculate dissTOM measurement
TOM = TOMsimilarity(weighted_adjacency_matrix)
dissTOM = 1-TOM
# hierarchical clustering of genes based on dissTOM
geneTree = hclust(as.dist(dissTOM), method = "average")

minModuleSize = 30;                                             # we like large modules, so we set the minimum module size relatively high
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)      # module identification using dynamic tree cut
dynamicColors = labels2colors(dynamicMods)                      # convert numeric lables into colors

table(dynamicMods)                                              # show number of genes for every detected module (" " is unassigned = grey)
table(dynamicColors)                                            # grey = unassigned


########################### 4. Merging of similar modules, based on eigengenes correlation ############################
MEs = (moduleEigengenes(dataExpr, colors = dynamicColors))$eigengenes     # calculate eigengenes
MEsDiss = 1-cor(MEs)                                                      # calculate dissimilarity of module eigengenes
MEsTree = hclust(as.dist(MEsDiss), method = "average")                    # cluster module eigengenes
MEsDissThres = 0.25                                                       # cut-off for dissimilarity merging
# merge similar modules, recalculate new colors and new MODULE EIGENGENES
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEsDissThres, verbose = 3)
mergedDynamicColors = merge$colors
mergedMEs = merge$newMEs
write.table(mergedMEs,paste(output_folder, "/network_mods.txt", sep=""), quote=FALSE, sep='\t', row.names=FALSE)

#print plots
pdf(paste(output_folder, "/network_mods.pdf", sep=""))

  # plot dendrogram of original modules
  plot(as.dendrogram(MEsTree, main = "Clustering of module eigengenes", xlab = "", sub = ""))
  abline(h=MEsDissThres, col = "red")

  # plot heatmap of correlations between modules (=correlation between eigengenes)
  heatmap.2(MEsDiss,
    main = "module eigengenes correlation",
    trace="none",
    dendrogram="row",
    Rowv=as.dendrogram(MEsTree),
    Colv=as.dendrogram(MEsTree),
    cexCol=1,
    cexRow=1,
    col=rev(brewer.pal(9, "RdYlGn")))
dev.off()

png(paste(output_folder, "/network_clust.png", sep=""),width=20,height=15,units="cm",res=1200,type="cairo")
  # check differences before and after merging modules NB: personally evaluate the results
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedDynamicColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste(output_folder, "/network_heatmap.png", sep=""),width=20,height=15,units="cm",res=1000,type="cairo")
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM^7
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  TOMplot(plotTOM, geneTree, mergedDynamicColors, main = "Network heatmap")

dev.off()

# rename to moduleColors and construct numerical labels corresponding to the colors
moduleColors = mergedDynamicColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;


########################### 5. Calculate membership of genes with theier respective modules ###########################

geneModuleMembership = as.data.frame(cor(dataExpr, mergedMEs, use = "p"))     # calculate the membership of a gene with his module
nSamples = nrow(dataExpr)
geneModuleMembershipPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) # test the significance of gene-module membership, return p-values


########################### 6. Data export ###########################
geneNames=colnames(dataExpr)
pValMember = geneModuleMembershipPvalue
colnames(pValMember) = gsub("ME","",colnames(pValMember))

#set up TOM matrix
myTOM = TOM
colnames(myTOM) = geneNames
rownames(myTOM) = geneNames

#select all modules (or create a subset of modules)
modules=unique(mergedDynamicColors)

#set up table
gmPvalue=data.frame( c(rep(0, length(geneNames))) )
rownames(gmPvalue) = geneNames


for (i in 1:length(geneNames)) {
  gmPvalue[i,1] = pValMember[i, mergedDynamicColors[i] ]
}

nodes = data.frame(geneNames, mergedDynamicColors, gmPvalue[,1])
#nodes is a data frame with 3 columns:
#   1. geneName
#   2. color of the module
#   3. p-value di appartenenza a tale modulo
typeModules = unique(nodes[,2])
for (mod in typeModules) {
  print(mod)
  thisMod = nodes[nodes[,2] == mod,]
  write.table(thisMod, paste(output_folder, "WGCNAnodes_", mod, ".txt", sep=""), quote=FALSE, sep='\t', row.names=FALSE)
}

inModule = is.finite(match(moduleColors, modules ))
modProbes = geneNames[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
weighted = TRUE,
threshold = 0.15,
nodeNames = modProbes,
nodeAttr = moduleColors[inModule]);

edges = cyt[[1]]
edges = edges[,1:3]

write.table(nodes, paste(output_folder, "WGCNAnodes.txt", sep=""), quote=FALSE, sep='\t', row.names=FALSE)
write.table(edges, paste(output_folder, "WGCNAedges.txt", sep=""), quote=FALSE, sep='\t', row.names=FALSE)

q()
