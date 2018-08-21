args<-commandArgs(TRUE)
input_file<-as.character(unlist(strsplit(args[1],"="))[2])
output_file<-as.character(unlist(strsplit(args[2],"="))[2])
min_var<-as.numeric(unlist(strsplit(args[3],"="))[2])               #variance filtering "James E. Peters, 2016" "Amber J Hackstadt, 2006" "Niklas Mahler, 2017"
min_median<-as.numeric(unlist(strsplit(args[4],"="))[2])

genes<-read.table(input_file)
exp_level <- rowSums(genes) / 92
var_level <- rowSums((genes - rowMeans(genes))^2) / (dim(genes)[2] - 1)
median <- apply(genes, 1, median, na.rm = TRUE)
exp_across_genotypes <- rowSums(genes>0)

#filtering
cond <- (var_level > min_var) & (median > min_median)
genes_filtered<-genes[cond,]
write.table(genes_filtered, output_file, sep="\t",quote=F, row.names=T, col.names=T)

#scrittura tabella
gene_summary <- data.frame(exp_level, exp_across_genotypes, var_level, exp_level, median, cond)
rownames(gene_summary) <- rownames(genes)
write.table(gene_summary, gsub("_merged_filtered.txt", "_summary.txt", output_file), sep="\t", quote=FALSE)
q()
