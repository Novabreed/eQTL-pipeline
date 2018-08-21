args<-commandArgs(TRUE)
snp_file<-as.character(unlist(strsplit(args[1],"="))[2])
exp_file<-as.character(unlist(strsplit(args[2],"="))[2])
snp_file_formatted<-as.character(unlist(strsplit(args[3],"="))[2])
exp_file_formatted<-as.character(unlist(strsplit(args[4],"="))[2])
gene_position<-as.character(unlist(strsplit(args[5],"="))[2])
snp_position<-as.character(unlist(strsplit(args[6],"="))[2])
gene_position_formatted<-as.character(unlist(strsplit(args[7],"="))[2])
snp_position_formatted<-as.character(unlist(strsplit(args[8],"="))[2])
covar<-as.character(unlist(strsplit(args[9],"="))[2])
covar_formatted<-as.character(unlist(strsplit(args[10],"="))[2])
nPC<-as.numeric(unlist(strsplit(args[11],"="))[2])
remove = as.character(unlist(strsplit(args[12],"="))[2])

library(data.table)
library(GenABEL)

#formattazione dati SNPs
snp<-fread(snp_file, data.table=FALSE, fill=TRUE)
colnames(snp)<-c("snp_id", colnames(snp)[1: (length(colnames(snp)) -1 ) ])
write.table(snp, snp_file_formatted, sep="\t",quote=F, row.names=F, col.names=T )
#formattazione dati di espressione
exp<-fread(exp_file, data.table=FALSE, fill=TRUE)
colnames(exp)<-c("gene_id", colnames(exp)[1:(length(colnames(exp)) -1 )])
write.table(exp, exp_file_formatted, sep="\t",quote=F, row.names=F, col.names=T )

#copia e filtering file di gene_pos e snp_position
SNP_pos = fread(snp_position, data.table=FALSE, fill=TRUE, stringsAsFactors = FALSE);
SNP_pos = SNP_pos[SNP_pos[,1] %in% snp[,1],]
write.table(SNP_pos, snp_position_formatted, sep="\t", quote=F, row.names=F, col.names=T )
GENE_pos = fread(gene_position, data.table=FALSE, fill=TRUE, stringsAsFactors = FALSE);
GENE_pos <- GENE_pos[GENE_pos[,1] %in% exp[,1],]
write.table(GENE_pos, gene_position_formatted, sep="\t", quote=F, row.names=F, col.names=T )

#formattazione covariate (PCs)
cv <- fread(covar, data.table=FALSE, fill=TRUE)
colnames(cv)<-c("id", colnames(cv)[1:91])
cv <- cv[1:nPC,]
cv[,remove] = NULL
write.table(cv, covar_formatted, sep="\t", quote=F, row.names=F, col.names=T )

q()
