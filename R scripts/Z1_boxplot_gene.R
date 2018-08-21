args<-commandArgs(TRUE)
gene_name<-as.character(unlist(strsplit(args[1],"="))[2])
snp_file<-as.character(unlist(strsplit(args[2],"="))[2])
exp_file<-as.character(unlist(strsplit(args[3],"="))[2])
results_file<-as.character(unlist(strsplit(args[4],"="))[2])
output_folder<-as.character(unlist(strsplit(args[5],"="))[2])

library(data.table)
library(ggplot2)
library(RColorBrewer)

snp<-fread(snp_file, fill=TRUE, data.table=FALSE)
colnames(snp) <- c("snp_id", colnames(snp)[1: (length(colnames(snp))-1) ])
row.names(snp)<-snp[,1]
snp[,1]<-NULL
snp<-as.data.frame(t(snp))

gene<-fread(exp_file, fill=TRUE, data.table=FALSE)
colnames(gene) <- c("gene_id", colnames(gene)[1: (length(colnames(gene))-1) ])
row.names(gene)<-gene[,1]
gene[,1]<-NULL
gene<-as.data.frame(t(gene))

#load results mEQTL
cis_file <- paste(results_file, "results_cis.txt", sep="")
trans_file <- paste(results_file, "results_trans.txt", sep="")
result_cis<-fread(cis_file, header=TRUE, data.table=FALSE)
result_trans<-fread(trans_file, header=TRUE, data.table=FALSE)
results <- rbind(result_cis, result_trans)


list<-results[results$gene==gene_name,]                     # create table with expression levels of gene_name and every SNP associated
check_table <- cbind(gene[,gene_name, drop=FALSE], lapply(snp[,list$SNP], factor) )

pdf(paste(output_folder, "Gene_", gene_name, ".pdf", sep="" ))        # create boxplots

for(i in 2:length(colnames(check_table)))
{
print(colnames(check_table)[i])
boxplot<-ggplot( data=check_table, aes(x=check_table[,i], y=check_table[,1], fill=check_table[,i]) ) +
   geom_boxplot(alpha=0.8, outlier.alpha = 0 ) +
   geom_jitter(position=position_jitter(0.1), alpha=0.5) +
   ggtitle(  paste(gene_name, colnames(check_table)[i], sep=" and ")  ) +
   xlab("Genotype") + ylab("Gene expression (normalized)") +
   scale_fill_brewer(palette="Reds")
print(boxplot)
}
dev.off()
q()
