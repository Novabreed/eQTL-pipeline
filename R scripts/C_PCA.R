args<-commandArgs(TRUE)
input_snp<-as.character(unlist(strsplit(args[1],"="))[2])
input_exp<-as.character(unlist(strsplit(args[2],"="))[2])
output_folder<-as.character(unlist(strsplit(args[3],"="))[2])

library(data.table)
library(ggplot2)

######## PCA genotype (you can use PCs as covariates in matrix_eQTL analysis)
data<-fread(input_snp, data.table=FALSE, fill=TRUE)
rownames(data)<-data[,1]
colnames(data)<-c("snp_id", colnames(data)[1:92])
data$snp_id<-NULL
data<-(data[complete.cases(data), ])
data<-t(data)

mypca<-prcomp(data, center = TRUE)
PCA_snp<-as.data.frame(mypca$x)
PCA_snp_plot <- ggplot(data = PCA_snp, aes(x = PC1, y = PC2, label = rownames(PCA_snp)) ) +
  geom_point(alpha=0.8) +
  geom_text(aes(label=rownames(PCA_snp)), cex=2.2, nudge_y=1) +
  ggtitle("PCA on genotype data") + xlab("PC1") + ylab("PC2")

######## PCA expression phenotype
data<-fread(input_exp, data.table=FALSE, fill=TRUE)
rownames(data)<-data[,1]
colnames(data)<-c("gene_id", colnames(data)[1: (length(colnames(data)) -1) ])
data$gene_id<-NULL
data<-t(data)

mypca<-prcomp(data, center = TRUE)
PCA_exp<-as.data.frame(mypca$x)

noRep<-c("berzamino_rep1", "cabernet-franc_rep1", "carignano_rep1", "chaouch-blanc_rep1",	"chasselas_rep1", "garganega_rep1" , "glera_rep1", "plechistik_rep1", "raboso-piave_rep1", "sahibi-safid_rep1", "traminer_rep1", "V278_rep1", "verduzzo_rep1", "vernaccia_rep1")
rep_col<- rep(1, length(rownames(PCA_exp)))
rep_col[rownames(PCA_exp) %in% noRep]<-2
rep_col<-as.integer(rep_col)


PCA_exp_plot <- ggplot( data=PCA_exp, aes(x=PC1, y=PC2, label = rownames(PCA_exp), color=rep_col)) +
      geom_point(alpha=0.8) +
      geom_text(aes(label=rownames(PCA_exp), colour=rep_col), , cex=2.2, nudge_y=1) +
      ggtitle("PCA on expression data") + xlab("PC1") + ylab("PC2")

# scrittura PC e grafici
write.table(t(PCA_snp), paste(output_folder, "PC_snp.txt", sep=""), sep="\t",quote=F,row.names=T, col.names=T   )
write.table(t(PCA_exp), paste(output_folder, "PC_exp.txt", sep=""), sep="\t",quote=F,row.names=T, col.names=T   )
pdf(paste(output_folder, "PCA.pdf", sep=""))
PCA_snp_plot
PCA_exp_plot
dev.off()
q()
