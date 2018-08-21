args<-commandArgs(TRUE)
inputFile = as.character(unlist(strsplit(args[1],"="))[2])
inputFileVg = as.character(unlist(strsplit(args[2],"="))[2])
outputFolder = as.character(unlist(strsplit(args[3],"="))[2])

library(data.table)
table1 = fread(inputFile, data.table=FALSE)
rownames(table1)= table1$gene_id_V2
table2= read.table(inputFileVg, header=TRUE)
table3 = merge (table1, table2, by="row.names", all=FALSE )

model<-lm(
	log10(sqmG.media+0.001)	~
	gBM +
	intergenic_TE_norm +
	bp_TE_intron_norm +
	SNP_nd_norm +
	REG_transcript_bp_norm +
	REG_SV_nd_norm,
	data=table3)

sink(paste(outputFolder, "summary_lm.txt", sep=""))
	print(summary(model))
sink()

pdf(paste(outputFolder, "plots_lm.pdf", sep=""))
	plot(model)
dev.off()

q()
