args=(commandArgs(TRUE))
if(length(args)>0)
input_dir<-as.character(unlist(strsplit(args[1],"="))[2])
output_file<-as.character(unlist(strsplit(args[2],"="))[2])


library(DESeq2)
library(data.table)
sample_files <- list.files(input_dir)
sample_condition <- rep("untreated",length(sample_files))
sample_table <- data.frame(sampleName = sample_files, fileName = sample_files, condition = sample_condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = input_dir, design= ~ 1)
norm_data<-counts(estimateSizeFactors(ddsHTSeq),normalized=T)

newnames<- gsub(".txt","", (colnames(norm_data)))                   #sistema i nomi_colonna
colnames(norm_data)<-newnames
write.table(norm_data, output_file, sep="\t",quote=F,row.names=T, col.names=T)
q()
