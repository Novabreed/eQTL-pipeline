args<-commandArgs(TRUE)
snp_file<-as.character(unlist(strsplit(args[1],"="))[2])
exp_file<-as.character(unlist(strsplit(args[2],"="))[2])
gene_position<-as.character(unlist(strsplit(args[3],"="))[2])
snp_position<-as.character(unlist(strsplit(args[4],"="))[2])
covar<-as.character(unlist(strsplit(args[5],"="))[2])
output_folder<-as.character(unlist(strsplit(args[6],"="))[2])
pvalue<-as.numeric(unlist(strsplit(args[7],"="))[2])

library(data.table)
library(GenABEL)
library(MatrixEQTL)

#imposto path delle tabelle input di SNP e EXP (NB: non si tratta delle originali ma di quelle formattate sopra e presenti in input_matrrixEQTL)
SNP_file_name = snp_file
EXP_file_name = exp_file
SNP_pos = fread(snp_position, data.table=FALSE, fill=TRUE, stringsAsFactors = FALSE)
GENE_pos = fread(gene_position, data.table=FALSE, fill=TRUE, stringsAsFactors = FALSE)
# load genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"        # the TAB character
snps$fileOmitCharacters = "9"   # denote missing values;
snps$fileSkipRows = 1            # one row of column labels
snps$fileSkipColumns = 1         # one column of row labels
snps$fileSliceSize = 2000        # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)
# load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"       # the TAB character
gene$fileOmitCharacters = "NA"  # denote missing values;
gene$fileSkipRows = 1           # one row of column labels
gene$fileSkipColumns = 1        # one column of row labels
gene$fileSliceSize = 2000       # read file in slices of 2,000 rows
gene$LoadFile(EXP_file_name)

#load covariates
covariates_file_name = covar
errorCovariance = numeric();
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"       # the TAB character
cvrt$fileOmitCharacters = "NA"  # denote missing values;
cvrt$fileSkipRows = 1           # one row of column labels
cvrt$fileSkipColumns = 1        # one column of row labels
cvrt$fileSliceSize = 2000
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

#cis-distance and min p-value
cisDist = 1e6;
pvOutputThreshold_cis = pvalue
pvOutputThreshold_tra = pvalue
useModel = modelLINEAR

#output paths
output_file_name_tra = paste(output_folder, "results_trans.txt", sep="")
output_file_name_cis = paste(output_folder, "results_cis.txt", sep="")


################ RUN THE MAIN ANALYSIS #####################
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = SNP_pos,
genepos = GENE_pos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

# plot results (cis + trans)
pdf(paste(output_folder, "qqplot.pdf", sep=""))
plot(me, pch = 16, cex = 0.7)
dev.off()




################ RUN THE LAMBDA ANALYSIS #####################
snp<-fread(snp_file, data.table=FALSE, fill=TRUE)
colnames(snp)<-c("snp_id", colnames(snp)[1:(length(colnames(snp)) -1 )])
snp <- snp[sample(1:nrow(snp), 2000, replace=FALSE),]
write.table(snp, gsub(".txt","_lambda.txt", snp_file), sep="\t",quote=F, row.names=F, col.names=T )

exp<-fread(exp_file, data.table=FALSE, fill=TRUE)
colnames(exp)<-c("gene_id", colnames(exp)[1:(length(colnames(exp)) -1 )])
exp <- exp[sample(1:nrow(exp), 2000, replace=FALSE),]
write.table(exp, gsub(".txt","_lambda.txt", exp_file), sep="\t",quote=F, row.names=F, col.names=T )

SNP_file_name = gsub(".txt","_lambda.txt", snp_file)
EXP_file_name = gsub(".txt","_lambda.txt", exp_file)
SNP_pos <- SNP_pos[SNP_pos[,1] %in% snp[,1],]
GENE_pos <- GENE_pos[GENE_pos[,1] %in% exp[,1],]
pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 1
output_file_name_tra = paste(output_folder, "lambda_trans.txt", sep="")
output_file_name_cis = paste(output_folder, "lambda_cis.txt", sep="")


# load genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"        # the TAB character
snps$fileOmitCharacters = "9"   # denote missing values;
snps$fileSkipRows = 1            # one row of column labels
snps$fileSkipColumns = 1         # one column of row labels
snps$fileSliceSize = 2000        # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

# load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"       # the TAB character
gene$fileOmitCharacters = "NA"  # denote missing values;
gene$fileSkipRows = 1           # one row of column labels
gene$fileSkipColumns = 1        # one column of row labels
gene$fileSliceSize = 2000       # read file in slices of 2,000 rows
gene$LoadFile(EXP_file_name)

# load covariates
covariates_file_name = covar
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"       # the TAB character
cvrt$fileOmitCharacters = "NA"  # denote missing values;
cvrt$fileSkipRows = 1           # one row of column labels
cvrt$fileSkipColumns = 1        # one column of row labels
cvrt$fileSliceSize = 2000
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

# run the analysis
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = SNP_pos,
genepos = GENE_pos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);


# plot results (cis + trans)
pdf(paste(output_folder, "lambda_qqplot.pdf", sep=""))
plot(me, pch = 16, cex = 0.7)
dev.off()


#Extract observed and expected distribution of p-values CIS
mycisp<-sort(me$cis$eqtls$pvalue)
#mycisq<-ppoints(mycisp)
#Transform pvalues in chisquare with 1df
mycischisq<-qchisq(mycisp,df=1,lower.tail=F)
#Estimate inflation of chisquare
cislambda<-estlambda(mycischisq, proportion=1, method="median")

#Extract observed and expected distribution of p-values TRANS
mytransp<-sort(me$trans$eqtls$pvalue)
#mytransq<-ppoints(mytransp)
#Transform pvalues in chisquare with 1df
mytranschisq<-qchisq(mytransp,df=1,lower.tail=F)
#Estimate inflation of chisquare
translambda<-estlambda(mytranschisq, proportion=1, method="median")

fileConn<-file( paste(output_folder , "lambda.txt", sep="" ))
writeLines(c("translambda = ", translambda$estimate, "cislambda = ", cislambda$estimate), fileConn)
close(fileConn)
q()
