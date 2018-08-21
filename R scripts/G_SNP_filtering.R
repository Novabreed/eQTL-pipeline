args<-commandArgs(TRUE)
input_file<-as.character(unlist(strsplit(args[1],"="))[2])
output_file<-as.character(unlist(strsplit(args[2],"="))[2])
min_gt<-as.numeric(unlist(strsplit(args[3],"="))[2])
min_gf<-as.numeric(unlist(strsplit(args[4],"="))[2])
remove<-as.character(unlist(strsplit(args[5],"="))[2])

library(data.table)
snp<-fread(input_file, data.table=FALSE, fill=TRUE)
colnames(snp)<-c("snp_id", colnames(snp)[1:92])
rownames(snp)<-snp[,1]
snp[,1]<-NULL
snp[,remove] = NULL

n_na <- rowSums(is.na(snp))
n_tot <- 92-n_na
n_gt0 <- rowSums((snp==0), na.rm=TRUE)
n_gt1 <- rowSums((snp==1), na.rm=TRUE)
n_gt2 <- rowSums((snp==2), na.rm=TRUE)
n_gt_min <- pmin(n_gt0, n_gt1, n_gt2)
f_al0 <- ( (n_gt0 * 2) + (1 * n_gt1) ) / (n_tot * 2)
f_al2 <- ( (n_gt2 * 2) + (1 * n_gt1) ) / (n_tot * 2)
maf <- pmin(f_al0, f_al2)
f_gt0 <- n_gt0 / n_tot
f_gt1 <- n_gt1 / n_tot
f_gt2 <- n_gt2 / n_tot
M_gf <- pmax(f_gt0, f_gt1, f_gt2)
m_gf <- pmin(f_gt0, f_gt1, f_gt2)

#ogni SNP deve soddisfare le condizioni
cond <- ( (n_tot/92) >= min_gt )   &            #ogni SNP deve essere presente in almeno >min_gt degli individui
        ( n_gt_min >= min_gf )                  #ogni classe genotipica deve avere almeno >min_gf% degli individui genotipizzati per quello SNP
snp_filtered<-snp[cond,]

#scrittura tabella SNP fitrati
write.table(snp_filtered, output_file, sep="\t")

#scrittura tabella summary
snp_summary = data.frame(n_na, n_tot, n_gt0, n_gt1, n_gt2, n_gt_min, f_al0, f_al2, f_gt0, f_gt1, f_gt2, maf, m_gf, M_gf, cond)
rownames(snp_summary) = rownames(snp)
write.table(snp_summary, gsub("filtered.txt", "summary.txt", output_file), sep="\t")

q()
