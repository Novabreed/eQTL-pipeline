args<-commandArgs(TRUE)
remove = as.character(unlist(strsplit(args[1],"="))[2])
input_file = as.character(unlist(strsplit(args[2],"="))[2])
output_file = as.character(unlist(strsplit(args[3],"="))[2])

expdata0 = read.table(input_file)
removeList = unlist(strsplit(remove,","))
expdata = expdata0[,!(colnames(expdata0) %in% removeList)]

names<-gsub( "_rep1", "", colnames(expdata) )
names<-gsub( "_rep2", "", names)
names<-unique(names)

#merge replicates
expdata_merged <- matrix(0, ncol = length(names), nrow = length(rownames(expdata)))
expdata_merged<-data.frame(expdata_merged)
colnames(expdata_merged)<-names
rownames(expdata_merged)<-rownames(expdata)

for (n in 1:length(names)) {
    print(names[n])
    pos<-grep(names[n], colnames(expdata))
    if (length(pos)==1) {
        expdata_merged[,n]<-expdata[,pos]
        }
    if (length(pos)==2) {
        a<-data.frame(expdata[,pos[1]], expdata[,pos[2]])
        mean<-rowMeans(a)
        expdata_merged[,n]<-mean
        }
}
write.table(expdata_merged, output_file, sep="\t",quote=F,row.names=T, col.names=T)

q()
