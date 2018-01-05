#use DESeq2 package to transform read count value to normalized and transformed expression data and calculate P value and fold change
#reference website:http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
#step1: install package
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2") #load package


#download RNA seq read count file from GEO database, for example, from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74729;
#raw txt file
#input data
input=read.delim("D:/demo/jinfeng/GSE74729_raw_counts.txt",head=T,sep="\t")
input2=input[,2:ncol(input)]#transform input to matrix
row.names(input2)=input[,1:1]#transform input to matrix
coldata=read.delim("D:/demo/jinfeng/coldata.txt",head=T,sep="\t") #this file is made by hand


#calculate Pvalue and fold change
dds=DESeqDataSetFromMatrix(countData=input2,colData=coldata,design=~condition)

dds$condition=relevel(dds$condition,ref="untreated")# let untreated group as reference

#output P value and fold change 
dds=DESeq(dds)
res=results(dds)
write.csv(res,"D:/demo/jinfeng/output2.csv")


######
#E2F1:ENSG00000101412
#IGF1R:ENSG00000140443
#plot specific gene's expression
plotCounts(dds,gene="ENSG00000101412",intgroup="condition",normalized=TRUE,transform=TRUE)
plotCounts(dds,gene="ENSG00000140443",intgroup="condition",normalized=TRUE,transform=TRUE)





