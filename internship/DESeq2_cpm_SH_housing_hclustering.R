library(data.table)
library(DESeq2)
library(dbplyr)
library(Glimma)
library(edgeR)

getwd()#Shows the default working directory 
setwd("C:/TaejoonLab") #to set the new working Directory

csvfile1 <- read.csv('DGRP_73.csv')
countData = csvfile1
rownames(countData) = countData$geneID
countData = countData[,c(2:ncol(countData))]
#countData = countData[,-6]
countData[1:10,1:6]

# If number of zeros in row greater than half, exclude the row
remove = apply(countData, 1, function(x) sum(x >= 1) >= 5 )
keep<- subset(countData, remove != FALSE)


first <- countData[,c(1,2,6,7)]
length(rownames(first[apply(first, 1, function(x) sum(x == 0) == 4) ,]))
#length(first[apply(first, 1, function(x) sum(x == 0) == 4) ,])
#remove1 = apply(first, 1, function(x) sum(x >= 1) >= 4 )
#keep1<- subset(first, remove1 != FALSE)
two <- countData[,c(3,4,5,8,9,10)]
length(rownames(two[apply(two, 1, function(x) sum(x == 0) == 6) ,]))
#length(two[apply(two, 1, function(x) sum(x == 0) == 6) ,])
#remove2 = apply(two, 1, function(x) sum(x >= 1) >= 6 )
#keep2<- subset(two, remove2 != FALSE)
two[1:10,1:6]

# table <- table[!(table$row == "condition")] #excluding the row
# table<- subset(table, condition)

coldata = data.frame(Housing = rep(c('group','isolate'), c(5,5)))
#coldata = data.frame(Housing =c('group','group','group','isolate','isolate','isolate','group','group','isolate','isolate','group','group','isolate','isolate','group','group','isolate','isolate','group','group','isolate','isolate','group','group','isolate','isolate'))
dds = DESeqDataSetFromMatrix(countData = countData, 
                             colData = coldata, design = ~ Housing)
dds = DESeq(dds)
cbind(resultsNames(dds))
res<-results(dds)
res<-res[order(res$padj),]
#res
#write.csv(res, "C:/TaejoonLab\\DESEQ_result_DSK_woZero.csv", row.names=TRUE)

#MDS plot
#glimmaMDS(dds)


#hclustering plot example distance=cor function
library(matrixStats)
vsd <- vst(dds, blind=FALSE)
rv <- rowVars(assay(vsd))
o <- order(rv,decreasing=TRUE)

dists <- dist(1-cor(assay(vsd), method = 'spearman'))
hc <- hclust(dists, method = 'average')
dend <- as.dendrogram(hc)
o.dend <- order.dendrogram(dend)
#labels(dend) <- vsd$sample[o.dend]
plot(dend)
plot(hc, labels=vsd$sample)

#round value
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}






