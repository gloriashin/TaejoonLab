library(data.table)
library(DESeq2)
library(dbplyr)
library(Glimma)
library(edgeR)
library(EnhancedVolcano)
library(matrixStats)

getwd()#Shows the default working directory 
setwd("C:/TaejoonLab") #to set the new working Directory

csvfile1 <- read.csv('DGRP_housing.csv')
countData = csvfile1
rownames(countData) = countData$geneID
countData = countData[,c(2:ncol(countData))]
countData = countData[,-c(3,4)]
countData[1:5,1:4]

#1to 1

rundata = countData[,c(5,6,7,8)]
rundata[1:5,]

deseq<- function(df){
  # If number of zeros in row greater than half, exclude the row
  remove = apply(rundata, 1, function(x) sum(x >= 1) >= 2)
  keep<- subset(rundata, remove != FALSE)
  coldata <- data.frame(Housing = rep(c('group','isolate'), each = 2))
  
  dds = DESeqDataSetFromMatrix(countData = keep, 
                               colData = coldata, design = ~ Housing)
  dds = DESeq(dds)
  cbind(resultsNames(dds))
  res<-results(dds)
  res<-res[order(res$padj),]
  res <- res[which(res$padj < 0.05 & abs( res$log2FoldChange) > 1),]
  return(res)
}

result<- deseq(rundata)
for (k in 1:12)#Deseq 조합 생각하기


write.csv(res, "C:/TaejoonLab\\DESEQ_result_DGRP_1_woZero.csv", row.names=TRUE)



EnhancedVolcano(res,
                lab =NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DGRP DEG Result',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 1.0)



# normalizing expression data
vsd <- vst(dds, blind=FALSE)
write.csv(assay(vsd), "C:/TaejoonLab\\vst_social_nonsocial_woZero.csv")
#saving file name with variables
for (i in x){
  readr::write_csv(filename, path = paste0('DEG_DGRP_' i, '.csv'))
}


#finding specific geneID
res[rownames(res) %in% c('FBgn0000500'),]

summary(res)

#The results res object contains the follow columns
mcols(res)$description
head(res)



Upregulated <- res[which(res$padj < 0.05 & res$log2FoldChange > 1),]
Upregulated
write.table(Upregulated, file="Upregulated.txt")


#ploting-lots of ways
# 1.kkj's method
reslfc = lfcShrink(dds, coef = "Age_2mM_vs_10mM", type ='apeglm')
write.table(as.data.frame(reslfc), file = "SH-SY5Y_IronTreat_2mMvs10mM.3pairs.DESeq2lfcShrink_output.txt",quote = FALSE, sep ='\t')
vsd = vst(dds, blind = FALSE)
norm.count = counts(dds, normalized = TRUE)
write.table(assay(vsd), file = "SH-SY5Y_IronTreat.3pairs.DESeq2normalized_GSEA_input_output.txt", quote = FALSE, sep ='\t')
tmp_cor = dist(1-cor(tbl.rpkm, method = 'spearman'))
tmp_clust = hclust(tmp_cor, method = 'average')
tmp_plot = as.dendrogram(tmp_clust)
plot(tmp_plot)

# 2. making other graphes-volcano plots(original)



pCutoff = 0.05
FCcutoff = 1.0

p = EnhancedVolcano(data.frame(res), x = 'log2FoldChange', y = 'padj',
                                        xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                                        pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                                        title = "Volcano plot", subtitle = "SSA/P vs. Normal",
                                        caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
                                        legend=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                                        legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
png("DGRP_so1nonso1_DGE_VolcanoPlots.STAR.png", width=7, height=7, units = "in", res = 300)
print(p)
dev.off()

#maplot

#MDS plot
glimmaMDS(dds)
