library(data.table)
library(DESeq2)
library(dbplyr)
library(Glimma)
library(edgeR)

getwd()#Shows the default working directory 
setwd("C:/TaejoonLab") #to set the new working Directory

csvfile1 <- read.csv('reads_total_merged.csv')
countData = csvfile1
rownames(countData) = countData$geneID
countData = countData[,c(3:ncol(countData))]
countData[1:10,1:6]

# If number of zeros in row greater than half, exclude the row
remove = apply(countData, 1, function(x) sum(x >= 1) >= 8)

keep<- subset(countData, remove != FALSE)


# table <- table[!(table$row == "condition")] #excluding the row
# table<- subset(table, condition)

coldata = data.frame(Housing = rep(c('housing','social'), c(4,12)))

dds = DESeqDataSetFromMatrix(countData = keep, 
                             colData = coldata, design = ~ Housing)
dds = DESeq(dds)
cbind(resultsNames(dds))
res<-results(dds)
res<-res[order(res$padj),]
#res
#write.csv(res, "C:/TaejoonLab\\DESEQ_result_DSK_woZero.csv", row.names=TRUE)

#MDS plot
glimmaMDS(dds)



# normalizing expression data
#vsd <- vst(dds, blind=FALSE)
#write.csv(assay(vsd), "C:/TaejoonLab\\vst_GH_SH.csv")

#drawing heatmap with vst result
#sampleDists <- dist(t(assay(vsd)))
#library("RColorBrewer")
#library(pheatmap)
#sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
#colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#pheatmap(sampleDistMatrix,
#         clustering_distance_rows=sampleDists,
#         clustering_distance_cols=sampleDists,
#         col=colors)






#finding specific geneID
#res[rownames(res) %in% c('FBgn0000500'),]

#summary(res)

#The results res object contains the follow columns
#mcols(res)$description
#head(res)



#Upregulated <- res[which(res$padj < 0.05 & res$log2FoldChange > 1),]
#Upregulated
#write.table(Upregulated, file="Upregulated.txt")


#ploting-lots of ways
# 1.kkj's method
#tmp_cor = dist(1-cor(tbl.rpkm, method = 'spearman'))
#tmp_clust = hclust(tmp_cor, method = 'average')
#tmp_plot = as.dendrogram(tmp_clust)
#plot(tmp_plot)

# 2. making other graphes-volcano plots(original)
#BiocManager::install('EnhancedVolcano')
#library(EnhancedVolcano)

#pCutoff = 0.05
#FCcutoff = 1.0

#p = EnhancedVolcano(data.frame(res), x = 'log2FoldChange', y = 'padj',
#                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
#                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
#                    title = "Volcano plot", subtitle = "SSA/P vs. Normal",
#                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
#                    legend=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
#                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)

#png("DGE_VolcanoPlots.STAR.png", width=7, height=7, units = "in", res = 300)
#print(p)
#dev.off()


# MA plot
#BiocManager::install('ggpubr')
library(ggpubr)
ggmaplot(res, main = expression("Housing(GH-SH)"),
         fdr = 0.05, fc = 2, size = 0.4,
         palette = c("#d1495b", "#edae49", "#66a182"),
         genenames = as.vector(rownames(res)),
         ylim=c(-10,20),
         legend = "top", top = 0,
         font.label = c("bold", 8),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())
