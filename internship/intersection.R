##social1
#social1 nonsocial1
csvfile_so1_no1 <- read.csv('DESEQ_result_so1_no1.csv')

rownames(csvfile_so1_no1) = csvfile_so1_no1[,1]
csvfile_so1_no1 = csvfile_so1_no1[,c(2:ncol(csvfile_so1_no1))]
#csvfile_so1_no1[1:5,1:4]
logp_so1_no1 <-csvfile_so1_no1[csvfile_so1_no1$log2FoldChange >0,]
logn_so1_no1 <-csvfile_so1_no1[csvfile_so1_no1$log2FoldChange <0,]

#social1 nonsocial2
csvfile_so1_no2 <- read.csv('DESEQ_result_so1_no2.csv')
rownames(csvfile_so1_no2) = csvfile_so1_no2[,1]
csvfile_so1_no2 = csvfile_so1_no2[,c(2:ncol(csvfile_so1_no2))]
#csvfile_so1_no2[1:5,1:4]
logp_so1_no2 <-csvfile_so1_no2[csvfile_so1_no2$log2FoldChange >0,]
logn_so1_no2 <-csvfile_so1_no2[csvfile_so1_no2$log2FoldChange <0,]

#social1 nonsocial3
csvfile_so1_no3 <- read.csv('DESEQ_result_so1_no3.csv')
rownames(csvfile_so1_no3) = csvfile_so1_no3[,1]
csvfile_so1_no3 = csvfile_so1_no3[,c(2:ncol(csvfile_so1_no3))]
#csvfile_so1_no3[1:5,1:4]
logp_so1_no3 <-csvfile_so1_no3[csvfile_so1_no3$log2FoldChange >0,]
logn_so1_no3 <-csvfile_so1_no3[csvfile_so1_no3$log2FoldChange <0,]

#intersections(overlapping)
inP_so1_no1no2 <- intersect(rownames(logp_so1_no1), rownames(logp_so1_no2))
inN_so1_no1no2 <- intersect(rownames(logn_so1_no1), rownames(logn_so1_no2))

inP_so1_no1no3 <- intersect(rownames(logp_so1_no1), rownames(logp_so1_no3))
inN_so1_no1no3 <- intersect(rownames(logn_so1_no1), rownames(logn_so1_no3))

inP_so1_no2no3 <- intersect(rownames(logp_so1_no3), rownames(logp_so1_no2))
inN_so1_no2no3 <- intersect(rownames(logn_so1_no3), rownames(logn_so1_no2))

inP_so1 <- intersect(rownames(logp_so1_no3), inP_so1_no1no2)
inN_so1 <- intersect(rownames(logn_so1_no3), inN_so1_no1no2)
write.csv(inP_so1, "C:/TaejoonLab\\inP_so1.csv")
write.csv(inN_so1, "C:/TaejoonLab\\inN_so1.csv")

##social2
#social2 nonsocial1
csvfile_so2_no1 <- read.csv('DESEQ_result_so2_no1.csv')

rownames(csvfile_so2_no1) = csvfile_so2_no1[,1]
csvfile_so2_no1 = csvfile_so2_no1[,c(2:ncol(csvfile_so2_no1))]
#csvfile_so2_no1[1:5,1:4]
logp_so2_no1 <-csvfile_so2_no1[csvfile_so2_no1$log2FoldChange >0,]
logn_so2_no1 <-csvfile_so2_no1[csvfile_so2_no1$log2FoldChange <0,]

#social2 nonsocial2
csvfile_so2_no2 <- read.csv('DESEQ_result_so2_no2.csv')
rownames(csvfile_so2_no2) = csvfile_so2_no2[,1]
csvfile_so2_no2 = csvfile_so2_no2[,c(2:ncol(csvfile_so2_no2))]
#csvfile_so2_no2[1:5,1:4]
logp_so2_no2 <-csvfile_so2_no2[csvfile_so2_no2$log2FoldChange >0,]
logn_so2_no2 <-csvfile_so2_no2[csvfile_so2_no2$log2FoldChange <0,]

#social2 nonsocial3
csvfile_so2_no3 <- read.csv('DESEQ_result_so2_no3.csv')
rownames(csvfile_so2_no3) = csvfile_so2_no3[,1]
csvfile_so2_no3 = csvfile_so2_no3[,c(2:ncol(csvfile_so2_no3))]
#csvfile_so2_no3[1:5,1:4]
logp_so2_no3 <-csvfile_so2_no3[csvfile_so2_no3$log2FoldChange >0,]
logn_so2_no3 <-csvfile_so2_no3[csvfile_so2_no3$log2FoldChange <0,]

#intersections(overlapping)
inP_so2_no1no2 <- intersect(rownames(logp_so2_no1), rownames(logp_so2_no2))
inN_so2_no1no2 <- intersect(rownames(logn_so2_no1), rownames(logn_so2_no2))

inP_so2_no1no3 <- intersect(rownames(logp_so2_no1), rownames(logp_so2_no3))
inN_so2_no1no3 <- intersect(rownames(logn_so2_no1), rownames(logn_so2_no3))

inP_so2_no2no3 <- intersect(rownames(logp_so2_no3), rownames(logp_so2_no2))
inN_so2_no2no3 <- intersect(rownames(logn_so2_no3), rownames(logn_so2_no2))

inP_so2 <- intersect(rownames(logp_so2_no3), inP_so2_no1no2)
inN_so2 <- intersect(rownames(logn_so2_no3), inN_so2_no1no2)
write.csv(inP_so2, "C:/TaejoonLab\\inP_so2.csv")
write.csv(inN_so2, "C:/TaejoonLab\\inN_so2.csv")



##social3
#social3 nonsocial1
csvfile_so3_no1 <- read.csv('DESEQ_result_so3_no1.csv')

rownames(csvfile_so3_no1) = csvfile_so3_no1[,1]
csvfile_so3_no1 = csvfile_so3_no1[,c(2:ncol(csvfile_so3_no1))]
#csvfile_so3_no1[1:5,1:4]
logp_so3_no1 <-csvfile_so3_no1[csvfile_so3_no1$log2FoldChange >0,]
logn_so3_no1 <-csvfile_so3_no1[csvfile_so3_no1$log2FoldChange <0,]

#social3 nonsocial2
csvfile_so3_no2 <- read.csv('DESEQ_result_so3_no2.csv')
rownames(csvfile_so3_no2) = csvfile_so3_no2[,1]
csvfile_so3_no2 = csvfile_so3_no2[,c(2:ncol(csvfile_so3_no2))]
#csvfile_so3_no2[1:5,1:4]
logp_so3_no2 <-csvfile_so3_no2[csvfile_so3_no2$log2FoldChange >0,]
logn_so3_no2 <-csvfile_so3_no2[csvfile_so3_no2$log2FoldChange <0,]

#social3 nonsocial3
csvfile_so3_no3 <- read.csv('DESEQ_result_so3_no3.csv')
rownames(csvfile_so3_no3) = csvfile_so3_no3[,1]
csvfile_so3_no3 = csvfile_so3_no3[,c(2:ncol(csvfile_so3_no3))]
#csvfile_so3_no3[1:5,1:4]
logp_so3_no3 <-csvfile_so3_no3[csvfile_so3_no3$log2FoldChange >0,]
logn_so3_no3 <-csvfile_so3_no3[csvfile_so3_no3$log2FoldChange <0,]

#intersections(overlapping)
inP_so3_no1no2 <- intersect(rownames(logp_so3_no1), rownames(logp_so3_no2))
inN_so3_no1no2 <- intersect(rownames(logn_so3_no1), rownames(logn_so3_no2))

inP_so3_no1no3 <- intersect(rownames(logp_so3_no1), rownames(logp_so3_no3))
inN_so3_no1no3 <- intersect(rownames(logn_so3_no1), rownames(logn_so3_no3))

inP_so3_no2no3 <- intersect(rownames(logp_so3_no3), rownames(logp_so3_no2))
inN_so3_no2no3 <- intersect(rownames(logn_so3_no3), rownames(logn_so3_no2))

inP_so3 <- intersect(rownames(logp_so3_no3), inP_so3_no1no2)
inN_so3 <- intersect(rownames(logn_so3_no3), inN_so3_no1no2)
write.csv(inP_so3, "C:/TaejoonLab\\inP_so3.csv")
write.csv(inN_so3, "C:/TaejoonLab\\inN_so3.csv")

inPso1so2<-intersect(inP_so1, inP_so2)
inPso2so2<-intersect(inP_so1, inP_so3)
inPso3so2<-intersect(inP_so3, inP_so2)
inPso1so2so3<-intersect(inPso1so2, inP_so3)


inNso1so2<-intersect(inN_so1, inN_so2)
inNso2so2<-intersect(inN_so1, inN_so3)
inNso3so2<-intersect(inN_so3, inN_so2)
inNso1so2so3<-intersect(inNso1so2, inN_so3)

##housing intersection
#nonsocial:317,360,707
#social: 73,563,370

csvfile_h317 <- read.csv('DESEQ_result_Housingno317g1.csv')

rownames(csvfile_h317) = csvfile_h317[,1]
csvfile_h317 = csvfile_h317[,c(2:ncol(csvfile_h317))]
#csvfile_so1_no1[1:5,1:4]
logp_h317 <-csvfile_h317[csvfile_h317$log2FoldChange >0,]
logn_h317 <-csvfile_h317[csvfile_h317$log2FoldChange <0,]

csvfile_h360 <- read.csv('DESEQ_result_Housingno360g1.csv')

rownames(csvfile_h360) = csvfile_h360[,1]
csvfile_h360 = csvfile_h360[,c(2:ncol(csvfile_h360))]
#csvfile_so1_no1[1:5,1:4]
logp_h360 <-csvfile_h360[csvfile_h360$log2FoldChange >0,]
logn_h360 <-csvfile_h360[csvfile_h360$log2FoldChange <0,]

csvfile_h707 <- read.csv('DESEQ_result_Housingno707g1.csv')

rownames(csvfile_h707) = csvfile_h707[,1]
csvfile_h707 = csvfile_h707[,c(2:ncol(csvfile_h707))]
#csvfile_so1_no1[1:5,1:4]
logp_h707 <-csvfile_h707[csvfile_h707$log2FoldChange >0,]
logn_h707 <-csvfile_h707[csvfile_h707$log2FoldChange <0,]

csvfile_h370 <- read.csv('DESEQ_result_Housingso370g1.csv')

rownames(csvfile_h370) = csvfile_h370[,1]
csvfile_h370 = csvfile_h370[,c(2:ncol(csvfile_h370))]
#csvfile_so1_no1[1:5,1:4]
logp_h370 <-csvfile_h370[csvfile_h370$log2FoldChange >0,]
logn_h370 <-csvfile_h370[csvfile_h370$log2FoldChange <0,]


csvfile_h563 <- read.csv('DESEQ_result_Housingso563g1.csv')

rownames(csvfile_h563) = csvfile_h563[,1]
csvfile_h3563 = csvfile_h563[,c(2:ncol(csvfile_h563))]
#csvfile_so1_no1[1:5,1:4]
logp_h563 <-csvfile_h563[csvfile_h563$log2FoldChange >0,]
logn_h563 <-csvfile_h563[csvfile_h563$log2FoldChange <0,]

csvfile_h73 <- read.csv('DESEQ_result_Housingso73g1.csv')

rownames(csvfile_h73) = csvfile_h73[,1]
csvfile_h73 = csvfile_h73[,c(2:ncol(csvfile_h73))]
#csvfile_so1_no1[1:5,1:4]
logp_h73 <-csvfile_h73[csvfile_h73$log2FoldChange >0,]
logn_h73 <-csvfile_h73[csvfile_h73$log2FoldChange <0,]

#nonsocial
inP_317_360 <- intersect(rownames(logp_h317), rownames(logp_h360))
inN_317_360 <- intersect(rownames(logn_h317), rownames(logn_h360))

inP_no <- intersect(inP_317_360, rownames(logp_h707))
inN_no <- intersect(inN_317_360, rownames(logn_h707))

#social
inP_370_563 <- intersect(rownames(logp_h370), rownames(logp_h563))
inN_370_563 <- intersect(rownames(logn_h370), rownames(logn_h563))

inP_so <- intersect(inP_370_563, rownames(logp_h73))
inN_so <- intersect(inN_370_563, rownames(logn_h73))

#intersect
inP_housing <-intersect(inP_no, inP_so)
inN_housing <-intersect(inN_no, inN_so)