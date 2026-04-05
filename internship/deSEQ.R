source("http://bioconductor.org/biocLite.R")
## Bioconductor version 3.1 (BiocInstaller 1.18.2), ?biocLite for help
biocLite("DESeq2")
## BioC_mirror: http://bioconductor.org
## Using Bioconductor version 3.1 (BiocInstaller 1.18.2), R version 3.2.0.
## Installing package(s) 'DESeq2'
## 
## The downloaded source packages are in
##  '/tmp/Rtmp2expZI/downloaded_packages'
## Old packages: 'BiocParallel', 'edgeR', 'manipulate', 'MASS', 'spatial'
library('DESeq2')
## Old packages: 'BiocParallel', 'edgeR', 'manipulate', 'MASS', 'spatial'
library('DESeq2')
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: IRanges
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
## Loading required package: Rcpp
## Loading required package: RcppArmadillo
matrix.1 <- read.table("~/kallisto/test/output.1/abundance.txt", header=T)
matrix.2 <- read.table("~/kallisto/test/output.2/abundance.txt", header=T)
matrix.3 <- read.table("~/kallisto/test/output.3/abundance.txt", header=T)
matrix.4 <- read.table("~/kallisto/test/output.4/abundance.txt", header=T)

matrix.1$length <- NULL
matrix.1$eff_length <- NULL
matrix.1$tpm <- NULL
matrix.1$est_counts <- round(matrix.1$est_counts)
write.table(matrix.1, file="untreated1.txt", row.names=FALSE, col.names=FALSE)

matrix.2$length <- NULL
matrix.2$eff_length <- NULL
matrix.2$tpm <- NULL
matrix.2$est_counts <- round(matrix.2$est_counts)
write.table(matrix.2, file="untreated2.txt", row.names=FALSE, col.names=FALSE)

matrix.3$length <- NULL
matrix.3$eff_length <- NULL
matrix.3$tpm <- NULL
matrix.3$est_counts <- round(matrix.3$est_counts)
write.table(matrix.3, file="treated1.txt", row.names=FALSE, col.names=FALSE)

matrix.4$length <- NULL
matrix.4$eff_length <- NULL
matrix.4$tpm <- NULL
matrix.4$est_counts <- round(matrix.4$est_counts)
write.table(matrix.4, file="treated2.txt", row.names=FALSE, col.names=FALSE)



directory<- "~"
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleFiles
## [1] "treated1.txt"   "treated2.txt"   "untreated1.txt" "untreated2.txt"
sampleCondition<-c("treated","treated", "untreated","untreated")
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
sampleTable

##       sampleName       fileName condition
## 1   treated1.txt   treated1.txt   treated
## 2   treated2.txt   treated2.txt   treated
## 3 untreated1.txt untreated1.txt untreated
## 4 untreated2.txt untreated2.txt untreated
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
ddsHTSeq
## class: DESeqDataSet 
## dim: 173259 4 
## exptData(0):
## assays(1): counts
## rownames(173259): ENST00000415118 ENST00000448914 ...
##   ENST00000630922 ENST00000630347
## rowRanges metadata column names(0):
## colnames(4): treated1.txt treated2.txt untreated1.txt
##   untreated2.txt
## colData names(1): condition
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))

dds<-DESeq(ddsHTSeq)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
res<-results(dds)
res<-res[order(res$padj),]
res
## log2 fold change (MAP): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 173259 rows and 6 columns
##                  baseMean log2FoldChange      lfcSE        stat
##                 <numeric>      <numeric>  <numeric>   <numeric>
## ENST00000341272  4335.776       2.886391 0.09226691    31.28306
## ENST00000268251 11247.083       1.919525 0.06374793    30.11117
## ENST00000244741 24918.736       4.482591 0.16090552    27.85853
## ENST00000274026  1465.294      -5.041049 0.19362874   -26.03461
## ENST00000508792  1598.958      -4.894165 0.19421428   -25.19982
## ...                   ...            ...        ...         ...
## ENST00000631343 4.5839868    -0.17319612  0.4173675 -0.41497271
## ENST00000630654 0.0000000             NA         NA          NA
## ENST00000625398 1.2315779    -0.05977933  0.3725500 -0.16045985
## ENST00000630922 2.5109605     0.02232550  0.4419170  0.05051968
## ENST00000630347 0.2362029    -0.12915518  0.2374691 -0.54388217
##                        pvalue          padj
##                     <numeric>     <numeric>
## ENST00000341272 7.933226e-215 4.409604e-210
## ENST00000268251 3.460703e-199 9.617985e-195
## ENST00000244741 8.490649e-171 1.573147e-166
## ENST00000274026 2.009711e-149 2.792694e-145
## ENST00000508792 4.023554e-140 4.472904e-136
## ...                       ...           ...
## ENST00000631343     0.6781619            NA
## ENST00000630654            NA            NA
## ENST00000625398     0.8725188            NA
## ENST00000630922     0.9597083            NA
## ENST00000630347     0.5865226            NA
Upregulated <- res[which(res$padj < 0.05 & res$log2FoldChange > 1),]
Upregulated
## log2 fold change (MAP): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 1460 rows and 6 columns
##                  baseMean log2FoldChange      lfcSE      stat
##                 <numeric>      <numeric>  <numeric> <numeric>
## ENST00000341272  4335.776       2.886391 0.09226691  31.28306
## ENST00000268251 11247.083       1.919525 0.06374793  30.11117
## ENST00000244741 24918.736       4.482591 0.16090552  27.85853
## ENST00000307365  4122.396       2.193406 0.09279652  23.63673
## ENST00000227507 27467.895       1.409140 0.06166382  22.85198
## ...                   ...            ...        ...       ...
## ENST00000619291  19.68681       1.640674  0.5807926  2.824888
## ENST00000488979 216.33648       1.039288  0.3679124  2.824825
## ENST00000526173  27.87218       1.632279  0.5780705  2.823668
## ENST00000484744 133.07497       1.162213  0.4116107  2.823574
## ENST00000559279  21.34107       1.576664  0.5589252  2.820885
##                        pvalue          padj
##                     <numeric>     <numeric>
## ENST00000341272 7.933226e-215 4.409604e-210
## ENST00000268251 3.460703e-199 9.617985e-195
## ENST00000244741 8.490649e-171 1.573147e-166
## ENST00000307365 1.616223e-123 1.122951e-119
## ENST00000227507 1.396554e-115 7.762608e-112
## ...                       ...           ...
## ENST00000619291   0.004729717    0.04876583
## ENST00000488979   0.004730648    0.04876639
## ENST00000526173   0.004747747    0.04889963
## ENST00000484744   0.004749151    0.04890270
## ENST00000559279   0.004789142    0.04921421
write.table(Upregulated, file="Upregulated.txt")