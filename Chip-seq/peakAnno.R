# Copyright © 2022 LVCS. All Rights Reserved
#导入相关包：
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#导入Chip-seq的bed格式文件
peak <- readPeakFile('test.bed')

#peak 在染色体上的分布
covplot(peak, chr = c("chr1", "chr2"))

#peak 在TSS位点附件的分布
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# 定义TSS上下游的距离
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peak, windows = promoter)
# 热图TSS分布
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = "red")
# 折线图TSS分布
plotAvgProf(
  tagMatrix,
  xlim = c(-3000, 3000),
  xlab = "Genomic Region (5'->3')",
  ylab = "Read Count Frequency")

#peak关联基因注释
peakAnno <- annotatePeak(
    peak,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db")

#绘制类别饼图
plotAnnoPie(peakAnno)

# csv表格文件导出
write.csv(
    as.data.frame(peakAnno),
    "peak.annotation.csv")
