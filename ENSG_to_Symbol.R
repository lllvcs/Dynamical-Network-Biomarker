# 加载相关包
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringi)
# 读取数据
data <- read.csv('FPKM.csv')
# 需要转换的Ensembl_ID
Ensembl_ID <- stri_sub(data$ENSG_ID, 1, 15)
# 查看org.Hs.eg.db 包提供的转换类型
keytypes(org.Hs.eg.db)
# 转换
gene_symbol <- mapIds(org.Hs.eg.db, keys = Ensembl_ID, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
# 查看转换的结果
head(gene_symbol)
data$symbol <- mapIds(org.Hs.eg.db, keys = Ensembl_ID, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
write.csv(data, file = "fine.csv")
