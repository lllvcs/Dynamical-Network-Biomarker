# Copyright © 2022 LVCS. All Rights Reserved
# 用于从KEGG的pathway中提取基因的symbol

# pathway的json文件获取网址 以阿兹海默症的pathway编号hsa05010为例
# http://togows.dbcls.jp/entry/pathway/hsa05010/genes.json

# 导入相关包
library(rjson)
library(tidyr)

# 在线读取pathway文件
r = fromJSON(file = "http://togows.dbcls.jp/entry/pathway/hsa05010/genes.json")

# 剪切方法二选一
gene_symbol <- sapply(1:length(r[[1]]), function(i) { strsplit(strsplit(r[[1]][i][[1]], split = "\\s+")[[1]][1], "[;]")[[1]] })
# gene_symbol <- as.character(sapply(r[[1]], function(x) sapply(strsplit(x[1], ";"), function(x) x[1])))

# csv表格文件导出
write.csv(gene_symbol, "kegg_symbol.csv")