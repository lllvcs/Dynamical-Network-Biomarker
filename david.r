f = read.csv("./david.csv", header = T)
enrich = f
enrich_signif = enrich[which(enrich$PValue < 0.05),]
enrich_signif = enrich_signif[, c(1:3, 5)]
GO_BP = enrich_signif[which(enrich_signif$Category == "GOTERM_BP_DIRECT"),]
GO_CC = enrich_signif[which(enrich_signif$Category == "GOTERM_CC_DIRECT"),]
GO_MF = enrich_signif[which(enrich_signif$Category == "GOTERM_MF_DIRECT"),]
KEGG = enrich_signif[which(enrich_signif$Category == "KEGG_PATHWAY"),]
##处理Term
library(stringi)
GO_BP$Term <- stri_sub(GO_BP$Term, 12, 100) ##去掉编号
GO_CC$Term <- stri_sub(GO_CC$Term, 12, 100)
GO_MF$Term <- stri_sub(GO_MF$Term, 12, 100)
KEGG$Term <- stri_sub(KEGG$Term, 10, 100)
##合并数据
library(ggplot2)
library(dplyr)
##排序，并取top5
GO_BP = arrange(GO_BP, GO_BP$PValue)[1:5,]
GO_CC = arrange(GO_CC, GO_CC$PValue)[1:5,]
GO_MF = arrange(GO_MF, GO_MF$PValue)[1:5,]
KEGG = arrange(KEGG, KEGG$PValue)[1:5,]

##合并数据
enrich_signif = rbind(GO_BP, rbind(GO_CC, rbind(GO_MF, KEGG)))

######
go = enrich_signif
go = arrange(go, go$Category, go$PValue)

##图例名称设置
m = go$Category
m = gsub("TERM", "", m)
m = gsub("_DIRECT", "", m)
go$Category = m

##
GO_term_order = factor(as.integer(rownames(go)), labels = go$Term)
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62", "red")
##
ggplot(data = go, aes(x = GO_term_order, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = COLS) +
  theme_bw() +
  xlab("Terms") +
  ylab("Gene_counts") +
  labs() +
  theme(axis.text.x = element_text(face = "bold", color = "gray50", angle = 70, vjust = 1, hjust = 1))