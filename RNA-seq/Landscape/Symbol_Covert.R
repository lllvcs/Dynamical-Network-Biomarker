library("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
gene <- read.csv("FPKM.csv", sep = ',', header = T)[, 1]
geneID <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "ensembl_peptide_id"), filters = "ensembl_gene_id", values = gene, mart = mart)
