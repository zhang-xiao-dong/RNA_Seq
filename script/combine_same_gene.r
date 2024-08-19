#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("合并matrix.counts.matrix和matrix.TMM.EXPR.matrix相同基因")

# Add command line arguments
p <- add_argument(p, "--counts.matrix", help="matrix.counts.matrix", type="character")
p <- add_argument(p, "--TMM.EXPR.matrix", help="matrix.TMM.EXPR.matrix", type="character")
p <- add_argument(p, "--gene_id_gene_name", help="gene_id_gene_name.txt", type="character")
p <- add_argument(p, "--output", help="output prefix", type="character")

# Parse the command line arguments
argv <- parse_args(p)

library(data.table)
library(tidyverse)

counts.matrix <- argv$counts.matrix
TMM.EXPR.matrix <- argv$TMM.EXPR.matrix
gene_id_gene_name <- argv$gene_id_gene_name
outFilePref <- argv$output


a1 <- fread(counts.matrix,header = T,data.table = F)
rownames(a1) <- a1$V1
a1 <- a1[,-1]
g2s <- fread(gene_id_gene_name,header = F,data.table = F, sep="\t") #载入从gencode的gtf文件中提取的信息文件
colnames(g2s)
colnames(g2s) <- c("geneid","symbol")

#row.names(a1)<-make.names(a1[,1],TRUE)

symbol <- g2s[match(rownames(a1),g2s$symbol),"symbol"] #匹配counts行名对应的symbol
table(duplicated(symbol))  #统计重复基因名
###使用aggregate根据symbol列中的相同基因进行合并 
counts <- aggregate(a1, by=list(symbol), FUN=sum)
counts <- column_to_rownames(counts,'Group.1')

a2 <- fread(TMM.EXPR.matrix,header = T,data.table = F)
rownames(a2) <- a2$V1
a2 <- a2[,-1]
tmm_expr <- aggregate(a2, by=list(symbol), FUN=sum)
tmm_expr <- column_to_rownames(tmm_expr,'Group.1')

outCountsFilePath <- paste(outFilePref, '_combine_counts.matrix', sep = '')
outTMM.EXPR.matrixFilePath <- paste(outFilePref, '_combine_TMM.EXPR.matrix', sep = '')
write.table(counts, outCountsFilePath, sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(tmm_expr, outTMM.EXPR.matrixFilePath, sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

