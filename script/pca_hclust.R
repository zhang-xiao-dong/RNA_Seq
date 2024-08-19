#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(tidyverse)
library(readr)
# setwd("d:/")
# Create a parser
p <- arg_parser("从matrix.TMM.EXPR.matrix和emapper.annotation获取gene_info、gene_exp")

# Add command line arguments
p <- add_argument(p, "--emapper_annotation", help="prefixes.emapper.annotations", type="character")
p <- add_argument(p, "--TMM.EXPR.matrix", help="prefixes.TMM.EXPR.matrix", type="character")
p <- add_argument(p, "--sample_info", help="sample_info.txt", type="character")
# Parse the command line arguments
argv <- parse_args(p)

emapper_annotation <- argv$emapper_annotation

TMM.EXPR.matrix <- argv$TMM.EXPR.matrix

sample_info <- argv$sample_info

gene_info <- read_delim(emapper_annotation,
                        delim = "\t", escape_double = FALSE, col_names = FALSE,
                        comment = "#", trim_ws = TRUE) %>%
  dplyr::select(ID = X1,
                GO = X10,
                KO = X12,
                Pathway = X13,
                Gene_Name = X1)
gene_info <- as.data.frame(gene_info)
#class(gene_info)
rownames(gene_info) <- make.unique(gene_info$ID)
gene_info <- gene_info[, -1]
gene_exp <- read.table(TMM.EXPR.matrix, header=T, row.names = 1)
#colnames(gene_exp)
sample_info <- read.table(file = sample_info, sep = "\t", header=T, row.names = 1)
#row.names(sample_info)
#计算距离
sample_cor <- cor(gene_exp)
#sample_cor1 <- round(sample_cor, digits = 2)
#画图
library(pheatmap)
#library(ggplot2)
#install.packages("svglite")
f1 <- pheatmap(sample_cor, display_numbers = T,fontsize = 10, angle_col = 45)
ggsave("相关性热图.png", f1, width = 10, height = 6, dpi = 300)
ggsave("相关性热图.svg", f1, width = 10, height = 6, dpi = 300)
# 聚类树状图
sample_dist <- dist(t(gene_exp))
sample_hc <- hclust(sample_dist)
f2 <- plot(sample_hc)
#ggsave("聚类树状图.png", f2)  # 导出图片不显示
#ggsave("聚类树状图.svg", f2)
# PCA绘图
#BiocManager::install("PCAtools")
library(PCAtools)
p <- pca(gene_exp, metadata = sample_info, removeVar = 0.1)
pca_loadings <- p$loadings #某基因对pc1\pc2\pc3\pc4的贡献
pca_rotated <- p$rotated #每个主成分与样本之间的关系
f3 <- screeplot(p)  #主成分对样本差异的解释度
ggsave("主成分对样本差异的解释度.png", f3, width = 10, height = 6, dpi = 300)
ggsave("主成分对样本差异的解释度.svg", f3, width = 10, height = 6, dpi = 300)
f4 <- biplot(p,
       x = 'PC1',
       y = 'PC2',
       colby = 'Group',
       shape = 'Group',
       legendPosition = 'right')
ggsave("PCA.png", f4, width = 10, height = 6, dpi = 300)
ggsave("PCA.svg", f4, width = 10, height = 6, dpi = 300)

