#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(tidyverse)
library(readr)
# setwd("d:/")
# Create a parser
p <- arg_parser("matrix.counts.matrix.control_vs_treated.DESeq2.DE_results绘制热图、火山图")

# Add command line arguments
p <- add_argument(p, "--de", help="matrix.counts.matrix.control_vs_treated.DESeq2.DE_results", type="character")
p <- add_argument(p, "--emapper_annotation", help="prefixes.emapper.annotations", type="character")
p <- add_argument(p, "--TMM.EXPR.matrix", help="prefixes.TMM.EXPR.matrix", type="character")
p <- add_argument(p, "--sample_info", help="sample_info.txt", type="character")
# Parse the command line arguments
argv <- parse_args(p)

de <- argv$de

emapper_annotation <- argv$emapper_annotation

TMM.EXPR.matrix <- argv$TMM.EXPR.matrix

sample_info <- argv$sample_info

DE <- read.delim(de) #导入差异分析结果
library(ggrepel)
my_palette <- c('#4DBBD5FF', '#999999','#E64B35FF' )
library(ggplot2)
DE$direction  <- as.factor(ifelse(DE$padj < 0.05 & abs(DE$log2FoldChange) > 1,
                                  ifelse(DE$log2FoldChange > 1 ,'up','down'),'ns'))
f1 <- ggplot(data = DE, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = direction, size = abs(log2FoldChange))) + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') + 
  geom_vline(xintercept = c(-1,1), linetype = 'dashed') +
  scale_color_manual(values = my_palette) + 
  scale_size(range = c(0.1,1.5)) + 
  guides(size = FALSE) +
  labs(x = 'log2 fold change',
       y = '-log10(pvalue)' ,
       title = 'Vocano plot',
       size = 'log2 fold change') +   
  theme_bw()
ggsave("火山图.png", f1, width = 10, height = 6, dpi = 300)
ggsave("火山图.svg", f1, width = 10, height = 6, dpi = 300)

gene_info <- read_delim(emapper_annotation,
                        delim = "\t", escape_double = FALSE, col_names = FALSE,
                        comment = "#", trim_ws = TRUE) %>%
  dplyr::select(ID = X1,
                GO = X10,
                KO = X12,
                Pathway = X13,
                Gene_Name = X1)
#gene_info <- as.data.frame(gene_info)
#rownames(gene_info) <- make.unique(gene_info$ID)
#gene_info <- gene_info[, -1]

gene_exp <- read.table(TMM.EXPR.matrix, header=T, row.names = 1)



#DE$ID<-rownames(DE)
#colnames(gene_info)[1] <- "GID"



library(tidyverse)
de_result <- filter(DE, abs(log2FoldChange) > 1 & padj < 0.05) %>%
  mutate(FC = 2 **log2FoldChange) %>%
  mutate(direction = if_else(
    padj > 0.05, 'ns', if_else(
      abs(log2FoldChange) < 1, 'ns', if_else(
        log2FoldChange >= 1, 'up', 'down')))) %>%
  left_join(gene_info, by = c('id' = 'ID')) %>%    # gene_info有GID,gene_exp无
  left_join(rownames_to_column(gene_exp, var = 'id'), by = 'id') %>%
  dplyr::select(-c(8:9))
#logFC>2,padi<0.05筛选差异基因，并关联基因信息表和表达矩阵
pheatmap_de <- select(de_result, -c(2:15))
row.names(pheatmap_de)<-make.names(pheatmap_de[,1],TRUE)
pheatmap_de <- pheatmap_de[,-1]
sample_info <- read.table(sample_info, header = T, row.names = 1)
library(pheatmap)
f2 <- pheatmap(pheatmap_de, 
         scale = 'row',
         color = colorRampPalette(c("green", "white", "red"))(200),
         annotation_col = dplyr::select(sample_info, Group),
         annotation_colors = list(
           stage = c(S1 = '#4DBBD5FF', 
                     S4 = '#E64B35FF')),
         cutree_rows = 2)

ggsave("基因表达热图.png", f2, width = 10, height = 6, dpi = 300)
ggsave("基因表达热图.svg", f2, width = 10, height = 6, dpi = 300)





