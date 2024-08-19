#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(tidyverse)
library(readr)
library(clusterProfiler)
# setwd("d:/")
# Create a parser
p <- arg_parser("matrix.counts.matrix.control_vs_treated.DESeq2.DE_results绘制热图、火山图")

# Add command line arguments
p <- add_argument(p, "--de", help="matrix.counts.matrix.control_vs_treated.DESeq2.DE_results", type="character")
p <- add_argument(p, "--emapper_annotation", help="prefixes.emapper.annotations", type="character")
p <- add_argument(p, "--TMM.EXPR.matrix", help="prefixes.TMM.EXPR.matrix", type="character")
p <- add_argument(p, "--org_db", help="org.Gspecies.eg.db_0.0.1.tar.gz", type="character")
# Parse the command line arguments
argv <- parse_args(p)

de <- argv$de

emapper_annotation <- argv$emapper_annotation

TMM.EXPR.matrix <- argv$TMM.EXPR.matrix

org_db <- argv$org_db
org_db_library <-  sub("_.*\\.tar\\.gz$", "", org_db)

dir.create('R_lib', recursive = T)
install.packages(org_db,repos = NULL,lib = 'R_lib')
do.call("library", list(org_db_library, lib.loc = 'R_lib'))
#library(org_db, lib = 'R_lib')
DE <- read.delim(de) #导入差异分析结果
library(ggrepel)
my_palette <- c('#4DBBD5FF', '#999999','#E64B35FF' )
library(ggplot2)
DE$direction  <- as.factor(ifelse(DE$padj < 0.05 & abs(DE$log2FoldChange) > 1,
                                  ifelse(DE$log2FoldChange > 1 ,'up','down'),'ns'))
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
gene <- filter(de_result, 
               abs(log2FoldChange) > 1 & padj < 0.05) %>%
  pull(id) 

geneList <- de_result$log2FoldChange
names(geneList) <- de_result$id
geneList <- sort(geneList, decreasing = T)

de_go <- enrichGO(gene = gene,
                  OrgDb = org_db_library,
                  keyType = 'GID',
                  ont = 'ALL',
                  qvalueCutoff = 0.05,
                  pvalueCutoff = 0.05,
                  minGSSize = 1)
de_go_df <- as.data.frame(de_go)
f1 <- barplot(de_go, showCategory = 10)
f2 <- dotplot(de_go, showCategory = 10)
#install.packages("ggnewscale")
library(ggnewscale)
f3 <- cnetplot(de_go, 
         foldChange = geneList, 
         showCategory = 5,
         node_label = "all", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)

ggsave("GO_barplot.png", f1, width = 10, height = 6, dpi = 300)
ggsave("GO_barplot.svg", f1, width = 10, height = 6, dpi = 300)
ggsave("GO_dotplot.png", f2, width = 10, height = 6, dpi = 300)
ggsave("GO_dotplot.svg", f2, width = 10, height = 6, dpi = 300)
ggsave("GO_cnetplot.png", f3, width = 10, height = 6, dpi = 300)
ggsave("GO_cnetplot.svg", f3, width = 10, height = 6, dpi = 300)

# KEGG分析

pathway2gene <- AnnotationDbi::select(get(org_db_library), 
                                      keys = keys(get(org_db_library)), 
                                      columns = c("Pathway","Ko")) %>% na.omit() %>% dplyr::select(Pathway, GID)
#columns(org_db)
load("kegg_info.RData")



library(clusterProfiler)
de_kegg <- enricher(gene,
                    TERM2GENE = pathway2gene,
                    TERM2NAME = pathway2name,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 1)

de_kegg_df <- as.data.frame(de_kegg)
f4 <- barplot(de_kegg, showCategory = 10)
f5 <- dotplot(de_kegg, showCategory = 10)
f6 <- cnetplot(de_kegg, 
         foldChange = geneList, 
         showCategory = 3,
         node_label = "category", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)
ggsave("kegg_barplot.png", f4, width = 10, height = 6, dpi = 300)
ggsave("kegg_barplot.svg", f4, width = 10, height = 6, dpi = 300)
ggsave("kegg_dotplot.png", f5, width = 10, height = 6, dpi = 300)
ggsave("kegg_dotplot.svg", f5, width = 10, height = 6, dpi = 300)
ggsave("kegg_cnetplot.png", f6, width = 10, height = 6, dpi = 300)
ggsave("kegg_cnetplot.svg", f6, width = 10, height = 6, dpi = 300)






