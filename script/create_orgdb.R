#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(tidyverse)
library(readr)
# setwd("d:/")
# Create a parser
p <- arg_parser("使用eggnog-mapper的注释结果构建OrgDB")

# Add command line arguments
p <- add_argument(p, "--emapper_annotation", help="prefixes.emapper.annotations", type="character")
p <- add_argument(p, "--tax_id", help="物种分类号", type="character")
p <- add_argument(p, "--genus", help="属名", type="character")
p <- add_argument(p, "--species", help="种名", type="character")
p <- add_argument(p, "--ko", help="ko00001.json(https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=)", type="character")

# Parse the command line arguments
argv <- parse_args(p)

emapper_annotation <- argv$emapper_annotation

tax_id <- argv$tax_id

genus <- argv$genus

species <- argv$species

ko <- argv$ko

library(tidyverse)
library(KEGGREST)
library(AnnotationForge) 
library(clusterProfiler) 
library(jsonlite) 
library(purrr) 
library(RCurl)
library(parallel)#用于后面的多线程

emapper <- read_delim(emapper_annotation,
                        delim = "\t", escape_double = FALSE, col_names = TRUE,
                        comment = "##", trim_ws = TRUE)
colnames(emapper)[1] <- 'query'
#将空值替换为NA，方便后续使用na.omit()函数提出没有注释到的行 
emapper[emapper==""]<-NA
###########################提取GO信息#############
#提取query列，Preferred_named列，GO列
gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit() 
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()
#构建一个空的数据框为后面填充数据
gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())
#填充gene2go数据框（这一步时间比较长，所以开启多线程）
# for (row in 1:nrow(gos)) { 
#   the_gid <- gos[row, "query"][[1]] 
#   the_gos <- str_split(gos[row,"GOs"], ",", simplify = FALSE)[[1]] 
#   df_temp <- data_frame(GID = rep(the_gid, length(the_gos)), 
#                         GO = the_gos, 
#                         EVIDENCE = rep("IEA", length(the_gos))) 
#   gene2go <- rbind(gene2go, df_temp)}

#####################多线程开始
#检测当前可用的cpu核数
cl.cores <- detectCores()
#使用指定核数8，开启集群任务
cl <- makeCluster(8)
###并行任务
#定义并行的函数
getgene2go <- function(row){
  the_gid <- gos[row, "query"][[1]] 
  the_gos <- str_split(gos[row,"GOs"], ",", simplify = FALSE)[[1]] 
  df_temp <- data_frame(GID = rep(the_gid, length(the_gos)), 
                        GO = the_gos, 
                        EVIDENCE = rep("IEA", length(the_gos))) 
  return(df_temp)
}
#开始并行分析
list_rows <- 1:nrow(gos)
clusterExport(cl,"gos")#导入外部变量gos
clusterEvalQ(cl, library(tidyverse))#加载外部环境的包tidyverse
gene2go <- parLapply(cl,list_rows,getgene2go) # lapply的并行版本
gene2go <- do.call('rbind',gene2go) # 整合结果
# 关闭集群任务
stopCluster(cl)
#####################多线程结束

#将“-”替换为NA，然后使用na.omit（）删除含有NA的行
gene2go$GO[gene2go$GO=="-"]<-NA 
gene2go<-na.omit(gene2go)
save(gene2go,file="gene2go.RData")

############提取KEGG信息##################
#将emapper中query列，KEGG_ko列提取出来
gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko) %>% na.omit() 
#将gene2ko的Ko列中"ko:"删除，不然后面找pathway会报错 
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)

#从https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir= 下载ko00001.json文件，主要是用来解析kegg的文件结构，并不是作为数据库，因此只下载这一个即可
update_kegg <- function(json = ko) { 
  pathway2name <- tibble(Pathway = character(), Name = character()) 
  ko2pathway <- tibble(Ko = character(), Pathway = character()) 
  kegg <- fromJSON(json) 
  for (a in seq_along(kegg[["children"]][["children"]])) { 
    A <- kegg[["children"]][["name"]][[a]] 
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) { 
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) { 
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]] 
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1] 
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "") 
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]] 
        kos <- str_match(kos_info, "K[0-9]*")[,1] 
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}} 
  save(pathway2name, ko2pathway, file = "kegg_info.RData")} 
# 调用函数后在本地创建kegg_info.RData文件 
update_kegg() 
# 载入kegg_info.RData文件 
load(file = "kegg_info.RData")

#根据gene2ko中的ko信息将gene2ko中的Pathway列提取出来
gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID, Pathway) %>% na.omit()

####去除重复
gene2go <- unique(gene2go) 
gene2go <- gene2go[!duplicated(gene2go),] 
gene2ko <- gene2ko[!duplicated(gene2ko),] 
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]

save(pathway2name, ko2pathway,gene2pathway,file="KEGG.RData")
#save.image(file="All.Rdata") #保存所有的数据到ALL.Rdata

############构建OrgDb(注意：邮箱(作者和维护者)中不能包括“_”)
#tax_id = "573"
#genus = "Klebsiella"
#species = "pneumoniae"

# 去除重复
gene_info <- dplyr::distinct(gene_info)


makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               maintainer='zhangxiaodong <zhangxiaodong8315@163.com>',#你自己的邮箱
               author='zhangxiaodong <zhangxiaodong8315@163.com>',#你自己的邮箱
               pathway=gene2pathway,
               version="0.0.1",
               outputDir =".",#输出路径
               tax_id=tax_id,#你的物种在NCBI的id
               genus=genus,#属名
               species=species,#种名
               goTable="go")
org <- paste0("org.", substr(genus, 1, 1), species, ".eg.db")
pkgbuild::build(org, dest_path = ".")
