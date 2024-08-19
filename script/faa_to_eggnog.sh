#!/bin/bash

### 用于替换faa文件的序列标识符，>开头到第一个空格之间的内容
###（保证用替换后的faa文件到eggnog-mapper注释后id和counts_cds.txt的id列一致）
# 检查参数数量
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <gene_id_gene_name_file> <faa_file>"
  exit 1
fi

# 获取传递的文件路径
gene_id_gene_name_file="$1"
fasta_file="$2"

# 检查文件是否存在
if [ ! -f "$gene_id_gene_name_file" ]; then
  echo "Error: File '$gene_id_gene_name_file' not found!"
  exit 1
fi

if [ ! -f "$fasta_file" ]; then
  echo "Error: File '$fasta_file' not found!"
  exit 1
fi

# 遍历 gene_id_gene_name_file 中的每一行
while IFS=$'\t' read -r fq1 fq2; do
  # 确保替换操作安全
  sed -i "s/${fq1}/${fq2}/g" "$fasta_file"
done < "$gene_id_gene_name_file"

echo "Replacement complete."

