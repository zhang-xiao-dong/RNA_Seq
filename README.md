# RNA_Seq
Updated on August 19, 2024

## preparation
awk '{if($3=="CDS") print $9}' ../reference/Mycobacterium_tuberculosis_H37Rv.gff3 | awk -F\; '{print $1"\t"$1}' | sed 's/ID=cds-//g' | sed 's/ID=cds-//g' | awk '{print "lcl|NC_000962.3_prot_"$1"_"NR"\tcds-"$2}' > ../inputfile/gene_id_gene_name.txt
cp ../reference/Mycobacterium_tuberculosis_H37Rv.faa ./
bash ../script/faa_to_eggnog.sh ../inputfile/gene_id_gene_name.txt Mycobacterium_tuberculosis_H37Rv.faa    # This replacement faa document is used for the eggnog annotation. 

## Quick start

### 1
cd output/
ls ../06_bam/*.bam > samples.txt
awk -F '/' '{print $3}' samples.txt | awk -F '_' '{print $1}' | awk '{print $1"\t../06_bam/"$1"_sorted.bam"}' > sample_path.txt
awk '{print "../script/run-featurecounts.R -b "$2" -g ../reference/Mycobacterium_tuberculosis_H37Rv.gff3 -o " $1" --isPairedEnd FALSE"}' sample_path.txt > run_Quantification.sh
bash run_Quantification.sh

### 2
../script/abundance_estimates_to_matrix.pl --est_method featureCounts *.count

### 3
Rscript ../script/combine_same_gene.R -c matrix.counts.matrix -g ../inputfile/gene_id_gene_name.txt -T matrix.TMM.EXPR.matrix -o out

### 4
../script/run_DE_analysis.pl --matrix out_combine_counts.matrix --method DESeq2 --samples_file ../inputfile/sampleinfo.txt --contrasts ../inputfile/contrasts.txt

### 5
Rscript ../script/pca_hclust.R -e ../inputfile/eggnog/out.emapper.annotations -T matrix.TMM.EXPR.matrix -s ../inputfile/sample_info.txt

### 6
sed -i '1s/^/id\t/' DESeq2.40467.dir/out_combine_counts.matrix.control_vs_10xMIC.DESeq2.DE_results
Rscript ../script/volcano_heatmap.R -d DESeq2.40467.dir/out_combine_counts.matrix.control_vs_10xMIC.DESeq2.DE_results -e ../inputfile/eggnog/out.emapper.annotations -T out_combine_TMM.EXPR.matrix -s ../inputfile/sample_info.txt

### 7
Rscript ../script/create_orgdb.R -e ../inputfile/eggnog/out.emapper.annotations -t 83332 -g Mycobacterium -s tuberculosisH37Rv -k ../inputfile/ko00001.json

### 8
Rscript ../script/enrich.R -d DESeq2.40467.dir/out_combine_counts.matrix.control_vs_10xMIC.DESeq2.DE_results -e ../inputfile/eggnog/out.emapper.annotations -T out_combine_TMM.EXPR.matrix -o org.MtuberculosisH37Rv.eg.db_0.0.1.tar.gz
