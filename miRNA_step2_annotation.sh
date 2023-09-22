#!/bin/bash -x
#PBS -N miRNA
#PBS -o miRNA.log
#PBS -e miRNA.err
#PBS -q workq
#PBS -j oe
#PBS -l nodes=1:ppn=10

cd $PBS_O_WORKDIR


source ./env.sh

# Data filtration
#echo "running Data filtration ..."
#echo "Starting run at $(date) on $(hostname)..."

#if [ -d ./14.fastqc ]
#then
#        echo "14.fastqc file existed."
#else
#       mkdir ./14.fastqc
#    	mkdir ./14.fastqc/raw_data
#    	mkdir ./14.fastqc/clean_data
#    	mkdir ./14.fastqc/uniq_data    
#       echo "14.fastqc file newly created."
#fi
#
#
## 将测序数据链接到raw_data目录
#ln -s $raw_data/*.gz 14.fastqc/raw_data/
#
## 使用TrimGalore进行small RNA数据过滤
#cat $sample | while read ref samp line
#do
#	trim_galore --small_rna --length 18 --max_length 30 --stringency 3 --phred33 --cores 4 --dont_gzip  \
#				-o 14.fastqc/clean_data 14.fastqc/raw_data/${samp}.fastq.gz
#done

#echo "
#clean_data/P1-1_trimmed.fq sequnces
#clean_data/P1-1.fq.gz_trimming_report.txt report
#"
##将small RNA 从 fq 转成无冗余的 fa 格式，从 line 部分反映每种序列出现的次数

#cat $sample | while read ref samp line
#do
#	mapper.pl 14.fastqc/clean_data/${samp}_trimmed.fq -e -g ${line} -h -m -s 14.fastqc/uniq_data/${line}.fa
#done

#echo "去重合并的 fasta 文件，用来做后续的 small RNA 分析。"

#echo "Data filtration complete at $(date) on $(hostname)."



# Stats
echo "running Stats ..."
echo "Starting run at $(date) on $(hostname)..."

if [ -d ./15.Stats ]
then
        echo "15.Stats file existed."
else
        mkdir ./15.Stats 
	mkdir -p ./15.Stats/Venn
        echo "15.Stats file newly created."
fi

## 数据量及长度统计

cd 15.Stats
cat $sample | while read ref samp line
do
	perl $script/stat_srna.pl ../14.fastqc/uniq_data/${line}.fa ${line}
done

cd ../
echo "## 生成文件如下：
					prefix.readstat  reads  :条数及碱基数统计
					prefix.len.total.txt 	:total reads 长度统计文件
					prefix.len.uniq.txt 	:unique reads 长度统计"


# 绘制长度分布图
cat $sample | while read ref samp line
do
	Rscript $script/draw_length_srna.R 15.Stats/${line}.len.total.txt 15.Stats/${line}.len.total
	Rscript $script/draw_length_srna.R 15.Stats/${line}.len.uniq.txt  15.Stats/${line}.len.uniq 
done

## 将处理好的fasta文件写成文件列表
cd 14.fastqc
ls uniq_data/*.fa |awk -F "/" '{print $2"\t"$0}' | sed 's/.fa//' > clean_srna.list

# 汇总每种序列在各个样品中的出现次数
perl $script/common_specific_reads.pl clean_srna.list > clean_srna.count.tsv
cd ../
## 绘制两两样品共有和特有序列的venn图
Rscript $script/draw_srna_pairwise_venn.R 14.fastqc/clean_srna.count.tsv 15.Stats/Venn/Venn

## 将单个样品长度统计表合并，进行绘图
cat 15.Stats/*.len.total.txt > 15.Stats/all.len.total.txt
cat 15.Stats/*.len.uniq.txt > 15.Stats/all.len.uniq.txt

Rscript $script/draw_length_srna.R 15.Stats/all.len.total.txt 15.Stats/all.len.total.txt
Rscript $script/draw_length_srna.R 15.Stats/all.len.uniq.txt 15.Stats/all.len.uniq.txt

echo "Stats complete at $(date) on $(hostname)."

# Class Notes
echo "running Class Notes ..."
echo "Starting run at $(date) on $(hostname)..."

if [ -d ./16.Class_Notes ]
then
        echo "16.Class_Notes file existed."
else
        mkdir ./16.Class_Notes
        mkdir ./16.Class_Notes/genome
        mkdir ./16.Class_Notes/mapped
        mkdir ./16.Class_Notes/miRNA
        mkdir ./16.Class_Notes/ncRNA
        mkdir ./16.Class_Notes/repeat
        mkdir ./16.Class_Notes/exon
        mkdir ./16.Class_Notes/intron
        echo "16.Class_Notes file newly created."
fi

cat $sample | while read ref samp line
do
        ## reads比对基因组
       bowtie -f -v 0 -p 10 -k 1 --al 16.Class_Notes/mapped/${line}.mapped.fa 13.genome_preparation/genome 14.fastqc/uniq_data/${line}.fa > 16.Class_Notes/genome/${line}.genome.bwt 2>16.Class_Notes/genome/${line}.genome.log
        ## 对比对结果进行统计
       perl $script/stat_bwt.pl 14.fastqc/uniq_data/${line}.fa 16.Class_Notes/genome/${line}.genome.bwt > 16.Class_Notes/genome/${line}.genome.stat.tsv
        ## 统计比对上的reads数目及长度
       perl $script/stat_srna.pl 16.Class_Notes/mapped/${line}.mapped.fa 16.Class_Notes/mapped/${line}.mapped
done

##进行分类注释
cat $sample | while read ref samp line
do
        ## 比对本物种miRNA前体数据库
       bowtie -f -v 0 -p 10 -a --best --strata 11.ref_DB/hairpin.$spc 	    16.Class_Notes/mapped/${line}.mapped.fa > 16.Class_Notes/miRNA/${line}.miRNA.bwt  2>16.Class_Notes/miRNA/${line}.miRNA.log
        ## 比对ncRNA数据库
       bowtie -f -v 2 -p 10 -a --best --strata 12.ncRNA/Rfam_rmMIR.fasta    16.Class_Notes/mapped/${line}.mapped.fa > 16.Class_Notes/ncRNA/${line}.ncRNA.bwt   2>16.Class_Notes/ncRNA/${line}.ncRNA.log
        ## 比对本物种repeat数据库
       bowtie -f -v 0 -p 10 -a --best --strata 13.genome_preparation/repeat 16.Class_Notes/mapped/${line}.mapped.fa > 16.Class_Notes/repeat/${line}.repeat.bwt 2>16.Class_Notes/repeat/${line}.repeat.log
        ## 比对本物种mRNA序列
       bowtie -f -v 0 -p 10 -a --best --strata 13.genome_preparation/exon   16.Class_Notes/mapped/${line}.mapped.fa > 16.Class_Notes/exon/${line}.exon.bwt      2>16.Class_Notes/exon/${line}.exon.log
        ## 比对本物种intron序列
       bowtie -f -v 0 -p 10 -a --best --strata 13.genome_preparation/intron 16.Class_Notes/mapped/${line}.mapped.fa > 16.Class_Notes/intron/${line}.intron.bwt 2>16.Class_Notes/intron/${line}.intron.log
done

## 对small rna进行分类注释
cat $sample | while read ref samp line
do
       perl $script/srna_anno.pl \
                 -fa 16.Class_Notes/mapped/${line}.mapped.fa\
                 -rfam 12.ncRNA/family.txt1\
                 -mirna 16.Class_Notes/miRNA/${line}.miRNA.bwt\
                 -ncrna 16.Class_Notes/ncRNA/${line}.ncRNA.bwt\
                 -intron 16.Class_Notes/intron/${line}.intron.bwt\
                 -repeat 16.Class_Notes/repeat/${line}.repeat.bwt\
                 -exon   16.Class_Notes/exon/${line}.exon.bwt\
                 -outpre 16.Class_Notes/${line}.out

       Rscript $script/draw_srna_annPie.R 16.Class_Notes/${line}.out.read_anno.stat 16.Class_Notes/${line}.out.read_anno
done

#echo "Class Notes complete at $(date) on $(hostname)."

