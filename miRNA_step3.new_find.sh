#!/bin/bash -x
#PBS -N miRNA
#PBS -o miRNA.log
#PBS -e miRNA.err
#PBS -q workq
#PBS -j oe
#PBS -l nodes=5:ppn=4

cd $PBS_O_WORKDIR


source ./env.sh

# New miRNA
echo "running New miRNA ..."
echo "Starting run at $(date) on $(hostname)..."

if [ -d ./17.New_miRNA ]
then
        echo "17.New_miRNA file existed."
else
        mkdir ./17.New_miRNA 
        echo "17.New_miRNA file newly created."
fi


cat $sample | while read ref samp line
do
	## 基于small rna 分类注释信息进行提取
	awk '$NF=="unknown" || $NF == "intron_sense" || $NF=="intron_antisense" ' 16.Class_Notes/${line}.out.read_anno.txt > 17.New_miRNA/${line}.reads.anno.txt
	## 将提取结果转成fa格式
	awk '{print ">"$1"\n"$2}' 17.New_miRNA/${line}.reads.anno.txt > 17.New_miRNA/${line}.reads.fa
done


cat $sample | while read ref samp line
do	
bash /share/home/off_wenhao/biosoft/miRDP2-v1.1.4/1.1.4/miRDP2-v1.1.4_pipeline.bash \
	--genome 13.genome_preparation/genome.fa\
	--index 13.genome_preparation/genome\
	--fasta \
	--input 16.Class_Notes/mapped/${line}.mapped.fa\
	--thread 20 \
	--output 17.New_miRNA

	perl $script/novel_name_mirdp2.pl -pbed 17.New_miRNA/${line}.mapped/${line}.mapped_filter_P_prediction -outpre 17.New_miRNA/${line}.novel_out
done


echo "
C01.novel_out.mature.fa
C01.novel_out.precursor.fa
C01.novel_out.prediction
"
echo "New miRNA complete at $(date) on $(hostname)."



#miRNA sequence analysis
echo "running miRNA sequence analysis ..."
echo "Starting run at $(date) on $(hostname)..."

if [ -d ./18.miRNA_seq_analysis ]
then
        echo "18.miRNA_seq_analysis file existed."
else
        mkdir ./18.miRNA_seq_analysis
        mkdir ./18.miRNA_seq_analysis/PDF
        echo "18.miRNA_seq_analysis file newly created."
fi
cd 18.miRNA_seq_analysis
ln -s ../16.Class_Notes/*.out.read_anno.txt ./
ln -s ../17.New_miRNA/*.novel_out.precursor.fa ./

### 从miRNA分类注释信息中提取已知miRNA并转成fasta格式
cat $sample | while read ref samp line
do
	awk '$NF=="miRNA"' ${line}.out.read_anno.txt | awk '{print ">"$1"\n"$2}' > ${line}.known_mirna.fa
done
### 统计已知miRNA的reads碱基分布
cat $sample | while read ref samp line
do
	perl $script/mirna_base_content.pl ${line}.known_mirna.fa ${line}.known
done

echo "
prefix.known.each_pos.total.txt : 各个位置碱基分布（total）
prefix.known.each_pos.uniq.txt : 各个位置碱基分布（unique）
prefix.known.first_base.total.txt : 首位碱基分布（total）
prefix.known.first_base.uniq.txt ：首位碱基分布（unique）
"
##对来源于新预测的 miRNA 的 srna 数据进行提取和汇总统计
## 新miRNA碱基分析
cat $sample | while read ref samp line
do
	awk '$NF=="unknown" || $NF == "intron_sense" || $NF=="intron_antisense"' ${line}.out.read_anno.txt | awk '{print ">"$1"\n"$2}' >  ${line}.reads_for_novel.fa


### 构建新miRNA前体bowtie index
bowtie-build ${line}.novel_out.precursor.fa ${line}.novel_out.precursor
done
### 将srna数据比对到新miRNA前体
cat $sample | while read ref samp line
do
	bowtie -f -v 0 -p 1 -a --best --strata\
	--al ${line}.novel_miRNA.fa ${line}.novel_out.precursor ${line}.reads_for_novel.fa > ${line}.novel_miRNA.bwt 2>${line}.novel_miRNA.log
done
### 统计新miRNA的reads碱基分布
cat $sample | while read ref samp line
do
	perl $script/mirna_base_content.pl ${line}.novel_miRNA.fa ./PDF/${line}.novel
done

### 绘制miRNA各个位置碱基分布图
cd ./PDF
ls ./*.each_pos.*.txt | while read aa ;
do
	Rscript $script/draw_base_content.R $aa $aa.pdf ;
done
cd ../../

