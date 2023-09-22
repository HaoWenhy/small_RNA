#!/bin/bash -x
#PBS -N miRNA
#PBS -o miRNA.log
#PBS -e miRNA.err
#PBS -q workq
#PBS -j oe
#PBS -l nodes=1:ppn=10

cd $PBS_O_WORKDIR

if [ -d ./19.miRNA_family ]
then
        echo "19.miRNA_family file existed."
else
        mkdir ./19.miRNA_family
        echo "19.miRNA_family file newly created."
fi

kn_hp=/share/home/off_wenhao/Work/small_RNA/11.ref_DB/hairpin.cpa.fa
nv_hp=/share/home/off_wenhao/Work/small_RNA/17.New_miRNA

cat $kn_hp $nv_hp/*.novel_out.precursor.fa > ./19.miRNA_family/all.hairpin.fa

blastn -db  ./12.ncRNA/Rfam_MIR.fasta -query ./19.miRNA_family/all.hairpin.fa -out ./19.miRNA_family/all.hairpin.Rfam_MIR.blast -outfmt 6 -evalue 0.01

perl $script/blastRfam2family.pl ./19.miRNA_family/all.hairpin.Rfam_MIR.blast ./12.ncRNA/family.txt.MIR > ./19.miRNA_family/all.hairpin.family.out

echo "miRNA_family complete at $(date) on $(hostname)."



if [ -d ./20.miRNA_experssion  ]
then
        echo "20.miRNA_experssion file existed."
else
        mkdir ./20.miRNA_experssion
        echo "20.miRNA_experssion file newly created."
fi

cd ./20.miRNA_experssion 
kn_mat=/share/home/off_wenhao/Work/small_RNA/11.ref_DB/mature.cpa.fa
kn_hp=/share/home/off_wenhao/Work/small_RNA/11.ref_DB/hairpin.cpa.fa
nv_hp=/share/home/off_wenhao/Work/small_RNA/17.New_miRNA

cat $kn_hp $nv_hp/*.novel_out.precursor.fa > all.hairpin.fa
cat $kn_hp $nv_hp/*.novel_out.mature.fa > all.mature.fa

cat ../14.fastqc/uniq_data/*.fa > all.reads.fa

/share/home/off_wenhao/biosoft/mirdeep2-0.1.3/bin/quantifier.pl -p all.hairpin.fa -m all.mature.fa -r all.reads.fa -g 0

less miRNAs_expressed_all_samples_*.csv | awk '{ for(i=5; i<=4+(NF-4)/2;i++){ c = c"\t"$i } ;print $1 c; c=""}'| awk '{if(NR==1 || !($1 in A) ){print $0}; A[$1] = 1}' > miRNAs_expressed.count.txt

less miRNAs_expressed_all_samples_*.csv | awk '{ for(i=5+(NF-4)/2 ; i<=NF;i++){ c = c"\t"$i }; print $1 c; c=""}'| awk '{if(NR==1 || !($1 in A) ){print $0} ; A[$1] = 1}' | sed 's/(norm)//g' > miRNAs_expressed.TPM.txt

echo "
15-4_XYh-F-6mm	C01
15-4_XYh-F-6mm	C02
15-4_XYh-F-6mm	C03
15-4_XYh-F-8mm	C04
15-4_XYh-F-8mm	C05
15-4_XYh-F-8mm	C06
15-4_XY-F-6mm	C07
15-4_XY-F-6mm	C08
15-4_XY-F-6mm	C09
15-4_XY-F-8mm	C10
15-4_XY-F-8mm	C11
15-4_XY-F-8mm	C12
YhYh-F-8mm	C13
YhYh-F-8mm	C14
YhYh-F-8mm	C15
15-4_YhYh_6mm	C16
15-4_YhYh_6mm	C17
15-4_YhYh_6mm	C18
20-1_XY_6mm	C19
20-1_XY_6mm	C20
20-1_XY_6mm	C21
20-1_XY_8mm	C22
20-1_XY_8mm	C23
20-1_XY_8mm	C24
20-1_XYh_6mm	C25
20-1_XYh_6mm	C26
20-1_XYh_6mm	C27
20-1_XYh_8mm	C28
20-1_XYh_8mm	C29
20-1_XYh_8mm	C30
AU9XY-F-6mm	C32
AU9XY-F-6mm	C33
AU9XY-F-6mm	C34
AU9XY-F-8mm	C35
AU9XY-F-9mm	C36
AU9XY-F-8mm	C37
SunupXYh-F-8mm	C38
SunupXYh-F-8mm	C39
SunupXYh-F-8mm	C40
Female_flower_replicate	C41
Female_flower_replicate	C42
Female_flower_replicate	C43
Male_flower_replicate	C46
Male_flower_replicate	C47
Male_flower_replicate	C48" > sample.txt

Rscript  $script/expression_draw.R miRNAs_expressed.count.txt miRNAs_expressed.TPM.txt sample.txt mirna_exp.TPM

echo "
15-4_XY-F-6mm	15-4_XY-F-8mm
15-4_XYh-F-6mm	15-4_XYh-F-8mm
AU9XY-F-6mm	AU9XY-F-8mm
15-4_YhYh_6mm	YhYh-F-8mm
20-1_XY_6mm	20-1_XY_8mm
20-1_XYh_6mm	20-1_XYh_8mm
15-4_XY-F-6mm	15-4_XYh-F-6mm
15-4_XY-F-8mm	15-4_XYh-F-8mm
20-1_XY_6mm	20-1_XYh_6mm
20-1_XY_8mm	20-1_XYh_8mm
AU9XY-F-6mm	15-4_YhYh_6mm
15-4_XYh-F-8mm	SunupXYh-F-8mm
20-1_XYh_8mm	SunupXYh-F-8mm
Female_flower_replicate	Male_flower_replicate
" >contrasts.txt

##进行差异表达分析
perl $script/run_DE_analysis.pl --matrix  miRNAs_expressed.count.txt --samples_file  sample.txt --contrasts contrasts.txt --output DE_result  --method DESeq2 

ls DE_result/*.DE_results | while read line
do
	awk 'NR>1 && ($5>1||$5<-1) && $9<0.05{print $1}' $line > $line.gene 
done

