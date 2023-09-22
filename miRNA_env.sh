#数据库准备
spc=cpa
classify=eudicotyledons
workdir=/share/home/off_wenhao/Work/small_RNA
mature=${workdir}/01.dataBase/data/miRBase/mature.fa
hairpin=${workdir}/01.dataBase/data/miRBase/hairpin.fa
organisms=${workdir}/01.dataBase/data/miRBase/organisms.txt
family=${workdir}/01.dataBase/data/Rfam/family.txt
fastadir=${workdir}/01.dataBase/data/Rfam/fasta_files/
#参考基因组文件准备
ref_genome=${workdir}/01.dataBase/data/ref_genome
script=${workdir}/script
genome=ref.20221224.fasta
gff3=ref.20221224.gff3
cds=ref.20221224.CDS.fasta
#miRNA测序准备
raw_data=${workdir}/00.data/02.cleandata
sample=${workdir}/sample.txt

#样本对应表
#$ cat sample.txt
#15-4F1XYh花6mm①         15-4_XYh-F-6mm_1        15XYh6-1
#15-4F1XYh花6mm②         15-4_XYh-F-6mm_2        15XYh6-2
#15-4F1XYh花6mm③         15-4_XYh-F-6mm_3        15XYh6-3
