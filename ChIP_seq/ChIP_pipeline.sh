#This pipeline is for processing ChIP-seq data.
#Usage: sh ChIP_pipeline.sh sample_name
#Chai Ruochen @ Tsinghua
#2024-01-29

sample=$1
mkdir ${sample}_logs/
# FastQC /////////////////////////////////////////////////////////////////////
mkdir fastqcLog
ls "${sample}"_*.fq.gz | time parallel --bar --results fastqcLog -j8 fastqc {}

# Quality trimming ///////////////////////////////////////////////////////////
export PATH=$PATH:~/.local/bin/
cutadapt -q 10,10 --minimum-length 100:100 --max-n 3 --pair-filter=any -o "${sample}"_QF_R1.fq.gz -p "${sample}"_QF_R2.fq.gz "${sample}"_R1.fq.gz "${sample}"_R2.fq.gz > "${sample}"_logs/QF.log


# Adaptor trimming N7 /////////////////////////////////////////////////////////
cutadapt --times 1 -a file:N7_barcodes.fasta -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o "${sample}"_fine_R1.fq.gz "${sample}"_QF_R1.fq.gz > "${sample}"_logs/adaptor_N7_trim.log

# Adaptor trimming I5 /////////////////////////////////////////////////////////
cutadapt --times 1 -a file:I5_barcodes.fasta -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -o "${sample}"_fine_R2.fq.gz "${sample}"_QF_R2.fq.gz > "${sample}"_logs/adaptor_I5_trim.log

# Mapping by paired-end //////////////////////////////////////////////////////
mkdir "${sample}"_logs/bowtie2_logs/

(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 "${sample}"_fine_R1.fq.gz -2 "${sample}"_fine_R2.fq.gz -S my_ChIP.sam)2>"${sample}"_logs/bowtie2_logs/ChIP_log.txt
samtools view -S -b my_ChIP.sam > my_ChIP.bam
rm -rf my_ChIP.sam

# Calculate depth
sh depth.sh -b my_ChIP.bam -o my_coverage.txt

#Normalize coverage and featureCounting
Rscript coverage_and_count.R