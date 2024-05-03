sample=$1
mkdir "${sample}"_logs/bowtie2_logs/

(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor5.R1.fastq.gz -2 end5-adaptor5.R2.fastq.gz -S lib5.sam)2>"${sample}"_logs/bowtie2_logs/adaptor5_log.txt
samtools view -S -b lib5.sam > lib5.bam
rm -rf lib5.sam

(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor6.R1.fastq.gz -2 end5-adaptor6.R2.fastq.gz -S lib6.sam)2>"${sample}"_logs/bowtie2_logs/adaptor6_log.txt
samtools view -S -b lib6.sam > lib6.bam
rm -rf lib6.sam

(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor7.R1.fastq.gz -2 end5-adaptor7.R2.fastq.gz -S lib7.sam)2>"${sample}"_logs/bowtie2_logs/adaptor7_log.txt
samtools view -S -b lib7.sam > lib7.bam
rm -rf lib7.sam

(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor8.R1.fastq.gz -2 end5-adaptor8.R2.fastq.gz -S lib8.sam)2>"${sample}"_logs/bowtie2_logs/adaptor8_log.txt
samtools view -S -b lib8.sam > lib8.bam
rm -rf lib8.sam

(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor9.R1.fastq.gz -2 end5-adaptor9.R2.fastq.gz -S lib9.sam)2>"${sample}"_logs/bowtie2_logs/adaptor9_log.txt
samtools view -S -b lib9.sam > lib9.bam
rm -rf lib9.sam

