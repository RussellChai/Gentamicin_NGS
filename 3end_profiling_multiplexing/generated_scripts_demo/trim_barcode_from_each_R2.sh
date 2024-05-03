sample=$1
mkdir "${sample}"_logs/trim_R2_logs/
# trim R2 for demultiplexed file with adaptor5
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -A file:RNA_barcode_R2_adaptor5.fasta -o end5-adaptor5.R1.fastq.gz -p end5-adaptor5.R2.fastq.gz end5-adaptor5_dem.R1.fastq.gz end5-adaptor5_dem.R2.fastq.gz > "${sample}"_logs/trim_R2_logs/adaptor5_log.txt

# trim R2 for demultiplexed file with adaptor7
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -A file:RNA_barcode_R2_adaptor7.fasta -o end5-adaptor7.R1.fastq.gz -p end5-adaptor7.R2.fastq.gz end5-adaptor7_dem.R1.fastq.gz end5-adaptor7_dem.R2.fastq.gz > "${sample}"_logs/trim_R2_logs/adaptor7_log.txt

# trim R2 for demultiplexed file with adaptor6
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -A file:RNA_barcode_R2_adaptor6.fasta -o end5-adaptor6.R1.fastq.gz -p end5-adaptor6.R2.fastq.gz end5-adaptor6_dem.R1.fastq.gz end5-adaptor6_dem.R2.fastq.gz > "${sample}"_logs/trim_R2_logs/adaptor6_log.txt

# trim R2 for demultiplexed file with adaptor9
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -A file:RNA_barcode_R2_adaptor9.fasta -o end5-adaptor9.R1.fastq.gz -p end5-adaptor9.R2.fastq.gz end5-adaptor9_dem.R1.fastq.gz end5-adaptor9_dem.R2.fastq.gz > "${sample}"_logs/trim_R2_logs/adaptor9_log.txt

# trim R2 for demultiplexed file with adaptor8
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -A file:RNA_barcode_R2_adaptor8.fasta -o end5-adaptor8.R1.fastq.gz -p end5-adaptor8.R2.fastq.gz end5-adaptor8_dem.R1.fastq.gz end5-adaptor8_dem.R2.fastq.gz > "${sample}"_logs/trim_R2_logs/adaptor8_log.txt

