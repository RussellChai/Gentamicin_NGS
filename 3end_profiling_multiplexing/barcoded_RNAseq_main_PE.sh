## Preprocessing NGS data (quality filter, adaptor removal and mapping) of a RNA-barcode-ligation based 3'-profiling protocol.
## In this protocol, link RNA barcode to the 3'-end of trasncripts (multiplexing if neccessary), RT with a relevant constant primer, 2nd strand synthesis (or template switch), tagmentation and PCR enrich for the RNA-barcode end.
## In this version, demultiplexing is based on barcodes anchored at 5' end of R1; after demultiplexing, 3' end of each R2 file is trimmed accordingly
## Written by Tianmin Wang
## Jintao Liu Lab, Tsinghua University
## Last updated Jun 17, 2022

# Accept arguments and create directory for log files ////////////////////////
sample=$1
mapping_version=$2
mkdir ${sample}_logs/

# Generate required scripts automatically ////////////////////////////////////
python generate_scripts.py adaptor_R2.fasta ${mapping_version} PE

# FastQC /////////////////////////////////////////////////////////////////////
mkdir fastqcLog
ls "${sample}"_*.fastq.gz | time parallel --bar --results fastqcLog -j8 fastqc {}

# Quality trimming ///////////////////////////////////////////////////////////
export PATH=$PATH:~/.local/bin/
cutadapt -q 10,10 --minimum-length 100:100 --max-n 3 --pair-filter=any -o "${sample}"_QF_R1.fastq.gz -p "${sample}"_QF_R2.fastq.gz "${sample}"_R1.fastq.gz "${sample}"_R2.fastq.gz > "${sample}"_logs/QF.log

# Adaptor trimming ///////////////////////////////////////////////////////////
cutadapt --times 1 -a CTGTCTCTTATACACATCTCCGAGCCCACG -a file:N7_barcodes.fasta -a CCCATTCACTCTGCGTTGATACCACTGCTT -G AAGCAGTGGTATCAACGCAGAGTGAATGGG -o "${sample}"_fine_R1.fastq.gz -p "${sample}"_fine_R2.fastq.gz "${sample}"_QF_R1.fastq.gz "${sample}"_QF_R2.fastq.gz > "${sample}"_logs/adaptor_trim.log

# Demultiplexing /////////////////////////////////////////////////////////////
# when the library quality is high, say, each adaptor can be found strictly anchored at the 5' of R1
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -e 0.12 --overlap 9 -g file:adaptor_R1.fasta -o end5-{name}_dem.R1.fastq.gz -p end5-{name}_dem.R2.fastq.gz "${sample}"_fine_R1.fastq.gz "${sample}"_fine_R2.fastq.gz > "${sample}"_logs/demultiplexing.log

# Trim relevant barcodes from 3' end of R2 for each demultiplexed file ///////
bash trim_barcode_from_each_R2.sh ${sample}

# Mapping to genome //////////////////////////////////////////////////////////
# get genome from NCBI RefSeq: NC000913.3: https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3 on Sep. 18, 2019
export PATH=$PATH:~/Documents/wtm/bowtie2-2.3.5.1/
export PATH=$PATH:/usr/local/bin/bin/
bash mapping_"${mapping_version}".sh ${sample}

# Summarize the result of demultiplexing and mapping /////////////////////////
python result_summary.py "${sample}"_logs/

# Count the features from alignments /////////////////////////////////////////
Rscript feature_counting_v2.R ./ PE > "${sample}"_logs/featureCounting.log

# Count depth on each nucleotide of both strands in genome ///////////////////
# Also count enrichment for significant ends across genome ///////////////////
bash count_depth_endenrich.sh

# Compress intermediate .bam files to save disk space ////////////////////////
bash compress_bam.sh

# Run preliminary QC to check the quality of each library ////////////////////
matlab -nodisplay -nosplash -nodesktop -r "run('rainbowseq_QC.m');exit;" | tail -n +11 > "${sample}"_logs/QC_analysis.log

# Calculate the activity of transcription modules ////////////////////////////
#python assess_module_activity.py -i TPM_rm_ncRNA.csv -t genomeFile/network_tf_gene.txt -s genomeFile/network_sigma_gene.txt -n log2 -f 1 -r 3 -o ./

# Plot the gene expression in pathways, defined by pathways_of_interest.csv //
#mkdir pathway_expression/
#Rscript plot_KEGG.R --expression TPM_rm_ncRNA.csv --poi genomeFile/pathways_of_interest.csv --output ./pathway_expression --cutoff 1 --zscore --log

# Relocate intermediate files ////////////////////////////////////////////////
mv "${sample}"_*fastqc* fastqcLog/
mkdir QF_files
mv "${sample}"_QF*.fastq.gz "${sample}"_fine*.fastq.gz end5-*_dem.R?.fastq.gz QF_files/
mkdir demultiplexed_fastq/
mv end5*.fastq.gz demultiplexed_fastq/
mkdir bam_files/
mv lib*.bam.gz lib*.bam.bai sample.csv bam_files/
mkdir rawcounts_tpm/
mv TPM*.csv rawcount.csv rawcounts_tpm/
mv "${sample}"_summary.txt QC_analysis/
mkdir generated_scripts/
mv trim_barcode_from_each_R2.sh mapping_*.sh count_depth_endenrich.sh compress_bam.sh generated_scripts/
mkdir RNA_barcodes_R2_to_be_trimmed/
mv RNA_barcode_R2_adaptor*.fasta RNA_barcodes_R2_to_be_trimmed/
rm -rf QF_files/
