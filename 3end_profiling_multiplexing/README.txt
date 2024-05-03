# Enable demultiplexing from barcodes anchored at 5' end of R1 (RNA barcodes ligated to 3' end of mRNA)
# Two PE files are expected, filePrefix_R1/2.fastq.gz
# v1: global alignment; v2: local alignment
cd /code/program/NGS/3end_profiling_multiplexing/
cp -r *.fasta *.py *.m *.R *.sh genomeFile/ MG1655_bowtie2Index/ README.txt target_path/
cd target_path
bash barcoded_RNAseq_main_P(S)E.sh filePrefix v1(2) 
