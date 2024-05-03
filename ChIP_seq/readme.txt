Please include those files in the same directory:
MG1655_bowtie2Index
  ¨NNC000913.3.1.bt2 
   NC000913.3.2.bt2 
   NC000913.3.3.bt2 
   NC000913.3.4.bt2 
   NC000913.3.fasta 
   NC000913.3.rev.1.bt2 
   NC000913.3.rev.2.bt2

genomeFile
  ¨Ngene_extension.bed 
   gene_stop.bed 
   NC000913.3.fasta 
   NC000913.3.gff3 
   ncRNA.txt 
   network_sigma_gene.txt 
   network_tf_gene.txt 
   pathways_of_interest.csv
   
ChIP_pipeline.sh 
coverage_and_count.R
I5_barcodes.fasta
N7_barcodes.fasta

sample_name_R1.fq.gz
sample_name_R2.fq.gz 

*Note that please change the names of fq files as abovementioned.

To run:

sh ChIP_pipeline.sh sample_name

 