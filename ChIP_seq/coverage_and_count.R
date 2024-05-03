#This script is to normalize coverage data and feature count for ChIP-seq data.
#Chai Ruochen @ Tsinghua
#2024-01-29

library(tidyverse)
count_cmd = str_c('samtools view -c my_ChIP.bam')
my_total = system(count_cmd,intern = T)
my_total = as.numeric(my_total)
my_cov = read.table('./my_coverage.txt',sep = '\t',header = F)

norm_cov = rbind(data.frame(position = c(my_cov$V2),
           depth = my_cov$V3*1000000/my_total,
           str = '+'),data.frame(position = c(my_cov$V2),
           depth = my_cov$V3*1000000/my_total,
           str = '-'))
write.csv(norm_cov,'./my_norm_coverage.csv',row.names = F)


Bamfilenames <- file.path(c('./my_ChIP.bam'))
file.exists(Bamfilenames)
library("Rsamtools")
Bamfiles <- BamFileList(Bamfilenames)
seqinfo(Bamfiles)

library("GenomicFeatures")
gff3file <- file.path("./genomeFile/NC000913.3.gff3")
txdb_gff <- makeTxDbFromGFF(gff3file, format = "gff3")

transc <- transcriptsBy(txdb_gff, by = "gene")

# Read counting step
library("GenomicAlignments")
library("BiocParallel")

se_ecoli <- summarizeOverlaps(feature = transc, 
                              reads = Bamfiles,                       
                              mode="Union",
                              ignore.strand=T,
                              singleEnd=F)

head(assay(se_ecoli))

write.csv(assay(se_ecoli),'./raw_count.csv',row.names = F)
ChIP_CPM = ((assay(se_ecoli)*1000000)/(apply(assay(se_ecoli),2,sum) %>% as.numeric())) %>% as.data.frame()

ChIP_CPM$gene = row.names(ChIP_CPM)

write.csv(ChIP_CPM,'./CPM.csv',row.names = F)


