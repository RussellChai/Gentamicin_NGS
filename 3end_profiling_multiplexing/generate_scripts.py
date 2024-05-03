# this script is used to generate bash scripts for trimming 3' ends of R2 to remove potential barcodes
# this script is used to generate bash scripts for mapping (bowtie2) multiple demultiplexed .fastq file
# this script generates another flat file "sample.csv", fed to the R scripts that perform feature counting
# this script generates bash scripts to count depth on both strands of the genome
# this script generates bash scripts to count 3' end enrichment 
# finally this script generates a bash scipt to compress the intermediate .bam file to save disk space
## Written by Tianmin Wang
## Jintao Liu Lab, Tsinghua University
## Last updated Jun 17, 2022

import os
import sys

# arguments processing
# ///////////////////////////////////////////////////////////////////////////////////////////
barcodes_file = sys.argv[1]
mapping_version = sys.argv[2]
mapping_type = sys.argv[3]
if mapping_version not in ['v1', 'v2']:
    print('incorrect mapping version, v1 or v2 is expected!')
    sys.exit(-1)

if mapping_type not in ['SE', 'PE']:
    print('incorrect mapping type, SE or PE is expected!')
    sys.exit(-1)

# find all barcodes
barcodes = []
R2_end3_barcodes = {}
f = open(barcodes_file, 'r')
for line in f:
    if line[0] == '>':
        barcode = line.rstrip()[1:]
        if barcode[:7] == 'adaptor':
            if barcode[7:] not in barcodes:
                barcodes.append(barcode[7:])
        else:
            print('incorrect barcode name for %s in %s'%(barcode, barcodes_file))
            sys.exit(-1)
    else:
        R2_end3_barcodes[barcode] = line.rstrip()
f.close()

# generate the shell script used to trim 3' ends of each demutiplexed R2 file
# ///////////////////////////////////////////////////////////////////////////////////////////
def generate_series_to_be_trimmed(barcode, barcode_seq):
    # function to generate a series of sequences to be trimmed, and store them in one file named after 'RNA_barcode_R2_adaptor?.fasta'
    os.system('cat /dev/null > RNA_barcode_R2_%s.fasta'%(barcode))
    g = open('RNA_barcode_R2_%s.fasta'%(barcode), 'r+')
    for i in range(len(barcode_seq)):
        truncated_adaptor = barcode_seq[:(len(barcode_seq) - i)] + '$'
        g.write('>RNA_%s_%d\n%s\n'%(barcode, i+1, truncated_adaptor))
    g.close()
    return 'RNA_barcode_R2_%s.fasta'%(barcode)

os.system('cat /dev/null > trim_barcode_from_each_R2.sh')
k = open('trim_barcode_from_each_R2.sh', 'r+')
k.write('sample=$1\n')
k.write('mkdir "${sample}"_logs/trim_R2_logs/\n')
for barcode in R2_end3_barcodes:
    barcode_seq = R2_end3_barcodes[barcode]
    specific_barcode_file = generate_series_to_be_trimmed(barcode, barcode_seq)
    k.write('# trim R2 for demultiplexed file with %s\n'%(barcode))
    k.write('cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -A file:%s -o end5-%s.R1.fastq.gz -p end5-%s.R2.fastq.gz end5-%s_dem.R1.fastq.gz end5-%s_dem.R2.fastq.gz > "${sample}"_logs/trim_R2_logs/%s_log.txt\n\n'%(specific_barcode_file, barcode, barcode, barcode, barcode, barcode))
k.close()

# generate the shell script used in mapping (bowtie2)
# ///////////////////////////////////////////////////////////////////////////////////////////
if mapping_type == 'SE':
    os.system('cat /dev/null > mapping_R1_%s.sh'%(mapping_version))
    g = open('mapping_R1_%s.sh'%(mapping_version), 'r+')
else:
    os.system('cat /dev/null > mapping_%s.sh'%(mapping_version))
    g = open('mapping_%s.sh'%(mapping_version), 'r+')

g.write('sample=$1\n')
g.write('mkdir "${sample}"_logs/bowtie2_logs/\n')
g.write('\n')

def mapping_command(mapping_version, mapping_type, barcode):
    # function to generate mapping command for each line
    this_command = ''
    if mapping_version == 'v1':
        if mapping_type == 'SE':
            this_command = this_command + '(bowtie2 -x MG1655_bowtie2Index/NC000913.3 -U end5-adaptor%s.R1.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode)            
        else: #PE
            this_command = this_command + '(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor%s.R1.fastq.gz -2 end5-adaptor%s.R2.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode, barcode)
    elif mapping_version == 'v2':
        if mapping_type == 'SE':
            this_command = this_command + '(bowtie2 --very-sensitive-local -N 1 -L 16 --no-1mm-upfront --score-min G,9,8 -p 1 -x MG1655_bowtie2Index/NC000913.3 -U end5-adaptor%s.R1.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode)
        else: # PE
            this_command = this_command + '(bowtie2 --very-sensitive-local -N 1 -L 16 --no-1mm-upfront --score-min G,9,8 -p 1 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor%s.R1.fastq.gz -2 end5-adaptor%s.R2.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode, barcode)
    this_command = this_command + 'samtools view -S -b lib%s.sam > lib%s.bam\n'%(barcode, barcode)
    this_command = this_command + 'rm -rf lib%s.sam\n'%(barcode)
    return this_command

for barcode in barcodes:
    this_command = mapping_command(mapping_version, mapping_type, barcode)
    g.write(this_command + '\n')
g.close()

# generate .csv file used in feature counting
# ///////////////////////////////////////////////////////////////////////////////////////////
os.system('cat /dev/null > sample.csv')
g = open('sample.csv', 'r+')
g.write('sample\n')
for barcode in barcodes: 
    g.write('lib%s\n'%(barcode))
g.close()

# generate shell script used in depth and end enrichment counting
# ///////////////////////////////////////////////////////////////////////////////////////////
os.system('cat /dev/null > count_depth_endenrich.sh')
g = open('count_depth_endenrich.sh', 'r+')
g.write('mkdir depth_end_enrichment/\n\n')
for barcode in barcodes:
    g.write('mkdir depth_end_enrichment/lib%s/\n'%(barcode))
    forward_depth_file = 'depth_end_enrichment/lib%s/lib%s_pos_depth.csv'%(barcode, barcode)
    reverse_depth_file = 'depth_end_enrichment/lib%s/lib%s_neg_depth.csv'%(barcode, barcode)
    depth_command = 'bash depth_by_strand.sh -b lib%s.bam -p %s -n %s'%(barcode, forward_depth_file, reverse_depth_file) 
    g.write(depth_command + '\n')
    end_enrichment_command = 'Rscript end_enrichment_from_bam.R -b lib%s.bam -p %s -n %s -o depth_end_enrichment/lib%s/'%(barcode, forward_depth_file, reverse_depth_file, barcode)
    g.write(end_enrichment_command + '\n')
    g.write('rm -rf %s %s\n\n'%(forward_depth_file, reverse_depth_file))
g.close()

# generate bash script to compress .bam files
# ///////////////////////////////////////////////////////////////////////////////////////////
os.system('cat /dev/null > compress_bam.sh')
g = open('compress_bam.sh', 'r+')
for barcode in barcodes: 
    g.write('gzip lib%s.bam\n'%(barcode))
g.close()
