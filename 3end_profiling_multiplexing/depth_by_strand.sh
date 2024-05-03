#ChaiRuochen
#2021-05-10
#Last modified: 2022-06-24
#This script is to calculate sequencing depth by strand from a PE mapped bam file.
#


#usage:
#sh depth_by_strand.sh -b <bam file path> -p <depth on positive strand> -n <depth on negative strand>
#Note that positive or negative depends on Read1

while getopts "b:p:n:" OPT;
do
  case $OPT in
      
      b)
        bam=$OPTARG
      ;;
      p)
        pos_depth=$OPTARG
      ;;
      n)
        neg_depth=$OPTARG
     
  esac
done

#sort bam and index
samtools sort $bam -o $bam
samtools index $bam $bam.bai

echo 'bam sorted and indexed'

#Reverse 
/usr/bin/bedtools2/bin/bedtools genomecov -strand - -ibam $bam -d -pc > $neg_depth

#Forward
/usr/bin/bedtools2/bin/bedtools genomecov -strand + -ibam $bam -d -pc > $pos_depth