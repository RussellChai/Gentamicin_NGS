while getopts "b:o:" OPT;
do
  case $OPT in
      
      b)
        bam=$OPTARG
      ;;
      o)
        out_put=$OPTARG
      
     
  esac
done

sort bam and index
samtools sort $bam -o $bam
samtools index $bam $bam.bai

echo 'bam sorted and indexed'

#depth 
/usr/bin/bedtools2/bin/bedtools genomecov -ibam $bam -d -pc > $out_put

