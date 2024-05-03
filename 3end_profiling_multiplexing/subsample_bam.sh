target="subsampled_depth_end_enrichment"
depth="4500000"
mkdir "${target}"_"${depth}"/

for sample in lib5 lib6 lib7 lib8 lib9
do
    mkdir "${target}"_"${depth}"/"${sample}"/

    gunzip bam_files/"${sample}".bam.gz
    Rscript /code/program/bam/sample_bam_PE.R -b bam_files/"${sample}".bam -s 1000 -o "${target}"_"${depth}"/"${sample}"/"${sample}"_sub -n "${depth}"
    gzip bam_files/"${sample}".bam

    bash depth_by_strand.sh -b "${target}"_"${depth}"/"${sample}"/"${sample}"_sub.bam -p "${target}"_"${depth}"/"${sample}"/"${sample}"_pos_depth.csv -n "${target}"_"${depth}"/"${sample}"/"${sample}"_neg_depth.csv
    Rscript end_enrichment_from_bam.R -b "${target}"_"${depth}"/"${sample}"/"${sample}"_sub.bam -p "${target}"_"${depth}"/"${sample}"/"${sample}"_pos_depth.csv -n "${target}"_"${depth}"/"${sample}"/"${sample}"_neg_depth.csv -o "${target}"_"${depth}"/"${sample}"/
    rm -rf "${target}"_"${depth}"/"${sample}"/*_depth.csv
    gzip "${target}"_"${depth}"/"${sample}"/"${sample}"_sub.bam
done
