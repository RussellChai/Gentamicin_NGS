mkdir depth_end_enrichment/

mkdir depth_end_enrichment/lib5/
bash depth_by_strand.sh -b lib5.bam -p depth_end_enrichment/lib5/lib5_pos_depth.csv -n depth_end_enrichment/lib5/lib5_neg_depth.csv
Rscript end_enrichment_from_bam.R -b lib5.bam -p depth_end_enrichment/lib5/lib5_pos_depth.csv -n depth_end_enrichment/lib5/lib5_neg_depth.csv -o depth_end_enrichment/lib5/
rm -rf depth_end_enrichment/lib5/lib5_pos_depth.csv depth_end_enrichment/lib5/lib5_neg_depth.csv

mkdir depth_end_enrichment/lib6/
bash depth_by_strand.sh -b lib6.bam -p depth_end_enrichment/lib6/lib6_pos_depth.csv -n depth_end_enrichment/lib6/lib6_neg_depth.csv
Rscript end_enrichment_from_bam.R -b lib6.bam -p depth_end_enrichment/lib6/lib6_pos_depth.csv -n depth_end_enrichment/lib6/lib6_neg_depth.csv -o depth_end_enrichment/lib6/
rm -rf depth_end_enrichment/lib6/lib6_pos_depth.csv depth_end_enrichment/lib6/lib6_neg_depth.csv

mkdir depth_end_enrichment/lib7/
bash depth_by_strand.sh -b lib7.bam -p depth_end_enrichment/lib7/lib7_pos_depth.csv -n depth_end_enrichment/lib7/lib7_neg_depth.csv
Rscript end_enrichment_from_bam.R -b lib7.bam -p depth_end_enrichment/lib7/lib7_pos_depth.csv -n depth_end_enrichment/lib7/lib7_neg_depth.csv -o depth_end_enrichment/lib7/
rm -rf depth_end_enrichment/lib7/lib7_pos_depth.csv depth_end_enrichment/lib7/lib7_neg_depth.csv

mkdir depth_end_enrichment/lib8/
bash depth_by_strand.sh -b lib8.bam -p depth_end_enrichment/lib8/lib8_pos_depth.csv -n depth_end_enrichment/lib8/lib8_neg_depth.csv
Rscript end_enrichment_from_bam.R -b lib8.bam -p depth_end_enrichment/lib8/lib8_pos_depth.csv -n depth_end_enrichment/lib8/lib8_neg_depth.csv -o depth_end_enrichment/lib8/
rm -rf depth_end_enrichment/lib8/lib8_pos_depth.csv depth_end_enrichment/lib8/lib8_neg_depth.csv

mkdir depth_end_enrichment/lib9/
bash depth_by_strand.sh -b lib9.bam -p depth_end_enrichment/lib9/lib9_pos_depth.csv -n depth_end_enrichment/lib9/lib9_neg_depth.csv
Rscript end_enrichment_from_bam.R -b lib9.bam -p depth_end_enrichment/lib9/lib9_pos_depth.csv -n depth_end_enrichment/lib9/lib9_neg_depth.csv -o depth_end_enrichment/lib9/
rm -rf depth_end_enrichment/lib9/lib9_pos_depth.csv depth_end_enrichment/lib9/lib9_neg_depth.csv

