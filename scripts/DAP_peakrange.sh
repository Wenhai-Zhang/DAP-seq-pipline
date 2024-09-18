wd="/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906"

mkdir -p $wd/result/peakranger
# samtools view -b -h $wd/result/alignment_BAM/N2-6_L2_filtered_sorted.bam HiC_scaffold_1 HiC_scaffold_2 \
# HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_8 > $wd/result/peakranger/N2-6_L2_filtered_sorted.bam
# samtools view -b -h $wd/result/alignment_BAM/N-W-5_L2_filtered_sorted.bam HiC_scaffold_1 HiC_scaffold_2 \
# HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_8 > $wd/result/peakranger/N-W-5_L2_filtered_sorted.bam
# samtools sort -n $wd/result/peakranger/N2-6_L2_filtered_sorted.bam -o $wd/result/peakranger/N2-6_L2_filtered_sorted_sortedn.bam
# samtools sort -n $wd/result/peakranger/N-W-5_L2_filtered_sorted.bam -o $wd/result/peakranger/N-W-5_L2_filtered_sorted_sortedn.bam
# peakranger nr \
# -d $wd/result/peakranger/N2-6_L2_filtered_sorted_sortedn.bam \
# -c $wd/result/peakranger/N-W-5_L2_filtered_sorted_sortedn.bam \
# -l 300 \
# --format bam 
# peakranger lc \
# -d $wd/result/peakranger/N2-6_L2_filtered_sorted_sortedn.bam 
# peakranger lc \
# -d $wd/result/peakranger/N-W-5_L2_filtered_sorted_sortedn.bam 

peakranger ranger \
-d $wd/result/peakranger/N2-6_L2_filtered_sorted_sortedn.bam \
-c $wd/result/peakranger/N-W-5_L2_filtered_sorted_sortedn.bam \
-o $wd/result/peakranger/ranger_result \
--format bam \
--report --gene_annot_file /home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/data/reference_genome/R108_HiC/20210817_MtrunR108_HiC_FunctionalAnnotation/annotation/MedtrR108_HiC.gtf