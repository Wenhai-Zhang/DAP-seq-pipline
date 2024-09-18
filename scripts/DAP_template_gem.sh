wd="/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906"

input_dir=$wd/data/raw_sequence_data/240902-A00133B
samples=("N-W-5_L2" "N2-6_L2")

# mkdir -p $wd/gem_result/reference/scaffold
# sed 's/HiC_scaffold_/chr/g' $wd/data/reference_genome/R108_HiC/MtrunR108_HiC.refseq.fasta > $wd/gem_result/reference/MtrunR108_HiC.refseq.fasta
# python size.py


# for sample_name in "${samples[@]}"
# do
#     read1=${sample_name}_1.fq.gz
#     read2=${sample_name}_2.fq.gz
#     # echo "============================================ adapter trimming ============================================"
#     # echo "----------------------------fastp----------------------------"
#     # mkdir -p $wd/gem_result/adapter_trimming/fastp
#     # fastp -i $input_dir/$read1 -I $input_dir/$read2 \
#     # -o $wd/gem_result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz} \
#     # -O $wd/gem_result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz} \
#     # --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     # --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     # -h $wd/gem_result/adapter_trimming/fastp/fastp.html -j $wd/gem_result/adapter_trimming/fastp/fastp.json -w 16

#     # echo "============================================ fastqc ============================================"
#     # # mkdir -p $wd/gem_result/fastqc/origin
#     # # fastqc -o $wd/gem_result/fastqc/origin -t 24 $input_dir/$read1
#     # # fastqc -o $wd/gem_result/fastqc/origin -t 24 $input_dir/$read2
#     # mkdir -p $wd/gem_result/fastqc/fastp
#     # fastqc -o $wd/gem_result/fastqc/fastp -t 24 $wd/gem_result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz}
#     # fastqc -o $wd/gem_result/fastqc/fastp -t 24 $wd/gem_result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz}

#     # echo "============================================ build bowtie index ============================================"

#     echo "============================================ bowtie2 read alignment ============================================"
#     # default we use fastp result
#     mkdir -p $wd/gem_result/alignment_SAM
#     mkdir -p $wd/gem_result/QC
#     bowtie2 --very-sensitive -p 64 \
#     -x $wd/gem_result/bowtie2/R108_HiC \
#     -1 $wd/result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz} \
#     -2 $wd/result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz} \
#     -S $wd/gem_result/alignment_SAM/${sample_name}.sam 2>&1 | tee $wd/gem_result/QC/${sample_name}_alignment.log

#     echo "convert SAM file to BAM file"
#     mkdir -p $wd/gem_result/alignment_BAM
#     samtools view -bS $wd/gem_result/alignment_SAM/${sample_name}.sam > $wd/gem_result/alignment_BAM/${sample_name}.bam
#     echo "sort bam file"
#     samtools sort -@ 16 -O BAM -o $wd/gem_result/alignment_BAM/${sample_name}_sorted.bam $wd/gem_result/alignment_BAM/${sample_name}.bam

#     echo "index bam file"
#     samtools index $wd/gem_result/alignment_BAM/${sample_name}_sorted.bam

#     echo "flagstat"
#     samtools flagstat $wd/gem_result/alignment_BAM/${sample_name}.bam > $wd/gem_result/alignment_BAM/${sample_name}_flagstat.txt

#     samtools addreplacerg -r "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:Illumina\tLB:${sample_name}" -o $wd/gem_result/alignment_BAM/${sample_name}_addrg.bam $wd/gem_result/alignment_BAM/${sample_name}_sorted.bam

#     echo "Mark duplicates"
#     export _JAVA_OPTIONS='-Xms10G -Xmx15G'
#     picard MarkDuplicates \
#     QUIET=true \
#     INPUT=$wd/gem_result/alignment_BAM/${sample_name}_addrg.bam \
#     OUTPUT=$wd/gem_result/alignment_BAM/${sample_name}_dupRemoved.bam \
#     METRICS_FILE=$wd/gem_result/QC/${sample_name}_dupMarked.sorted.metrics \
#     REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
#     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#     echo "Remove multi-mapped reads"
#     samtools view -h -q 20 $wd/gem_result/alignment_BAM/${sample_name}_dupRemoved.bam > $wd/gem_result/alignment_BAM/${sample_name}_mapQ20.bam

#     echo "Remove reads unmapped 4, mate unmapped 8, not primary alignment 256, reads failing platform 512"
#     samtools view -h -b -F 780 -f 2 $wd/gem_result/alignment_BAM/${sample_name}_mapQ20.bam > $wd/gem_result/alignment_BAM/${sample_name}_filtered.bam

#     echo "sort bam file"
#     samtools sort -@ 16 -O BAM -o $wd/gem_result/alignment_BAM/${sample_name}_filtered_sorted.bam $wd/gem_result/alignment_BAM/${sample_name}_filtered.bam

#     echo "index bam file"
#     samtools index $wd/gem_result/alignment_BAM/${sample_name}_filtered_sorted.bam

#     # echo "bamCoverage"
#     # mkdir -p $wd/gem_result/bigWig
#     # bamCoverage -b $wd/gem_result/alignment_BAM/${sample_name}_filtered_sorted.bam \
#     # --normalizeUsing RPKM \
#     # --extendReads -p max/2 \
#     # -o $wd/gem_result/bigWig/${sample_name}.bw
    
# done

# bamCompare  -b1 $wd/gem_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
#         -b2 $wd/gem_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
#         -o $wd/gem_result/bigWig/bamCompare.peaks.bw \
#         --binSize 80 --operation ratio --scaleFactorsMethod SES -n 1000


# mkdir -p $wd/gem_result/peaks
# macs3 callpeak -t $wd/gem_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
# -c $wd/gem_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
# -f BAMPE --keep-dup all --cutoff-analysis -g 399348944 -n N2-6_L2_treatment -B --outdir $wd/gem_result/peaks 2> $wd/gem_result/peaks/N2-6_L2_treatment_macs3.log 

export _JAVA_OPTIONS='-Xms10G -Xmx15G'
mkdir -p $wd/gem_result/GEM
gem --t 32 \
--f BAM \
--d $wd/result/GEM/Read_Distribution_default.txt \
--g $wd/gem_result/reference/size.txt \
--genome $wd/gem_result/reference/scaffold \
--expt $wd/gem_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
--ctrl $wd/gem_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
--outBED \
--out $wd/gem_result/GEM/N2-6_L2 \
--k_min 6 --kmax 20 --k_seqs 600 --k_neg_dinu_shuffle
