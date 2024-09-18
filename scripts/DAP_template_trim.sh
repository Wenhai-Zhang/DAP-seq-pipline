

wd="/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906"

# check md5
# cd $wd/data/raw_sequence_data/240902-A00133B
# md5sum -c md5.md5
# cd $wd

input_dir=$wd/data/raw_sequence_data/240902-A00133B
# samples=("N-W-5_L2" "N2-6_L2")
# for sample_name in "${samples[@]}"
# do
#     read1=${sample_name}_1.fq.gz
#     read2=${sample_name}_2.fq.gz
#     # echo "============================================ adapter trimming ============================================"
#     # # echo "----------------------------trim_galore----------------------------"
#     # # mkdir -p $wd/trim_result/adapter_trimming/trim_galore
#     # # mkdir -p $wd/trim_result/fastqc/trim_galore
#     # # trim_galore -j 4 --fastqc --fastqc_args "--outdir $wd/trim_result/fastqc/trim_galore" \
#     # # --gzip --output_dir $wd/trim_result/adapter_trimming/trim_galore \
#     # # --paired $input_dir/$read1 $input_dir/$read2 \
#     # # -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
#     # # echo "----------------------------cutadapt----------------------------"
#     # # mkdir -p $wd/trim_result/adapter_trimming/cutadapt
#     # # cutadapt -j 24 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     # # -o $wd/trim_result/adapter_trimming/cutadapt/${read1//_1.fq.gz/_cutadapt_1.fq.gz} \
#     # # -p $wd/trim_result/adapter_trimming/cutadapt/${read2//_2.fq.gz/_cutadapt_2.fq.gz} \
#     # # $input_dir/$read1 $input_dir/$read2
#     # echo "----------------------------fastp----------------------------"
#     # mkdir -p $wd/trim_result/adapter_trimming/fastp
#     # fastp -i $input_dir/$read1 -I $input_dir/$read2 \
#     # -o $wd/trim_result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz} \
#     # -O $wd/trim_result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz} \
#     # --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     # --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     # -f 2 -F 2 \
#     # -h $wd/trim_result/adapter_trimming/fastp/fastp.html -j $wd/trim_result/adapter_trimming/fastp/fastp.json -w 16

#     # echo "============================================ fastqc ============================================"
#     # # mkdir -p $wd/trim_result/fastqc/origin
#     # # fastqc -o $wd/trim_result/fastqc/origin -t 24 $input_dir/$read1
#     # # fastqc -o $wd/trim_result/fastqc/origin -t 24 $input_dir/$read2
#     # # mkdir -p $wd/trim_result/fastqc/cutadapt
#     # # fastqc -o $wd/trim_result/fastqc/cutadapt -t 24 $wd/trim_result/adapter_trimming/cutadapt/${read1//_1.fq.gz/_cutadapt_1.fq.gz}
#     # # fastqc -o $wd/trim_result/fastqc/cutadapt -t 24 $wd/trim_result/adapter_trimming/cutadapt/${read2//_2.fq.gz/_cutadapt_2.fq.gz}
#     # mkdir -p $wd/trim_result/fastqc/fastp
#     # fastqc -o $wd/trim_result/fastqc/fastp -t 24 $wd/trim_result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz}
#     # fastqc -o $wd/trim_result/fastqc/fastp -t 24 $wd/trim_result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz}

#     # echo "============================================ build bowtie index ============================================"

#     # echo "============================================ bowtie2 read alignment ============================================"
#     # # default we use fastp result
#     # mkdir -p $wd/trim_result/alignment_SAM
#     # mkdir -p $wd/trim_result/QC
#     # bowtie2 --very-sensitive -p 64 \
#     # -x $wd/result/bowtie2/R108_HiC \
#     # -1 $wd/trim_result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz} \
#     # -2 $wd/trim_result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz} \
#     # -S $wd/trim_result/alignment_SAM/${sample_name}.sam 2>&1 | tee $wd/trim_result/QC/${sample_name}_alignment.log

#     # echo "convert SAM file to BAM file"
#     # mkdir -p $wd/trim_result/alignment_BAM
#     # samtools view -@ 16 -bS $wd/trim_result/alignment_SAM/${sample_name}.sam > $wd/trim_result/alignment_BAM/${sample_name}.bam
#     # echo "sort bam file"
#     # samtools sort -@ 16 -O BAM -o $wd/trim_result/alignment_BAM/${sample_name}_sorted.bam $wd/trim_result/alignment_BAM/${sample_name}.bam

#     # echo "index bam file"
#     # samtools index $wd/trim_result/alignment_BAM/${sample_name}_sorted.bam

#     # echo "flagstat"
#     # samtools flagstat $wd/trim_result/alignment_BAM/${sample_name}.bam > $wd/trim_result/alignment_BAM/${sample_name}_flagstat.txt

#     # samtools addreplacerg -r "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:Illumina\tLB:${sample_name}" -o $wd/trim_result/alignment_BAM/${sample_name}_addrg.bam $wd/trim_result/alignment_BAM/${sample_name}_sorted.bam

#     echo "Mark duplicates"
#     export _JAVA_OPTIONS='-Xms10G -Xmx15G'
#     picard MarkDuplicates \
#     QUIET=true \
#     INPUT=$wd/trim_result/alignment_BAM/${sample_name}_addrg.bam \
#     OUTPUT=$wd/trim_result/alignment_BAM/${sample_name}_dupRemoved.bam \
#     METRICS_FILE=$wd/trim_result/QC/${sample_name}_dupMarked.sorted.metrics \
#     REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
#     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#     echo "Remove multi-mapped reads"
#     samtools view -@ 16 -h -q 20 $wd/trim_result/alignment_BAM/${sample_name}_dupRemoved.bam > $wd/trim_result/alignment_BAM/${sample_name}_mapQ20.bam

#     echo "Remove reads unmapped 4, mate unmapped 8, not primary alignment 256, reads failing platform 512"
#     samtools view -@ 16 -h -b -F 780 -f 2 $wd/trim_result/alignment_BAM/${sample_name}_mapQ20.bam > $wd/trim_result/alignment_BAM/${sample_name}_filtered.bam

#     echo "sort bam file"
#     samtools sort -@ 16 -O BAM -o $wd/trim_result/alignment_BAM/${sample_name}_filtered_sorted.bam $wd/trim_result/alignment_BAM/${sample_name}_filtered.bam

#     echo "index bam file"
#     samtools index $wd/trim_result/alignment_BAM/${sample_name}_filtered_sorted.bam

#     # echo "bamCoverage"
#     # mkdir -p $wd/trim_result/bigWig
#     # bamCoverage -b $wd/trim_result/alignment_BAM/${sample_name}_filtered_sorted.bam \
#     # --normalizeUsing RPKM \
#     # --extendReads -p max/2 \
#     # -o $wd/trim_result/bigWig/${sample_name}.bw
    
# done

# bamCompare  -b1 $wd/trim_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
#         -b2 $wd/trim_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
#         -o $wd/trim_result/bigWig/bamCompare.peaks.bw \
#         --binSize 80 --operation ratio --scaleFactorsMethod SES -n 1000


mkdir -p $wd/trim_result/peaks
macs3 callpeak -t $wd/trim_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
-c $wd/trim_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
-f BAMPE --keep-dup all --cutoff-analysis -g 399348944 -n N2-6_L2_treatment -B --outdir $wd/trim_result/peaks 2> $wd/trim_result/peaks/N2-6_L2_treatment_macs3.log 


# mkdir -p $wd/trim_result/peaks
# macs3 callpeak -t $wd/trim_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
# -c $wd/trim_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
# -f BAMPE -g 399348944 -n N2-6_L2_treatment_q3 -p 0.001 -B --outdir $wd/trim_result/peaks 2> $wd/trim_result/peaks/N2-6_L2_treatment_macs3_q3.log

# mkdir -p $wd/trim_result/peaks
# macs3 callpeak -t $wd/trim_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
# -c $wd/trim_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
# -f BAMPE -g 399348944 -n N2-6_L2_treatment_q2 -q 0.01 -B --outdir $wd/trim_result/peaks 2> $wd/trim_result/peaks/N2-6_L2_treatment_macs3_q2.log 

# mkdir -p $wd/trim_result/peaks
# macs3 callpeak -t $wd/trim_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
# -c $wd/trim_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
# -f BAMPE -g 399348944 -n N2-6_L2_treatment__q3_real -q 0.001 -B --outdir $wd/trim_result/peaks 2> $wd/trim_result/peaks/N2-6_L2_treatment_macs3_q3_real.log 

# mkdir -p $wd/trim_result/peaks
# macs3 callpeak -t $wd/trim_result/alignment_BAM/N2-6_L2_sorted.bam \
# -c $wd/trim_result/alignment_BAM/N-W-5_L2_sorted.bam \
# -f BAMPE -g 399348944 -n N2-6_L2_treatment_undup -B --outdir $wd/trim_result/peaks 2> $wd/trim_result/peaks/N2-6_L2_treatment_macs3_undup.log 


# mkdir -p $wd/trim_result/macs2/peaks
# macs2 callpeak -t $wd/trim_result/alignment_BAM/N2-6_L2_filtered_sorted.bam \
# -c $wd/trim_result/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
# -f BAM --keep-dup all --cutoff-analysis -g 399348944 -n N2-6_L2_treatment_macs2_bam -B --outdir $wd/trim_result/macs2/peaks 2> $wd/trim_result/macs2/peaks/N2-6_L2_treatment_macs2_bam.log 
