

wd="/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906"

# check md5
# cd $wd/data/raw_sequence_data/240902-A00133B
# md5sum -c md5.md5
# cd $wd

input_dir=$wd/data/raw_sequence_data/240902-A00133B
samples=("N-W-5_L2" "N2-6_L2")
# for sample_name in "${samples[@]}"
# do
#     read1=${sample_name}_1.fq.gz
#     read2=${sample_name}_2.fq.gz
#     # echo "============================================ adapter trimming ============================================"
#     # # echo "----------------------------trim_galore----------------------------"
#     # # mkdir -p $wd/result_macs2/adapter_trimming/trim_galore
#     # # mkdir -p $wd/result_macs2/fastqc/trim_galore
#     # # trim_galore -j 4 --fastqc --fastqc_args "--outdir $wd/result_macs2/fastqc/trim_galore" \
#     # # --gzip --output_dir $wd/result_macs2/adapter_trimming/trim_galore \
#     # # --paired $input_dir/$read1 $input_dir/$read2 \
#     # # -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
#     # # echo "----------------------------cutadapt----------------------------"
#     # # mkdir -p $wd/result_macs2/adapter_trimming/cutadapt
#     # # cutadapt -j 24 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     # # -o $wd/result_macs2/adapter_trimming/cutadapt/${read1//_1.fq.gz/_cutadapt_1.fq.gz} \
#     # # -p $wd/result_macs2/adapter_trimming/cutadapt/${read2//_2.fq.gz/_cutadapt_2.fq.gz} \
#     # # $input_dir/$read1 $input_dir/$read2
#     # echo "----------------------------fastp----------------------------"
#     # mkdir -p $wd/result_macs2/adapter_trimming/fastp
#     # fastp -i $input_dir/$read1 -I $input_dir/$read2 \
#     # -o $wd/result_macs2/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz} \
#     # -O $wd/result_macs2/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz} \
#     # --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     # --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     # -h $wd/result_macs2/adapter_trimming/fastp/fastp.html -j $wd/result_macs2/adapter_trimming/fastp/fastp.json -w 16

#     # echo "============================================ fastqc ============================================"
#     # # mkdir -p $wd/result_macs2/fastqc/origin
#     # # fastqc -o $wd/result_macs2/fastqc/origin -t 24 $input_dir/$read1
#     # # fastqc -o $wd/result_macs2/fastqc/origin -t 24 $input_dir/$read2
#     # # mkdir -p $wd/result_macs2/fastqc/cutadapt
#     # # fastqc -o $wd/result_macs2/fastqc/cutadapt -t 24 $wd/result_macs2/adapter_trimming/cutadapt/${read1//_1.fq.gz/_cutadapt_1.fq.gz}
#     # # fastqc -o $wd/result_macs2/fastqc/cutadapt -t 24 $wd/result_macs2/adapter_trimming/cutadapt/${read2//_2.fq.gz/_cutadapt_2.fq.gz}
#     # mkdir -p $wd/result_macs2/fastqc/fastp
#     # fastqc -o $wd/result_macs2/fastqc/fastp -t 24 $wd/result_macs2/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz}
#     # fastqc -o $wd/result_macs2/fastqc/fastp -t 24 $wd/result_macs2/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz}

#     # echo "============================================ build bowtie index ============================================"

#     # echo "============================================ bowtie2 read alignment ============================================"
#     # # default we use fastp result
#     # mkdir -p $wd/result_macs2/alignment_SAM
#     # mkdir -p $wd/result_macs2/QC
#     # bowtie2 --very-sensitive -p 64 \
#     # -x $wd/result/bowtie2/R108_HiC \
#     # -1 $wd/result/adapter_trimming/fastp/${read1//_1.fq.gz/_fastp_1.fq.gz} \
#     # -2 $wd/result/adapter_trimming/fastp/${read2//_2.fq.gz/_fastp_2.fq.gz} \
#     # -S $wd/result_macs2/alignment_SAM/${sample_name}.sam 2>&1 | tee $wd/result_macs2/QC/${sample_name}_alignment.log

#     # echo "convert SAM file to BAM file"
#     # mkdir -p $wd/result_macs2/alignment_BAM
#     # samtools view -bS $wd/result_macs2/alignment_SAM/${sample_name}.sam > $wd/result_macs2/alignment_BAM/${sample_name}.bam
#     # echo "sort bam file"
#     # samtools sort -@ 16 -O BAM -o $wd/result_macs2/alignment_BAM/${sample_name}_sorted.bam $wd/result_macs2/alignment_BAM/${sample_name}.bam

#     # echo "index bam file"
#     # samtools index $wd/result_macs2/alignment_BAM/${sample_name}_sorted.bam

#     # echo "flagstat"
#     # samtools flagstat $wd/result_macs2/alignment_BAM/${sample_name}.bam > $wd/result_macs2/alignment_BAM/${sample_name}_flagstat.txt

#     # samtools addreplacerg -r "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:Illumina\tLB:${sample_name}" -o $wd/result_macs2/alignment_BAM/${sample_name}_addrg.bam $wd/result_macs2/alignment_BAM/${sample_name}_sorted.bam

#     # echo "Mark duplicates"
#     # export _JAVA_OPTIONS='-Xms10G -Xmx15G'
#     # picard MarkDuplicates \
#     # QUIET=true \
#     # INPUT=$wd/result_macs2/alignment_BAM/${sample_name}_addrg.bam \
#     # OUTPUT=$wd/result_macs2/alignment_BAM/${sample_name}_dupRemoved.bam \
#     # METRICS_FILE=$wd/result_macs2/QC/${sample_name}_dupMarked.sorted.metrics \
#     # REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
#     # MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#     echo "Remove multi-mapped reads"
#     samtools view -h -q 20 $wd/result_macs2/alignment_BAM/${sample_name}_dupRemoved.bam > $wd/result_macs2/alignment_BAM/${sample_name}_mapQ20.bam

#     # echo "Remove reads unmapped 4, mate unmapped 8, not primary alignment 256, reads failing platform 512"
#     # samtools view -h -b -F 780 -f 2 $wd/result_macs2/alignment_BAM/${sample_name}_mapQ20.bam > $wd/result_macs2/alignment_BAM/${sample_name}_filtered.bam

#     echo "sort bam file"
#     samtools sort -@ 16 -O BAM -o $wd/result_macs2/alignment_BAM/${sample_name}_filtered_sorted.bam $wd/result_macs2/alignment_BAM/${sample_name}_mapQ20.bam

#     echo "index bam file"
#     samtools index $wd/result_macs2/alignment_BAM/${sample_name}_filtered_sorted.bam

#     # echo "bamCoverage"
#     # mkdir -p $wd/result_macs2/bigWig
#     # bamCoverage -b $wd/result_macs2/alignment_BAM/${sample_name}_filtered_sorted.bam \
#     # --normalizeUsing RPKM \
#     # --extendReads -p max/2 \
#     # -o $wd/result_macs2/bigWig/${sample_name}.bw
    
# done

# bamCompare  -b1 $wd/result_macs2/alignment_BAM/N2-6_L2_filtered_sorted.bam \
#         -b2 $wd/result_macs2/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
#         -o $wd/result_macs2/bigWig/bamCompare.peaks.bw \
#         --binSize 80 --operation ratio --scaleFactorsMethod SES -n 1000


mkdir -p $wd/result_macs2/peaks
macs2 callpeak -t $wd/result_macs2/alignment_BAM/N2-6_L2_filtered_sorted.bam \
-c $wd/result_macs2/alignment_BAM/N-W-5_L2_filtered_sorted.bam \
-f BAMPE --keep-dup all --cutoff-analysis -g 399348944 -n N2-6_L2_treatment -B --outdir $wd/result_macs2/peaks 2> $wd/result_macs2/peaks/N2-6_L2_treatment_macs3.log 



