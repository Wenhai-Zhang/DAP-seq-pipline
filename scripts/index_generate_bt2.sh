
wd="/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906"

# mkdir -p $wd/result/bowtie2/
# bowtie2-build $wd/data/reference_genome/R108_HiC/MtrunR108_HiC.refseq.fasta $wd/result/bowtie2/R108_HiC --threads 64

# index for GEM
mkdir -p $wd/gem_result/bowtie2/
bowtie2-build $wd/gem_result/reference/MtrunR108_HiC.refseq.fasta $wd/gem_result/bowtie2/R108_HiC --threads 64