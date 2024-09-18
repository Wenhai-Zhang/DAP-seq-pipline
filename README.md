# DAP-seq-pipline

## Note
I cannot guarantee that the pipeline is completely correct. I ran this pipeline on a dataset but did not get good results. I am not sure if it is a problem with the pipeline, but I believe it is more of a problem with the data.

## Pipline
The pipline includes removing adapter with fastp, quality control with fastqc, alignment with bowtie, removing duplication with picard, filtering reads(such as quality below 20, removing unmapped reads and so on) with samtools, peak calling, peak annotation with ChIPseeker.

## peak calling
Peak calling uses three tools: [macs3](https://github.com/macs3-project/MACS), [GEM](https://github.com/ZJU-Robotics-Lab/GEM), and [peakranger](https://github.com/drestion/peakranger)

## reference
1. https://github.com/ndu-bioinfo/Dap-Seq-pipeline
2. https://github.com/SinghLabUCSF/Diapause-multiomics