library(ChIPseeker)
library(GenomicFeatures)
library(txdbmaker)

MtrunR108_txdb <-makeTxDbFromGFF("/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/data/reference_genome/R108_HiC/20210817_MtrunR108_HiC_FunctionalAnnotation/annotation/MedtrR108_HiC.gtf", format = "gtf")
# file <- "/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/trim_result/peaks/N2-6_L2_treatment_peaks.xls"
# file <- "/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/result/peaks/N2-6_L2_treatment_peaks.xls"
file <- "/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/result/peakranger/ranger_result_region.bed"
peak <- readPeakFile(file)
peakAnno <- annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=MtrunR108_txdb)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
#查看TSS附近peak的分布频率
plotAvgProf2(peak, TxDb=MtrunR108_txdb, upstream=2500, downstream=2500,             
             xlab="Genomic Region (5'->3')",              
             ylab = "Peak Count Frequency",conf = 0.95, resample = 500)
peak_df <- as.data.frame(peakAnno)
filtered_df <- peak_df[peak_df$annotation == "Promoter", ]
write.table(filtered_df, file = "output.tsv", row.names = FALSE, sep = "\t")



library(dplyr)
library(stringr)

# 读取 GFF3 文件
gff_data <- read.table("/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/data/reference_genome/R108_HiC/20210817_MtrunR108_HiC_FunctionalAnnotation/annotation/MedtrR108_HiC.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE, skip = 2)

# 确保 GFF3 第9列存在并从中提取 `ID` 和 `locus_tag`
# 假设第9列是注释字段
gff_data %>%
  mutate(ID = str_extract(V9, "ID=[^;]+"),  # 提取ID
         locus_tag = str_extract(V9, "locus_tag=[^;]+")) %>%  # 提取locus_tag
  select(ID, locus_tag) -> extracted_data

# 清理ID和locus_tag列的前缀
extracted_data$geneId <- str_replace(extracted_data$ID, "ID=", "")
extracted_data$LOCUS_TAG <- str_replace(extracted_data$locus_tag, "locus_tag=", "")
extracted_data <- subset(extracted_data, select = -c(ID, locus_tag))
merged_df <- merge(filtered_df, extracted_data, by = "geneId", all = FALSE)
summary_file <- "/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/data/reference_genome/R108_HiC/20210817_MtrunR108_HiC_FunctionalAnnotation/annotation/MedtrR108_HiC.summary.tsv"
summary_file_data <- read.csv(summary_file, sep = "\t", stringsAsFactors = FALSE)
merged_df <- merge(merged_df, summary_file_data, by = "LOCUS_TAG", all = FALSE)
write.table(merged_df, file = "output.tsv", row.names = FALSE, sep = "\t")
