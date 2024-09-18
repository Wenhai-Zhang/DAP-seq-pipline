# GEM need this format
# each chromosome matched to a fasta file and must be chr*
from Bio import SeqIO

def calculate_chromosome_lengths(fasta_file):
    chromosome_lengths = {}
    
    # 读取 FASTA 文件
    for record in SeqIO.parse(fasta_file, "fasta"):
        chromosome = record.id  # 获取染色体名称
        length = len(record.seq)  # 获取染色体长度
        chromosome_lengths[chromosome] = length
        with open(f"/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/gem_result/reference/scaffold/{chromosome}.fa", "w") as f:
            f.write(f">{chromosome}\n")
            f.write(str(record.seq))
    
    return chromosome_lengths

def print_chromosome_lengths(chromosome_lengths):
    # print(f"{'Chromosome':<20} {'Length (bp)':<15}")
    # print("-" * 35)
    # for chromosome, length in chromosome_lengths.items():
    #     print(f"{chromosome:<20} {length:<15}")
    with open("/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/gem_result/reference/size.txt", "w") as f:
        for chromosome, length in chromosome_lengths.items():
            f.write(f"{chromosome}\t{length}\n")


if __name__ == "__main__":
    fasta_file = "/home/work/wenhai/bioinformatics_analysis/DAPseq_20240906/gem_result/reference/MtrunR108_HiC.refseq.fasta"  # 替换为你的基因组FASTA文件路径
    chromosome_lengths = calculate_chromosome_lengths(fasta_file)
    print_chromosome_lengths(chromosome_lengths)
