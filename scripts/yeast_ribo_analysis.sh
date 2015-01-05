#!/bin/bash
work_dir="$(pwd)"
work_dir=${work_dir%/*}/
data_dir=${work_dir}data/
sra_dir=${data_dir}sra/
fasta_dir=${data_dir}fasta/
#=============================
# make directories
#=============================
mkdir -p ${sra_dir}
mkdir -p ${fasta_dir}
#=============================
# step 1: download sra
#=============================
echo "downloading data..."
wget -P ${sra_dir} ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476344/SRR1177156/SRR1177156.sra
wget -P ${sra_dir} ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476345/SRR1177157/SRR1177157.sra
#==============================
# step 2: convert sra to fastq
#==============================
echo "converting sra to fastq.gz..."
for file in ${sra_dir}*.sra; do
    fastq-dump $file -O ${fasta_dir} --gzip
done
#==============================
# step 3: rename files 
#==============================
echo "renaming fastq files to be more informative..."
mv ${fasta_dir}SRR1177156.fastq.gz ${fasta_dir}BY_mRNA.fastq.gz
mv ${fasta_dir}SRR1177157.fastq.gz ${fasta_dir}BY_FP.fastq.gz
#==============================
# step 4: run ribomap
#==============================
rnaseq_fq=${fasta_dir}BY_mRNA.fastq.gz
riboseq_fq=${fasta_dir}BY_FP.fastq.gz
transcript_fa=/home/hw1/scratch/ribojamdetector/transcriptome/protein_coding_33_filtered.fasta
contaminant_fa=/home/hw1/scratch/ribojamdetector/transcriptome/rRNA.fa
cds_range=/home/hw1/scratch/ribojamdetector/transcriptome/cds_range.txt
star_idx_dir=${work_dir}/StarIndex/yeast/
./run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --work_dir ${work_dir} --star_idx_dir ${star_idx_dir}
