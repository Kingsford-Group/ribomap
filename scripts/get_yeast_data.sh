#!/bin/bash
data_dir=/home/hw1/scratch/ribomap/data/
sra_dir="${data_dir}"sra
fasta_dir="${data_dir}"fasta
#=============================
# make directories
#=============================
mkdir -p ${data_dir}
mkdir -p ${sra_dir}
mkdir -p ${fasta_dir}
#=============================
# step 1: download sra
#=============================
cd ${sra_dir}
echo "downloading data..."
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476344/SRR1177156/SRR1177156.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476345/SRR1177157/SRR1177157.sra
#==============================
# step 2: convert sra to fastq
#==============================
echo "converting sra to fastq.gz..."
fastq-dump SRR1177156.sra -O ${fasta_dir} --gzip
fastq-dump SRR1177157.sra -O ${fasta_dir} --gzip
#==============================
# step 3: rename files 
#==============================
echo "renaming fastq files to be more informative..."
cd ${fasta_dir}
mv SRR1177156.fastq.gz BY_mRNA.fastq.gz
mv SRR1177157.fastq.gz BY_FP.fastq.gz

