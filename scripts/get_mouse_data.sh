#!/bin/bash
work_dir=/home/hw1/scratch/ribomap/
data_dir=${work_dir}data/
sra_dir=${data_dir}sra
fasta_dir=${data_dir}fasta
ref_dir=${work_dir}ref/mouse/
#=============================
# make directories
#=============================
mkdir -p ${data_dir}
mkdir -p ${sra_dir}
mkdir -p ${fasta_dir}
mkdir -p ${ref_dir}
#=============================
# step 1: download sra & refs
#=============================
# echo "downloading data..."
# cd ${sra_dir}
# # ES cell feeder-free, w/ LIF 60 s CYH (100 ug/ml) mrna_mesc_yeslif Illumina GAII
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084812/SRR315595/SRR315595.sra
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084812/SRR315596/SRR315596.sra
# # ES cell feeder-free, w/ LIF 60 s CYH (100 ug/ml) ribo_mesc_chx Illumina GAII
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084815/SRR315601/SRR315601.sra
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084815/SRR315602/SRR315602.sra
# # ES cell feeder-free, w/ LIF 60 s CYH (100 ug/ml) ribo_mesc_yeslif Illumina GAII
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084824/SRR315624/SRR315624.sra
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084824/SRR315625/SRR315625.sra
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084824/SRR315626/SRR315626.sra
cd ${ref_dir}
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.pc_transcripts.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.pc_translations.fa.gz
gunzip *.gz
exit
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

