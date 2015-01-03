#!/bin/bash
work_dir=/home/hw1/scratch/ribomap/
fasta_dir=${work_dir}data/fasta
ref_dir=${work_dir}ref/human/
#=============================
# make directories
#=============================
mkdir -p ${fasta_dir}
mkdir -p ${ref_dir}
#=============================
# step 1: download sra & refs
#=============================
echo "downloading data..."
cd ${ref_dir}
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.pc_transcripts.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.pc_translations.fa.gz
echo "unzipping data..."
gunzip *.gz
cd ${fasta_dir}
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz
#=============================
# step 2: process reference
# build cds range file
#=============================
echo "filtering transcripts...."
gtf=${ref_dir}gencode.v21.annotation.gtf
tfa=${ref_dir}gencode.v21.pc_transcripts.fa
pfa=${ref_dir}gencode.v21.pc_translations.fa
python filter_gencode_transcript.py ${gtf} ${tfa} ${pfa}
# echo "getting rid of comment lines in gtf..."
# mv ${gtf} ${gtf}.bak
# sed '/^#/ d' < ${gtf}.bak > ${gtf}
# rm ${gtf}.bak
#=============================
# step 3: run ribomap
#=============================

