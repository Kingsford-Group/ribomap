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
gtf_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz
tfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.pc_transcripts.fa.gz
pfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.pc_translations.fa.gz
# wget -P ${ref_dir} ${gtf_url}
# wget -P ${ref_dir} ${tfa_url}
# wget -P ${ref_dir} ${pfa_url}
# echo "unzipping data..."
# gunzip -f ${ref_dir}*.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz
#=============================
# step 2: process reference
# build cds range file
#=============================
echo "filtering transcripts...."
gtf=${ref_dir}${gtf_url##*/}
gtf=${gtf%.gz}
tfa=${ref_dir}${tfa_url##*/}
tfa=${tfa%.gz}
pfa=${ref_dir}${pfa_url##*/}
pfa=${pfa%.gz}
python filter_gencode_transcript.py ${gtf} ${tfa} ${pfa}
# echo "getting rid of comment lines in gtf..."
# mv ${gtf} ${gtf}.bak
# sed '/^#/ d' < ${gtf}.bak > ${gtf}
# rm ${gtf}.bak
#=============================
# step 3: run ribomap
#=============================

