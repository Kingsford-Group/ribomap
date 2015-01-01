#!/bin/bash
#=============================
# default parameters
#=============================
# please fill this line out by yourself
work_dir=/home/hw1/scratch/ribomap/
# references -- required!!! please fill out these lines
transcript_fa=/home/hw1/scratch/ribojamdetector/transcriptome/protein_coding_33_filtered.fasta
rrna_fa=/home/hw1/scratch/ribojamdetector/transcriptome/rRNA.fa
# reads --required!!! please provide fastq reads file names
rnaseq_fq=/home/hw1/scratch/ribojamdetector/data/fasta/BY_mRNA.fastq.gz
riboseq_fq=/home/hw1/scratch/ribojamdetector/data/fasta/BY_FP.fastq.gz
# default parameters: you can leave them alone
adapter=CTGTAGGCACCATCAAT
min_fplen=27
max_fplen=33
nproc=15 # threads
nmismatch=1
offset=15 # P-site offset
#==============================
# preprocess
fasta_dir=${work_dir}data/fasta/
# star index
star_idx_dir=${work_dir}StarIndex/yeast/
rrna_idx=${star_idx_dir}contaminant/
transcript_idx=${star_idx_dir}transcript/
# star outputs
tmp_dir=${work_dir}tmp/
# star params
align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax ${nmismatch} --outFilterIntronMotifs RemoveNoncanonical"
SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMprimaryFlag AllBestScore"
# salmon
sm_odir=${work_dir}sm_quant
# ribomap
output_dir=${work_dir}outputs
#=============================
# functions
#=============================
# check whether file generated successfully at each single step
# quit program if file failed to be generated
# check_file $file_name $error_msg
check_file ()
{
    if [ ! -f "$1" ]; then
	echo "$2"
	exit
    fi
}
#=============================
# make directories
#=============================
mkdir -p ${fasta_dir}
mkdir -p ${tmp_dir}
mkdir -p ${output_dir}
rna_core=${rnaseq_fq##*/}
rna_core=${rna_core%%.*}
ribo_core=${riboseq_fq##*/}
ribo_core=${ribo_core%%.*}
#=============================
# step 1: preprocess reads
#=============================
echo "preprocessing reads (quality control + trim adapter + trim first base + collapse duplicate reads + fastq to fasta"
fastx_pipe="fastx_clipper -Q33 -a ${adapter} -l ${min_fplen} -n -v | fastq_to_fasta -v"
rna_fa=${fasta_dir}${rna_core}.fa
ribo_fa=${fasta_dir}${ribo_core}.fa
ribo_size_fa=${fasta_dir}${ribo_core}-size.fa
if [ ! -f ${rna_fa} ]; then
    zcat ${rnaseq_fq} | ${fastx_pipe} -o ${rna_fa}
    check_file ${rna_fa} "pipeline failed at preprocessing rnaseq_fq: ${rnaseq_fq}!"
fi
if [ ! -f ${ribo_fa} ]; then
    zcat ${riboseq_fq} | ${fastx_pipe} | fastx_collapser -v -o ${ribo_fa}
    check_file ${ribo_fa} "pipeline failed at preprocessing riboseq_fq: ${riboseq_fq}!"
fi
if [ ! -f ${ribo_size_fa} ]; then
    python filter_reads_by_size.py ${ribo_fa} ${ribo_size_fa} ${min_fplen} ${max_fplen}
    check_file ${ribo_size_fa} "pipeline failed at filtering riboseq with the wrong size!"
fi
#=============================
# step 2: filter rrna
#=============================
ornaprefix=${tmp_dir}${rna_core}_rrna_
oriboprefix=${tmp_dir}${ribo_core}_rrna_
rna_nrrna_fa=${ornaprefix}Unmapped.out.mate1
ribo_nrrna_fa=${oriboprefix}Unmapped.out.mate1
echo "filtering contaminated reads"
if  [ ! -d ${rrna_idx} ];  then
    echo "building contaminant index..."
    mkdir -p ${rrna_idx}
    STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${rrna_idx} --genomeFastaFiles ${rrna_fa} --genomeSAindexNbases 5 --genomeChrBinNbits 11
fi
if  [ ! -f ${rna_nrrna_fa} ];  then
    echo "filtering contaminants in RNA_seq..."
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${rna_fa} --outFileNamePrefix ${ornaprefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${align_params} > /dev/null
    check_file ${rna_nrrna_fa} "pipeline failed at filtering rrna in RNA_seq!"
fi
if  [ ! -f ${ribo_nrrna_fa} ];  then
    echo "filtering contaminants in ribo_seq..."
    STAR --runThreadN $nproc --genomeDir ${rrna_idx} --readFilesIn ${ribo_fa} --outFileNamePrefix ${oriboprefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS ${align_params} > /dev/null
    check_file ${ribo_nrrna_fa} "pipeline failed at filtering rrna in ribo_seq!"
fi
#========================================
# step 3: map to transcriptome
#========================================
ornaprefix=${tmp_dir}${rna_core}_transcriptb_
oriboprefix=${tmp_dir}${ribo_core}_transcriptb_
rna_bam=${ornaprefix}Aligned.out.bam
ribo_bam=${oriboprefix}Aligned.out.bam
echo "aligning reads to transcriptome"
if  [ ! -d ${transcript_idx} ];  then
    echo "building transcriptome index..."
    mkdir -p ${transcript_idx}
    STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${transcript_idx} --genomeFastaFiles ${transcript_fa} --genomeSAindexNbases 11 --genomeChrBinNbits 12
fi
if [ ! -f ${rna_bam} ]; then
    echo "aligning RNA_seq to transcriptome..."
    STAR --runThreadN $nproc --genomeDir ${transcript_idx} --readFilesIn ${rna_nrrna_fa} --outFileNamePrefix ${ornaprefix} ${SAM_params} ${align_params}
    check_file ${rna_bam} "pipeline failed at mapping RNA_seq to transcriptome!"
fi
if [ ! -f ${ribo_bam} ]; then
    echo "aligning ribo_seq to transcriptome..."
    STAR --runThreadN $nproc --genomeDir ${transcript_idx} --readFilesIn ${ribo_nrrna_fa} --outFileNamePrefix ${oriboprefix} ${SAM_params} ${align_params}
    check_file ${ribo_bam} "pipeline failed at mapping ribo_seq to transcriptome!"
fi
#============================================
# step 4: salmon expression quantification
#============================================
sm_out=${sm_odir}/quant.sf
if [ ! -f ${sm_out} ]; then
    echo "running salmon quant..."
    salmon quant -t ${transcript_fa} -l U -a ${rna_bam} -o ${sm_odir} -p $nproc --bias_correct
    check_file ${sm_out} "pipeline failed at expression quantification!"
fi
#=============================
# step 5: run ribomap
#=============================
ribomap_out=${output_dir}/${ribo_core}_norm.profile
if [ ! -f ${ribomap_out} ]; then
    echo "running riboprof..."
    ./riboprof --mrnabam ${rna_bam} --ribobam ${ribo_bam} --fasta ${transcript_fa} --sf ${sm_out} --out ${ribomap_out} --offset ${offset}
    check_file ${ribomap_out} "pipeline failed at ribosome profile generation!"
fi
