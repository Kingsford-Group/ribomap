#!/bin/bash
#=============================
# default parameters
#=============================
adapter=CTGTAGGCACCATCAAT
min_fplen=27
max_fplen=33
nproc=15 # threads
nmismatch=1
offset=12 # P-site offset
#=============================
# pre-filled parameters
#=============================
src_dir=`dirname $0`
bin_dir=${src_dir}/../bin/
lib_dir=${src_dir}/../lib/
export PATH=${bin_dir}:$PATH
export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${lib_dir}:$DYLD_LIBRARY_PATH
work_dir=${src_dir}/../
star_idx_dir=${work_dir}StarIndex/
# star index
rrna_idx=${star_idx_dir}contaminant/
transcript_idx=${star_idx_dir}transcript/
# star params
align_params="--seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax ${nmismatch} --outFilterIntronMotifs RemoveNoncanonical"
SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS" # --outSAMprimaryFlag AllBestScore"
#=============================
# functions
#=============================
# print error message and usage message
# quit program if indicated
# error_msg $error_msg $quit
error_msg ()
{
    echo "$1"
    if [ "$2" = true ]; then
	echo "Usage: ./run_ribomap.sh --rnaseq_fq rnaseq.fq.gz --riboseq_fq riboseq.fq.gz --contaminant_fa contaminant.fa --transcript_fa transcript.fa --cds_range cds_range.txt"
	exit
    fi
}

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
# read in command line args
# space separated
#=============================
while [[ $# > 1 ]]
do
    key="$1"
    shift
    case $key in
	--rnaseq_fq)
	    rnaseq_fq="$1"
	    shift
	    ;;
	--riboseq_fq)
	    riboseq_fq="$1"
	    shift
	    ;;
	--transcript_fa)
	    transcript_fa="$1"
	    shift
	    ;;
	--contaminant_fa)
	    contaminant_fa="$1"
	    shift
	    ;;
	--cds_range)
	    cds_range="$1"
	    shift
	    ;;
	--work_dir)
	    work_dir="$1"
	    shift
	    ;;
	--adapter)
	    adapter="$1"
	    shift
	    ;;
	--min_fplen)
	    min_fplen="$1"
	    shift
	    ;;
	--max_fplen)
	    max_fplen="$1"
	    shift
	    ;;
	--nproc)
	    nproc="$1"
	    shift
	    ;;
	--nmismatch)
	    nmismatch="$1"
	    shift
	    ;;
	--offset)
	    offset="$1"
	    shift
	    ;;
	--fasta_dir)
	    fasta_dir="$1"
	    shift
	    ;;
	--star_idx_dir)
	    star_idx_dir="$1"
	    rrna_idx=${star_idx_dir}contaminant/
	    transcript_idx=${star_idx_dir}transcript/
	    shift
	    ;;
	--alignment_dir)
	    tmp_dir="$1"
	    shift
	    ;;
	--sailfish_dir)
	    sm_odir="$1"
	    shift
	    ;;
	--output_dir)
	    output_dir="$1"
	    shift
	    ;;
	*)
            # unknown option
	    ;;
    esac
done

if [ -z "${riboseq_fq}" ]; then 
    error_msg "ribo-seq reads not provided!" true
elif [ ! -f ${riboseq_fq} ]; then
    error_msg "ribo-seq file not exist! ${riboseq_fq}" true
elif [ -z "${rnaseq_fq}" ]; then
    error_msg "RNA-seq reads not provided!" true
elif [ ! -f ${rnaseq_fq} ]; then
    error_msg "RNA-seq file not exist! ${rnaseq_fq}" true
elif [ -z "${contaminant_fa}" ]; then
    error_msg "contaminant fasta not provided! Filter step skipped." false
elif [ ! -f ${contaminant_fa} ]; then
    error_msg "contaminant fasta not exist! ${contaminant_fa}" false
elif [ -z "${cds_range}" ]; then
    error_msg "cds range not provided! assume transcript fasta only contain cds regions." false
elif [ ! -f ${cds_range} ]; then
    error_msg "cds range file not exist! ${cds_range}" false
fi

#=============================
# make directories
#=============================
# preprocess
fasta_dir=${work_dir}data/fasta/
# star outputs
tmp_dir=${work_dir}alignment/
# salmon
sm_odir=${work_dir}sm_quant
# ribomap
output_dir=${work_dir}outputs
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
fastx_pipe="fastx_clipper -Q33 -a ${adapter} -l ${min_fplen} -c -n -v | fastq_to_fasta -v"
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
    python ${src_dir}/filter_reads_by_size.py ${ribo_fa} ${ribo_size_fa} ${min_fplen} ${max_fplen}
    check_file ${ribo_size_fa} "pipeline failed at filtering riboseq with the wrong size!"
fi
#=============================
# step 2: filter rrna
#=============================
ornaprefix=${tmp_dir}${rna_core}_rrna_
oriboprefix=${tmp_dir}${ribo_core}_rrna_
rna_nrrna_fa=${ornaprefix}Unmapped.out.mate1
ribo_nrrna_fa=${oriboprefix}Unmapped.out.mate1
if [ ! -z "${contaminant_fa}" ] && [ -f ${contaminant_fa} ]; then
    echo "filtering contaminated reads"
    if  [ ! -d ${rrna_idx} ];  then
	echo "building contaminant index..."
	mkdir -p ${rrna_idx}
	STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${rrna_idx} --genomeFastaFiles ${contaminant_fa} --genomeSAindexNbases 5 --genomeChrBinNbits 11
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
else
    echo "skipped filter read step."
    # TODO change rna_nrna_fa file name here
fi
#========================================
# step 3: map to transcriptome
#========================================
ornaprefix=${tmp_dir}${rna_core}_transcript_
oriboprefix=${tmp_dir}${ribo_core}_transcript_
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
ribomap_out=${output_dir}/${ribo_core}
options="--mrnabam ${rna_bam} --ribobam ${ribo_bam} --fasta ${transcript_fa} --sf ${sm_out} --offset ${offset} --out ${ribomap_out}"
if [ ! -z "${cds_range}" ] && [ -f ${cds_range} ]; then
    options+=" --cds_range ${cds_range}"
fi
if [ ! -f ${ribomap_out}.base ]; then
    echo "running riboprof..."
    riboprof ${options}
    check_file ${ribomap_out}.base "pipeline failed at ribosome profile generation!"
fi
