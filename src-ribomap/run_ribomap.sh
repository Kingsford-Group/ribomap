#!/bin/bash
#=============================
# default values of variables
#=============================
transcript_fa=/data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa
transcript_gtf=/data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf
rrna_fa=/home/hw1/riboseq/data/human_mouse_rp_Guo/fasta/rrna_human.fasta
rnaseq_fq=/home/hw1/riboseq/data/human_mouse_rp_Guo/fasta/GSM546921_filtered_sequence.fq
riboseq_fq=/home/hw1/riboseq/data/human_mouse_rp_Guo/fasta/GSM546920_filtered_sequence.fq
# current working directory
ribo_dir=/home/hw1/ribomap/
# bowtie
#bowtie_idx_dir=/home/hw1/riboseq/data/BowtieIndex/
bowtie_idx_dir=../test/
nproc=8
# trim footprint reads to this length
seedlen=25
# offset of reads where P site maps to the codon location
offset=15
#=============================
# read in command line args
# equals separated
#=============================
for i in "$@"
do
case $i in
    --gtf=*)
    transcript_gtf="${i#*=}"
    shift
    ;;
    --ref_fa=*)
    transcript_fa="${i#*=}"
    shift
    ;;
    --rnaseq=*)
    rnaseq_fq="${i#*=}"
    shift
    ;;
    --riboseq=*)
    riboseq_fq="${i#*=}"
    shift
    ;;
    --rrna_fa=*)
    rrna_fa="${i#*=}"
    shift
    ;;
    --seedlen=*)
    seedlen="${i#*=}"
    shift
    ;;
    --offset=*)
    offset="${i#*=}"
    shift
    ;;
    --nproc=*)
    nproc="${i#*=}"
    shift
    ;;
    *)
            # unknown option
    ;;
esac
done
if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi
#=============================
# fill out other variables
#=============================
ribo_core=${riboseq_fq##*/}
ribo_core=${ribo_core%.*}
# sailfish
sf_idx_dir="${ribo_dir}"sf_idx/
sf_odir="${ribo_dir}"sf_quant/
#sf_odir=../test/
# folder for teporarily holding intermediate result
tmp_dir="${ribo_dir}"tmp
tmp_dir=../test/
# ribomap
src_dir="${ribo_dir}"src
output_dir="${ribo_dir}"outputs
output_dir=../test/
cache_dir="${ribo_dir}"cache 
#=============================
# make directories
#=============================
mkdir -p ${tmp_dir}
mkdir -p ${output_dir}
mkdir -p ${cache_dir}
mkdir -p ${bowtie_idx_dir}
#=============================
# step 1: sailfish
#=============================
# create sailfish index if not exist
if [ ! -d "$sf_idx_dir" ]; then
    sailfish index -t ${transcript_fa} -o ${sf_idx_dir} -k 20 -p $nproc -f
fi
#sailfish quant -l "T=SE:S=U" -i ${sf_idx_dir} -o ${sf_odir} -r ${rnaseq_fq} -p $nproc -a -f
#=============================
# step 2: trim sequences
#=============================
riboseq_core=${riboseq_fq##*/}
riboseq_core=${riboseq_core%.*}
nodup_fa="${tmp_dir}${riboseq_core}_${seedlen}_nodup.fa"
#./merge_fq_to_fa -i ${riboseq_fq} -o ${nodup_fa} -l ${seedlen}
#=============================
# step 3: filter rrna
#=============================
rrna_core=${rrna_fa##*/}
rrna_core=${rrna_core%.*}
if  [[ ! $(ls "${bowtie_idx_dir}${rrna_core}"*) ]];  then
    echo "bowtie index for rrna not exist, build it"
    bowtie-build -f ${rrna_fa} ${bowtie_idx_dir}${rrna_core}
fi
ndup_nrrna_fa="${tmp_dir}${riboseq_core}_${seedlen}_nodup_norrna.fa"
#bowtie -p $nproc ${bowtie_idx_dir}${rrna_core} -f ${nodup_fa} --un=${ndup_nrrna_fa} > /dev/null
#=================================
# step 4: map riboseq with bowtie
#=================================
transcript_core=${transcript_fa##*/}
transcript_core=${transcript_core%.*}
if  [[ ! $(ls "${bowtie_idx_dir}${transcript_core}"*) ]];  then
    echo "bowtie index for transcriptome not exist, build it"
    bowtie-build -f ${transcript_fa} ${bowtie_idx_dir}${transcript_core}
fi
bam_out="${tmp_dir}${riboseq_core}_nodup.bam"
#bowtie -p $nproc --chunkmbs 300 -a --best --strata -m 255 -n 1 ${bowtie_idx_dir}${transcript_core} -f ${ndup_nrrna_fa} -S | samtools view -bS -o ${bam_out} -
#=============================
# step 5: run ribomap
#=============================
./ribomap --bam ${bam_out} --fasta ${transcript_fa} --gtf ${transcript_gtf} --sf ${sf_odir}quant_bias_corrected.sf --out "${tmp_dir}${riboseq_core}.profile" --offset 15
