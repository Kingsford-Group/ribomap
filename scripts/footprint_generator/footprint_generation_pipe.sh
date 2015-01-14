#!/bin/bash
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/hela/
sm_fn=/home/hw1/scratch/scratch2/ribomap-playground/hela/sm_quant/quant_bias_corrected.sf
ref_old="${work_dir}"ref/gencode.v18.pc_transcripts_filter.fa
cds_range="${work_dir}"ref/gencode.v18.pc_transcripts_cds.txt
rlsim_exe=/home/hw1/code_store/rlsim/rlsim
output_dir="${work_dir}"synthetic_data/
ref_to_rlsim="${output_dir}"ref_with_abd.fa
frags_from_rlsim="${output_dir}"synth_rnaseq_frags.fa
reads_from_rlsim="${output_dir}"synth_rnaseq.fq
read_cnt=100000000
nproc=15
read_len=30
offset=12
footprint_generator=./footprint_generator
erate_fn=elongation_rate_human.txt
ilow=0.03
ihigh=0.3
mkdir -p ${output_dir}
#===========================================
# RNA-seq generation
#===========================================
#===========================================
# Step 1: append abdance to reference fasta
#===========================================
#python make_ref_fa_with_abd.py ${sm_old_dir}quant_bias_corrected.sf ${ref_old} ${ref_to_rlsim}
#===========================================
# Step 2: synthetically generate fragments
#===========================================
#${rlsim_exe} -n ${read_cnt} -t ${nproc} ${ref_to_rlsim} > ${frags_from_rlsim}
#===========================================
# Step 3: trim fragments to read length
#===========================================
#python read_trimmer.py ${frags_from_rlsim} ${reads_from_rlsim} ${read_len}
#===========================================
# Ribo-seq generation
#===========================================
#===========================================
# Step 1: generate footprint reads
#===========================================
#${footprint_generator} ${ref_old} ${cds_range} ${sm_fn} ${erate_fn} ${ilow} ${ihigh} ${read_cnt} ${read_len} ${offset} ${output_dir}synth_riboseq.profile ${output_dir}synth_riboseq.fq ${nproc}
#===========================================
# Step 2: add errors to reads
#===========================================
for rate in 0.005 0.01 0.02; do
    python mutate_fq.py --in=${output_dir}synth_riboseq.fq --out=${output_dir}synth_riboseq_${rate}.fq --rate=${rate}
done

