#!/usr/bin/env python
"""
input fastq file
pipeline: 
1) transcript abundance on RNA-seq data using sailfish
2) ribo-seq preprocess: trim reads, merge duplicate count, and convert to fasta
3) trim reads and save to fasta (for bowtie-best)
4) filter out rrna reads
5) bowtie align
"""
import sys
import ConfigParser
import os
import glob
import read_trimmer

def get_core(s): return s.split("/")[-1].split(".")[0]

if len(sys.argv)!=2:
    print "Usage: python ribomap.py config_file"
    exit(1)

cfg_fn = sys.argv[1]
config = ConfigParser.SafeConfigParser()
config.read(cfg_fn)
transcript_fa = config.get("Data", "transcript_fa")
transcript_gtf = config.get("Data", "transcript_gtf")
rrna_fa = config.get("Data", "rrna_fa")
rnaseq_fq = config.get("Data", "rnaseq_fq")
riboseq_fq = config.get("Data", "riboseq_fq")
sf_idx_dir = config.get("Data", "sf_idx_dir")
sf_odir = config.get("Data", "sf_odir")
bowtie_idx_dir = config.get("Data","bowtie_idx_dir")
nproc = config.getint("Data", "nproc")
tmp_dir = config.get("Data", "tmp_dir")
seedlen = config.getint("Data","seedlen")
src_dir = config.get("Data","src_dir")
cpp_cmp = config.get("Data","cpp_cmp")
cpp_exe = config.get("Data", "cpp_exe")
output_dir = config.get("Data", "output_dir")
offset = config.getint("Data", "offset")

fasta_core = get_core(riboseq_fq)
rrna_core = get_core(rrna_fa)
ref_core = transcript_fa.split("/")[-1].rstrip(".fa")

def transcript_abundance(forced=False):
    """ transcript abundance estimation from RNA-seq with Sailfish"""
    if forced or len(glob.glob('{0}/*'.format(sf_idx_dir)))==0:
        cmd = "sailfish index -t {0} -o {1} -k 20 -p 15 -f".format(transcript_fa, sf_idx_dir)
        print cmd
        os.system(cmd)
    cmd = 'sailfish quant -l "T=SE:S=U" -i {0} -o {1} -r {2} -p {3} --no_bias_correct -a -f'.format(sf_idx_dir, sf_odir, rnaseq_fq, nproc)
    print cmd
    os.system(cmd)

def mkdirs():
    # make directories
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # change working directory
    os.chdir(tmp_dir)

#=========================================
# merge duplicate reads then bowtie align
#=========================================
def trim_fq_merge_dup():
    outfn = '{0}_{1}_nodup.fa'.format(fasta_core, seedlen)
    read_trimmer.trim_fastq_to_fasta_merge_dup(riboseq_fq,outfn,seedlen)

def filter_rrna_fa(fn):
    # Step 3: filter out rRNA footprint
    # create a fasta list that bowtie will read
    # build bowtie index
    if len(glob.glob('{0}{1}*'.format(bowtie_idx_dir,rrna_core)))==0:
        cmd = "bowtie-build -f {0} {1}{2}".format(rrna_fa,bowtie_idx_dir,rrna_core)
        print cmd
        os.system(cmd)
    # bowtie align
    cmd = "bowtie -p {0} {1}{2} -f {3}.fa --un={3}_norrna.fa >/dev/null".format(nproc, bowtie_idx_dir, rrna_core, fn)
    print cmd
    os.system(cmd)

def align2Bowtie1():
    if len(glob.glob("{0}{1}*".format(bowtie_idx_dir, ref_core)))==0:
        cmd = "bowtie-build -f {0} {1}{2}".format(transcript_fa, bowtie_idx_dir, ref_core)
        print cmd
        os.system(cmd)
    cmd = "bowtie -p {0} --chunkmbs 300 -a --best --strata -m 255 -n 1 {1}{2} -f {3}_{4}_nodup_norrna.fa -S | samtools view -bS -o {3}_nodup.bam -".format(nproc, bowtie_idx_dir, ref_core, fasta_core, seedlen)
    print cmd
    os.system(cmd)

#=========================================
# no merging for bowtie-best
#=========================================
def trim_fq_to_fa():
    outfn = '{0}_{1}.fa'.format(fasta_core, seedlen)
    read_trimmer.trim_fastq_to_fasta(riboseq_fq,outfn,seedlen)

def align2Bowtie1best():
    if len(glob.glob("{0}{1}*".format(bowtie_idx_dir, ref_core)))==0:
        cmd = "bowtie-build -f {0} {1}{2}".format(transcript_fa,bowtie_idx_dir,ref_core)
        print cmd
        os.system(cmd)
    cmd = "bowtie -p {0} --chunkmbs 300 --best -m 255 -n 1 {1}{2} -f {3}_{4}_norrna.fa -S | samtools view -bS -o {3}_best.bam -".format(nproc, bowtie_idx_dir, ref_core, fasta_core, seedlen)
    print cmd
    os.system(cmd)

def run_ribomap():
    # Step 5: read count summary
    os.chdir(src_dir)
    if not os.path.exists(cpp_exe):
        print cpp_cmp
        os.system(cpp_cmp)
    # Usage: ./ribomap input_dir bam_core transcript_fasta gtf_fn sailfish_result output_dir footprint_offset
    cpp_cmd = "./{0} {1} {2} {3} {4} {5}quant.sf {6} {7}".format(cpp_exe, tmp_dir, fasta_core, transcript_fa, transcript_gtf, sf_odir,output_dir,offset)
    print cpp_cmd
    os.system(cpp_cmd)

#=========================================
# not in current pipeline
#=========================================
def sra2fasta(fastq_exe, sra_dir, sra_fn):
    # not used in human data
    # Step 1: convert sra into fasta
    for sra_fn in sra_fn_list:
        cmd = "{0} -O {1} --fasta {2}.sra".format(fastq_exe, fasta_dir,sra_dir+sra_fn)
        os.system(cmd)

def trim_seqs(sra_fn):
    # not used in human data
    # Step 2: trim the sequence before aligning it to the reference genome
    for sra_fn in sra_fn_list:
        cmd = "python {0}trim_fasta.py {1}.fasta {1}_{2}.fasta {2}".format(src_dir,sra_fn,seedlen)
        os.system(cmd)

def filter_rrna_fq():
    # Step 3: filter out rRNA footprint
    # create a fasta list that bowtie will read
    fn = "{0}_{1}.fq".format(fasta_fn,seedlen)
    # build bowtie index
    if len(glob.glob('{0}{1}*'.format(bowtie_idx_dir,rrna_fn)))==0:
        cmd = "bowtie-build -f {0}.fasta {1}{0}".format(rrna_fn,bowtie_idx_dir)
        print cmd
        os.system(cmd)
    # bowtie align
    cmd = "bowtie -p {0} {1}{2} -q {3} --un={4}_norrna.fq >/dev/null".format(nproc, bowtie_idx_dir, rrna_fn, fn,fasta_fn)
    print cmd
    os.system(cmd)

def plot_mp():
    # Step 6: look at histogram of mappabilities
    os.chdir(src_dir)
    cmd = "python {0}{1} {2}".format(src_dir,plot_script,cfg_fn)
    os.system(cmd)

if __name__ == "__main__":
    filter_rrna_nodup = lambda: filter_rrna_fa("{0}_{1}_nodup".format(fasta_core,seedlen))
    filter_rrna = lambda: filter_rrna_fa("{0}_{1}".format(fasta_core,seedlen))
    task_list = [ transcript_abundance, trim_fq_merge_dup, filter_rrna_nodup, align2Bowtie1, trim_fq_to_fa, filter_rrna, align2Bowtie1best, run_ribomap]
    start  = 1
    end = 8
    mkdirs()
    for task in task_list[start:end]: task()

