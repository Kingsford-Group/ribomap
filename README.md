Overview
------
Ribomap is a package that generates isoform-level ribosome profiles from ribosome profiling data. Ribosome profiling is a recently developed high-throughput sequencing technique that captures approximately 30 bp long ribosome-protected mRNA fragments during translation. Because of alternative splicing and genomic repetitive sequences, a ribosome-protected read may map to many places in the transcriptome, leading to discarded or arbitrary mappings when standard approaches are used. Ribomap addresses this problem by assigning reads to potential origins in the transcriptome proportional to the estimated transcript abundance. This results in a more accurate estimation of the ribosome pileup compared to naive read assignment.

Prerequisites packages
------
* [__Sailfish__](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html) for transcript abundance estimation
* [__Bowtie__](http://bowtie-bio.sourceforge.net/index.shtml) for read mapping

Compile from Source code
------
### Prerequisites
* [boost](http://www.boost.org/)
* [cereal (v1.0.0)](http://uscilab.github.io/cereal/)
* [seqan (v1.4.1)](http://www.seqan.de/)

### compile
    cd src
    make all

This will generate two executables: 
* `merge_fq_to_fa`: preprocess the ribosome profiling reads
* `ribomap`: assign ribosome profiling reads to transcript locations

Run Ribomap
------
under the `src` directory, run:

      ./run_ribomap.sh [options]
The list of options are as follows:
* __--rnaseq__ Input RNA-seq read fastq file for transcript abundance estimation
* __--riboseq__	     Input ribosome profiling (riboseq) read fastq file
* __--ref_fa__ Input trascriptome reference fasta file
* __--gtf__ Input transcriptome	  annotation gtf file
* __--rrna_fa__	(default ribomap/data/rrna_human.fasta) Input ribosome RNA sequence fasta file (human ribosome RNA sequences are included in directory `data`)
* __--seedlen__ (default 25) Seed length to trim down the riboseq reads
* __--offset__ (default 15) Offset location in a read that the ribosome P-site maps to
* __--nproc__ (default 15) Number of threads can be used by ribomap

One example of using the shell script:
    ./run_ribomap.sh --rnaseq=../data/GSM546921_filtered_sequence.fq --riboseq=../data/GSM546920_filtered_sequence.fq

Test case
------
### run test case
under the `src` directory, run:

      ./get_data.sh
      ./run_ribomap.sh

`get_data.sh` automatically downloads the transcriptome fasta, gtf, a RNA-seq data and a riboseq data. The transcriptome fasta file is preprocessed with a _python_ script `transcript_filter.py` to excludes the following transcript

1. transcripts without verified start codon
2. transcripts with stop codon in the middle
3. transcripts with duplicated sequences
4. peptide sequence length less than 3 after getting rid of start and end of the seq

### test case data sets
* __RNA-seq__ [GSM546921](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz)
* __riboseq__ [GSM546920](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz)
* __human transcriptome reference fasta__ [gencode.v18.pc_transcripts.fa](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_transcripts.fa.gz)
* __human transcriptome annotation gtf__ [gencode.v18.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz)

