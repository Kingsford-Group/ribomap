Overview
------
Ribomap is a package that generates isoform-level ribosome profiles from ribosome profiling data. Ribosome profiling is a recently developed high-throughput sequencing technique that captures approximately 30 bp long ribosome-protected mRNA fragments during translation. Because of alternative splicing and genomic repetitive sequences, a ribosome-protected read may map to many places in the transcriptome, leading to discarded or arbitrary mappings when standard approaches are used. Ribomap addresses this problem by assigning reads to potential origins in the transcriptome proportional to the estimated transcript abundance. This results in a more accurate estimation of the ribosome pileup compared to naive read assignment.

Prerequisites for Ribomap
------
* [__Sailfish__ (v0.6.3)](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html) for transcript abundance estimation
* [__Bowtie__ (v1.1.0)](http://bowtie-bio.sourceforge.net/index.shtml) for read mapping

Compile from Source code
------
### Prerequisites
* [boost](http://www.boost.org/)
* [cereal (v1.0.0)](http://uscilab.github.io/cereal/)
* [seqan (v1.4.1)](http://www.seqan.de/)

### Compile
    cd src
    make all INC="-I/opt/local/include"

This will generate two executables: 
* `merge_fq_to_fa`: preprocess the ribosome profiling reads
* `ribomap`: assign ribosome profiling reads to transcript locations
Please add the path for the prerequistite headers with flag `INC="-I<path/to/include/>`

Run Ribomap
------
### Run Ribomap with transcript abundance estimation
`run_ribomap.sh` is a pipeline that takes in the riboseq data and the RNA-seq data and automatically estimates the transcript abundance, then assigns riboseq reads to transcript locations based on the estimated transcript abundance. 

Under the `src` directory, run:

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

    ./run_ribomap.sh  \
    --rnaseq=../data/GSM546921_filtered_sequence.fq \
    --riboseq=../data/GSM546920_filtered_sequence.fq

Please connect the parameter flags and the parameters with `=` (no space in between).

### Run Ribomap without transcript abundance estimation
Ribomap supports transcript abundance estimation files from [*Sailfish*](http://www.cs.cmu.edu/~ckingsf/software/sailfish/), [*Cufflinks*](http://cufflinks.cbcb.umd.edu/index.html) and [*eXpress*](http://bio.math.berkeley.edu/eXpress/overview.html). Mapping the ribosome footprint can be performed providing any of the three transcript abundance esitmation files listed above.

Under the `src` directory, run:

      ./ribomap [options]

The list of options are as follows:

* __-b | --bam__ Bam file of ribosome footprint read mapping to the transcriptome
* __-f | --fasta__ Transcriptome reference fasta file
* __-g | --gtf__ Transcriptome annotation gtf file
* __-s | --sf__ Transcript abundance estimation produced by Sailfish
* __-c | --cl__ Transcript abundance estimation produced by Cufflinks
* __-e | --ep__ Transcript abundance estimation produced by eXpress
* __-o | --out__ Output file name of the estimated ribosome profile
* __-p | --offset__ Offset location in a read that the ribosome P-site maps to

One example of using the executable:

    ./ribomap --bam ../tmp/GSM546920_filtered_sequence_nodup.bam \
    --fasta ../data/gencode.v18.pc_transcripts_filter.fa \
    --gtf ../data/gencode.v18.annotation.gtf \
    --sf ../sf_quant/quant_bias_corrected.sf \
    --out ../outputs/GSM546920_filtered_sequence.profile \
    --offset 15

Please connect the parameter flags with the parameters with space.

Ribomap output files
------
Ribomap produces three output files:
#### _XXX.profile_
The ribosome profiles for the expressed transcripts with other analysis statistics. Each entry of a specific transcript looks like this:
~~~~~~
  refID: 425
  tid: ENST00000234875.4
  rabd:  1358.65
  tabd:  3.88878e-07
  te: 3.49377e+09
  rprofile: 7 0 0 0 54 0 0 0 1 4 ...
~~~~~~  
* __refID__ The transcript fai index in the transcriptome fasta file.
* __tid__ Transcript ensemble ID.
* __rabd__ Total ribosome loads, which is the sum of the __rprofile__ vector.
* __tabd__ Relative transcript abundance per codon (alpha_m) derived from Sailfish’s result.
* __te__ Relative translational efficiency, which is the ratio between __rabd__ and __tabd__.	
* __rprofile__ Ribosome profile vector of the CDS regions of the transcript. Each number in the vector is the number of ribosome footprints that are estimated to be from the corresponding codon location.

#### _XXX_abundant.list_ 
A list of transcripts whose total ribosome abundance is more than expected given the transcript abundance. 
There is one transcript record per row. The columns are defined as follows:

| Column number | Description |
|---------------|-------------|
| 1 | transcript Ensembl ID | 
| 2 | relative transcript abundance |
| 3 | total ribosome footprint count |
| 4 | pencentile ranking of the transcript abundance |
| 5 | percentile ranking of the total ribosome footprint count |
| 6 | difference between the transcript abundance rank and the total ribosome footprint count rank

#### _XXX_scarce.list_
A list of transcripts whose total ribosome abundance is less than expected given the transcript abundance.
The file format is the same as _XXX_abundant.list_.

Test case
------
### Run test case
under the `src` directory, run:

      ./get_data.sh
      ./run_ribomap.sh

`get_data.sh` automatically downloads the transcriptome fasta, gtf, a RNA-seq data and a riboseq data. The transcriptome fasta file is preprocessed with a _python_ script `transcript_filter.py` to excludes the following transcript

1. transcripts without verified start codon
2. transcripts with stop codon in the middle
3. transcripts with duplicated sequences
4. peptide sequence length less than 3 after getting rid of start and end of the seq

### Test case data sets
* __RNA-seq__ [GSM546921](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz)
* __riboseq__ [GSM546920](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz)
* __human transcriptome reference fasta__ [gencode.v18.pc_transcripts.fa](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_transcripts.fa.gz)
* __human transcriptome annotation gtf__ [gencode.v18.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz)

