Overview
------
Ribomap is a package that generates isoform-level ribosome profiles from ribosome profiling data. Ribosome profiling is a recently developed high-throughput sequencing technique that captures approximately 30 bp long ribosome-protected mRNA fragments during translation. Because of alternative splicing and genomic repetitive sequences, a ribosome-protected read may map to many places in the transcriptome, leading to discarded or arbitrary mappings when standard approaches are used. Ribomap addresses this problem by assigning reads to potential origins in the transcriptome proportional to the estimated transcript abundance. This results in a more accurate estimation of the ribosome pileup compared to naive read assignments.

Prerequisites for Ribomap
------
<!---
* [__FASTX-Toolkit__] (http://hannonlab.cshl.edu/fastx_toolkit/index.html) for preprocessing reads
-->
* [__Sailfish__ (latest improved version: Salmon v0.2.3)](https://github.com/kingsfordgroup/sailfish/releases/tag/v0.2.3) for transcript abundance estimation
* [__STAR__ (v2.4.0h)](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.0h1) for read mapping

Compile from Source code
------
### Prerequisites
* [boost](http://www.boost.org/)
* [seqan (v1.4.1)](http://www.seqan.de/)

### Compile
a c++ compiler that support c++11 features (for instance g++ >= 4.7) is required.

    cd src
    make riboprof INC="-I/opt/local/include"
    make install

This will generate a c++ executable `riboprof` that assign ribosome profiling reads to transcript locations and copy the executable to the `bin` directory.

Please add the path for the prerequistite headers with flag `INC="-I<path/to/include/>"`

Run Ribomap
------
#### Run Ribomap with automatic transcript abundance estimation
`run_ribomap.sh` is an automatic pipeline for ribosome profiling data. It takes in the riboseq data and the RNA-seq data and automatically estimates the transcript abundance, then assigns riboseq reads to transcript locations based on the estimated transcript abundance. 

Under the `scripts` directory, run:

      ./run_ribomap.sh [options]
The list of options are as follows:
* __--rnaseq_fq__ (required) Input RNA-seq read fastq.gz file for transcript abundance estimation.
* __--riboseq_fq__ (required) Input ribosome profiling (riboseq) read fastq.gz file.
* __--transcript_fa__ (required) Input trascriptome reference fasta file.
* __--contaminant_fa__ Input contaminant sequence fasta file (human ribosome RNA sequences are included in directory `data`).
* __--cds_range__ A text file that includes the coding sequence (CDS) range for all transcripts (see description below). If such an option is not provided, the transcript fasta file is assume to only include the CDS regions.
* __--work_dir__ (default the parent directory of `scripts`) The working directory where all intermediate and final results will write to.
* __--nproc__ (default 15) Number of threads can be used by ribomap.
* __--adapter__ (default `CTGTAGGCACCATCAAT`) The linker sequence attached to the 5' end of the ribo-seq reads.
* __--nmismatch__ (default 1) Number of mismatches allowed in the read alignments.
* __softClipping__ (default `true`) Whether reads are allowed to be soft-clipped by STAR when aligning to the transcriptome.
* __--min_fplen__ (default 27) Minimun read length to keep for downstream analysis.
* __--max_fplen__ (default 33) Maximum riboseq read length to keep for downstream analysis.
* __--offset__ (default 12) Offset location in a read that the ribosome P-site maps to, or a text file name that defines the P-site offset based on read length (see description below).
*__--rnaUnstranded__ (default `false`) Whether the RNA-seq protocol is stranded. If the RNA-seq protocol is unstranded, the `librarytype` to run Salifish is set to `-l U`; otherwise the `librarytype` is set to `-l SF`, and alignments with the RC flag set in the RNA-seq data are discarded.
* __--tabd_cutoff__ (default 0) Transcript abundance threshold to be considered expressed..
*__--useSecondary__ (default `true`) Whether multi-mapping alignments are used when assigning footprints to candidate loci.
* __--star_idx_dir__ (default `$work_dir/StarIndex/`) Directory to store Star index.
* __--alignment_dir__ (default  `$work_dir/alignment/`) Directory to store alignment results output by STAR.
* __--sailfish_dir__ (default `$work_dir/sm_quant/`) Directory to store sailfish result.
* __--output_dir__ (default `$work_dir/outputs/`) Directory to store ribomap's outputs.
* __--force__ Force ribomap to regenerate all intermediate steps.
One example of using the shell script:
~~~~~~
    ./run_ribomap.sh \
    --rnaseq_fq rnaseq.fq.gz \
    --riboseq_fq riboseq.fq.gz \
    --contaminant_fa contaminant.fa \
    --transcript_fa transcript.fa \
    --cds_range cds_range.txt
~~~~~~
Please connect the parameter flags and the parameters with a white space.

* __CDS range file__ A plain text file that includes the CDS regions of transcriptome. Each line in the file should be in the following format:

       `transcript_id start stop`

`transcript_id` should be consistent with the transcript ID in the fasta file, `start` is the start base of the coding region in the transcript fasta file (zero-based), and `stop` is one base pass the stop position of the coding sequence. 

* __offset file__ A plain text file that describes the read length and the P site offset for that read length. Each line in the file should be in the following format:

  	   `read_length P-site_offset`

Only a proper range of read length should be included, reads with a length not specified in this file will be discarded for downstream analysis.

#### Run Ribomap by providing the transcript abundance estimation file
Ribomap supports transcript abundance estimation files from [*Sailfish*](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html), [*Cufflinks*](http://cufflinks.cbcb.umd.edu/index.html) and [*eXpress*](http://bio.math.berkeley.edu/eXpress/overview.html). Mapping the ribosome footprint can be performed providing any of the three transcript abundance esitmation files listed above.

Under the `bin` directory, run:

      ./riboprof [options]

The list of options are as follows:

* __-r | --ribobam__ Bam file of ribo-seq read mappings to the transcriptome
* __-m | --mrnabam__ Bam file of RNA-seq read mappings to the transcriptome
* __-f | --fasta__ Transcriptome reference fasta file
* __-cds | --cds_range__ CDS range file
* __-s | --sf__ Transcript abundance estimation produced by Sailfish
* __-c | --cl__ Transcript abundance estimation produced by Cufflinks
* __-e | --ep__ Transcript abundance estimation produced by eXpress
* __-o | --out__ Output file prefix of Ribomap's result
* __-p | --offset__ Offset location in a read that the ribosome P-site maps to, or the name of a offset file that specifies P-site offset for different read length

One example of using the executable:
~~~~~~
    ./riboprof  \
    --mrnabam mRNA.bam --ribobam ribo.bam \
    --fasta transcript.fa --cds_range cds_range.txt \
    --sf quant.sf \
    --offset 12 \
    --out ../outputs/ribomap\n
~~~~~~
Please connect the parameter flags and the parameters with a space.

Ribomap output files
------
Ribomap produces five output files:
#### _XXX.codon_
The ribosome profiles within the CDS of reach transcript. Each entry of a specific transcript looks like this:
~~~~~~
	refID: 0
	tid: YAL001C
	ribo profile: 0 0 0 74 68 ...
	mRNA profile: 31 35 50 73 87 96 104 ...
	normalized ribo profile: 0 0 0 1.0137 0.781609 0.0208333 0.125 ...
~~~~~~  
* __refID__ The transcript fai index in the transcriptome fasta file.
* __tid__ Transcript ensemble ID.
* __ribo profile__ Ribosome profile vector of the CDS regions of the transcript. Each number in the vector is the number of ribosome footprints that are estimated to be from the corresponding codon location.
* __mRNA profile__ RNA-seq profile vector of the CDS regions of the transcript. Each number in the vector is the read coverage count that are esimated on the corresponding codon location.
* __normalized ribo profile__ Ribosome profile vector of the CDS regions of the transcript after bias correction. Each number in the vetor is the ratio between the ribo profile count and the mRNA profile count

#### _XXX.base_
The nucleotide-level ribosome profiles including the UTR regions. The file format is exactly the same as the the _XXX.codon_ file.

#### _XXX.stats_
The summarized statistics for each transcripts. Each entro of a specific transcript looks like this: 
~~~~~
	refID: 0
	tid: YAL001C
	rabd: 3959
	tabd: 0.000209384
	te: 1.89078e+07
~~~~~~
* __refID__ The transcript fai index in the transcriptome fasta file.
* __tid__ Transcript ensemble ID.
* __rabd__ Ribosome loads, which is the total number of ribosome reads that are esimated from this trascript.
* __tabd__ Relative transcript abundance from Sailfish’s result.
* __te__ Relative translational efficiency, which is the ratio between __rabd__ and __tabd__.


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
### Run test case (TODO)
under the `script` directory, run:

      ./get_data.sh
      ./run_ribomap.sh

`get_data.sh` automatically downloads the transcriptome fasta, gtf, a RNA-seq data and a riboseq data. The transcriptome fasta file is preprocessed with a _python_ script `transcript_filter.py` to excludes the following transcripts:

1. transcripts without verified start codon
2. transcripts with stop codon in the middle of the CDS region
3. transcripts with duplicated sequences
4. peptide sequence length less than 3 when the start and stop codons are not included

### Test case data sets
* __RNA-seq__ [GSM546921](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz)
* __riboseq__ [GSM546920](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz)
* __human transcriptome reference fasta__ [gencode.v18.pc_transcripts.fa](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_transcripts.fa.gz)
* __human transcriptome annotation gtf__ [gencode.v18.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz)

