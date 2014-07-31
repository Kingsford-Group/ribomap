Overview
------
Ribomap is a package that generates isoform-level ribosome profiles from ribosome profiling data. Ribosome profiling is a recently developed high-throughput sequencing technique that captures approximately 30 bp long ribosome-protected mRNA fragments during translation. Because of alternative splicing and genomic repetitive sequences, a ribosome-protected read may map to many places in the transcriptome, leading to discarded or arbitrary mappings when standard approaches are used. Ribomap addresses this problem by assigning reads to potential origins in the transcriptome proportional to the estimated transcript abundance. This results in a more accurate estimation of the ribosome pileup compared to naive read assignment.

Compile from Source code
------
### Prerequisites
* [boost]:(http://www.boost.org/)
* [cereal (v1.0.0)]:(http://uscilab.github.io/cereal/)
* [seqan (v1.4.1)]:(http://www.seqan.de/)

### compile
> cd src
> make all
This will generate two executables: 
* `merge_fq_to_fa`: preprocess the ribosome profiling reads
* `ribomap`: assign ribosome profiling reads to transcript locations

Run Ribomap
------
under the `src` directory, run:
> ./run_ribomap [options]

# dataset
  paper: Guo, Huili, et al. "Mammalian microRNAs predominantly act to decrease target mRNA levels." Nature 466.7308 (2010): 835-840.
  GEO: GSE22004
  RNA-seq: GSM546921_filtered_sequence.txt.gz
  Ribo-seq: GSM546920_filtered_sequence.txt.gz


# run code:
./run_ribomap [options]
 