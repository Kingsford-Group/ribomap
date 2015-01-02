#!/usr/bin/env python
from transcript_codon_check import get_transcript_frame
"""
sequence not included in fasta if
(1) transcripts without verified start codon
(2) transcripts with stop codon in the middle
(3) transcripts with duplicated sequences
(4) peptide sequence length less than 3 after getting rid of start and end of the seq
"""

def check_mid_stop(pseq,stop_codon="U"):
    """ return true if there are stop codon in the middle of the peptide sequence"""
    idx_stop = pseq.find(stop_codon)
    return (idx_stop != -1 and idx_stop != len(pseq)-1)

def check_tseq_dup(tseq,tid,tseq2tids):
    """ return true if transcript sequence is seen before"""
    tseq2tids.setdefault(tseq,[]).append(tid)
    return len(tseq2tids[tseq])>1

def peptide_len(twords,frame):
    # get start, stop codon from transcript header
    for w in twords[7:]:
        if w.startswith("CDS:"):
            start,stop = map(int,w.lstrip("CDS:").split("-"))
    # get the right in-frame sequence
    start += frame
    tlen = stop - start + 1
    if tlen%3 != 0:
        stop -= tlen%3;
    # ignore the first and the last codon
    start += 3
    stop -= 3
    return (stop-start+1)/3


# global varialbles: duplicate sequence transcript IDs
tseq2tids = {}

print "getting frame for all transcripts"
tid2frame = get_transcript_frame("../data/gencode.v18.annotation.gtf")
print "filtering fasta file..."
itfile = open("../data/gencode.v18.pc_transcripts.fa","r")
"""
transcript.fa header:
0 transcript-id|
1 gene-id|
2 Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
3 Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
4 transcript-name|
5 gene-name|
6 sequence-length|
7 5'-UTR (3'-UTR if reverse strand) location in the transcript|
8 CDS location in the transcript|
9 3'-UTR (5'-UTR if reverse strand) location in the transcript
"""
ipfile = open("../data/gencode.v18.pc_translations.fa","r")
"""
peptide.fa header
0 transcript-id|
1 gene-id|
2 Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
3 Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
4 transcript-name|
5 gene-name|
6 sequence-length
"""
otfile = open("../data/gencode.v18.pc_transcripts_filter.fa","w")
opfile = open("../data/gencode.v18.pc_translations_filter.fa","w")
while True:
    tline = itfile.readline()
    if not tline: break
    pline = ipfile.readline()
    if not pline: break
    if tline.startswith(">") and pline.startswith(">"):
        twords = tline.lstrip(">").split("|")
        pwords = pline.lstrip(">").split("|")
        # make sure transcript id between two files match
        assert twords[0] == pwords[0]
        # make sure sequence length in header match the real length
        tseq = itfile.readline()
        pseq = ipfile.readline()
        # skip if
        # (1) transcripts without verified start codon
        # (2) transcripts with stop codon in the middle
        # (3) transcripts with duplicated sequences
        # (4) peptide sequence length less than 3 after getting rid of start and end of the seq
        if len(twords)<8 or check_mid_stop(pseq) or check_tseq_dup(tseq, twords[0], tseq2tids) or peptide_len(twords,tid2frame[twords[0]])<3:
            continue
        otfile.write(tline)
        otfile.write(tseq)
        opfile.write(pline)
        opfile.write(pseq)
itfile.close()
ipfile.close()
otfile.close()
opfile.close()
