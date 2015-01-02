#!/usr/bin/env python
def convert_codon_to_aa(aa2codon):
    return {codon:aa for aa in aa2codon for codon in aa2codon[aa]}

aa2codon = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    'U': ('TAA', 'TAG', 'TGA'),
}

codon2aa = convert_codon_to_aa(aa2codon)

def get_transcript_frame(gtfname):
    """
    gtf fields:
    0: seqname
    1: source
    2: feature
    3: start
    4: end
    5: score
    6: strand
    7: frame
    8: attribute
      gene_id ENSGXXXXXXXXXXX *
      transcript_id ENSTXXXXXXXXXXX *
      gene_type list of biotypes
      gene_status {KNOWN, NOVEL, PUTATIVE}
      gene_name string
      transcript_type list of biotypes
      transcript_status {KNOWN, NOVEL, PUTATIVE}
      transcript_name string
      level 1 (verified loci), 2 (manually annotated loci), 3 (automatically annotated loci)
    """
    gtfile = open(gtfname)
    tid2frame= {}
    for line in gtfile:
        if line.startswith("#"): continue
        words = line.rstrip('\n').split('\t')
        if words[2] == "CDS":
            frame = int(words[7])
            assert words[8].find("exon_number") != -1
            attributes = words[8].split("; ")
            for a in attributes:
                aname, aval = a.lstrip().split(" ")
                aval = aval.lstrip('"').rstrip('"')
                if aname == "transcript_id":
                    tid = aval
                if aname == "exon_number":
                    enum = int(aval)
            if tid in tid2frame:
                e_pre, f_pre = tid2frame[tid]
                if enum < e_pre:
                    tid2frame[tid] = (enum, frame)
            else:
                tid2frame[tid] = (enum, frame)
    gtfile.close()
    return {tid:frame for tid,(enum,frame) in tid2frame.iteritems()}

def encode_peptide(tseq):
    return "".join([ codon2aa[ tseq[i:i+3] ] for i in xrange(0,len(tseq),3) ])
    
def hamming_dist(a,b):
    d = 0
    slen = min(len(a),len(b))
    for i in xrange(slen):
        if a[i]!=b[i]:
            d += 1
    return d

def main():
    print "getting frame for all transcripts"
    tid2frame = get_transcript_frame("gencode.v18.annotation.gtf")
    print "attempting to encoding peptide and validating result..."
    tfile = open("gencode.v18.pc_transcripts.fa")
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
    pfile = open("gencode.v18.pc_translations.fa")
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
    diff = set([])
    tids = []
    tids_cpp = ["ENST00000361227.2", "ENST00000361335.1", "ENST00000361381.2", "ENST00000361390.2", "ENST00000361453.3", "ENST00000361567.2", "ENST00000361624.2", "ENST00000361681.2", "ENST00000361739.1", "ENST00000361789.2", "ENST00000361851.1", "ENST00000361899.2", "ENST00000362079.2"]
    Mcnt = 0
    Ucnt = 0
    small_pep_cnt = 0
    tcnt = 0

    while True:
        tline = tfile.readline()
        if not tline: break
        pline = pfile.readline()
        if not pline: break
        if tline.startswith(">") and pline.startswith(">"):
            twords = tline.lstrip(">").split("|")
            pwords = pline.lstrip(">").split("|")
            # make sure transcript id between two files match
            assert twords[0] == pwords[0]
            # make sure sequence length in header match the real length
            tseq = tfile.readline().strip()
            assert len(tseq) == int(twords[6]) 
            pseq = pfile.readline().strip()
            assert len(pseq) == int(pwords[6])
            if len(twords)<8: continue
            # get start, stop codon from transcript header
            for w in twords[7:]:
                if w.startswith("CDS:"):
                    start,stop = map(int,w.lstrip("CDS:").split("-"))
            # get the right in-frame sequence
            frame = tid2frame[twords[0]]
            plen = len(pseq)
            if frame: plen -= 1
            start += frame
            tlen = stop - start + 1
            if tlen%3 != 0: 
                stop -= tlen%3;
            stop_plen = start + plen*3 -1
            if stop - stop_plen !=3 and stop_plen != stop:
                print twords[0], stop, stop_plen, 
                stop = stop_plen

            pconvert = encode_peptide(tseq[start-1:stop])

            # re-adjust start and stop
            start += 3
            stop -= 3


            if (stop-start+1)/3 < 3:
                small_pep_cnt += 1

            tcnt += 1
            if pconvert.startswith('M'): Mcnt += 1
            if pconvert.endswith('U'): Ucnt += 1
            if tid2frame[twords[0]] != 0:
                pconvert = "X" + pconvert
            pconvert = pconvert.rstrip('U')
            if len(pconvert) != len(pseq):
                print 'len', twords[0], len(pconvert), len(pseq)

            print "{0} {1}-{2}".format(twords[0], start, stop)
            if tcnt>20: break

            # # only compare sequence after the first codon (skip start codon)
            # # don't care the encoding part that's longer than the peptide
            # if pconvert != pseq:
            #     tids.append(twords[0])
            #     print twords[0]
            #     print pconvert
            #     print pseq
            #     if twords[0] in tids_cpp:
            #         for i in xrange(len(pseq)):
            #             if pconvert[i] != pseq[i]:
            #                 diff.add((pconvert[i],pseq[i]))

    tfile.close()
    pfile.close()

    # print "transcripts whose encoding not agree with their peptides: "
    # print "total: ", len(tids)
    # #for t in tids:
    # #    print t,
    # print "\nmisinterpreted codons:"
    # print "total: ", len(diff)
    # for aa_convert, aa_original in diff:
    #     print aa_convert, aa_original, "|",
    # print "\n"

    print "# peptides starts with 'ATG': {0} ({1})".format(Mcnt, Mcnt/float(tcnt))
    print "# peptides ends with 'U': {0} ({1})".format(Ucnt, Ucnt/float(tcnt))
    print "tiny peptides (len<3) : {0} ({1})".format(small_pep_cnt, small_pep_cnt/float(tcnt))

if __name__ == "__main__": main()
