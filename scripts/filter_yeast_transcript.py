#!/usr/bin/env python
"""
prefilter yeast transcriptome
prerequest python module: Biopython
"""
from translation import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def get_cds_boundries(cds_fa, utr_fa):
    # get cds range on the genomic coordinates
    print "getting cds ranges..."
    cds_dna = {}
    tf = open(cds_fa)
    for line in tf:
        if line.startswith('>'):
            tid = line.lstrip('>').split(' ')[0]
            istart = line.find("from ")+len("from ")
            ispan = line[istart: ].find(", ")
            exon_str = line[istart: istart+ispan].split(',')
            exon_str = [ e.split('-') for e in exon_str ]
            cds_start = int(exon_str[0][0])
            cds_end = int(exon_str[-1][1])
            cds_dna[tid] = (cds_start, cds_end)
    tf.close()
    # get utr range on the genomic coordinates
    utr_dna = {}
    tf = open(utr_fa)
    for line in tf:
        if line.startswith('>'):
            tid = line.lstrip('>').split(' ')[0]
            istart = line.find("from ")+len("from ")
            ispan = line[istart: ].find(", ")
            exon_str = line[istart: istart+ispan].split('-')
            utr_start = int(exon_str[0])
            utr_end = int(exon_str[1])
            utr_dna[tid] = (utr_start, utr_end)
    tf.close()
    # compute cds range in the fasta coordinates
    cds_range = {}
    for tid, (cds_start, cds_end) in cds_dna.iteritems():
        (utr_start, utr_end) = utr_dna[tid]
        start = cds_start - utr_start
        end = cds_end - utr_end
        if start<0:
            cds_range[tid] = (-start, -end)
        else:
            cds_range[tid] = (start, end)
    return cds_range

def exclude_gene_with_introns(fname, glist):
    tf = open(fname)
    tf.readline()
    for line in tf:
        gid = line.strip().split('\t')[1]
        if gid in glist:
            glist[gid] = False
    tf.close()
    return glist

def exclude_overlapping_genes(fname, glist):
    print "excluding overlapping genes..."
    tf = open(fname)
    tf.readline()
    for line in tf:
        words = line.strip().split('\t')
        tid1 = words[1]
        tid2 = words[9]
        if tid1 in glist and tid2 in glist:
            glist[tid1] = False
            glist[tid2] = False
    tf.close()
    return glist

def compare_pseq(pconvert, pseq, stop_codon):
    if pconvert[-1] != stop_codon:
        pconvert += stop_codon;
    l = min(len(pconvert), len(pseq))
    cset = set([])
    for i in xrange(l):
        if pconvert[i] != pseq[i]:
            cset.add((pconvert[i], pseq[i]))
        return cset

def transcript_codon_check(tfname, pfname, glist, cds_range):
    """
    exclude :
    1) transcripts that I cannot correctly encode to peptides
    2) transcripts with stop codon in the middle
    3) transcripts with duplicated sequences
    4) peptide sequence length less than 3 after getting rid of start and end of the seq
    """
    Mcnt = 0
    Ucnt = 0
    tcnt = 0
    ms_tid = []
    short_tid = []
    dup_tid = []
    failed_tid = []
    ccset = set([])
    tseq2tids = {}
    print "checking transcripts..."
    tfile = open(tfname, "rU")
    pfile = open(pfname, "rU")
    for trec, prec in zip(SeqIO.parse(tfile, "fasta"), SeqIO.parse(pfile, "fasta")):
        tcnt += 1
        tid = trec.id
        assert(tid == prec.id)
        if not glist[tid]: continue
        if cds_range:
            start, stop = cds_range[tid]
            tseq = str(trec.seq[start:stop]) 
        else:
            tseq = str(trec.seq)
        pconvert = encode_peptide(tseq, codon2aa)
        if is_start_codon(tseq[:3]): Mcnt += 1
        if pconvert.endswith(stop_codon): Ucnt += 1
        if check_mid_stop(pconvert, stop_codon):
            glist[tid] = False
            ms_tid.append(tid)
        elif check_tseq_dup(tseq, tid, tseq2tids):
            glist[tid] = False
            dup_tid.append(tid)
        elif peptide_len(tseq, codon2aa, stop_codon)<3:
            glist[tid] = False
            short_tid.append(tid)
        else:
            pseq = str(prec.seq)
            cset = compare_pseq(pconvert, pseq, stop_codon)
            if cset:
                glist[tid] = False
                failed_tid.append(tid)
                ccset |= cset
                print  "{0} pconvert_len: {1}, pseq_len: {2}".format(tid, len(pconvert), len(pseq))
                print pconvert
                print pseq
                print cset
    tfile.close()
    pfile.close()
    print "# peptides with start codon: {0} ({1})".format(Mcnt, Mcnt/float(tcnt))
    print "# peptides ends with stop codon: {0} ({1})".format(Ucnt, Ucnt/float(tcnt))
    print "# middle stop codon: {0} ({1})".format(len(ms_tid), len(ms_tid)/float(tcnt))
    print "# duplicate transcripts: {0} ({1})".format(len(dup_tid), len(dup_tid)/float(tcnt))
    print "# tiny peptides (len<3) : {0} ({1})".format(len(short_tid), len(short_tid)/float(tcnt))
    print "{0} transcripts failed to encode correctly".format(len(failed_tid))
    print "codons that are encoded wrong: ", ccset
    print "total transcripts: {0}".format(tcnt)
    return glist

def build_filtered_cds_with_boundry(cds_fn, utr_fn, glist, cds_range, otfa, boundry):
    # exclude transcripts that I cannot correctly encode to peptides
    print "building filtered transcript fasta with boundries"
    cds_file = open(cds_fn, "rU")
    utr_file = open(utr_fn, "rU")
    otfile = open(otfa, "w")
    for cds_rec, utr_rec in zip(SeqIO.parse(cds_file, "fasta"), SeqIO.parse(utr_file, "fasta")):
        tid = cds_rec.id
        assert(tid == utr_rec.id)
        if not glist[tid]: continue
        start, end = cds_range[tid]
        utr_seq = str(utr_rec.seq)
        if start >= boundry: 
            utr5 = utr_seq[start-boundry : start]
        else: 
            utr5 = utr_seq[:start]
        if len(utr_seq[end: ]) >= boundry:
            utr3 = utr_seq[end : end+boundry]
        else: 
            utr3 = utr_seq[end:]
        tseq = utr5+str(cds_rec.seq)+utr3
        start = len(utr5)
        end = len(tseq) - len(utr3)
        cds_range[tid] = (start, end)
        record = SeqRecord(Seq(tseq,IUPAC.unambiguous_dna),id=tid, description="length: {0} | CDS: {1}-{2}".format(len(tseq), start, end))
        SeqIO.write(record, otfile, "fasta")
    cds_file.close()
    utr_file.close()
    otfile.close()
    return cds_range

def record_cds_range(ofn, tid2cds, glist):
    print "writing cds file..."
    ofile = open(ofn, 'w')
    for tid in tid2cds:
        if glist[tid]:
            start, stop = tid2cds[tid]
            ofile.write("{0}\t{1}\t{2}\n".format(tid, start, stop))
    ofile.close()

def main():
    dirc = "/home/hw1/scratch/ribojamdetector/transcriptome/"
    cds_fa = dirc+"orf_coding.fasta"
    transcript_fa = dirc+"orf_genomic_1000.fasta"
    peptide_fa = dirc+"orf_trans.fasta"
    overlapping_gene_fn = dirc+"overlapping_genes.tsv"
    #intron_gene_fn = dirc+"gene_with_introns.tsv"
    boundry = 33
    filtered_tfa = dirc+"protein_coding_{0}_filtered.fasta".format(boundry)
    filtered_pfa = dirc+"protein_coding_trans_{0}_filtered.fasta".format(boundry)
    cds_fn = dirc+"cds_range.txt"

    global stop_codon, codon2aa
    stop_codon = '*'
    codon2aa = build_codon_to_aa(stop_codon)
    cds_range = get_cds_boundries(cds_fa, transcript_fa)
    glist = get_gene_list(transcript_fa)
    glist = exclude_overlapping_genes(overlapping_gene_fn,glist)
    glist = transcript_codon_check(cds_fa, peptide_fa, glist, None)
    print "transcripts included in study: {0}".format(sum(glist.values()))
    cds_range = build_filtered_cds_with_boundry(cds_fa, transcript_fa, glist, cds_range, filtered_tfa, boundry)
    record_cds_range(cds_fn, cds_range, glist)
    filter_fasta(peptide_fa, filtered_pfa, glist)
    transcript_codon_check(filtered_tfa, filtered_pfa, glist, cds_range)

if __name__ == "__main__": main()
