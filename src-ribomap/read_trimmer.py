#!/usr/bin/env python
"""
functions to trim fasta/fastq reads and saved the trimmed reads into fasta/fastq
"""
import sys

def trim_fasta(infn,outfn,trimlen):
    in_file = open(infn,"r")
    out_file = open(outfn, "w")
    for line in in_file:
        if line[0] == '>':
            line_buf = [line]
        else:
            line_buf.append(line[:trimlen]+'\n')
            out_file.writelines(line_buf)
    in_file.close()
    out_file.close()

def trim_fastq_zin(infn,outfn,trimlen):
    import gzip
    in_file = gzip.open(infn,"rb")
    out_file = open(outfn, "wb")
    line = in_file.readline()
    while(line):
        if line[0] == '@':
            line_buf = [line]
            line_buf.append(in_file.readline()[:trimlen]+'\n')
            line_buf.append(in_file.readline())
            line_buf.append(in_file.readline()[:trimlen]+'\n')
            out_file.writelines(line_buf)
            line = in_file.readline()
    in_file.close()
    out_file.close()

def trim_fastq(infn,outfn,trimlen):
    in_file = open(infn,"rb")
    out_file = open(outfn, "wb")
    line = in_file.readline()
    while(line):
        if line[0] == '@':
            line_buf = [line]
            line_buf.append(in_file.readline()[:trimlen]+'\n')
            line_buf.append(in_file.readline())
            line_buf.append(in_file.readline()[:trimlen]+'\n')
            out_file.writelines(line_buf)
            line = in_file.readline()
    in_file.close()
    out_file.close()

def trim_fastq_to_fasta(infn,outfn,trimlen):
    in_file = open(infn,"rb")
    out_file = open(outfn, "wb")
    line = in_file.readline()
    while(line):
        if line.startswith('@'):
            line_buf = ['>'+line[1:]]
            line_buf.append(in_file.readline()[:trimlen]+'\n')
            in_file.readline()
            in_file.readline()
            out_file.writelines(line_buf)
            line = in_file.readline()
    in_file.close()
    out_file.close()

def trim_fastq_to_fasta_merge_dup(infn,outfn,trimlen):
    seq2ids = {}
    print "reading in fastq file..."
    in_file = open(infn,'r')
    cnt = 0
    line = in_file.readline()
    while(line):
        if line.startswith('@'):
            header = line.lstrip('@').split()[0]
            seq = in_file.readline()[:trimlen]
            seq2ids.setdefault(seq,[]).append(header)
            cnt += 1
            if (cnt%1000==0): 
                sys.stdout.write("processed %d reads\r"%cnt)
                sys.stdout.flush()
        line = in_file.readline()
    in_file.close()
    print "\nwriting to fasta file..."
    out_file = open(outfn, 'w')
    for seq, headers in seq2ids.iteritems():
        header = ">{0}_{1}\n".format(headers[0],len(headers))
        out_file.writelines([header, seq+'\n'])
    out_file.close()

if __name__ == "__main__": print("main currently blank")

