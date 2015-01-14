#!/usr/bin/env python
"""
functions to trim fasta/fastq reads and saved the trimmed reads into fasta/fastq
"""
import sys

def trim_fasta(infn,outfn,trimlen):
    ifile = open(infn,"r")
    ofile = open(outfn, "w")
    line = ifile.readline()
    while line:
        if line.startswith('>'):
            header = line
            seq = ifile.readline()
            if len(seq)>trimlen+1:
                seq = seq[:trimlen]+'\n'
                ofile.writelines([header, seq])
        line = ifile.readline()
    ifile.close()
    ofile.close()

def trim_fasta_to_fastq(infn,outfn,trimlen):
    ifile = open(infn,"r")
    ofile = open(outfn, "w")
    line = ifile.readline()
    while line:
        if line.startswith('>'):
            header = line.lstrip('>')
            seq = ifile.readline()
            if len(seq)>trimlen+1:
                seq = seq[:trimlen]+'\n'
                ofile.writelines(['@'+header, seq, '+'+header, ''.join(['a']*trimlen)+'\n' ])
        line = ifile.readline()
    ifile.close()
    ofile.close()

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

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print "Usage: python read_trimmer.py input_fa output_fa trimlen"
        exit(1)

    infn = sys.argv[1]
    outfn = sys.argv[2]
    trimlen = int(sys.argv[3])
    trim_fasta_to_fastq(infn,outfn,trimlen)

