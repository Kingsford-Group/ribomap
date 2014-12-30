#!/usr/bin/env python
import sys
if len(sys.argv)!=5:
    print "Usage: python filter_reads_by_size.py ifname ofname min_len max_len"
    exit(1)

ifname = sys.argv[1]
ofname = sys.argv[2]
min_len = int(sys.argv[3])
max_len = int(sys.argv[4])

ifile = open(ifname, 'r')
ofile = open(ofname, 'w')

for line in ifile:
    if line.startswith('>'):
        header = line
    else:
        readlen = len(line)-1
        if readlen <= max_len and readlen >= min_len:
            ofile.write(header)
            ofile.write(line)
ifile.close()
ofile.close()
