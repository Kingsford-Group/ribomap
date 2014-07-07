import itertools
import random
import copy
import math
import sys

mutList = {'A' : ['T', 'C', 'G'],
           'T' : ['A', 'C', 'G'],
           'C' : ['A', 'G', 'T'],
           'G' : ['A', 'C', 'T'],
           'a' : ['c', 'g', 't'],
           't' : ['a', 'c', 'g'],
           'c' : ['a', 'g', 't'],
           'g' : ['a', 'c', 't'],
           'N' : ['A', 'T', 'C', 'G'],
           'n' : ['a', 't', 'c', 'g']}

def mutate(base):
    ml = mutList[base]
    return ml[random.randint(0, len(ml)-1)]

def randCeilFloor(x):
    if x < 1.0 or random.random() < 0.5:
        return int(math.ceil(x))
    else:
        return int(math.floor(x))

def main(args):
    mutRate = float(args['--rate'])
    refIn = args['--in']
    refOut = args['--out']
    ifile = open(refIn, 'r')
    ofile = open(refOut,'w')
    currBase = 0
    nMut = 0
    totLen = 0
    i = 0
    nextMut = currBase + randCeilFloor(random.expovariate(mutRate))
    iline = ifile.readline()
    while iline:
        # read in one entry
        if iline.startswith("@"): 
            # description and readid seperated by space
            idx = iline.find(" ")
            # skip first char '@'
            rid = iline[1:idx]
            # skip last char '\n'
            description = iline[idx+1:-1] if idx>0 else ''
            seq_orig = ifile.readline().strip()
            ifile.readline()
            quality = ifile.readline()
            # get ready for processing the next read
            iline = ifile.readline()
            if len(seq_orig)==0: continue
            i += 1
        lenS = len(seq_orig)
        firstBase = currBase
        lastBase = currBase + lenS
        totLen += lenS
        seq_mut = list(seq_orig)
        while nextMut < lastBase:
            offset = nextMut - firstBase
            orig = seq_mut[offset]
            seq_mut[offset] = mutate(orig)
            rid += "_{0}{1}".format(offset, orig)
            nMut += 1
            currBase = nextMut
            nextMut = currBase + randCeilFloor(random.expovariate(mutRate))
        currBase = lastBase
        seq_mut = "".join(seq_mut)
        text = [ "@"+rid+" "+description,
                 seq_mut,
                 "+"+rid+" "+description,
                 quality ]
        ofile.write("\n".join(text))
        if totLen != i*lenS:
            print "firstBase: {0} lastBase: {1} lenS: {2} totLen: {3} currBase: {4} nextMut: {5}".format(firstBase, lastBase, lenS, totLen, currBase, nextMut)
            print rid, description, seq_orig
            break
        if (i % 1000 == 0):
            sys.stderr.write("processed {} records; performed {} mutations; rate = {:.2f}%\r\r".format(i, nMut, (100.0 * nMut) / totLen))
    ifile.close()
    ofile.close()
    print "firstBase: {0} lastBase: {1} lenS: {2} totLen: {3} currBase: {4} nextMut: {5}".format(firstBase, lastBase, lenS, totLen, currBase, nextMut)
    print i,
    print("\ndone.")

if __name__ == "__main__":
    if len(sys.argv) != 4 :
        print "Usage: mutate_fq.py --in=<input> --out=<output> --rate=<rate>"
        exit(1)
    arguments = {}
    for argv in sys.argv[1:]:
        k,v = argv.split("=")
        arguments[k] = v
    main(arguments)
