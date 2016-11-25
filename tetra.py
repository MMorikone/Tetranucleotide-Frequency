#Current developers are: Michael Morikone, Kathleen Conner
#Original pipeline uses R, conversion to entirely Python is underway
## Any line with '##' are comments from Senthil

##!/usr/bin/python -tt

## Senthil/UNLV October 15, 2012
## Given a multifasta file calculate tetranucleotide frequency
## for each of the sequence and output a 257 column, that contains
## the name of the sequence and the 256 possible tetranucleotide frequency
## Requires python version > 2.7 (for collections)

## July 14 2015: Modified since python3 version gives error

import sys
from random import choice
from Bio import SeqIO
from itertools import product, islice
from collections import Counter
import numpy
import pandas as pd
from pandas import DataFrame
import string
import subprocess
import os
import csv
from sklearn import manifold
from sklearn.metrics import euclidean_distances

def main():
    if len(sys.argv) != 2:
        sys.exit('Simple Tetra. Usage: python '
                 + sys.argv[0] + ' <in_fasta>')
##   iupac_N = ('A', 'T', 'G', 'C')
    infile = sys.argv[1]
    f = open(infile, 'rU')
    origlist = make_tetra()
    ## print tetramers once
    ## old print 'Contigs'
    print '\t',
    for x in origlist:
      print x +"\t",
    for i in SeqIO.parse(f,'fasta'):
        ## do something with
        ## i.id and i.seq
        ## rcseq reverse complement sequence
        ## fandrseq forward+reverse complement seq
        if len(i.seq) >= 2000:
            rcseq = i.seq.reverse_complement()
            fandrseq = (i.seq + rcseq)
            fandrseq = fandrseq.upper()
        ## print "\n", fandrseq
        ## bases = list(fandrseq)
        ## for base in range(len(bases)):
        ##     if bases[base] == 'N':
        ##         bases[base] = choice(iupac_N)
        ## fandrseq = ''.join(bases)
        ## print fandrseq
            movwindow = ["".join(x) for x in window(fandrseq, 4)]
            cfreq = calc_freq(movwindow, origlist)
            print_result(i.id, origlist, cfreq)
    print "\n"

    DataFrame = pd.read_csv(filepath_or_buffer='/home/pogo/Desktop/TNF_WIP/results/op9_mw_tetra.csv', sep='\t', names=None)	#input for pandas csv

#    seed = numpy.random.RandomState(seed=3)	#Test with default numpy global as well as fixed seed

    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity='precomputed', n_jobs=1) #from sklearn.manifold, now that everything else works work on implementation

    DataFrame.to_csv(path_or_buf='/home/pogo/Desktop/TNF_WIP/results/op9_mw_tetra2.csv', sep='\t', header=False, index = False) #output for pandas csv

def make_tetra():
    a = []
    for x in product('ATGC', repeat=4):
        a.append(''.join(x))
    return a

def window(seq, n=4):
    """ Returns a sliding window (of width n) over data from the
    iterable s -> (s0, s1, ... s[n-1]), (s1, s2, ..., sn), ... """
    it = iter(seq)
    result = tuple(islice(it, n))
## This does not make sense so commented out April 7 2015.
##    if len(result) == n:
##        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        

def calc_freq(mw, ol):
    cnt = Counter()
    for quad in mw:
        if quad in ol:
            cnt[quad] += 1
    return cnt

def print_result(cname, ol, cf):
    print "\n" +  cname + "\t",
    ## totaltet means the sum of occurrences of all tetranucleotide in
    ## the sequence
    totaltet = sum(cf.values())
    ## Nobel et al 1998
    for x in ol:
        print str(cf[x]/float(totaltet)) + "\t",

if __name__ == '__main__':
    main()
