import sys
import pysam
import scipy.stats
import random
import getopt
import os.path
from Bio.Seq import Seq
import re
import argparse

## NOTE: this may require python 2.7.5...

# Argument parsing:
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, metavar='<inFile>',
                        help='input bam file')
    parser.add_argument('-o', type=str, required=True, metavar='<outFile>',
                        help='output bam file')
    parser.add_argument('-b', type=str, required=True, metavar='<bcCountFile>',
                        help='barcode counts file')
    parser.add_argument('-n', type=int, required=True, metavar='<countThresh>',
                        help='minimum read threshold')

    args = parser.parse_args()
    return args

def run(args):
    inFile = args.i    # input BAM file
    outFile = args.o   # output BAM file
    bcFile = args.b    # barcode counts file (output by extractValidR2reads_wCorrection.py)
    minReads = args.n  # minimum numer of reads per barcode

    # load the valid barcode dictionary:
    fIn = open(bcFile,'r')
    bcSet = {}
    while 1:
        line = fIn.readline()
        if not line:
            break
        # skip the header line:
        if line.startswith('BC'):
            continue
        if line.endswith('\n'):
            line=line[:-1]

        fields = line.split('\t')
        bc = fields[0]
        count = int(fields[1])
        if count>=minReads:
            bcSet.setdefault(bc,count)

    fIn.close()
    print 'Number of useable barcodes: %d' % len(bcSet.keys())

    # open the input file:
    inBam = pysam.AlignmentFile(inFile, 'rb')

    # create the output file:
    outBam = pysam.AlignmentFile(outFile, 'wb', template=inBam) 

    countMod = 100000
    readSum = 0
    keepCount = 0
    discardCount = 0

    # loop over all reads:
    for rd in inBam.fetch():
        # get the read name:
        rName = rd.query_name
        fields = rName.split(':')   # split the query name on ':" to separate read name, BC and UMI
        readSum+=1
	if readSum%countMod == 0:
	    print readSum
        if len(fields)<4:
            print 'Error in read name format: %s' % rName
            continue
        # build the full barcode:
        bc = fields[-2]
        if bcSet.has_key(bc):
            # write the read to the output file
            outBam.write(rd)
            keepCount+=1
        else:
            discardCount+=1

    # close the input files:
    inBam.close()
    outBam.close()

    # Print counts:
    print 'Total input reads: %d, Kept %d reads, discarded %d reads, total reads: %d' % (readSum, keepCount, discardCount, keepCount+discardCount)

    return

    
if __name__ == "__main__":

    args = get_args()  # process the input line
    run(args)
