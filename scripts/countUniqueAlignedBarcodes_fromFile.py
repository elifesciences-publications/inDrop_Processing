import sys
import pysam
import getopt
import os.path
import glob
import argparse
from Bio.Seq import Seq

MINCOUNT = 500      # barcodes with fewer than MINCOUNT counts will not be saved in the output file

## NOTE: this may require python 2.7.5...
# Argument parsing:
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, metavar='<alignments list file>',
                        help='input alignments list file')
    parser.add_argument('-o', type=str, required=True, metavar='<outFile>',
                        help='output counts file')
    parser.add_argument('-m', type=int, required=False, metavar='<min BC count>',default=MINCOUNT,
                        help='minimum barcode read count')
    args = parser.parse_args()
    return args

def run(args):
    inFP = args.i       # input alignment file list file
    outFile = args.o    # output file name
    minCount = args.m   # minimum read count for a valid barcode

    ## open and read the list of input files:
    fIn = open(inFP, 'r')
    inFiles = []
    while 1:
        line = fIn.readline()
        if not line:
            break
        # skip commented lines:
        if line.startswith('#'):
            continue
        if line.endswith('\n'):
            line=line[:-1]
        
        fields = line.split('\t')
        inFiles.append(fields[0])
    fIn.close()

    print '%d input files...' % len(inFiles)

    bcDict = {}
    nReadsIn = 0
    nReadsOut = 0
    bcCountOut = 0
    readSum = 0
    nUnique = 0
    nhMissing = 0

    for f in inFiles:
        print "processing %s..." % f
        # open the input file:
        inBam = pysam.AlignmentFile(f, 'rb')

        for rd in inBam.fetch():
            # get the read name:
            rName = rd.query_name
            fields = rName.split(':')   # split the query name on ':" to separate read name, BC and UMI
            readSum+=1
            if len(fields)<3:
                print 'Error in read name format: %s' % rName
                continue

            # check for unique alignments:
            try:
                nAligns = rd.get_tag('NH')     # NH=number of lignments for this read
                if nAligns==1:
                    isUnique = True
                    nUnique+=1
                else:
                    isUnique = False
            except:
                # if the NH tag isn't present assume it's unique
                isUnique = True
                nhMissing+=1
            
            if isUnique:
                # get the barcode:
                bc = fields[-2]    # new barcode format (<read name>:<BC>:<UMI>)
                bcDict.setdefault(bc,0)
                bcDict[bc]+=1

        # close the input files:
        inBam.close()

    print 'Total reads in: %d\t unique barcodes: %d' % (readSum, len(bcDict.keys()))
    print 'Unique alignments: %d\tNH tag missing: %d\n' % (nUnique, nhMissing)

    # Write the new barcode counts file:
    fOut = open(outFile, 'w')
    fOut.write('BC\tcount\n')
    for bc in bcDict.keys():
        count = bcDict[bc]
        if count >= minCount:
            nReadsOut+=count
            bcCountOut+=1
            fOut.write('%s\t%d\n' % (bc,count))
            
    fOut.close()
    
    # Print the statistics:
    print 'Total reads out: %d\tbarcodes out (>%d reads): %d\n' % (nReadsOut, minCount, bcCountOut)
    
    return

    
if __name__ == "__main__":
    args = get_args()
    run(args)


