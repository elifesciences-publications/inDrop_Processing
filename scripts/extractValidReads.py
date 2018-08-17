import sys
import cPickle
import pysam
import scipy.stats
import random
import getopt
import os.path
from Bio.Seq import Seq
import re
import argparse

## NOTE: this may require python 2.7.5...

BC1_LEN = 8    # BC1 length
BC2_LEN = 8    # BC2 length
UMI_LEN = 6    # UMI length
MAX_BCMISMATCH = 1   # number of bases that can be error-corrected
TOPN = 10     # top number of unassigned barcodes to print out

# Argument parsing:
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, metavar='<R1 file>',
                        help='input R1 fastq file')
    # NOTE: output files will be named <outfile>_R1_valid.fastq
    parser.add_argument('-o', type=str, required=True, metavar='<outFile>',
                        help='output file basename')
    parser.add_argument('-d', type=str, required=False,default=None, metavar='<outDir>',
                        help='output file directory')
    parser.add_argument('-b', type=str, required=True, metavar='<validBcFile>',
                        help='valid barcode list file')
    parser.add_argument('-u', type=int, required=False,default=6, metavar='UMI length',
                        help='UMI length')

    args = parser.parse_args()
    return args

def patMatch(x, y, n):
    # default match location:
    loc = None

    # attempt to find pattern 'x' in query sequence 'y', allowing a maximum of 'n' mismatches:
    lx = len(x)
    ly = len(y)
    
    for i in range(ly-lx):   # loop over each possible start location
        nm = 0    # number of mismatches found
        for j in range(lx):   # loop over each base of the target pattern
            if x[j]!=y[i+j]:
                nm+=1        # mismatch
                if nm>n:
                    break    # too many mismatches
        if nm<=n:
            loc = i   # match location
            break

    return loc

### Barcode error-correction:
def correctBC(bc, bcDict):
    bcNew = None    # default (uncorrectable)
    nMatch = 0
    for k in bcDict.keys():
        if hammingDist(bc, k)<=MAX_BCMISMATCH:
            bcNew=k   # corrected barcode
            nMatch+=1 # count the number of possible corrections
            ##break     # don't bother looking any further
    
    if nMatch>1:
        bcNew = None  # barcode can't be unambiguously corrected

    return bcNew   # return corrected barcode

def hammingDist(x, y):
    # find the Hamming distance between two input strings:
    if len(x)!=len(y):
        hd = len(x)
    else:
        hd = 0
        for i in range(len(x)):
            if x[i]!=y[i]:
                hd+=1    # count mismatches
    return hd

def parseBarcodeAndUmiV3(bc1, bc2, umi, bcDict):
    # returned flags:
    bc1Valid = 0
    bc2Valid = 0
    umiValid = 0
    bc1Corr = 0
    bc2Corr = 0
    bailOut = 0

    # test BC1:
    if bcDict.has_key(bc1):
        bc1Valid = 1
    else:
        bc1 = correctBC(bc1, bcDict)  # corrected barcode or None if not correctable
        if bc1!=None: 
            bc1Corr = 1
        else:
            bailOut = 1    # exit without testing BC2 and UMI

    if not bailOut:
        # test BC2:
        bc2 = str(Seq(bc2).reverse_complement())
        if bcDict.has_key(bc2):
            bc2Valid = 1
        else:
            bc2 = correctBC(bc2, bcDict)  # corrected barcode or None if not correctable
            if bc2!=None: 
                bc2Corr = 1
            bailOut = 1    # exit without testing UMI

    # check UMI:
    if not bailOut:
        if umi.count('N')==0:
            umiValid=1
        else:
            umi=None

    return (bc1, bc2, umi, bc1Valid, bc1Corr, bc2Valid, bc2Corr, umiValid)
        
def fastqWrite(f, r, rName):
    # write each field of the read:
    f.write('@%s %s\n' % (rName, r.comment))     # read name and comment
    f.write('%s\n' % r.sequence)                  # the sequence
    f.write('+%s %s\n' % (rName, r.comment))     # read name and comment (filler?)
    f.write('%s\n' % r.quality)                   # the quality string
    return

#def run(basename, outdir):
def run(args):
    R1file = args.i       # R1 (sequence read) file 
    R2file = args.i.replace('R1','R2')      # R2 (BC1) file 
    R3file = args.i.replace('R1','R3')      # R3 (BC2 + UMI) file 
    outdir = args.d   # output directory
    if outdir==None:
        outdir = os.path.dirname(R1file)  # if not provided, put the output file in the same directory as the input
    outbase = args.o  # output file basename
    outFile = os.path.join(outdir,'%s_R1_valid.fastq' % outbase)
    bcFile = args.b   # valid barcode file
    umi_len = args.u  # umi length

    # load the valid barcode dictionary:
    fIn = open(bcFile,'r')
    bcSet = {}
    while 1:
        line = fIn.readline()
        if not line:
            break
        # skip the header line:
        if line.startswith('well'):
            continue
        if line.endswith('\n'):
            line=line[:-1]

        fields = line.split('\t')

        # forward barcode (trimmed to BC1_LEN bases:
        bcFwd = fields[0][:BC1_LEN]
        bcSet.setdefault(bcFwd,0)

    fIn.close()

    ## storage for the output file pointer and statistics counters:
    samp = {}
    samp['name'] = outbase
    oFile = open(outFile,'w')  # open outut file for this sample
    samp['file'] = oFile
    # initialize counters:
    samp['total'] = 0   # total reads
    samp['SBC'] = 0     # sample barcode corrected
    samp['valid'] = 0   # total valid reads
    samp['BC1v'] = 0    # valid BC1
    samp['BC2v'] = 0    # valid BC2
    samp['UMIv'] = 0    # valid UMI
    samp['BC1c'] = 0    # corrected BC1
    samp['BC2c'] = 0    # corrected BC2

    # open the input files:
    fq1 = pysam.Fastqfile(R1file)
    fq2 = pysam.Fastqfile(R2file)
    #fq3 = pysam.Fastqfile(R3file)
    fq3 = pysam.Fastqfile(R3file)

    # counters:
    nBc1Valid = 0
    nBc2Valid = 0
    nBc1Corr = 0
    nBc2Corr = 0
    nUmiValid = 0

    countMod = 100000
    unassigned = 0
    rCount = 0

    unassignedBC = {}   # collect counts on unassigned barcodes

    # loop over all reads:
    while 1:
        try:
            r1 = fq1.next()     # mRNA sequence read
            r2 = fq2.next()     # BC1
            #r3 = fq3.next()     # sample index
            r3 = fq3.next()     # BC2 + UMI
            rCount+=1           # read counter
            if not rCount%countMod:
                print 'read %d' % rCount
        except StopIteration:
            break      # last item
        except:
            print 'pysam.Fastqfile iterator error.'
            eFlag = True
            break

        # parse out the two halves of the cell barcode, and the UMI:
        bc1 = r2.sequence                          # first half of the cell barcode
        bc2 = r3.sequence[:BC2_LEN]                 # second half of the cell barcode
        umi = r3.sequence[BC2_LEN:(BC2_LEN+umi_len)]  # UMI sequence

        # check the barcodes and UMI and update counts:
        (bc1, bc2, umi, bc1Valid, bc1Corr, bc2Valid, bc2Corr, umiValid) = parseBarcodeAndUmiV3(bc1, bc2, umi, bcSet)
        samp['total'] += 1    # total reads
        samp['BC1v'] += bc1Valid    # valid BC1
        samp['BC2v'] += bc2Valid    # valid BC2
        samp['BC1c'] += bc1Corr     # corrected BC1
        samp['BC2c'] += bc2Corr     # corrected BC2
        samp['UMIv'] += umiValid    # valid UMI

        # write out the sequence read if bc1, bc2 and umi are all valid:
        if bc1!=None and bc2!=None and umiValid:
            samp['valid'] += 1   # total valid reads for this sample
            ## create the new read name:
            rName = '%s:%s%s:%s' % (r1.name,bc1,bc2,umi)
            fastqWrite(samp['file'], r1, rName)
        
    # close the input files:
    fq1.close()
    fq2.close()
    #fq3.close()
    fq3.close()

    # print counts:
    print 'Total reads: %d' % rCount
    # close the output file:
    samp['file'].close()

    # print sample-by-sample stats:
    print 'sample\ttotal\tvalid\tBC1valid\tBC1corr\tBC2valid\tBC2corr\tUMIvalid'
    x = samp
    sOut = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (x['name'],x['total'],x['valid'],x['BC1v'],x['BC1c'],x['BC2v'],x['BC2c'],x['UMIv'])
    print sOut

    return
    
if __name__ == "__main__":
    args = get_args()
    run(args)

