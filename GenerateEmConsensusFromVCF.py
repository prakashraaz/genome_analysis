#!/usr/local/bin/python

# This script should take the all-sites vcf file 
# for a set of scaffolds of the non-reference sample,
# reduce polymorphic sites to just the more abundant base pair,
# and replace base pairs with coverage <3 with N's.
# It should produce a consensus fasta similar to that output by SOAP.

# Pseudo-code:
# Read in VCF file. Iterate over records. For each record, replace:
# 1) low quality (lowercase) bases with coverage >=3 with uppercase
# 2) ambiguous bases with the more prevalent base pair (or one at random if equal)
# 3) any bases with coverage <3 with N's.
# If any records are converted entirely to N's in the process, remove them.

import sys
sys.path.append("/mnt/lustre/home/wynn/bin/Python-2.7.3")
sys.path.append("/mnt/lustre/home/wynn/bin/lib/python2.7/site-packages/PyVCF-0.6.3-py2.7-linux-x86_64.egg")
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re
import gzip
import numpy
from time import gmtime, strftime
import random
import csv

#qcut = 30
depthlo1 = 3
#depthhi1 = 104
#depthlo2 = 10
#depthhi2 = 42
#mapcut = 30

vcfone = sys.argv[1]
samplename = sys.argv[2] #Use 'EmHarmonia'
outbase = sys.argv[3]

try:
    import vcf
except:
    print 'VCF import failed for '+vcfone+' at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())
    sys.exit[1]

print 'Running conversion to consensus fasta for '+vcfone+' at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())

emvcf_reader = vcf.Reader(open(vcfone, 'r'))
scaffoldlist = []
seqs = []
rawseq=''
tot = 0
for record1 in emvcf_reader:
    tot = tot + 1
    if tot % 5000000 is 0:
        print str(tot)+' records of '+vcfone+' parsed at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())
        sys.stdout.flush()
    if record1.CHROM not in scaffoldlist:
    #Append the old record and start a new one.
        scaffoldlist.append(record1.CHROM)
        if len(scaffoldlist)>1:
            ns = re.findall('[Nn]', rawseq)
            if len(ns) == len(rawseq):
            # No high quality base pairs exist in this scaffold.
                print 'There are no high quality bp for '+samplename+' in '+newid
            else:
                rec = SeqRecord(Seq(rawseq), id=newid, description='')
                seqs.append(rec)
            rawseq=''
        newid=record1.CHROM
    if record1.genotype(samplename)['GT']!=None:
        meetsdepth = record1.genotype(samplename)['DP']>= depthlo1
    else:
    #This occurs when there is no data at this position in the sample.
        meetsdepth = False
    refem = record1.genotype(samplename)['GT'] == '0/0'
    emsnp = record1.genotype(samplename)['GT'] == '0/1'
    newem = record1.genotype(samplename)['GT'] == '1/1'
    newsnp = record1.genotype(samplename)['GT'] == '1/2'
    if meetsdepth:
    #If there are at least 3 reads, output a nucleotide.
        if refem:
            rawseq+=record1.REF
        elif newem:
            rawseq+=str(record1.ALT[0])
        elif emsnp:
            if record1.genotype(samplename)['AD'][0] > record1.genotype(samplename)['AD'][1]:
            # Incorporate reference allele.
                rawseq+=record1.REF
            elif record1.genotype(samplename)['AD'][1] > record1.genotype(samplename)['AD'][0]:
            # Incorporate alternate allele.
                rawseq+=str(record1.ALT[0])
            elif record1.genotype(samplename)['AD'][0] == record1.genotype(samplename)['AD'][1]:
            # Flip a coin.
                coinflip = random.sample([0,1],1)
                if coinflip == [0]:
                    rawseq+=record1.REF
                else:
                    rawseq+=str(record1.ALT[0])
            else:
            # Something weird is going on.
                print 'There is an issue with read depth in line '+tot
                print record1.genotype(samplename)
                rawseq+='n'
       	elif newsnp:
       	    if record1.genotype(samplename)['AD'][1] > record1.genotype(samplename)['AD'][2]:
       	    # Incorporate first alternate allele.
       	       	rawseq+=str(record1.ALT[0])
       	    elif record1.genotype(samplename)['AD'][2] > record1.genotype(samplename)['AD'][1]:
       	    # Incorporate second alternate allele.
       	       	rawseq+=str(record1.ALT[1])
       	    elif record1.genotype(samplename)['AD'][1] == record1.genotype(samplename)['AD'][2]:
       	    # Flip a coin.
                coinflip = random.sample([0,1],1)
                if coinflip == [0]:
                    rawseq+=str(record1.ALT[0])
       	       	else:
       	       	    rawseq+=str(record1.ALT[1])
            else:
            # Something weird is going on.
                print 'There is an issue with read depth in line '+tot
                print record1.genotype(samplename)
        else:
        # Something weird is going on.
            print 'There is an issue with genotype in line '+tot
            print record1.genotype(samplename)
            rawseq+='n' 
    else:
    # If there are not at least 3 high quality reads, output 'N'
        rawseq+='N'

#Add the final record:
ns = re.findall('[Nn]', rawseq)
if len(ns) == len(rawseq):
# No high quality base pairs exist in this scaffold.
    print 'There are no high quality bp for '+samplename+' in '+newid
else:
    rec = SeqRecord(Seq(rawseq), id=newid, description='')
    seqs.append(rec)

#Write the records to an output file. These can be concatenated.
with open(outbase+".fasta",'w') as f_out:
    SeqIO.write(seqs, f_out, "fasta")

#f_out.close()
print 'Finished running conversion to consensus fasta for '+vcfone+' at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())

