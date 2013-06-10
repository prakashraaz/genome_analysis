#!/usr/local/bin/python

# This script aims to calculate nucleotide diversity within species and
# divergence between species by counting the number of polymorphic or divergent sites per base pair.
# It will report these numbers for EACH scaffold within the file.
# Pseudo-code:
# Generate counters for:
# 1) total high quality bases in each genome, and in both.
# 2) number of high quality bases where E. flavifrons genotype is heterozygous (0/1)
# 3) number of high quality bases where E. macaco genotype is heterozygous (0/1)
# 4) within the above, subset at which the other genotype is
#      a) callable
#      b) heterozygous
#      c) heterozygous with the same alternate allele
# 5) number of high quality bases where E. flavifrons and E. macaco genotype are different homozygotes (0/0 and 1/1).
# 6) total bases
# Ideally, this program should take as input the desired quality score cutoff and depth requirements.
# Also ideally, this program should output what it calculates.

import sys
from Bio import SeqIO
import os
import re
import gzip
import numpy
from time import gmtime, strftime
import vcf
import random
import csv

hiqef = 0
hiqem = 0
hiqboth = 0
snpefmonoem = 0
snpemmonoef = 0
snpboth = 0
samesnpboth = 0
snpefnem = 0
snpemnef = 0
div = 0
tot = 0
unc = 0

qcut = 30
depthlo1 = 17
depthhi1 = 104
depthlo2 = 10
depthhi2 = 42
mapcut = 30
nonrefsamplename = 'EmHarmonia'

vcfone = sys.argv[1]
vcftwo = sys.argv[2]
outbase = sys.argv[3]

print 'Getting diversity and divergence estimates for files '+vcfone+' and '+vcftwo+' at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())

efvcf_reader = vcf.Reader(open(vcfone, 'r'))
emvcf_reader = vcf.Reader(open(vcftwo, 'r'))
csvfile = open(outbase+".txt",'w') 
outwriter = csv.writer(csvfile, delimiter='\t')
scaffoldlist = []

#missingscaf = []
#missingbp = []
for record1 in efvcf_reader:
    tot = tot + 1
    if tot % 5000000 is 0:
        print str(tot)+' records of '+vcfone+' parsed at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())
        sys.stdout.flush()
    if record1.CHROM not in scaffoldlist:
        if len(scaffoldlist) > 0:
            outwriter.writerow([str(scaffoldlist[len(scaffoldlist)-1]),str(tot),str(hiqef+hiqboth),str(hiqem+hiqboth),str(hiqboth),str(snpefmonoem),str(snpemmonoef),str(snpboth),str(samesnpboth),str(snpefmonoem+snpboth+snpefnem),str(snpemmonoef+snpboth+snpemnef),str(div),str(unc)])
            hiqef = 0
            hiqem = 0
            hiqboth = 0
            snpefmonoem = 0
            snpemmonoef = 0
            snpboth = 0
            samesnpboth = 0
            snpefnem = 0
            snpemnef = 0
            div = 0
            tot = 0
            unc = 0
        scaffoldlist.append(record1.CHROM)
    try:
        record2 = emvcf_reader.next() #corresponding row (hopefully) for black lemur
    except:
#This means there are no more E.m. records, but there are more E.f. records.
#In this case, just do E.f. stuff and skip the rest of the loop. Try using 'continue.'
#        missingscaf.append(record1.CHROM)
#        missingbp.append(record1.POS)
        efhiq = record1.QUAL >= qcut and record1.INFO['DP'] >= depthlo1 and record1.INFO['DP'] <= depthhi1 and record1.INFO['MQ'] >= mapcut
        efsnp = record1.ALT[0] is not None
        if efhiq:
           hiqef = hiqef + 1
           if efsnp:
                snpefnem = snpefnem + 1
        continue
    while (record2.POS > record1.POS or record2.CHROM not in scaffoldlist):
#There is no data for this position in the E. macaco vcf file.
#Try assuming that it means there's no data for the E. macaco sequence at this position.
#        missingscaf.append(record1.CHROM)
#        missingbp.append(record1.POS)
        efhiq = record1.QUAL >= qcut and record1.INFO['DP'] >= depthlo1 and record1.INFO['DP'] <= depthhi1 and record1.INFO['MQ'] >= mapcut
        efsnp = record1.ALT[0] is not None
        if efhiq:
           hiqef = hiqef + 1
           if efsnp:
                snpefnem = snpefnem + 1
        record1 = efvcf_reader.next()
        tot = tot + 1
        if record1.CHROM not in scaffoldlist:
            scaffoldlist.append(record1.CHROM)
#Try doing a bunch of logicals and then updating counters based on their intersection.
    efhiq = record1.QUAL >= qcut and record1.INFO['DP'] >= depthlo1 and record1.INFO['DP'] <= depthhi1 and record1.INFO['MQ'] >= mapcut
    emhiq = record2.QUAL >= qcut and record2.INFO['DP'] >= depthlo2 and record2.INFO['DP'] <= depthhi2 and record2.INFO['MQ'] >= mapcut
    efsnp = record1.ALT[0] is not None
    emsnp = record2.genotype(nonrefsamplename)['GT'] == '0/1'
    newem = record2.genotype(nonrefsamplename)['GT'] == '1/1'
    newemsnp = record2.genotype(nonrefsamplename)['GT'] == '1/2'
#That should be enough to categorize everything.
    if efhiq and not emhiq:
        hiqef = hiqef + 1
        if efsnp:
            snpefnem = snpefnem + 1
    elif emhiq and not efhiq:
        hiqem = hiqem + 1
        if (emsnp or newemsnp):
            snpemnef = snpemnef + 1
    elif efhiq and emhiq:
        hiqboth = hiqboth + 1
        if not efsnp and not (emsnp or newemsnp or newem):
            continue
        elif efsnp and not (emsnp or newemsnp):
            snpefmonoem = snpefmonoem + 1
            # Pick one of the E.f. alleles. If it's not the one E. m. has, call it divergent.
            # This is tricky because of weird formatting (the REF is a '_Substitution object). 
            # If E. m. is homozygous for either the REF or the ALT allele, flip a coin.
            # Otherwise, it's definitely divergent.
            if record2.ALT[0] is None or (newem and record1.ALT == record2.ALT):
                coinflip = random.sample([0,1],1)
                if coinflip == [0]:
                    div = div + 1
            else:
                div = div + 1
        elif emsnp and not efsnp:
            snpemmonoef = snpemmonoef + 1
            # This is just a coin flip; there's no way for the alternate allele to be present in E. f.
            coinflip = random.sample([0,1],1)
            if coinflip == [0]:
                div = div + 1
        elif newemsnp and not efsnp:
        #Neither allele in E.m. is present in E.f.
            snpemmonoef = snpemmonoef + 1
            div = div + 1
        elif efsnp and emsnp:
            snpboth = snpboth + 1
            if record1.ALT[0] == record2.ALT[0]:
                samesnpboth = samesnpboth + 1
                coinflip = random.sample([0,1],1)
                if coinflip == [0]:
                    div = div + 1
            else:
            # The reference is the same but the alternate allele is different.
                oneinfour = random.sample([0,0,0,1],1)
                if oneinfour == [0]:
                    div = div + 1
        elif efsnp and newemsnp:
            snpboth = snpboth + 1
            # Em has a 2nd possible alternate allele.
            if record1.ALT[0] in record2.ALT[0:1]:
            # The genotypes have one allele in common.
                oneinfour = random.sample([0,0,0,1],1)
                if oneinfour == [0]:
                    div = div + 1
            else:
                div = div + 1
        elif newem and not efsnp:
            div = div + 1
        else:
            unc = unc + 1
            print 'The record for position '+str(tot)+' is high quality but does not fall into a standard category.'
            print record1
            print record2

#Add the final values
outwriter.writerow([str(scaffoldlist[len(scaffoldlist)-1]),str(tot),str(hiqef+hiqboth),str(hiqem+hiqboth),str(hiqboth),str(snpefmonoem),str(snpemmonoef),str(snpboth),str(samesnpboth),str(snpefmonoem+snpboth+snpefnem),str(snpemmonoef+snpboth+snpemnef),str(div),str(unc)])
csvfile.close()

print 'Finished getting diversity and divergence estimates for files '+vcfone+' and '+vcftwo+' at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())
