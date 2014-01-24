##Wynn Meyer
##14 July 2013
##Modified from:
##John Blischak
##28 Oct 2011

fasta_file1='EfOCA2RegionsSurroundingFDSites.fasta'  ## Do not include path to the file. Just make sure file is in same folder as python script.
fasta_file2='EmOCA2RegionsSurroundingFDSites.fasta'

## Goals for steps to add:
## Generate fasta files (or strings) with 1 kb to each side of "fixed difference" sites for both species.
## Only find primers that do not contain sites that are polymorphic in either species or "fixed differences."
## Might be easier to do in bedtools: PyBedtools
## Find the set of amplicons that covers the most "fixed difference" sites overall.

##  Earlier steps.
##  Make a file with fasta sequences that contain 2 kb upstream and 100 bp downstream of a TSS
##  New: Can these be smaller regions? Does one of the other scripts do this?
##  Maybe: Download BLAST+ from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
##  Follow directions: http://www.ncbi.nlm.nih.gov/books/NBK52637/
##  Download mEMBOSS tools (http://www.interactive-biosoftware.com/embosswin/embosswin.html) from ftp://emboss.open-bio.org/pub/EMBOSS/windows/
##  Make sure mEMBOSS is in path environment variable (it should do this automatically)

import os
from Bio import SeqIO
from Bio import Entrez
#Entrez.email = 'jdblischak@uchicago.edu'
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from pprint import pprint  #able to print dictionary output nicely
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
####Needed to run Primer3 
import sys
import subprocess
import re
from Bio.Emboss.Applications import Primer3Commandline ##Constructs a command line to submit to eprimer3 (EMBOSS tool) which then calls primer3_core
from Bio.Emboss import Primer3 ##Has the functions to parse Primer3 output
####

###################
##Definitions######
###################


def parse_fasta_descr(id):
    '''
    Input: The description that follows the > in a fasta sequence (as a string).
    Output: The desired name for the sequence to be used as a key in a dictionary or other purpose.
    Current function: The description is SYM_prom, and the function returns SYM.
    '''
    return id.split('_')[0]

def read_in_fasta_multiple(filename, species):
    '''
    Input: The name of the file that contains the fasta sequences.
    Output: A dictionary containing the Seq records for each fasta sequence
    Note: The keys of the dictionary are determined by the function parse_fasta_descr.
    Currently works only for DNA sequences.
    '''
    dic={}
    for seq_rec in SeqIO.parse(filename,'fasta',IUPAC.IUPACAmbiguousDNA()):
        key=parse_fasta_descr(seq_rec.id)
        dic[key]=seq_rec
        filename='fasta_lemur_seqs/'+key+'_'+species+'.fa'
        if not os.path.isfile(filename):
            save_file = open(filename, "w")
            save_file.write('>'+key+'\n')
            save_file.write(str(seq_rec.seq))
            save_file.close()           
    return dic

def ensure_dir(f):
    '''
    Input: Name of a folder (directory) as a string
    Output: None.
    Function: If the folder exists, nothing happens. If the folder does not exist, the folder is created.
    Inspired by: http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    '''
    if not os.path.exists(f):
        os.makedirs(f)
        
def make_consensus(seq1,seq2,filename):
    ''' '''
    #Modified to allow seq1 and seq2 to be Seq objects.
    align1 = MultipleSeqAlignment([
#             SeqRecord(Seq(seq1, generic_dna), id="seq1"),
#             SeqRecord(Seq(seq2, generic_dna), id="seq2")
             SeqRecord(seq1, id="seq1"),
             SeqRecord(seq2, id="seq2")
         ])
    consensus=AlignInfo.SummaryInfo.dumb_consensus(AlignInfo.SummaryInfo(align1),ambiguous='N')
    #Over-write!
#    if not os.path.isfile(filename):
    save_file = open(filename, "w")
    save_file.write(str(consensus))
    save_file.close()
    return consensus   

def filter_out(keys,dic1,dic2):
    ''' '''
    flagged=[]
    for key in keys:
        if len(str(dic1[key])) < 2000:
            flagged.append(key)
        elif str(dic2[key].seq)[1995:2006] not in str(dic1[key]):
            flagged.append(key)
##        elif str(dic[key]) > 2201:
##            flagged.append(key)
##        elif str(dic[key]).count('N') > 50:
##            flagged.append(key)
##        elif str(dic[key]).count('N') / len(str(dic[key])) > .2:
##            flagged.append(key)

    for x in flagged:
        keys.remove(x)

    return keys
                
def call_primer3(keys,dic):
    ''' '''
    cline=Primer3Commandline()
    cline.auto=True ##Turns off the prompts
    cline.maxdifftm='3'
    cline.otm='60'
    cline.numreturn='2'
    cline.explainflag=True
    cline.prange='900-1400'
    cline.psizeopt='1250'
##    cline.includedregion='1,1000,2000,2100'
    cline.maxtm='70'
    cline.minsize='17'
        
    for key in keys:
        cline.target=str(len(d_consensus[key])-100)+'-'+str(len(d_consensus[key])-98)
        cline.sequence='consensus_seqs/'+key+'_cons.txt'
        cline.outfile='primer3_results/'+key+'_primer3.txt'
        # print cline
        cline()       

def parse_primer3(keys):
    ''' '''
    flagged=[]
    dic={}
    regexp = re.compile('[RrYyKkMmSsWwNn]')
    for key in keys:
        filename='primer3_results/'+key+'_primer3.txt'
        open_handle=open(filename,'r')
        dic[key]=Primer3.read(open_handle)
        open_handle.close()
        if dic[key].primers==[]:
            flagged.append(key)
            del dic[key]
        #Added: Also remove if there are polymorphic/divergent sites w/in primer:
        else:
            try: 
                if regexp.search(str(dic[key].primers[0].forward_seq)) is not None:
                    flagged.append(key)
                    del dic[key]
                elif regexp.search(str(dic[key].primers[0].reverse_seq)) is not None:
                    flagged.append(key)
                    del dic[key]
            except:
                print dic[key].primers[0]
                exit(1)
    for x in flagged:
        keys.remove(x)

    return keys, dic
        
###################
##Script###########
###################

## Generate fasta files with 1 kb to each side of "fixed difference" sites for both species.
fdfile = "/scratch/wynn/FilesForCluster/scaffold2503_100kb_div.txt"
seqfile = "/mnt/lustre/home/wynn/LemurFiles/Pigmentation/BothLemursFastaFromFastq071813.fasta" #contains Ef and Em sequences, plus another "scaffold_2503"
with open(fdfile) as f:
    diffsites = f.read().splitlines()
effa = []
emfa = []
seqhandle = open(seqfile, "r")
for record in SeqIO.parse(seqhandle, "fasta"):
    if record.id == "Ef_scaffold2503":
         for diff in diffsites:
             efsubrec = record.seq[max(0,int(diff)-1500):min(len(record),int(diff)+1500)]
             effa.append(SeqRecord(efsubrec, id=str(diff)))
    elif record.id == "Em_scaffold2503":
         for diff in diffsites:
             emsubrec = record.seq[max(0,int(diff)-1500):min(len(record),int(diff)+1500)]
       	     emfa.append(SeqRecord(emsubrec, id=str(diff)))
efout = open(fasta_file1, "w")
emout = open(fasta_file2, "w")
SeqIO.write(effa, efout, "fasta")
SeqIO.write(emfa, emout, "fasta")
seqhandle.close()
efout.close()
emout.close()

## Make a folder to hold invdividual fasta files
ensure_dir('fasta_lemur_seqs')
## Read in the fasta sequences, store in dictionary, and also write to file.
## key = gene name ; value = Seq record object
ef_seqs=read_in_fasta_multiple(fasta_file1, "Ef")
em_seqs=read_in_fasta_multiple(fasta_file2, "Em")
## Create a list that will keep track of the genes that can still have primers designed
## This should be a list of fixed difference sites or regions.
##viable_proms=['PRDM1']
viable_sites=ef_seqs.keys()

##Create a folder to store blast results if there is not already one
##Not necessary for this project.
#ensure_dir('blast_results')
## Blast the human promoter to find orthologous chimp sequence
#blast_recs=blast_chimp_local(viable_proms,d_seqs)
##
##Create a folder to store blast results if there is not already one
ensure_dir('consensus_seqs')
## Make the consensus sequences and store in dictionary
## Do this from the two fastas, rather than a MSA.
d_consensus={}
for key in viable_sites:
    filename='consensus_seqs/'+key+'_cons.txt'
#    d_consensus[key]=make_consensus(blast_recs[key].alignments[0].hsps[0].query,blast_recs[key].alignments[0].hsps[0].sbjct,filename)
    d_consensus[key]=make_consensus(ef_seqs[key].seq,em_seqs[key].seq,filename)

## Remove bad promoters
## Can the "filter_out" sub-process be used to remove other undesirable features?
#viable_proms=filter_out(viable_proms,d_consensus,d_seqs)

##Create a folder to store blast results if there is not already one
ensure_dir('primer3_results')
## Call Primer3 with consensus sequences
call_primer3(viable_sites,d_consensus)

## Parse through primer3 output
d_p3recs={} 
viable_sites,d_p3recs=parse_primer3(viable_sites)

save_handle=open('primers.txt','w')
save_handle.write('Name\tSequence\tGC\tTm\tLength\tProduct\tStart\n')
for key in viable_sites:
    save_handle.write(key+'_F'+'\t'+str(d_p3recs[key].primers[0].forward_seq)+'\t'+str(d_p3recs[key].primers[0].forward_gc)+'\t'+
                      str(d_p3recs[key].primers[0].forward_tm)+'\t'+str(d_p3recs[key].primers[0].forward_length)+'\t'+str(d_p3recs[key].primers[0].size)+
                      '\t'+str(d_p3recs[key].primers[0].forward_start)+'\n'+
                      key+'_R'+'\t'+str(d_p3recs[key].primers[0].reverse_seq)+'\t'+str(d_p3recs[key].primers[0].reverse_gc)+'\t'+
                      str(d_p3recs[key].primers[0].reverse_tm)+'\t'+str(d_p3recs[key].primers[0].reverse_length)+'\t'+str(d_p3recs[key].primers[0].size)+
                      '\t'+str(d_p3recs[key].primers[0].reverse_start)+'\n')
save_handle.close()

## Get the predicted amplicons for Harlow and Harmonia for each primer pair.
ensure_dir('predicted_amplicons')
#From R code:
for key in viable_sites:
    inf1 = open('fasta_lemur_seqs/'+key+'_Ef.fa', "fasta")
    harl=SeqIO.read(inf1, "fasta")
    subharl = harl[int(d_p3recs[key].primers[0].forward_start)+35:int(d_p3recs[key].primers[0].reverse_start)-15]
    outf1 = open('predicted_amplicons/'+key+'_Harlow.fasta', 'w')
    SeqIO.write([subharl],outf1, "fasta")
    inf1.close()
    outf1.close()
    inf2 = open('fasta_lemur_seqs/'+key+'_Em.fa', "fasta")
    harm=SeqIO.read(inf2, "fasta")
    subharm = harm[int(d_p3recs[key].primers[0].forward_start)+35:int(d_p3recs[key].primers[0].reverse_start)-15]
    outf2 = open('predicted_amplicons/'+key+'_Harmonia.fasta', 'w')
    SeqIO.write([subharl],outf2, "fasta")
    inf2.close()
    outf2.close()
    
#What is the following for?
#key='PRDM1'

#Find the combination of primers (one set within each 30 kb interval) that will amplify the most "fixed difference" sites.
#Do this in R? See what Primer3 output looks like.

##Primer3 code that I didn't use for my script

##for key in viable_genes:
##    file_name=seq_recs_1[key].id
##    file_extension=r'C:\Users\John\Desktop\LAB\TF_tagging_project\primer_design\tf_primer3'
##    open_sequence=open(file_extension+'\\'+file_name+'_sequence.txt','w')
##    open_sequence.write(str(primer_F_GW[key].upper()[:34])+str(cds_human[key][:-1])+str(primer_R_GW[key].reverse_complement().upper()[21:]))
##    open_sequence.close()
##    ##Call primer3
##    cline = Primer3Commandline()
##    cline.sequence=file_extension+'\\'+file_name+'_sequence.txt' ##The file with sequence
##    cline.outfile=file_extension+'\\'+file_name+'_primer3.txt'   ##Where to store results
##    cline.auto=True ##Turns off prompts. Does not seem to be necessary, but other programmers changed it to True.
##    cline.explainflag=True
##    cline.forwardinput=str(primer_F_GW[key][-35:]) ##forward primer
##    cline.outfile=file_extension+'\\'+file_name+'_primer3.txt'
##    ##Make parameters not stringent
##    cline.maxgc=100.0
##    cline.maxpolyx=10
##    cline.maxsize=35 ##Primer3 can't calculate Tm for primer over this length
##    cline.maxtm=150
##    cline.mingc=0
##    cline.mintm=0
##    cline.numreturn=1
####    cline.prange=90000
##    
##    child = subprocess.Popen(str(cline), \
##                        stdout=subprocess.PIPE, \
##                        stderr=subprocess.PIPE, \
##                        shell=(sys.platform!="win32"))
##
##    MeltingTemp.Tm_staluc(str(primer_R_GW[key]))
    
##open_handle=open(r'C:\Users\John\Desktop\LAB\TF_tagging_project\primer_design\myresults.out','r')
##results=Primer3.read(open_handle)
##open_handle.close()

##child.stderr.read() ##Read error message
