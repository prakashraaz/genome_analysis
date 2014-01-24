#Use the calculated cutoff for the tail of the PWM score distribution
#to get a list of sites within the tail and for which the blue-eyed allele
#has a score < {folddiff} times the brown-eyed allele score.
#This is the same as CalculatePWMScores.py, but it only produces output for SNPs.

import gzip, sys, math, numpy
args = sys.argv
tfname = str(args[1])
quant = float(args[2])
folddiff = float(args[3])
indir="/mnt/lustre/home/wynn/LemurFiles/FilesForBryce"
outdir="/scratch/wynn/FilesForCluster/TFBSOut"
ref_file=open(indir+"/EfReference.fasta")
vcfbrown=gzip.open(indir+"/EmVCFForScaff2503Only.vcf.gz")
vcfblue=gzip.open(indir+"/NewEfVCFForScaff2503Only.vcf.gz")
result_file=open(outdir+"/"+tfname+"Q"+str(quant)+"FD"+str(folddiff)+".results.txt","w")
#newpwm_file=open(outdir+"/"+tfname+".newpwmscores.txt","w")

### function to get PWM score given sequence and matrix
def get_PWM_score(seq,pwm):
    lscore=0
    try:
        for i in range(len(seq)):
            if seq[i]=="A":
                lscore+=pwm[i][0]
            elif seq[i]=="C":
                lscore+=pwm[i][1]
            elif seq[i]=="G":
                lscore+=pwm[i][2]
            elif seq[i]=="T":
                lscore+=pwm[i][3]
    except:
        sys.stderr.write(seq)
        exit()
    return lscore

#First, determine the cutoff for PWM score.
cutoffs=open(outdir+"/PWMScoreCutoffsFor"+str(quant)+"TailOfRefBr.txt","r")
for line in cutoffs:
    tfcut=line.strip().split()
    if tfcut[0]==tfname:
        relcutoff=float(tfcut[1])
        break
cutoffs.close()
print "Cutoff for "+tfname+" is "+str(relcutoff)
print "Fold difference is "+str(folddiff)

### specify the PWM you want
###     A   C   G   T ####

# E-box
eboxpwm=[[16 ,59,     21   ,   23],
[0    , 118   ,  0     ,  1],
[106  ,   7   ,    1    ,   5],
[4    ,   94  ,    15   ,   6],
[23   ,   44  ,    48   ,   4],
[0    ,   0   ,    1    ,   118],
[1    ,   4   ,    112  ,   2],
[47   ,   21  ,    21   ,   30],
[5    ,   71  ,    14   ,   29],
[16  ,    47  ,    26  ,    30]]


# LEF1 M01022
lef1pwm = [[2 ,      8 ,      8 ,      4],
[12  ,    3   ,    1   ,    7],
[7   ,    1   ,    1   ,    14],
[1   ,    22  ,    0   ,    0],
[23  ,    0   ,    0   ,    0],
[22  ,    1   ,    0   ,    0],
[21  ,    0   ,    1   ,    1],
[0   ,    2   ,    21  ,    0],
[4   ,    3   ,    9   ,    7],
[5   ,    6   ,    10  ,    2]]

#RUSH
rushpwm=[[16 ,     10 ,     12 ,     13],
[20     , 12  ,    5   ,    16],
[26     , 33  ,    0   ,    0],
[0      , 59  ,    0   ,    0],
[25     , 0   ,    0   ,    34],
[0      , 0   ,    0   ,    59],
[24     , 5   ,    12  ,    18],
[0      , 0   ,    22  ,    37],
[15     , 14  ,    14  ,    8],
[12     , 12  ,    10  ,    15]]

if tfname == "E-box":
    pwm = eboxpwm
elif tfname == "LEF1":
    pwm = lef1pwm
elif tfname == "RUSH":
    pwm = rushpwm
else:
    print "Provided TF name ("+tfname+") does not match a known TF \
          (E-box, HLTF, LEF1, or RUSH). In order to use a different TF, \
          you must provide a file with a PWM, as well as a cutoff \
          for inclusion in the top hits file."
    exit

for i in range(len(pwm)):
    pwm[i]=numpy.log((numpy.array(pwm[i])+1)/sum(numpy.array(pwm[i])+1.0))-numpy.log(.25)

rpwm=[]
for i in range(len(pwm)):
    rpwm.insert(0,[pwm[i][3],pwm[i][2],pwm[i][1],pwm[i][0]])

refblue=""
altblue=""
refbrown=""
altbrown=""

ref_line=ref_file.readline()
ref_line=ref_file.readline()

vcfblue_line=vcfblue.readline()
vcfblue_line=vcfblue.readline()

vcfbrown_line=vcfbrown.readline()
vcfbrown_line=vcfbrown.readline()
i=1
j=0


#### This chunk loads the first PWM-lengthed sequence
while len(refblue)<len(pwm):
    vcfblue_split=vcfblue_line.strip().split()
    vcfbrown_split=vcfbrown_line.strip().split()

    if j>=len(ref_line):
        ref_line=ref_file.readline()
        j=0
    if int(vcfblue_split[1])==i:
        refblue+=vcfblue_split[3]
        refbrown+=vcfbrown_split[3]

        if vcfblue_split[4] != ".":
            altblue+=vcfblue_split[4]
            altbrown+=vcfbrown_split[4]
        else:
            altblue+=vcfblue_split[3]
            altbrown+=vcfbrown_split[3]

        vcfblue_line=vcfblue.readline()
        vcfbrown_line=vcfbrown.readline()
    else:
        refblue+=ref_line[j]
        refbrown+=ref_line[j]
        altblue+=ref_line[j]
        altbrown+=ref_line[j]

    i+=1
    j+=1

while vcfblue_line:

    #get all the pwm scores
    refbl=max(get_PWM_score(refblue,pwm),get_PWM_score(refblue,rpwm))
    refbr=max(get_PWM_score(refbrown,pwm),get_PWM_score(refbrown,rpwm))
    altbl=max(get_PWM_score(altblue,pwm),get_PWM_score(altblue,rpwm))
    altbr=max(get_PWM_score(altbrown,pwm),get_PWM_score(altbrown,rpwm))    


    #here you can specify under which conditions to write
    if refbr>=refbl+math.log(folddiff) and refbr>=relcutoff and refbr==altbr and refbl==altbl: 
        snps=[]
        for k in range(len(refbrown)):
            if refbrown[k]!=refblue[k]:
                snps.append(str(k+i-len(pwm)))
        for s in range(len(snps)):
            result_file.write("%d\t%f\t%f\t%f\t%f\t%s\n"%(i-len(pwm),refbl,altbl,refbr,altbr,snps[s]))
            #This writes each SNP separately; previously these were written to the same line:
            #result_file.write("%d\t%f\t%f\t%f\t%f\t%s\n"%(i-len(pwm),refbl,altbl,refbr,altbr,",".join(snps)))
    #newpwm_file.write("%d\t%f\t%f\t%f\t%f\n"%(i,refbl,altbl,refbr,altbr))
    
    
    ### Slide the window over one
    vcfblue_split=vcfblue_line.strip().split()
    vcfbrown_split=vcfbrown_line.strip().split()

    if j>=len(ref_line):
        ref_line=ref_file.readline()
        j=0
    if int(vcfblue_split[1])==i:
        refblue=refblue[1:]+vcfblue_split[3]
        refbrown=refbrown[1:]+vcfbrown_split[3]

        if vcfblue_split[4] != ".":
            if len(vcfblue_split[4])==3:  #if there at 2 "alt" alleles
                refblue=refblue[:-1]+vcfblue_split[4].split(",")[0]
                altblue=altblue[1:]+vcfblue_split[4].split(",")[1]
            elif vcfblue_split[9].split(":")[0]=="1/1":
                refblue=refblue[:-1]+vcfblue_split[4]
                altblue=altblue[1:]+vcfblue_split[4]
            else:
                altblue=altblue[1:]+vcfblue_split[4]  
        else:
            altblue=altblue[1:]+vcfblue_split[3]

        if vcfbrown_split[4] != ".":
            if len(vcfbrown_split[4])==3:  #if there at 2 "alt" alleles
                refbrown=refbrown[:-1]+vcfbrown_split[4].split(",")[0]
                altbrown=altbrown[1:]+vcfbrown_split[4].split(",")[1]
            elif vcfbrown_split[9].split(":")[0]=="1/1":
                refbrown=refbrown[:-1]+vcfbrown_split[4]
                altbrown=altbrown[1:]+vcfbrown_split[4]
            else:
                altbrown=altbrown[1:]+vcfbrown_split[4]  
        else:
            altbrown=altbrown[1:]+vcfbrown_split[3]
            

        vcfblue_line=vcfblue.readline()
        vcfbrown_line=vcfbrown.readline()
    else:
        refblue=refblue[1:]+ref_line[j]
        refbrown=refbrown[1:]+ref_line[j]
        altblue=altblue[1:]+ref_line[j]
        altbrown=altbrown[1:]+ref_line[j]

    i+=1
    j+=1

result_file.close()
ref_file.close()
vcfbrown.close()
vcfblue.close()
