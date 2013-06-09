genome_analysis
===============

Custom scripts developed for the processing of genome-wide data for Eulemur macaco and Eulemur flavifrons.

File descriptions:

GenerateEmConsensusFromVCF.py: used to generate a consensus sequence with no polymorphic sites and only 
                               low-coverage sites masked from an all-sites VCF file generated using the 
                               alignment (.bam) of reads from one species to the reference sequence of 
                               another species, produced using GATK
