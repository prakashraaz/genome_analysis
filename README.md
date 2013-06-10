genome_analysis
===============

Custom scripts developed for the processing of genome-wide data for Eulemur macaco and Eulemur flavifrons.

File descriptions:

GenerateEmConsensusFromVCF.py: used to generate a consensus sequence with no polymorphic sites and only 
                               low-coverage sites masked from an all-sites VCF file generated using the 
                               alignment (.bam) of reads from one species to the reference sequence of 
                               another species, produced using GATK

vcfutils_mod.pl: contains modified vcf2fq and vcf2fqnonref functions used to generate SNP-inclusive
                 consensuses for E. flavifrons and E. macaco, respectively, from all-sites VCF files
                 generated using GATK

CalculateDiversityAndDivergence.py: used to calculate per-scaffold diversity within and divergence
                                    between species, using high quality base pairs from the all-sites
                                    VCF files from GATK
