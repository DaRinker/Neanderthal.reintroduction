#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon May 14 13:48:07 2018

@author: drinker
this script follows the running of "LDsnps_annotate.py" and requires that there be a pickle binary ojject named "addSNPs2haploLDpopulations"

the required pickel object for this script is a dictionary keyed on position for all additional SNPs that Vernot reported as bieng in LD with Neand tag SNPs

this script will 1000 genomes vcf files and extract information for every position in the addSNPs2haploLDpopulations dictionary

"""
usage = """
USAGE:
-h     This message.
"""

import os
import sys
import pickle
import os.path
import pandas
import numpy
import gzip


#####_FUNCTIONS__#########

def save_obj(obj, savefilename ):
        with open ('obj/' + savefilename + '.pkl', 'w+b') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL) #pickle.HIGHEST_PROTOCOL is a binary format; protocol 0 zero is a text format

####################

addSNPs2oneKG = {}
addSNPs2haploLDpopulations = {}
addSNPs2haploLDpopulations = pickle.load(open('obj/addSNPs2haploLDpopulations.pkl', 'rb')) 


CHROMOSOMES = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

#CHROMOSOMES = ['21','22']

for CHROM in CHROMOSOMES:
     with gzip.open('../data/1kgph3_chr' + CHROM + '.gz', 'rt') as f:
    
        for line in f:
            
            if line.startswith('#'):
            		continue
            if "VT=SNP" not in line:
                continue
            
            seg = line.split('\t')
            loc = int(seg[1])
            
            if (CHROM,loc) not in addSNPs2haploLDpopulations:
                continue
                    
            ref = str(seg[3])
            alt = str(seg[4])
            avgAF = seg[7].split(';')[1].split('=')[1]
            easAF = seg[7].split(';')[5].split('=')[1]
            amrAF = seg[7].split(';')[6].split('=')[1]
            afrAF = seg[7].split(';')[7].split('=')[1]
            eurAF = seg[7].split(';')[8].split('=')[1]
            sasAF = seg[7].split(';')[9].split('=')[1]
            AA = seg[7].split(';')[10].split('=')[1].split('|||')[0]
    #        print(CHROM, loc, ref, alt, avgAF, easAF, amrAF, afrAF, eurAF, sasAF, AA)
            addSNPs2oneKG[(CHROM, loc)] = [(ref, alt, avgAF, easAF, amrAF, afrAF, eurAF, sasAF, AA)]

save_obj(addSNPs2oneKG, 'addSNPs2oneKG')
