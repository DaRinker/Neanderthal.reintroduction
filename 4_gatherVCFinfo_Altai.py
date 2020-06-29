#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 09:48:49 2018

@author: rinkerd
"""

#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon May 14 13:48:07 2018

@author: drinker
this script follows the running of "LDsnps_annotate.py" and requires that there be a pickle binary ojject named "addSNPs2haploLDpopulations"

the required pickel object for this script is a dictionary keyed on position for all additional SNPs that Vernot reported as bieng in LD with Neand tag SNPs

this script will query the Altai genoem (Pruffer 2013) vcf files and extract information for every position in the addSNPs2haploLDpopulations dictionary


"""
usage = """
USAGE:
-h     This message.
"""

import os
import sys
import pickle
import os.path
import gzip


#####_FUNCTIONS__#########

def save_obj(obj, savefilename ):
        with open ('obj/' + savefilename + '.pkl', 'w+b') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL) #pickle.HIGHEST_PROTOCOL is a binary format; protocol 0 zero is a text format

####################

addSNPs2Altai = {}
addSNPs2haploLDpopulations = {}
addSNPs2haploLDpopulations = pickle.load(open('obj/addSNPs2haploLDpopulations.pkl', 'rb')) 


CHROMOSOMES = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

for CHROM in CHROMOSOMES:

    with gzip.open('<path to Altai Neanderthal vcf files>/' + CHROM + '.mod.vcf.gz', 'rt') as f:
        for line in f:
            
            if line.startswith('#'):
            		continue
            if "AC=" not in line:
                continue
            if "LowQual"  in line:
                continue
            
            seg = line.split('\t')
            loc = int(seg[1])    
                
            if (CHROM, loc) not in addSNPs2haploLDpopulations:
                continue
            
            ref = str(seg[3])
            alt = str(seg[4])
            altAC = seg[7].split('AC=')[1].split(';')[0]	# note ';AF=' also captures the same lines...
            
            addSNPs2Altai[(CHROM, loc)]=[(ref, alt, altAC)]

save_obj(addSNPs2Altai, 'addSNPs2Altai')
                
