#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 13 12:36:42 2018

@author: rinkerd
this script follow the running of "V16_to_1KG_LD_extract.py" and requires that there be a pickle binary of the dictionary 

The output of this script is a new pickle-ed dictionary containing the summarized LD information for all populations; 

"""
usage = """
USAGE:
-h	 This message.
"""

import os
import sys
import pymysql
import pickle
import os.path
import pathlib
from pathlib import Path
import pandas
import numpy

hap2rsq = {}
addSNPs2haploLDpopulations = {}
tagSNP2AF = {}

POP = ['EAS', 'EUR', 'SAS']

#####_FUNCTIONS__#########

def save_obj(obj, savefilename ):
		with open ('obj/' + savefilename + '.pkl', 'w+b') as f:
			pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL) #pickle.HIGHEST_PROTOCOL is a binary format; protocol 0 zero is a text format

####################

tagSNP2AF = pickle.load(open('obj/tagSNP2AF.pkl', 'rb'))

for pop in POP:
	if os.path.isfile('obj/' + pop + 'hap2rsq.dict.pkl'):
		hap2rsq = pickle.load(open('obj/' + pop + 'hap2rsq.dict.pkl', 'rb'))
	else:
		print('File: obj/' + pop + 'hap2rsq.dict.pkl does not exist. Omitting from output...')
		continue
		
##_GENERATE LD SUMMARY dict keyed on (CHR, LOC) of addSNP with values for (POP, maxLD value w/ tagSNP, avg LD w/ tagSNPs, total number tgSNPs in LD, total number of tagSNPs in perfect LD, nead_haplo, tagSNP AF for first tagSNP in maxLD)
#################################
	print("Beginning data for " + pop)	
	addSNPs2haploLD = {}
	for hap in hap2rsq:
		for chrom, loc, r2, tag_loc in hap2rsq[hap]:
			if not r2:
				continue
			else:
				if (chrom, loc) not in addSNPs2haploLD:
					addSNPs2haploLD[(chrom, loc)] = [(r2, hap, tag_loc)]
				else:
					addSNPs2haploLD[(chrom, loc)].append((r2, hap, tag_loc))
	 
	for (chr, loc), v in addSNPs2haploLD.items():
		sum = 0
		maxld = 0
		neandhaplo = (v[0])[1]
		perfld = 0
		
		for i in v:
			sum += i[0]
			if i[0] == 1.0:
				perfld += 1
			if i[0] > maxld:
				maxld = i[0]
				maxldtagSNPaf = (tagSNP2AF[chr, i[2], pop])[0]	#ie, assign the AF of the Neand tagSNP in the current population
		avg = sum / len(v)
		if (chr, loc) not in addSNPs2haploLDpopulations:
			addSNPs2haploLDpopulations[(chr,loc)] = [(pop, maxld, avg , len(v), perfld, neandhaplo, maxldtagSNPaf)]
		else:
			addSNPs2haploLDpopulations[(chr,loc)].append((pop, maxld, avg , len(v), perfld, neandhaplo, maxldtagSNPaf))

save_obj(addSNPs2haploLDpopulations, 'addSNPs2haploLDpopulations')
print("Dictionary object addSNPs2haploLDpopulations written to binary .pkl")

#print('{}\t{}\t{:0.3f}\t{:0.3f}\t{}\t{}\t{}'.format (chr, loc, pop, maxld, avg , len(v), perfld, neandhaplo))


