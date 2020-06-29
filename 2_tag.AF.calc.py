# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:31:36 2018

@author: drinker

this script will create a dictionary keyed on Vernot 2016 tag SNPs from each populaiton; the values are the AF of the tags in the respective populaitons

this dictionary will feed into 'LDsnps_summarize.py', providing the tagSNP AF for the maxLD tagSNP for each addSNP; this will be used to infer the "introgressed" variant at each addSNP position...

"""
usage = """
USAGE:
-h	 This message.
"""

import os
import sys
import pickle
import numpy as np

tagSNP2AF = {}

POPULATIONS = ['EAS', 'EUR','SAS'] 
# POPULATIONS = ['EUR'] 

for POP in POPULATIONS:
	print("Beginning " + POP)
	for line in open("../data/" + "V16_"+ POP + "_tagSNPs_pos_AF.txt"):
		seg = line.rstrip().split('\t')
		chr = seg[0]
		loc = int(seg[1])
		AFinPOP = float(seg[2])
		tagSNP2AF[(chr, loc, POP)] = [AFinPOP]

with open ('obj/tagSNP2AF.pkl', 'w+b') as f:
	pickle.dump(tagSNP2AF, f, pickle.HIGHEST_PROTOCOL) #pickle.HIGHEST_PROTOCOL is a binary format; protocol 0 zero is a text format

###################

