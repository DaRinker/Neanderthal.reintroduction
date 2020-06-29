#!/usr/bin/env python3
"""
@author: drinker

this script will calculate allele frequencies of variants of interest in sub Saharan subpopulations from 1000 Genomes phase 3
"""

import sys
import gzip
import pickle

SAMPLE_ID_FILE = "<path to 1000 Genomes phase 3 sample panel>"


POPS_OF_INTEREST = ['ESN', 'GWD', 'LWK', 'MSL', 'YRI']
SAMPLES_OF_INTEREST = []

#####_FUNCTIONS__#########

def save_obj(obj, savefilename ):
        with open ('obj/' + savefilename + '.pkl', 'w+b') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL) #pickle.HIGHEST_PROTOCOL is a binary format; protocol 0 zero is a text format

##########################

addSNPs2haploLDpopulations = {}
addSNPs2haploLDpopulations = pickle.load(open('obj/addSNPs2haploLDpopulations.pkl', 'rb'))

print("addSNPs2haploLDpopulations loaded")


with open(SAMPLE_ID_FILE, 'r') as f:
	for line in f:
		if line.startswith('#'):
			continue
		elif not line.startswith('#'):
			indiv_id, family_id, pop, popdescription, gender = line.strip().split('\t')
			if pop in POPS_OF_INTEREST:
				SAMPLES_OF_INTEREST.append(indiv_id)

addSNPs2AFinAFR = {}
CHROMOSOMES = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

for CHROM in CHROMOSOMES:
	INDICES = []
	with gzip.open('../data/1kgph3_chr' + CHROM + '.gz', 'rt') as f:
		for line in f:
			t = line.strip().split()
			if line.startswith('#CHR'):
				print(t[0:7])
				for sample in SAMPLES_OF_INTEREST:
					if sample not in t:
						continue
					else:
						INDICES.append(t.index(sample))

				print(len(INDICES))

			elif not line.startswith('#'):
				
				chrom, loc, rsid, ref, alt, qual, filt, data, alleles = t[:9]
				loc = int(loc)
				if (CHROM, loc) not in addSNPs2haploLDpopulations:
					# print("{} {} not in dict......".format(CHROM,loc))
					continue
				else:
					genotype = [t[i] for i in INDICES]
					afr_allele_freq = "".join(genotype).count('1') / (2. * len(genotype))
					# print(CHROM, loc, rsid, afr_allele_freq)
					addSNPs2AFinAFR[(CHROM, loc)] = [(afr_allele_freq)]

save_obj(addSNPs2AFinAFR, 'addSNPs2AFinAFR')
