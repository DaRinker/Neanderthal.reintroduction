#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 11:39:27 2018
@author: drinker

Extracts all LD (1KGp3 r2) for Neand tag SNPS (as defined in Vernot) with additional SNPS (Vernot 2016) in given population

Outputs tab separated file containing (chromosome, V16_additional_snp_position, r2 in given population, Neand tag SNP position)


"""
usage = """
USAGE:
-h     This message.
"""

import os
import sys
import pymysql
import pickle

POP = sys.argv[1]

## INPUT FILES.

DATAPATH = "/dors/capra_lab/data/archaic_dna/vernot2016/introgressed_tag_snp_frequencies"

if POP == "EAS":
    POP_V16 = "ASN" #Vernot uses ASN for East Asians
else:
    POP_V16 = POP
    
TAG_FILE = DATAPATH + "/all_tag_snps." + POP_V16 + ".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed"
SNP_FILE = DATAPATH + "/all_tag_snps." + POP_V16 + ".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed.extended_LD.sorted"
haplotype_list_FILE = DATAPATH + "/all_tag_snps." + POP_V16 + ".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.median_af.bed.extended_LD"

#####_FUNCTIONS__#########

def bed_2_ChrPos (in_file):
    out_list=[]
    for line in open(in_file):
        seg = line.rstrip().split('\t')
        CHR = seg[0].split('r')[1]
        POS1 = int(seg[2])
        HAPLO = seg[-1]
        out_list += [(CHR,POS1,HAPLO)]
      
    return out_list

def msqlquery (MYSQLTABLE, CHR, POS1, POS2): 
    rsq = []
    cursor = cnx.cursor()
    query = ("SELECT Rsquared FROM "+ MYSQLTABLE + " WHERE (CHROM = %s) AND (POS1 = %s) AND (POS2 = %s)")
    cursor.execute(query, (CHR,POS1,POS2))
    for row in cursor:
        rsq = (float(row['Rsquared']))
        
    return rsq

def save_obj(obj, savefilename ):
        with open ('obj/' + savefilename + '.pkl', 'w+b') as f:
#        with open (savefilename + '.pkl', 'w+b') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL) #pickle.HIGHEST_PROTOCOL is a binary format; protocol 0 zero is a text form

####################
tagsnp_list = []    #all Neand. tag SNP positions from Vernot 2016
ld_snp_list = []    #all Neand. tag SNP positions along with all 1KG SNPS in LD with r2>0.8 from Vernot 2016
additional_ld_snps = []  #only SNPs in LD r2>0.8 with Neand. tag SNP
haplotypes= []             #set of all Neand haplotypes that correlate with introgressed SNP (LD r2>0.8 according to Vernot)
hap2tagsnp={}     #dict of tag SNPs keyed on linked haplotype_listlotype
hap2snp={}          #dict of non-tag SNPs keyed on linked haplotype_listlotype
hap2rsq = {}    #dic3

tagsnp_list = bed_2_ChrPos(TAG_FILE)
ld_snp_list = bed_2_ChrPos(SNP_FILE)
additional_ld_snps = set(ld_snp_list) - set(tagsnp_list)

haplotypes=set([x[2] for x in additional_ld_snps])


    
for chrom, pos, haplotype in tagsnp_list: 
    
    if haplotype not in hap2tagsnp:
        hap2tagsnp[haplotype] = [(chrom,pos)]
    else:
        hap2tagsnp[haplotype].append((chrom,pos))

for chrom, pos, haplotype in additional_ld_snps:
    if haplotype not in hap2snp:
        hap2snp[haplotype] = [(chrom, pos)]
    else:
        hap2snp[haplotype].append((chrom,pos))

## Get LD partners.

if POP=='EUR':
    DATABASETABLE='SNP_LD_1KG_PHASE3v5'
else:
    DATABASETABLE='SNP_LD_1KG_PHASE3v5_' + POP

cnx = pymysql.connect(host='chgr2.accre.vanderbilt.edu',
                             user='colbrall',
                             password='wretched-calculator',
                             db='1kg_snpld_phase3',
                             charset='utf8mb4',
                             cursorclass=pymysql.cursors.DictCursor)

for haplotype in haplotypes: 
    tagsnps = hap2tagsnp[haplotype]
    snps = hap2snp[haplotype]
    for tagsnp in tagsnps: # for tag_snp in tag_snps:
        for snp in snps: # for cand_snp in cand_snps:
            CHR=snp[0]
            if tagsnp[1] < snp[1]:
                POS1 = tagsnp[1]
                POS2 = snp[1]
            else:
                POS1 = snp[1]
                POS2 = tagsnp[1]
            rsq = msqlquery(DATABASETABLE, CHR, POS1, POS2)
    
            if haplotype not in hap2rsq:
                hap2rsq[haplotype] = [(CHR,snp[1],rsq,tagsnp[1])]
            else:
                hap2rsq[haplotype].append((CHR,snp[1],rsq,tagsnp[1]))
        
cnx.close()

save_obj(hap2rsq, POP + "hap2rsq.dict")