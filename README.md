# Neanderthal.reintroduction
Companion data files and scripts to Rinker DC, et al (2020). Neanderthal introgression reintroduced functional ancestral alleles lost in Eurasian populations. Nature Ecology and Evolution  

RESULTS
* Rinker_Suppl.File1.tar.gz -- tab delimited (.bed) files containing introgressed alleles in each of the three 1000 Genomes Eurasian populations (EAS,EUR,SAS).

SCRIPTS
* 1_Vernot16tags_to_1KG_LD_extract.py -- extracts populations specific LD (r2) for all introgressed, Neanderthal tag SNPs previously identified (Vernot 2016)
* 2_tag.AF.calc.py creates --extracts allele frequencies for tag SNPs from 1000 Genomes superpopulations
* 3_LDsnps_summarize_v2.py --summarized LD information for all identified tagSNPs in each population
* 4_gatherVCFinfo_Altai.py --extract Alti genotype informaiton all positions identified 
* 5_gatherVCFinfo_1kg.py --associates tag SNPs, high LD SNPs, and allele frequencies
* 6_ancestry_AF_calc.py --calculates sub-Saharan allele frequencies for all SNPs
* 7_ID_RAcandidates.py --compiles all gathered informaiton into summary files

The scripts are used to implement workflow illustrated in Extended Data Figure 1. The individual scripts must be run to completion sequentially before the next script may begin (however most of the individual scripts are intended to be run in parallel on a per chromosome basis).NOTE To run the scripts, the user must first have locally available:
  1) 1000 genomes vcf files (phase3, hg19)
  2) a local LD database (mysql) constructed from 1000 Genomes variants
  3) Altai Neanderthal vcf files (hg19)
  4) Neanderthal tag SNPs reported by Vernot 2016 (bed format)

