# create and activate the conda environment necessary to run everything
conda env create -f environment.yml
source activate tempus

# build vep cache for GRCh37
vep_install -a cf -s homo_sapiens -y GRCh37 # can use ctrl+C to prevent cache conversion, will save time, though could make gnomAD annotation faster

# VEP annotation step
# we should be using gnomAD 2.1 instead of ExAC... (more up-to-date, but there is no REST API for gnomAD as of yet)
# also grabs sequence ontology consequence, gene name, transcript, allele number as the VCF is not decomposed 
vep -i Challenge_data\ \(1\).vcf --cache --symbol --allele_number -o challenge-vep.vcf --vcf --fields Consequence,SYMBOL,Feature,ALLELE_NUM --offline --fork 12 --force_overwrite

# thoroughly annotates each allele since we were not permitted to use vt decompose/normalize; I also suspect many more variants would have exac AFs (that is, INDELs) if they were properly normalized
python annotate.py challenge-vep.vcf > final.vcf

#################################################################################################################################
# FINAL FILE ("final.vcf") has "Custom" column at the end of the columns.  This contains:
# DELCSQ - the most deleterious consequence for the variant for each alternate allele
# exacAF - the ExAC allele frequency for each alternate allele
# normaldepth - the "normal" sample read depth
# vaf5depth - the "vaf5" sample read depth
# normalaltreads - the "normal" sample's alternate read count support for each alternate allele
# vaf5altreads - the "vaf5" sample's alternate read count support for each alternate allele
# normalaltpctsupport - the "normal" sample's percentage of read support for each alternate allele relative to reference reads
# vaf5altpctsupport - the "vaf5" sample's percentage of read support for each alternate allele relative to reference reads
################################################################################################################################
