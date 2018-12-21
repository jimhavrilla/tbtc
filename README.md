Overview
========
This repository is for the Tempus Bioinformatics Technical Challenge as described in the TempusBioinformaticsChallenge.pdf.

The instructions were:

For this challenge, you are asked to prototype a variant annotation tool. We will provide you with a VCF file, and you will create a small software program to output a table annotating each variant in the file. Each variant must be annotated with the following pieces of information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.

Additional instructions were "It is acceptable to use variant annotation tools  (snpEff, Oncotator, VEP, etc.) for addition of gene, effect and transcript information. All I/O, transformations, and external source interactions must be done with standard libraries (numerical libraries or data frame libraries such as pandas are also acceptable).

Please use this opportunity to demonstrate your knowledge of software development and your abilities to think critically about how a toy program like this should be engineered."

How to Run
===========
Run `build.sh` to create the environment if you like, but the `environment.yml` file uploaded simplifies the process and all you need to run is simply :
```
bash annotate.sh
```
to create the environment from the YML file, it installs VEP, annotates with VEP, then runs `annotate.py` and creates the final resulting VCF file, `final.vcf`.


#### FINAL FILE ("final.vcf") has "Custom" column at the end of the columns.  This contains:
Column              | Description |
--------            | ----------- |
DELCSQ              | the most deleterious consequence for the variant for each alternate allele
exacAF              | the ExAC allele frequency for each alternate allele
normaldepth         | the "normal" sample read depth
vaf5depth           | the "vaf5" sample read depth
normalaltreads      | the "normal" sample's alternate read count support for each alternate allele
vaf5altreads        | the "vaf5" sample's alternate read count support for each alternate allele
normalaltpctsupport | the "normal" sample's percentage of read support for each alternate allele relative to reference reads
vaf5altpctsupport   | the "vaf5" sample's percentage of read support for each alternate allele relative to reference reads
