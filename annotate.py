"""
This program will annotate the VCF with variant impact information and read support calculations
Also note that it would be much easier and optimal if permitted to vt decompose and normalize and then use cyvcf2
CSQ ranking taken from: https://gemini.readthedocs.io/en/latest/content/database_schema.html#details-of-the-impact-and-impact-severity-columns
"""
import sys
import requests
import json
csqrank=["exon_loss_variant",
"frameshift_variant",
"splice_acceptor_variant",
"splice_donor_variant",
"start_lost",
"stop_gained",
"stop_lost",
"initiator_codon_variant",
"initiator_codon_variant",
"rare_amino_acid_variant",
"chromosomal_deletion",
"missense_variant",
"protein_altering_variant",
"inframe_insertion",
"inframe_deletion",
"coding_sequence_variant",
"disruptive_inframe_deletion",
"disruptive_inframe_insertion",
"5_prime_UTR_truncation + exon_loss_variant",
"3_prime_UTR_truncation + exon_loss_variant",
"splice_region_variant",
"mature_miRNA_variant",
"regulatory_region_variant",
"TF_binding_site_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
"TFBS_ablation",
"TFBS_amplification",
"stop_retained_variant",
"synonymous_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"intron_variant",
"coding_sequence_variant",
"upstream_gene_variant",
"downstream_gene_variant",
"intergenic_variant",
"intragenic_variant",
"gene_variant",
"transcript_variant",
"exon_variant",
"5_prime_UTR_premature_start_codon_gain_variant",
"start_retained_variant",
"conserved_intron_variant",
"nc_transcript_variant",
"non_coding_transcript_variant",
"non_coding_transcript_exon_variant",
"NMD_transcript_variant",
"incomplete_terminal_codon_variant",
"non_coding_exon_variant",
"transcript_ablation",
"transcript_amplification",
"feature_elongation",
"feature_truncation"]
VCF = open(sys.argv[1], "r") # no need to use argparse for one specific task for one toy file right?
j = 0 # to grab the format as a dict (Could vary depending on file)
url = "http://exac.hms.harvard.edu/rest/variant/variant/"
for variant in VCF:
# this is all already done by the tool cyvcf2 which is permitted for commercial use, but was told not to use external sources here
    if variant.startswith("#") and not variant.startswith("#CHROM"):
        print (variant, end="")
    if variant.startswith("##INFO=<ID=CSQ"):
        # parses the CSQ string in the header into keys for a dict for parsing variant impact and allele frequency
        kcsq = variant.strip('\n ">').split(":")[1].strip(' "').split("|")
    elif variant.startswith("#CHROM"):
        # parses the column header into keys for a dict for the fields
        vcffields = variant.strip("# \n").split("\t")
        print (variant.strip() + "\t" + "Custom")
    elif not variant.startswith("#"):
        # creates the dict for columns
        columns = variant.strip().split("\t")
        fields = dict(zip(vcffields, columns))
        # creates the dict for CSQ strings
        rawcsq = fields["INFO"].split("CSQ=")[1].split(";")[0]
        csqs = [dict(zip(kcsq, c.split("|"))) for c in rawcsq.split(",")]
        """
        Creates the fields for:
           - Depth of sequence coverage at the site of variation (grabbed from DP or DPR).
           - Number of reads supporting the variant (grabbed from AO or DPR).
           - Percentage of reads supporting the variant versus those supporting reference reads (grabbed from DPR).
        for each allele, for normal and vaf5.
        Using the format
        GT:GQ:DP:DPR:RO:QR:AO:QA
        DPR (1 REF, 3 ALT), AO (3 ALT), QA (3 ALT)
        """
        if j == 0:
            formats = fields["FORMAT"].split(":")
            j += 1
        normal = dict(zip(formats,fields["normal"].split(":")))
        ndpr = normal["DPR"].split(",")
        vaf5 = dict(zip(formats,fields["vaf5"].split(":")))
        vdpr = vaf5["DPR"].split(",")
        # preps for multi allelics
        allelecsq = []; alleleaf = [];
        # preps for read data
        naltreads = []; valtreads = [];
        npctsupport = []; vpctsupport = [];
        # depth is based on reference, which doesn't change from allele to allele
        ndepth = float(ndpr[0])
        vdepth = float(vdpr[0])
        # for getting allele frequency
        alts = fields["ALT"].split(",")
        for i in range(0,len(alts)): # unfortunately, not decomposed so we have to iterate through alleles
            # getting the read information for each allele
            nalt=float(ndpr[i+1]); valt=float(vdpr[i+1])
            # normal and vaf5 alt reads 
            naltreads.append(str(nalt))
            valtreads.append(str(valt))
            # normal and vaf5 percent read support for alternate vs reference
            npctsupport.append("%.3f" % (nalt/ndepth))
            vpctsupport.append("%.3f" % (valt/vdepth))
            # sets maximum possible index
            maxcsq = len(csqrank)-1
            # grabs exac_AF for the allele (same for all consequences)
            response = requests.get(url+fields["CHROM"]+"-"+fields["POS"]+"-"+fields["REF"]+"-"+alts[i])
            jdata = json.loads(response.content)
            try: 
                exacaf = "%.3f" % jdata['allele_freq']
            except KeyError:
                exacaf = "None"
            for csq in csqs:
                if int(csq["ALLELE_NUM"]) != i+1:
                    continue
                else:
                    # grabs index of maximum worst consequence of variant for the allele in csqrank
                    for c in csq["Consequence"].split("&"):
                        maxcsq = min(maxcsq, csqrank.index(c))
            allelecsq.append(csqrank[maxcsq])
            alleleaf.append(exacaf)
        print(variant.strip() + "\t" + \
             ";".join(["DELCSQ="+",".join(allelecsq), \
                        "exacAF="+",".join(alleleaf), \
                        "normaldepth="+str(ndepth), \
                        "vaf5depth="+str(vdepth), \
                        "normalaltreads="+",".join(naltreads), \
                        "vaf5altreads="+",".join(valtreads), \
                        "normalaltpctsupport="+",".join(npctsupport), \
                        "vaf5altpctsupport="+",".join(vpctsupport)]))
