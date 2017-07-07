from sys import argv
script, vcf, output =argv
import numpy as np
from collections import defaultdict
import gzip

'''
script to get allele counts from VCF binned by race
requires a file that maps VCF IDs to race
here getting AC for AFR (african), NFE( non-Finniah European), and EAS (East Asian) from TCGA germline VCF
'''

# make a dictionary that maps each TCGA ID to the PCA predicted ancestry
ethnicity_file="/cellar/users/abuckley/ref/TCGA.pred.ancestry.dedup"

mapping=defaultdict(lambda : "NA")
with open(ethnicity_file) as E:
    for line in E:
        ID=line.strip().split('\t')[0]
        ETH=line.strip().split('\t')[1]
        mapping[ID]=ETH
 
OUT=open(output, 'w')

def parse_geno(vcf_field, alt_num=1):
    line=vcf_field.split(":")
    geno=line[0]
    if geno == "./.":
        return 0
    else:
        allele_count=geno.count(str(alt_num))
        return allele_count

with gzip.open(vcf) as v:
    for v_line in v:
# find and write sample file names to output file
        if "#CHROM" in v_line:
            sample_name=v_line.strip().split("\t")[9:]
# map sample names to ethnicity
            eth_array=np.array([mapping[i] for i in sample_name])
# create np masks for each ethnicity
            AFR=np.where(eth_array == "BLACK")
            NFE=np.where(eth_array == "WHITE")
            EAS=np.where(eth_array == "ASIAN")
# write a header line
            OUT.write("CHR\tPOS\tREF\tALT\tAC_AFR\tAC_NFE\tAC_EAS\n")
	    continue
        elif "#" in v_line:
	    continue
        parsed=v_line.strip().split("\t")
        chr=parsed[0]
        pos=parsed[1]
        ref=parsed[3]
        alt=parsed[4]
        alt_list=alt.split(",")
        samples=parsed[9:]

# If multilple alternate aleles, parse each alt allele separately
        if len(alt_list)>1:
            for i in range(len(alt_list)):
                alti=alt_list[i]
                altnum=int(i)+1
# remove sites with GATK * notation (sites not present due to a nearby deletion in that sample)
                if alti == "<*:DEL>" or alti == "*":
                    continue
# parse snps coded as pseudo indels at multiallelic sites
# ex: ATGC/TTGC to A/T
                if (len(ref )>1) & (len(ref) == len(alti)):
                    if ref[1:] == alti[1:]:
                        ref=ref[0]
                        alt=alti[0]
                        parsed_fields=np.array([parse_geno(j,alt_num=altnum) for j in samples])
		        AC_AFR=sum(parsed_fields[AFR])
		        AC_NFE=sum(parsed_fields[NFE])
                        AC_EAS=sum(parsed_fields[EAS])
                        format="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alt, AC_AFR, AC_NFE, AC_EAS)
                        OUT.write(format)
                else:
                    parsed_fields=np.array([parse_geno(j,alt_num=altnum) for j in samples])        
		    AC_AFR=sum(parsed_fields[AFR])
                    AC_NFE=sum(parsed_fields[NFE])
                    AC_EAS=sum(parsed_fields[EAS])
                    format="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alti, AC_AFR, AC_NFE, AC_EAS)
	            OUT.write(format)
        else:
            parsed_fields=np.array([parse_geno(j) for j in samples])
	    AC_AFR=sum(parsed_fields[AFR])
            AC_NFE=sum(parsed_fields[NFE])
            AC_EAS=sum(parsed_fields[EAS])
            format="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alt, AC_AFR, AC_NFE, AC_EAS)
            OUT.write(format)
