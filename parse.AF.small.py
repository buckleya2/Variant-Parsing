from sys import argv

script, vcf, output =argv

'''
script that converts vcf coordinates with multiallelic sites to output with a single alt allele per line
input:tsv with chromosome, position, ref, alt, as first 4 columns
output:tsv with one line for each alt allele and a column indicating what number alt allele for each line
ex: A T/TTC/TTCCC  T = 1, TTC = 2, TTCCC = 3
'''

OUT=open(output, 'w')

v = open(vcf)
# read through file line by line
# make list from alt alleles
# convert alt genotypes to alleles by using alt number as list index
# if genotype is ref, fill with "NA"
for v_line in v.readlines():
    parsed=v_line.strip().split("\t")
    chr=parsed[0]
    pos=parsed[1]
    ref=parsed[2]
    alt=parsed[3]
    alt_list=alt.split(",")
#accomodate files with more than 4 columns
    extracol=[]
    if len (parsed) > 4:
        extracol=parsed[4:]
#If multilple alternate aleles, parse each alt allele separately
    if len(alt_list)>1:
	for i in range(len(alt_list)):
        	alti=alt_list[i]
#remove sites with GATK * notation (sites not present due to a nearby deletion in that sample)
		if alti == "<*:DEL>" or alti == "*":
                    continue
#parse snps coded as pseudo indels at multiallelic sites
#ex: ATGC/TTGC to A/T
                if (len(ref )>1) & (len(ref) == len(alti)):
                    if ref[1:]== alti[1:]:
                        refn=ref[0]
			alti=alti[0]
			format="%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,refn,alti,i+1,'\t'.join(map(str,extracol)))
			OUT.write(format)
                else:
		    format="%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alti,i+1,'\t'.join(map(str,extracol)))
                    OUT.write(format)
    else:
    	format = "%s\t%s\t%s\t%s\t1\t%s\n" % (chr,pos,ref,alt,'\t'.join(map(str,extracol)))
    	OUT.write(format)
