
'''
script that takes in variant positions with one alt allele per line and parses down to smallest indel notation ( a problem at multialleleic sites)
ex:  CG/C,TG 
     CG/C   remains 
     CG/TG  becomes  C/T 
input: tsv with chromsome, position, ref, alt as first 4 columns
output: tsv with new coordinates
'''

import fileinput

for line in fileinput.input():
    parsed=line.strip().split('\t')
    chr=parsed[0]
    pos=parsed[1]
    ref=parsed[2]
    alt=parsed[3]
#accomodate files with more than 4 columns
    extracol=[]
    if len (parsed) > 4:
        extracol=parsed[4:]
#remove sites with GATK * notation (sites not present due to a nearby deletion in that sample)
    if alt == "<*:DEL>" or alt == "*":
        continue
    if ref[-1]== alt[-1] and len(ref) > 1 and len(alt) >1:
        refn=ref
        alti=alt
        while refn[-1] == alti[-1]  and len(refn) > 1 and len(alti) >1:
	    refn=refn[0:-1]
            alti=alti[0:-1]
	format="%s\t%s\t%s\t%s\t%s" % (chr,pos,refn,alti,'\t'.join(map(str,extracol)))
        print format
    else:
	format="%s\t%s\t%s\t%s\t%s" % (chr,pos,ref,alt,'\t'.join(map(str,extracol)))
        print format
