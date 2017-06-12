'''
takes stdout from samtools mpileup -uv -t INFO/AD 
outputs variable sites only, splits multiallaleic sites to different lines
outputs ref 
'''

def parse_info(x):
    info_dict={}
    for i in x:
        tag=i.split("=")[0]
        val=i.split("=")[1]
        info_dict[tag]=val
    return info_dict['AD']  

'''
if ref and alt end in the same sequence, trim off base pairs from the end one at a time to get smallest indel notation
always leave at least one ref/alt base per GATK VCF notation
ex:     AGATGAT AGAT
becomes AGAT    A
'''

def simplify_indel(ref,alt):
    if ref[-1]== alt[-1] and len(ref) > 1 and len(alt) >1:
        while ref[-1] == alt[-1]  and len(ref) > 1 and len(alt) >1:
            ref=ref[0:-1]
            alt=alt[0:-1]
        return(ref,alt)
# if indel can't be shortened, output as-is
    else:
        return(ref,alt)

import fileinput

#print header line
format= "CHR\tPOS\tREF\tALT\tREF_COUNT\tALT_COUNT"
print format

for line in fileinput.input():
#skip vcf header
    if "#" in line:
        continue
    parsed=line.strip().split('\t')
    chr=parsed[0]
    pos=parsed[1]
    ref=parsed[3]
    alt=parsed[4]
    info=parsed[7]

    read_depth=parse_info(info.split(";")[1:])
    ref_depth=read_depth.split(",")[0]
    alt_depth=read_depth.split(",")[1:]

#skip sites with no variation
    if alt == "<*>":
        continue
#for multiallelic sites, output each allele on a separate line
    if "," in alt:
        alt_list=alt.split(',')
        for i in range(len(alt_list)):
#skip sites with no variation
	    alti=alt_list[i]
            alt_depthi=alt_depth[i]
            if alti == "<*>":
                continue
	    elif len(ref) > 1 or len(alt) > 1:
		ref,alti=simplify_indel(ref,alti) 
            format= "%s\t%s\t%s\t%s\t%s\t%s" % (chr,pos,ref,alti,ref_depth,alt_depthi )
            print format
        
    elif len(ref) > 1 or len(alt) > 1:
        ref,alt=simplify_indel(ref,alt)
	format= "%s\t%s\t%s\t%s\t%s\t%s" % (chr,pos,ref,alt,ref_depth,alt_depth[0] )
        print format
    else:
	format= "%s\t%s\t%s\t%s\t%s\t%s" % (chr,pos,ref,alt,ref_depth,alt_depth[0] )
        print format
