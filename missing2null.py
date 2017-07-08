import numpy as np
from sys import argv

script, vcf, output=argv

# here I made a example dict with the format position : IDs
coord_dict={861330: ['H_LR-4E-A92E-10A-01D-A37N-09','H_LR-2E-A9G8-10A-01D-A403-09','H_LR-5S-A9Q8-10A-01D-A403-09']}

# baby function to find and replace ./. to 0/0
def change_geno(x):
    return x.replace("./.", "0/0")

# need to vectorize to use on np.vector objects
change_geno_v=np.vectorize(change_geno)

#open output file
OUT=open(output, 'w')

# read through VCF line by line
with open(vcf) as v:
    for v_line in v:
# write header lines to file
        if "##" in v_line:
            OUT.write(v_line)
            continue
#  parse VCF IDs into np.aray           
        if "#CHROM" in v_line:
            sample_name=np.array(v_line.strip().split("\t")[9:])
            OUT.write(v_line)
            continue
        parsed=v_line.strip().split("\t")
        pos=int(parsed[1])
        info=parsed[:9]
        parsed_fields=np.array(parsed[9:])
# check if position is one to modify
	if pos in coord_dict: 
# here find the sample IDs for the samples you want to change to 0/0
            to_fix=coord_dict[pos]
    
# this makes a bool array indicating which samples to change    
            geno_mask=np.in1d(sample_name,to_fix)

# only change the fields of interest
            parsed_fields[geno_mask] = change_geno_v(parsed_fields[geno_mask])

# write new line to file
            format="%s\t%s\n" % ('\t'.join(map(str,info)), '\t'.join(map(str,parsed_fields)))
            OUT.write(format)
# if not, just write the line to file
        else:
	    OUT.write(v_line)

