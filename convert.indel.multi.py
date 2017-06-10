from sys import argv

script, infile, output =argv

OUT=open(output, 'w')

from pyfaidx import Fasta

'''
script that converts indels from "-" format to VCF format where preceding base is reported
ex: AT/- vs. CAT/C;  -/TT vs. A/ATT
input: tsv with chromosome,position,ref,alt of sites to be converted (GRCh37 reference)
output: tsv with translated chromsome, position, ref, alt, as well as original position, ref, alt
'''
 
#read in grch37 fasta as dict object
grch=Fasta('/cellar/users/abuckley/ref/hs37d5.fa')

coord= open(infile)
prevpos=0
#for each line in input file, if line is an indel, find preceding base in reference sequence
for coord_line in coord.readlines():
    parsed=coord_line.strip().split("\t")
    chr=parsed[0]
    pos=parsed[1]
    ref=parsed[2]
    alt=parsed[3]
#accomodate files with more than 4 columns
    extracol=[]
    if len (parsed) > 4:
        extracol=parsed[4:]
#for insertions, position remains the same

    	if ref == "-":
		startloc=int(pos)-1
		startbase=grch[chr][startloc].seq
		format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,startbase,startbase+alt, pos,ref,alt, '\t'.join(map(str,extracol)))
		OUT.write(format)
        	prevpos=pos
        	prevbase=startbase
#for deletions, start position is one base prior to that reported
    	elif alt == "-":
		startloc=int(pos)-2
		startbase=grch[chr][startloc].seq
		format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,int(pos)-1,startbase+ref,startbase,pos,ref,alt,'\t'.join(map(str,extracol)))	    
		OUT.write(format)
        	prevpos=pos
        	prevbase=startbase

	elif pos == prevpos:
                multiref=prevbase+ref
                multialt=prevbase+alt
                format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,int(pos)-1,multiref,multialt,pos,ref,alt,'\t'.join(map(str,extracol)))
                OUT.write(format)

	elif len(alt) > len(ref):
		startloc=int(pos)-2
                startbase=grch[chr][startloc].seq
                multiref=startbase+ref
                multialt=startbase+alt
                format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,int(pos)-1,multiref,multialt,pos,ref,alt,'\t'.join(map(str,extracol)))
                OUT.write(format)
		prevpos=pos
		prevbase=startbase
    	else:
		format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alt,pos,ref,alt,'\t'.join(map(str,extracol)))
		OUT.write(format)
