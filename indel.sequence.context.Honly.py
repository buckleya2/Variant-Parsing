'''
Script that looks for homopolymer tracts in sequence surrounding variants
defines homopolymer as 4 conseutive repeats of any base pair
takes in a TSV with chrom,pos,ref,alt
outputs a file with varaint coordinates, hompolymer start pos, end pos, base content, hompolymer call (see check_hp function), GC content, and sequence +/- 10 bp from variant start pos
'''

from sys import argv
from pyfaidx import Fasta
import re

#parse command line args
script, infile, output =argv

#open file to write results
OUT=open(output, 'w')

# function to check for A,T,C,G homopolymer tracts
# take in string of sequence, outputs start and end position of the longest homopolymer tract in input sequence
# start and end position are relative to indel start position (ex: 10 = 10 bp upstream, -10 10 bp downstream)
# use start pos = -1  and homopolymer base = inserted base as indication of high confidence hompolymer indel
# use start pos = -1  and homopolymer base != inserted base as indication of low confidence hompolymer indel
# ex: A/ATTT   TAAAGCGTCAATTTTTATGCT = HC homopolymer indel

def check_hp( genome, ref, alt ):
    startpos="NA"
    endpos="NA"
    base="NA"
    maxval=0
    reobj="A{4,}|T{4,}|C{4,}|G{4,}"
    match=re.finditer(reobj, genome)
    for i in match:
        size=len(i.group())
        if size > maxval:
            startpos=10 - i.start()
            endpos=10 - i.end()
            maxval=size
            base=i.string[i.start()]
    if startpos == -1 and (ref[-1] == base or alt[-1] == base):
        call="HC"
    elif startpos == -1 and (ref[-1] != base and alt[-1] != base):
        call="LC"
    else:
        call="NON"
    return {'START':startpos, 'END':endpos, 'BASE' : base, 'CALL' : call}


#Function to get GC fraction
#takes in string of sequence, outputs GC fraction
def get_GC(x):
    GC=(x.count('G')+x.count('C'))/ float(len(x))
    return(GC)

#create fasta object 
grch=Fasta('/cellar/users/abuckley/ref/hs37d5.fa')


# write header to file
OUT.write('CHR\tPOS\tREF\tALT\tSTART\tEND\tBASE\tCALL\tGC\tSEQUENCE\n')

#read theough input file line by line
#use pyfaidx to pull out sequence +/- 10 bp from indel start position
#look for homopolymers in flanking sequence
with open(infile) as coord:
	for coord_line in coord:
	    if "POS" in coord_line:
		continue
	    parsed=coord_line.strip().split("\t")
    	    chr=str(parsed[0])
    	    pos=int(parsed[1])
    	    ref=parsed[2]
    	    alt=parsed[3]
    	    genome=grch[chr][pos-11:pos+10].seq
    	    GC=get_GC(genome)
    	    out=check_hp(genome, ref, alt)
    	    format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alt,out['START'],out['END'],out['BASE'], out['CALL'], GC,genome )
    	    OUT.write(format)
