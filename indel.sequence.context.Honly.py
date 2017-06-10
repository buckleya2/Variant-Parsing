'''
Script that looks for homopolymer tracts in sequence surrounding variants
defines homopolymer as 4 conseutive repeats of any base pair
takes in a TSV with chrom,pos,ref,alt
outouts a file with coordinates, hompolymer start and end pos, GC content, and sequence +/- 10 bp from variant start pos
'''

from sys import argv
from pyfaidx import Fasta
import re

#parse command line args
script, infile, output =argv

#open file to write results
OUT=open(output, 'w')

#Function to check for A,T,C,G homopolymer tracts
#Take in string of sequence, outputs start and end position of the longest homopolymer tract in input sequence
def check_hp( genome ):
    startpos="NA"
    endpos="NA"
    maxval=0
    reobj="A{4,}|T{4,}|C{4,}|G{4,}"
    match=re.finditer(reobj, genome)
    for i in match:
        size=len(i.group())
        if size > maxval:
            startpos=10 - i.start()
            endpos=10 - i.end()
            maxval=size
    return {'START':startpos, 'END':endpos }


#Function to get GC fraction
#takes in string of sequence, outputs GC fraction
def get_GC(x):
    GC=(x.count('G')+x.count('C'))/ float(len(x))
    return(GC)

#create fasta object 
grch=Fasta('/cellar/users/abuckley/ref/hs37d5.fa')

#read theough input file line by line
#use pyfaidx to pull out sequence +/- 10 bp from indel start position
#look for homopolymers in flanking sequence
coord= open(infile)
for coord_line in coord.readlines():
    parsed=coord_line.strip().split("\t")
    chr=str(parsed[0])
    pos=int(parsed[1])
    ref=parsed[2]
    alt=parsed[3]
    genome=grch[chr][pos-11:pos+10].seq
    GC=get_GC(genome)
    out=check_hp(genome)
    format="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alt,out['START'],out['END'],GC,genome )
    OUT.write(format)
