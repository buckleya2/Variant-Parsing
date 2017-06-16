'''
script that uses STRING protein interaction data to look for double hits (germline/somatic)
expand double hit analysis from two hits in the same gene to two hits in interacting genes
takes in two variant indicator files in the format ID GENE COUNT 
	ex: AAVF   TP53   1 
	    AAVF   ABL1   1
outputs long format list with both interacting proteins and sample ID
by defualt the first column indicates the gene with a germline LOF, the second somatic LOF
'''

from collections import defaultdict
from sys import argv

script, output= argv

germfile="/nrnb/users/abuckley/regroup_gvcf/LOF_matrix/germline.QDfilt.racefile.LOF.indicator"
#somfile="/nrnb/users/abuckley/somatic_vcf/LOF_matrix/somatic.LOF.indicator"
somfile="/nrnb/users/abuckley/methylation/firehose/meth_expr/scaled_files/methyl.for.doublehit"
STRING="/nrnb/users/abuckley/two_hit/STRING/STRING.HUGO.no.self"

OUT=open(output, "w")

# make dict of all germline LOF - ID pairs
GERM=defaultdict(list)
with open(germfile) as g:
    for line in g:
        parse=line.strip().split('\t')
	GERM[parse[1]].append(parse[0])

# make dict of all somatic LOF - ID pairs
SOM=defaultdict(list)
with open(somfile) as s:
    for line in s:
        parse=line.strip().split('\t')
        SOM[parse[1]].append(parse[0])

# read through STRING interaction pairs
with open(STRING) as t:
    for line in t:
	parse=line.strip().split('\t')
	gene1=parse[0]
	gene2=parse[1]
# get indivdiuals with germlie hits in gene 1, somatic hits in gene 2
	G=GERM[gene1]
	S=SOM[gene2]
# for each indivdiual with a germline and somatic hit in the interaction pair 
# write indivdiual ID and interaction pair
	double_hit=set(G).intersection(S)
	for i in double_hit:
	    format= "%s\t%s\t%s\n" % (gene1, gene2, i)
	    OUT.write(format)
