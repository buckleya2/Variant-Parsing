'''
script to translate EMSEMBL protein IDs to HUGO gene symbol
specifically translate STRING protein-protein interaction data file to HUGO symbol
'''

from collections import defaultdict
from sys import argv

script, output = argv

# file with ENSEMBL protein ID as first column, HUGO symbol as second
translation_file="/nrnb/users/abuckley/two_hit/STRING/ENSEMBL.protein.map"
# STRING protein protein interaction file with EMSEBL protein IDs as both columns
to_translate="/nrnb/users/abuckley/two_hit/STRING/STRING.parsed"

# open file to write trasnalted values
OUT=open(output, 'w')

# initialize a defaultdict dictionary
# when translating ENSEMBL ID to HUGO symbol, if the ENSEMBL symbol is not in dict, will return "NA" 
translate=defaultdict(lambda: "NA")
with open(translation_file) as f:
    for line in f:
	translate[line.strip().split("\t")[0]] = line.strip().split("\t")[1]

# for each interaction pair, find mapping to HUGO ID using translation dictionary
with open(to_translate) as f:
    for line in f:
        new_1=translate[line.strip().split("\t")[0]]
        new_2=translate[line.strip().split("\t")[1]]
	format= "%s\t%s\n" % (new_1, new_2)
	OUT.write(format)
