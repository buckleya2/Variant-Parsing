'''
script that takes in binary mutation indicator files for two different mutation types and finds co-occurances
ex: germline and somaic mutations in the same gene
takes in files produced by agg.bygene.py
format of input files SHORT_ID, GENE, GENO, TYPE (ind is a binary mutation indicator and type is a string that distingusihes mutation types)
	ex: A704   TP53    1    SOMATIC
'''

import pandas as pd
from sys import argv

script, df1, df2, output = argv

# read in mutation indicator files
DF1=pd.read_csv(df1, sep='\t')
DF2=pd.read_csv(df2, sep='\t')

# combine
CONCAT=pd.concat([DF1,DF2])

# if GENO column is count and not an indicator, change to 1/0 indicator values 
CONCAT.GENO[CONCAT.GENO > 1]=1

# aggregate mutation types by gene, ID
AGG=CONCAT.groupby(['SHORT_ID','GENE']).agg({'GENO' : sum}).reset_index()

# pull out double hit gene-ID pairs
double=AGG[AGG.GENO == 2]

# write to file
double.to_csv(output, sep='\t', index=False)
