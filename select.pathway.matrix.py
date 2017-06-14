import pandas as pd
from sys import argv

script, genes, col_file, matrix,  output =argv

# read in list of genes of interest, list of colnames in the natrix
select=pd.read_csv(genes, sep='\t', header=None, names=['GENE'])
cols=pd.read_csv(col_file, sep='\t', header=None, names=['GENE'])

# find cols of interest
select=cols[cols.GENE.isin(select.GENE)].index.tolist()

# only read in cols of interest from large matrix
DF=pd.read_csv(matrix, sep='\t', usecols=select)

print DF.shape


