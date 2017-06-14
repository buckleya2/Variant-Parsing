
'''
script that converts individual TCGA VCF files into one wide dataframe (similar to parse.vcf.readcounts.py)
takes in one long file with chr, pos,ref,alt,ref count, alt count, and ID for all TCGA samples
converts to wide format with individual entries with ref and alt counts 
ex: 54:23 (ref reads : alt reads )
'''

from sys import argv
script, table, output =argv
import pandas as pd


DF=pd.read_csv(table, sep="\t")

# create identifier column 
DF['IDEN']=DF.CHR.map(str) + "-" + DF.POS.map(str) + "-" + DF.REF.map(str) + "-" + DF.ALT.map(str)
DF['GENO']=DF.REF_READS.map(str) + ":" + DF.ALT_READS.map(str)

# reshape from wide to long
DFsmall=DF[['IDEN','GENO','ID']]

# because of how I combined VCFs of repeated samples, some duplicate emtries w/ different read counts
# choose first
WIDE=pd.pivot_table(DFsmall, aggfunc="first", index='IDEN', columns="ID", values="GENO").reset_index()

# fill non variant sites to 0
WIDE.fillna(0, inplace=True)

WIDE['CHR']=WIDE.IDEN.apply(lambda x : str(x).split("-")[0])
WIDE['POS']=WIDE.IDEN.apply(lambda x : str(x).split("-")[1])
WIDE['REF']=WIDE.IDEN.apply(lambda x : str(x).split("-")[2])
WIDE['ALT']=WIDE.IDEN.apply(lambda x : str(x).split("-")[3])

WIDE.drop('IDEN', axis=1, inplace=True)

WIDE.to_csv(output, index=False, sep='\t')
