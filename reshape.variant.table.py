from sys import argv
script, table, output =argv
import pandas as pd

'''
script that converts parsed vcf files generated by parse.vcf.readcounts.py from wide to long format
removes inviduals with non variant sites
each row represents an individual/ variant site pair
'''

DF=pd.read_csv(table, sep="\t")

# create identifier column 
DF['IDEN']=DF.CHR.map(str) + "-" + DF.POS.map(str) + "-" + DF.REF.map(str) + "_" + DF.ALT.map(str)

# reshape from wide to long
LONG=pd.melt(DF.iloc[:,5:], id_vars='IDEN', var_name="ID", value_name="GENO")

# remove non-variant sites
LONG=LONG[LONG.GENO.map(lambda x : str(x).split(":")[0]) != "0"]

LONG.to_csv(output, index=False, sep='\t')