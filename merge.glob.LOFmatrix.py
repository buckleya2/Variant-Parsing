'''
generic script to merge together files in a folder based on a common prefix
takes in a path, regex search string, and an output filname
'''

import pandas as pd
import glob
from sys import argv

script, output=argv

path="/nrnb/users/abuckley/somatic_vcf/LOF_matrix/"
search="RESHAPE.*"

merge_files=glob.glob(path + search)

DFlist=[]
for i in merge_files:
	tmp=pd.read_csv(i, sep='\t')
	DFlist.append(tmp)

#fin=reduce(lambda x, y: pd.merge(x, y ), DFlist)
fin=pd.concat(DFlist)
fin.fillna(0, inplace=True)

fin.to_csv(output, sep='\t', index=False)
