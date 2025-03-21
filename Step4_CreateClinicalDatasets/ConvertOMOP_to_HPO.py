import pandas as pd 
import tqdm
import numpy as np

#This script takes the output from 'BuildOMOPCodeMatrix.py' and converts the OMOP concept codes into HPO terms using the alignment constructed by this manuscript. The alignment is included in the Supplementary Materials and in the Auxillary_Data directory

FPR = 20

omop_map=pd.read_pickle('/Path/to/Auxillary_Data/OMOP_HPO_Map.pth')
omop_map.set_index('OMOP',inplace=True)
disease_dx_table=pd.read_pickle('/Path/to/Auxillary_Data/CollapsedDxTable.pth')
subject_dx_data=pd.read_pickle('/Path/to/Clinical/Data/OMOP_FirstDxOnly.pth')

has_dx_data=pd.isna(subject_dx_data['DxData'])==False

# doublt check to make sure symptoms are not included
disallowed_omop=set([str(x) for x in set().union(*disease_dx_table['OMOP'].values)])



omop_conversion_table={}
for idx,data in omop_map[omop_map['{0:d}_PERC_FPR'.format(FPR)]==True].iterrows():
	if str(idx) not in disallowed_omop:
		try:
			omop_conversion_table[str(idx)].update([data['HPO']])
		except KeyError:
			omop_conversion_table[str(idx)]=set([data['HPO']])
omop_conv_set=set(omop_conversion_table.keys())


hpo_table=pd.DataFrame([],index=subject_dx_data.index)
hpo_table['HPO']=pd.Series([set() for x in range(subject_dx_data.shape[0])],index=subject_dx_data.index)
hpo_table.loc[subject_dx_data.index[has_dx_data],'HPO']=subject_dx_data.loc[has_dx_data]['DxData'].apply(lambda x: set().union(*[omop_conversion_table[str(x)] for x in set().union(*[list(map(str,y['FirstDx'])) for y in x]).intersection(omop_conv_set)]))
hpo_table.to_pickle('/Path/to/Clinical/Data/HPODx.pth')

