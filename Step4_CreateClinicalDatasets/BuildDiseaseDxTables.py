import pandas as pd 
import numpy as np
import tqdm

#This script produces tables of disease diagnoses using OMOP diagnostic codes. One table is produced for each disease. It requires the following data files:
#1) The table of diseases and OMOP Dx codes (in Auxillary_Data)
#2) A pandas data table of unique OMOP Dx codes for the biobank subjects. This table is created by BuildOMOPCodeMatrix.py


disease_dx_table=pd.read_pickle('/Path/to/Auxillary_Data/CollapsedDxTable.pth')
subject_dx_data=pd.read_pickle('/Path/to/OMOP/Data/OMOP_FirstDxOnly.pth')

def _parse_codes(codes):
    try:
        return set().union(*[y['FirstDx'] for y in codes])
    except TypeError:
        return set()

has_emr_data=pd.isna(subject_dx_data['DxData'])==False
allowed_omop=set().union(*subject_dx_data.loc[has_emr_data]['DxData'].apply(lambda x: _parse_codes(x)))

def return_first_dx_age(dx_info,target_codes):
	not_found=True
	onset_age=pd.NA
	for dx_set in dx_info:
		if len(set(target_codes).intersection(dx_set['FirstDx']))>0:
			return dx_set['Age']
	if not_found==True:
		raise ValueError('{0:s} not found'.format(','.join(list(target_codes))))

	return onset_age



for dname in tqdm.tqdm(disease_dx_table.index):
	dx_table=pd.DataFrame([],index=subject_dx_data.index,columns=['HasDisDx','LeftObsInterval','RightObsInterval'])
	dx_codes=disease_dx_table.loc[dname]['OMOP']
	if len(dx_codes.intersection(allowed_omop))==0:
		print('OMOP codes for {0:s} are not found in the UKBB. Setting to NaN.'.format(dname))
		print('DX Codes: {0:s}'.format(','.join(list(dx_codes))))
	else:
		dx_codes=dx_codes.intersection(allowed_omop)
		dx_table.loc[has_emr_data,'HasDisDx']=subject_dx_data.loc[has_emr_data]['DxData'].apply(lambda x: 1 if len(dx_codes.intersection(set().union(*[y['FirstDx'] for y in x])))>0 else 0)
		dx_table.loc[dx_table.HasDisDx==1,'RightObsInterval']=subject_dx_data.loc[dx_table.HasDisDx==1]['DxData'].apply(lambda x: return_first_dx_age(x,dx_codes))
		dx_table.loc[dx_table.HasDisDx==0,'RightObsInterval']=subject_dx_data.loc[dx_table.HasDisDx==0]['DxData'].apply(lambda x: max([y['Age'] for y in x]))
		dx_table.loc[has_emr_data,'LeftObsInterval']=subject_dx_data.loc[has_emr_data]['DxData'].apply(lambda x: min([y['Age'] for y in x]))
	dx_table.to_pickle('/Path/to/Clinical/Datasets/DisDxDatasets/{0:s}_DxDataset.pth'.format(dname))
	dx_table.to_csv('/Path/to/Clinical/Datasets/DisDxDatasets/{0:s}_DxDataset.txt'.format(dname),sep='\t')