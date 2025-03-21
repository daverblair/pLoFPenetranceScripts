import numpy as np
import pandas as pd
from datetime import datetime
from tqdm import tqdm

#This script takes a tab-delmitted text file (AllSubjects_OMOP_FirstDx.txt) of all unique, patient specific OMOP diagnoses and converts this into a custom, sparse matrix reprentation. The text file has the following structure and can easily be created using biobank-specific SQL queries:
#Column 1: Subject ID
#Column 2: OMOP Concept ID
#Column 3: First date of concept ID diagnosis.
#The file must be sorted by subject id, and dates must have the following format: year-month-day.
#The output is a pandas data table with one entry per subject. There is one column, which contains all of the dates in which a new diagnosis was given for this subject. Multiple new diagnoses per date are possible. The covariate data table must have a column called 'BirthDate'. This sparse data structure is quite condensed can be loaded into the RAM of most modern laptops.


covariate_table = pd.read_csv('/Path/to/Covariate/Data/Covariates.txt',sep='\t')
covariate_table.set_index('Subject_ID',inplace=True)



dx_table=pd.DataFrame([],columns=['DxData'],index=covariate_table.index)
counter=0
with open("Path/to/OMOP/Data/AllSubjects_OMOP_FirstDx.txt",'r') as f:
	table_header=f.readline().strip('\n').split('\t')
	current_patient_data=[f.readline().strip('\n').split('\t')]
	pbar = tqdm(total=dx_table.shape[0])
	while len(current_patient_data[-1])>1:
		next_file_line=f.readline().strip('\n').split('\t')
		if next_file_line[0]==current_patient_data[0][0]:
			current_patient_data+=[next_file_line]
		else:
			person_id=current_patient_data[0][0]
			all_subject_dx = pd.DataFrame(current_patient_data,columns=table_header,index=list(range(len(current_patient_data))))
			all_subject_dx.replace('',pd.NA,inplace=True)
			if len(all_subject_dx)>0:
			    birth_date=covariate_table.loc[person_id,'BirthDate'].date()
			    all_subject_dx['date']=pd.Series([datetime.strptime(x, '%Y-%m-%d').date() for x in all_subject_dx['min(condition_start_datetime)']],index=all_subject_dx.index)
			    all_subject_dx['dx_age']=pd.Series([int(np.round((x-birth_date).days/365)) for x in all_subject_dx['date']],index=all_subject_dx.index,dtype=np.int32)
			    all_subject_dx['dx_year']=pd.Series([x.year for x in all_subject_dx['date']],index=all_subject_dx.index,dtype=np.int32)
			    all_subject_dx['dx_age_year']=pd.Series(zip(all_subject_dx['dx_age'],all_subject_dx['dx_year']),index=all_subject_dx.index)
			    for current_age_year in sorted(all_subject_dx['dx_age_year'].unique()):
			        current_dx=all_subject_dx.loc[all_subject_dx.dx_age_year==current_age_year]
			        all_dxs=set(current_dx['condition_concept_id'].unique())
			        try:
			            dx_table.loc[person_id,'DxData']+=[{'Age':current_age_year[0],'Year':current_age_year[1],'FirstDx':all_dxs}]
			        except TypeError:
			            dx_table.loc[person_id,'DxData']=[{'Age':current_age_year[0],'Year':current_age_year[1],'FirstDx':all_dxs}]


			pbar.update(1)
			current_patient_data=[next_file_line]


	pbar.close()

dx_table.to_pickle('/Path/to/Clinical/Data/OMOP_FirstDxOnly.pth')