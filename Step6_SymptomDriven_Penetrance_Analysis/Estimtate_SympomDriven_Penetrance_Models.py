import pandas as pd 
import numpy as np 
from tqdm import tqdm 
import pickle
import os
from SymptomSetModel.SymptomSetModel import SymptomaticDiseaseModel

#This script estimates the symptom-driven penetrance models using the SymptomSetModel software package. It also estimates the carrier-level expression probabilities. 


disease_dx_table = pd.read_pickle('Path/to/Auxillary_Data/CollapsedDxTable.pth')
covariate_table = pd.read_csv('/Path/to/Covariate/Data/Covariates.pth',sep='\t')
covariate_table.set_index('Subject_ID',inplace=True)
carrier_info_path='Path/to/CarrierDatasets/CarrierInfoFiles/'

with open('/Path/to/Clinical/Data/HPODx_SparseCSR.pth','rb') as f:
	symptom_dx_dataset=pickle.load(f)

symptom_table=pd.read_pickle('/Path/to/Auxillary_Data/Disease_to_Symptoms.pth')

lof_penetrance_table={}

symptom_model=SymptomaticDiseaseModel(symptom_dx_dataset['SubjectIndex'].index,symptom_dx_dataset['HPOColumns'].index,symptom_dx_dataset['SparseSymptomMatrix'])
for dis_abbrev in tqdm(disease_dx_table.index):
	with open(carrier_info_path+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		lof_info_table=pickle.load(f)['lof_table']
	if len(symptom_table.loc[dis_abbrev]['HPO'].intersection(symptom_dx_dataset['HPOColumns'].index))>0 and (len(lof_info_table.index.intersection(covariate_table.index))>0):

		symptom_model.BuildBackgroundModel(list(symptom_table.loc[dis_abbrev]['HPO']),'Independent',TrainingIndex=covariate_table.index)
		output = symptom_model.FitPenetranceModel(lof_info_table.index.intersection(covariate_table.index))
		lof_penetrance_table[dis_abbrev]=output[0].copy()
		symptom_model.SavePenetranceModel('/Path/to/Penetrance/Models/{0:s}_SymptomSetPosterior'.format(dis_abbrev))
		
		penetrance_posterior = pd.DataFrame({'Alpha':[output[1][0]],'Beta':[output[1][1]]})
		penetrance_posterior.to_csv('/Path/to/Penetrance/Models/{0:s}_PenetrancePosterior.txt'.format(dis_abbrev),sep='\t',index=False)
		penetrance_posterior.to_pickle('/Path/to/Penetrance/Models/{0:s}_PenetrancePosterior.pth'.format(dis_abbrev))

	else:
		lof_penetrance_table[dis_abbrev]=pd.Series([np.nan]*len(lof_info_table.index),index=lof_info_table.index)
with open('/Path/to/Penetrance/Models/LoFSymptomDrivenPenetrance.pth','wb') as f:
	pickle.dump(lof_penetrance_table,f)

