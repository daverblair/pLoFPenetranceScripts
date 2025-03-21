import pandas as pd 
import numpy as np 
from tqdm import tqdm 
import pickle
from SymptomSetModel.SymptomSetModel import SymptomaticDiseaseModel

#This table produces one matrix of PheRS (outlier scores), with each column corresponding to a unique disease and each row a unique biobank partipant. It requires the multiple Auxillary_Datasets plus the Sparse CSR matrix of HPO diagnoses created by ConstructSparseHPOMatrices.py


disease_table = pd.read_pickle('/Path/to/Auxillary_Data/CollapsedDxTable.pth')
covariate_table = pd.read_csv('/Path/to/Covariate/Data/Covariates.pth',sep='\t')
covariate_table.set_index('Subject_ID',inplace=True)



with open('/Path/to/Clinical/Data/HPODx_SparseCSR.pth','rb') as f:
	symptom_dx_dataset=pickle.load(f)

symptom_table=pd.read_pickle('/Path/to/Auxillary_Data/Disease_to_Symptoms.pth')

outlier_score_matrix=pd.DataFrame(np.zeros((covariate_table.shape[0],disease_dx_table.shape[0]),dtype=np.float64),index=covariate_table.index,columns=disease_dx_table.index)

symptom_model=SymptomaticDiseaseModel(symptom_dx_dataset['SubjectIndex'].index,symptom_dx_dataset['HPOColumns'].index,symptom_dx_dataset['SparseSymptomMatrix'])
for dis_abbrev in tqdm(disease_dx_table.index):
	if len(symptom_table.loc[dis_abbrev]['HPO'].intersection(symptom_dx_dataset['HPOColumns'].index))>0:
		symptom_model.BuildBackgroundModel(list(symptom_table.loc[dis_abbrev]['HPO']),'Independent',TrainingIndex=covariate_table.index)
		full_background_scores = symptom_model.ComputeSurprisal(covariate_table.index)
		outlier_score_matrix.loc[covariate_table.index,dis_abbrev]=full_background_scores.loc[covariate_table.index]

	else:
		outlier_score_matrix.loc[covariate_table.index,dis_abbrev]=np.nan
outlier_score_matrix.to_pickle('Path/to/OutlierMatrices/OutlierScoreMatrix.pth')