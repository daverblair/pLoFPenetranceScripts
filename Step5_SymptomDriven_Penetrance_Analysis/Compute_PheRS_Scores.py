import pandas as pd 
import numpy as np 
from tqdm import tqdm 
import pickle
from SymptomSetModel.SymptomSetModel import SymptomaticDiseaseModel

hpo_alignmnet_fpr=20
annot_rate_cutoff=100


disease_dx_table = pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/CollapsedDxTable.pth')
covariate_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_Covariates.pth')
symptom_table_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/7_SymptomTables/UKBBSymptomTables/'
symptom_matrix_direc = '/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/HPODatasets/'

print('#'*20)
print('Building Outlier Score Matrix w/HPO FPR {0:d}, Annot Rate Threshold {1:d}'.format(hpo_fpr,annot_rate))
print('#'*20)

with open(symptom_matrix_direc+'SparseHPOMatrix_{0:d}.pth'.format(hpo_fpr),'rb') as f:
	symptom_dx_dataset=pickle.load(f)

with open(symptom_table_direc+'SymptomDatasets_{0:d}.pth'.format(annot_rate),'rb') as f:
	symptom_tables=pickle.load(f)

outlier_score_matrix=pd.DataFrame(np.zeros((covariate_table.shape[0],disease_dx_table.shape[0]),dtype=np.float64),index=covariate_table.index,columns=disease_dx_table.index)

symptom_model=SymptomaticDiseaseModel(symptom_dx_dataset['SubjectIndex'].index,symptom_dx_dataset['HPOColumns'].index,symptom_dx_dataset['SparseSymptomMatrix'])
for dis_abbrev in tqdm(disease_dx_table.index):
	if len(symptom_tables[dis_abbrev].index.intersection(symptom_dx_dataset['HPOColumns'].index))>0:
		symptom_model.BuildBackgroundModel(symptom_tables[dis_abbrev].index,'Independent',TrainingIndex=covariate_table.index)
		full_background_scores = symptom_model.ComputeBackgroundScores(covariate_table.index)
		outlier_score_matrix.loc[covariate_table.index,dis_abbrev]=full_background_scores.loc[covariate_table.index]

	else:
		outlier_score_matrix.loc[covariate_table.index,dis_abbrev]=np.nan
outlier_score_matrix.to_pickle('../OutlierMatrices/OutlierScoreMatrix_HPO_FPR-{0:d}_ANNOT_RATE-{1:d}.pth'.format(hpo_fpr,annot_rate))