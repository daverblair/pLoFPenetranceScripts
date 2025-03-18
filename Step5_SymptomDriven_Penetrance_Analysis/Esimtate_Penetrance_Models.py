import pandas as pd 
import numpy as np 
from tqdm import tqdm 
import pickle
import os
from SymptomSetModel.SymptomSetModel import SymptomaticDiseaseModel

hpo_alignmnet_fprs=[1,5,10,20,25]
annot_rate_cutoffs=[100,50,25,10,5]


disease_dx_table = pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/CollapsedDxTable.pth')
covariate_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_Covariates.pth')
symptom_table_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/7_SymptomTables/UKBBSymptomTables/'
symptom_matrix_direc = '/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/HPODatasets/'
carrier_info_path='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CarrierInfoFiles/'

for annot_rate in annot_rate_cutoffs:
	if not os.path.exists('../PenetranceModels/AnnotRate_{0:d}'.format(annot_rate)):
	    os.makedirs('../PenetranceModels/AnnotRate_{0:d}'.format(annot_rate))

	for hpo_fpr in hpo_alignmnet_fprs:

		if not os.path.exists('../PenetranceModels/AnnotRate_{0:d}/HPO_FPR_{1:d}'.format(annot_rate,hpo_fpr)):
		    os.makedirs('../PenetranceModels/AnnotRate_{0:d}/HPO_FPR_{1:d}'.format(annot_rate,hpo_fpr))

		print('#'*20)
		print('Building Outlier Score Matrix w/HPO FPR {0:d}, Annot Rate Threshold {1:d}'.format(hpo_fpr,annot_rate))
		print('#'*20)

		with open(symptom_matrix_direc+'SparseHPOMatrix_{0:d}.pth'.format(hpo_fpr),'rb') as f:
			symptom_dx_dataset=pickle.load(f)

		with open(symptom_table_direc+'SymptomDatasets_{0:d}.pth'.format(annot_rate),'rb') as f:
			symptom_tables=pickle.load(f)

		lof_penetrance_table={}

		symptom_model=SymptomaticDiseaseModel(symptom_dx_dataset['SubjectIndex'].index,symptom_dx_dataset['HPOColumns'].index,symptom_dx_dataset['SparseSymptomMatrix'])
		for dis_abbrev in tqdm(disease_dx_table.index):
			with open(carrier_info_path+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
				lof_info_table=pickle.load(f)['lof_table']
			if (len(symptom_tables[dis_abbrev].index.intersection(symptom_dx_dataset['HPOColumns'].index))>0) and (len(lof_info_table.index.intersection(covariate_table.index))>0):

				symptom_model.BuildBackgroundModel(symptom_tables[dis_abbrev].index,'Independent',TrainingIndex=covariate_table.index)
				output = symptom_model.FitPenetranceModel(lof_info_table.index.intersection(covariate_table.index))
				lof_penetrance_table[dis_abbrev]=output[0].copy()
				symptom_model.SavePenetranceModel('../PenetranceModels/AnnotRate_{0:d}/HPO_FPR_{1:d}/{2:s}_SymptomSetPosterior'.format(annot_rate,hpo_fpr,dis_abbrev))
				
				penetrance_posterior = pd.DataFrame({'Alpha':[output[1][0]],'Beta':[output[1][1]]})
				penetrance_posterior.to_csv('../PenetranceModels/AnnotRate_{0:d}/HPO_FPR_{1:d}/{2:s}_PenetrancePosterior.txt'.format(annot_rate,hpo_fpr,dis_abbrev),sep='\t',index=False)
				penetrance_posterior.to_pickle('../PenetranceModels/AnnotRate_{0:d}/HPO_FPR_{1:d}/{2:s}_PenetrancePosterior.pth'.format(annot_rate,hpo_fpr,dis_abbrev))

			else:
				lof_penetrance_table[dis_abbrev]=pd.Series([np.nan]*len(lof_info_table.index),index=lof_info_table.index)
		with open('../LoFPenetranceTables/LoFPenetrance_HPO_FPR-{0:d}_ANNOT_RATE-{1:d}.pth'.format(hpo_fpr,annot_rate),'wb') as f:
			pickle.dump(lof_penetrance_table,f)

