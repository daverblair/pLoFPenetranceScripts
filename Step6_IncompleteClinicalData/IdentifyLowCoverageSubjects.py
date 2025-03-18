import pandas as pd 
from scipy import sparse
import numpy as np 
import pickle
import copy
import tqdm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import auc,roc_curve
from sklearn.utils import shuffle
from joblib import dump, load


combined_penetrance_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/13_BuildCombinedPenetranceDataset/UKBB/Results/GlobalPenetranceTable.pth')

ehr_coverage_design_matrix = pd.DataFrame([],columns=['FirstObsAge','RecruitmentAge','LastObsAge','NumTotalVisits','Onset'],index = combined_penetrance_table.index)
ehr_coverage_design_matrix['FirstObsAge'] =combined_penetrance_table.FIRST_OBS_AGE
ehr_coverage_design_matrix['RecruitmentAge']=combined_penetrance_table.RECRUITMENT_AGE
ehr_coverage_design_matrix['LastObsAge']=combined_penetrance_table.LAST_OBS_AGE
ehr_coverage_design_matrix['NumTotalVisits']=combined_penetrance_table.NUM_TOTAL_VISITS
ehr_coverage_design_matrix['Onset']=combined_penetrance_table.ONSET

ehr_coverage_design_matrix=ehr_coverage_design_matrix[(pd.isna(ehr_coverage_design_matrix)).sum(axis=1)==0]
no_disease_symptoms = 1-np.ceil(combined_penetrance_table[['PENETRANCE_PROB','HAS_DX']].max(axis=1))

num_splits=5
np.random.seed(102381807)
was_flagged = pd.DataFrame([],columns=['Subject_ID','DisId','Onset','Score','LowEHRDataFlag'],index=no_disease_symptoms.index)
for onset_type in ehr_coverage_design_matrix.Onset.unique():

	target_rows = ehr_coverage_design_matrix.index[ehr_coverage_design_matrix.Onset==onset_type]

	results =pd.DataFrame([],columns=['AUC','TPR_at_FivePercFPR','Threshold_at_FivePercFPR'],index=['CV_{0:d}'.format(i) for i in range(1,num_splits+1)])

	kf = KFold(n_splits=num_splits,shuffle=True)
	for i, (train_index, test_index) in tqdm.tqdm(enumerate(kf.split(target_rows)),total=num_splits):
		log_mod=LogisticRegression().fit(ehr_coverage_design_matrix[ehr_coverage_design_matrix.columns[:-1]].loc[target_rows[train_index]].values,no_disease_symptoms.loc[target_rows[train_index]].values)
		pred_no_symptoms = pd.Series(log_mod.predict_proba(ehr_coverage_design_matrix[ehr_coverage_design_matrix.columns[:-1]].loc[target_rows[test_index]].values)[:,1],index=target_rows[test_index])
		fpr,tpr,thresholds = roc_curve(no_disease_symptoms.loc[target_rows[test_index]],pred_no_symptoms)

		five_perc_fpr_thresh = thresholds[np.where(fpr<=0.05)[0][-1]]
		true_pos_rate_at_five_perc = tpr[np.where(fpr<=0.05)[0][-1]]
		auc_results=auc(fpr,tpr)
		results.loc['CV_{0:d}'.format(i+1),'AUC']=auc_results
		results.loc['CV_{0:d}'.format(i+1),'TPR_at_FivePercFPR']=true_pos_rate_at_five_perc
		results.loc['CV_{0:d}'.format(i+1),'Threshold_at_FivePercFPR']=five_perc_fpr_thresh

		was_flagged.loc[target_rows[test_index],'Subject_ID']=combined_penetrance_table.loc[target_rows[test_index]]['S_ID']
		was_flagged.loc[target_rows[test_index],'DisId']=combined_penetrance_table.loc[target_rows[test_index]]['DIS']
		was_flagged.loc[target_rows[test_index],'Onset']=combined_penetrance_table.loc[target_rows[test_index]]['ONSET']
		was_flagged.loc[target_rows[test_index],'LowEHRDataFlag']=(pred_no_symptoms>five_perc_fpr_thresh)
		was_flagged.loc[target_rows[test_index],'Score']=pred_no_symptoms
	full_log_mod=LogisticRegression().fit(ehr_coverage_design_matrix[ehr_coverage_design_matrix.columns[:-1]].loc[target_rows].values,no_disease_symptoms.loc[target_rows].values)
	dump(full_log_mod, '../Models/EHRCoverageModel_{0:s}.pth'.format(onset_type)) 

	results.to_csv('../Results/LowCoveragePerformanceTable_{0:s}.txt'.format(onset_type),sep='\t')
was_flagged.to_pickle('../Results/UKBB_LowEHRDataFlags.pth')
was_flagged.to_csv('../Results/UKBB_LowEHRDataFlags.txt',sep='\t')

