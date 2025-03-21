#!/usr/bin/env python
# coding: utf-8

import pandas as pd 
import numpy as np 
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import roc_curve,auc,precision_recall_curve,average_precision_score,f1_score
from sklearn.utils import shuffle,resample
from sklearn.preprocessing import OneHotEncoder
import tqdm
from joblib import dump, load
import pickle

#This script applies the splice change model to the validation dataset. 


splice_table=pd.read_pickle('/Path/to/ML/Model/Data/SpliceChangeExpressionData_ValidationSet.pth')



design_matrix = splice_table[['CADD','SPLICE_AI_SCORE','FRAC_AA_IMPACTED']]
loftee_splice_type = OneHotEncoder(sparse_output=False,categories=[['DONOR_LOSS','DONOR_GAIN','ACC_LOSS','ACC_GAIN','INDETERMINATE']]).fit_transform(splice_table['SPLICE_TYPE'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(loftee_splice_type,columns = ['DONOR_LOSS','DONOR_GAIN','ACC_LOSS','ACC_GAIN','INDETERMINATE'],index=splice_table.index)],axis=1)
tx_annot_one_hot  = OneHotEncoder(sparse_output=False,categories=[['NONE','MANE_PLUS_CLINICAL','MANE_SELECT']]).fit_transform(splice_table['TX_ANNOT'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(tx_annot_one_hot,columns = ['NONE','MANE_PLUS_CLINICAL','MANE_SELECT'],index=splice_table.index)],axis=1)
loftee_one_hot = OneHotEncoder(sparse_output=False,categories=[['LC','HC']]).fit_transform(splice_table['LOFTEE_CLASS'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(loftee_one_hot,columns = ['LC','HC'],index=splice_table.index)],axis=1)
design_matrix['OUTSIDE_CODING_REGION']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: 1 if 'OUTSIDE_CODING_REGION' in set([y[0] for y in x]) else 0)

def _return_persistent_orig_score(splice_rescue_list):
	for rescue_flag in splice_rescue_list:
		if rescue_flag[0]=='ORIG_SPLICESITE_INTACT':
			return rescue_flag[1]
	return 0.0
design_matrix['PERSISTENT_ORIG_SPLICESITE_SCORE']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: _return_persistent_orig_score(x))
design_matrix['LAST_CODING_EXON']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: 1 if 'LAST_CODING_EXON' in set([y[0] for y in x]) else 0)
design_matrix['INFRAME_EXON_SKIP']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: 1 if 'INFRAME_EXON_SKIP' in set([y[0] for y in x]) else 0)
design_matrix['POSSIBLE_MET_RESCUE']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: 1 if 'POSSIBLE_MET_RESCUE' in set([y[0] for y in x]) else 0)
design_matrix['INFRAME_INTRON_RETENTION']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: 1 if 'INFRAME_INTRON_RETENTION' in set([y[0] for y in x]) else 0)
def _return_cryptic_rescue_score(splice_rescue_list):
	for rescue_flag in splice_rescue_list:
		if rescue_flag[0]=='CRYPTIC_RESCUE':
			return rescue_flag[1]
	return 0.0
design_matrix['CRYPTIC_RESCUE_SCORE']=splice_table.SPLICE_RESCUE_FLAGS.apply(lambda x: _return_cryptic_rescue_score(x))


design_matrix=design_matrix.loc[pd.isna(design_matrix).sum(axis=1)==0]


rf_mod = load('/Path/to/ML/Model/Models/SpliceChangeRandomForestModel.pth')
lin_mod = load('/Path/to/ML/Model/Models/SpliceChangeLogisticRegressionModel.pth')


prob_symptomatic_rf = rf_mod.predict_proba(design_matrix.values)
prob_symptomatic_logreg=lin_mod.predict_proba(design_matrix.values)


penetrance_score_table = pd.DataFrame([],columns=['RF_Score','LogReg_Score','RF_Expressed_Flag','LogReg_Expressed_Flag','RF_Asymptomatic_Flag','LogReg_Asymptomatic_Flag'],index=design_matrix.index)
penetrance_score_table.loc[design_matrix.index,'RF_Score']=prob_symptomatic_rf[:,1]
penetrance_score_table.loc[design_matrix.index,'LogReg_Score']=prob_symptomatic_logreg[:,1]


penetrance_score_table.to_csv('/Path/to/ML/Model/Results/SpliceChangePenetrance_ValidationSubset.txt',sep='\t')
penetrance_score_table.to_pickle('/Path/to/ML/Model/Results/SpliceChangePenetrance_ValidationSubset.pth')

