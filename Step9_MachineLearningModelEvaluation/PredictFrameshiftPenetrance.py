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


#This script predicts the frameshift penetrance scores in a new dataset.


frameshift_table=pd.read_pickle('/Path/to/ML/Model/Data/FrameshiftExpressionData_ValidationSet.pth')

design_matrix = frameshift_table[['CADD','IN_LAST_CODING_EXON','FRAC_AA_IMPACTED','POSSIBLE_MET_RESCUE']]
loftee_one_hot = OneHotEncoder(sparse_output=False,categories=[['LC','HC']]).fit_transform(frameshift_table['LOFTEE_CLASS'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(loftee_one_hot,columns = ['LC','HC'],index=frameshift_table.index)],axis=1)
tx_annot_one_hot  = OneHotEncoder(sparse_output=False,categories=[['NONE','MANE_PLUS_CLINICAL','MANE_SELECT']]).fit_transform(frameshift_table['TX_ANNOT'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(tx_annot_one_hot,columns = ['NONE','MANE_PLUS_CLINICAL','MANE_SELECT'],index=frameshift_table.index)],axis=1)


design_matrix=design_matrix.loc[pd.isna(design_matrix).sum(axis=1)==0]


rf_mod = load('Path/to/ML/Model/Models/FrameshiftRandomForestModel.pth')
lin_mod = load('Path/to/ML/Model/Models/FrameshiftLogisticRegressionModel.pth')


prob_symptomatic_rf = rf_mod.predict_proba(design_matrix.values)
prob_symptomatic_logreg=lin_mod.predict_proba(design_matrix.values)



penetrance_score_table = pd.DataFrame([],columns=['RF_Score','LogReg_Score','RF_Expressed_Flag','LogReg_Expressed_Flag','RF_Asymptomatic_Flag','LogReg_Asymptomatic_Flag'],index=design_matrix.index)
penetrance_score_table.loc[design_matrix.index,'RF_Score']=prob_symptomatic_rf[:,1]
penetrance_score_table.loc[design_matrix.index,'LogReg_Score']=prob_symptomatic_logreg[:,1]

penetrance_score_table.loc[design_matrix.index,'RF_Expressed_Flag']=prob_symptomatic_rf[:,1]>=score_max_f1_threshold_rf
penetrance_score_table.loc[design_matrix.index,'LogReg_Expressed_Flag']=prob_symptomatic_logreg[:,1]>=score_max_f1_threshold_logreg

penetrance_score_table.loc[design_matrix.index,'RF_Aymptomatic_Flag']=prob_symptomatic_rf[:,0]>=five_perc_fpr_thresh_rf
penetrance_score_table.loc[design_matrix.index,'LogReg_Aymptomatic_Flag']=prob_symptomatic_logreg[:,0]>=five_perc_fpr_thresh_logreg
penetrance_score_table.to_csv('/Path/to/ML/Model/Results/FrameshiftPenetranceScores_ValidationSet.txt',sep='\t')
penetrance_score_table.to_pickle('/Path/to/ML/Model/Results/FrameshiftPenetranceScores_ValidationSet.pth')

