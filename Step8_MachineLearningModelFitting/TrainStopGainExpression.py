import pandas as pd 
import numpy as np 
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import roc_curve,auc,precision_recall_curve,average_precision_score,f1_score
from sklearn.utils import shuffle
from sklearn.preprocessing import OneHotEncoder
import tqdm
from joblib import dump, load
import pickle

#Trains the machine learning penetrance prediction models for stop gain variants. The final output is a table of penetrance scores and a set of trained models.


num_splits=5

stop_gain_table=pd.read_pickle('../Results/StopGainExpressionData.pth')
stop_gain_table=stop_gain_table.loc[stop_gain_table.IN_HQ_SET==True]

#predictors
design_matrix = stop_gain_table[['CADD','FRAC_AA_IMPACTED','POSSIBLE_MET_RESCUE']]
loftee_one_hot = OneHotEncoder(sparse_output=False,categories=[['LC','HC']]).fit_transform(stop_gain_table['LOFTEE_CLASS'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(loftee_one_hot,columns = ['LC','HC'],index=stop_gain_table.index)],axis=1)
tx_annot_one_hot  = OneHotEncoder(sparse_output=False,categories=[['NONE','MANE_PLUS_CLINICAL','MANE_SELECT']]).fit_transform(stop_gain_table['TX_ANNOT'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(tx_annot_one_hot,columns = ['NONE','MANE_PLUS_CLINICAL','MANE_SELECT'],index=stop_gain_table.index)],axis=1)
nmd_one_hot=OneHotEncoder(sparse_output=False,categories=[['FALSE','LAST_CODING_EXON', 'LARGE_EXON','FIRST_EXON_LEQ_150NT_FROM_START', 'LEQ_50NT_FROM_LAST_EJ']]).fit_transform(stop_gain_table['PRED_NMD_ESCAPE'].values.reshape(-1,1))
design_matrix=pd.concat([design_matrix,pd.DataFrame(nmd_one_hot,columns = ['FALSE','LAST_CODING_EXON', 'LARGE_EXON','FIRST_EXON_LEQ_150NT_FROM_START', 'LEQ_50NT_FROM_LAST_EJ'],index=stop_gain_table.index)],axis=1)

#outcomes 
is_expressed = stop_gain_table.BINARIZED_EXPRESSION

logreg_output_table = pd.DataFrame([],columns=['RandomClassifierAvgPrecScore','ExpressionAvgPrecScore','ExpressionMaxF1','ExpressionMaxF1Thresh','ExpressionMaxF1Precision','ExpressionMaxF1Recall','AsymptomaticAUC','Asymptomatic5PercFPRThresh','Asymptomatic5PercFPR_TPR'],index=['CV_{0:d}'.format(i+1) for i in range(num_splits)])

rf_output_table = pd.DataFrame([],columns=['RandomClassifierAvgPrecScore','ExpressionAvgPrecScore','ExpressionMaxF1','ExpressionMaxF1Thresh','ExpressionMaxF1Precision','ExpressionMaxF1Recall','AsymptomaticAUC','Asymptomatic5PercFPRThresh','Asymptomatic5PercFPR_TPR'],index=['CV_{0:d}'.format(i+1) for i in range(num_splits)])

penetrance_score_table = pd.DataFrame([],columns=['RF_Score','LogReg_Score','RF_Expressed_Flag','LogReg_Expressed_Flag','RF_Aymptomatic_Flag','LogReg_Aymptomatic_Flag'],index=stop_gain_table.index)

np.random.seed(10231981)
kf = KFold(n_splits=num_splits,shuffle=True)
for i, (train_index, test_index) in tqdm.tqdm(enumerate(kf.split(design_matrix)),total=num_splits):
	rf_mod = RandomForestClassifier(min_samples_leaf=5,min_samples_split=10,criterion='log_loss',n_estimators=500)
	fitted_model=rf_mod.fit(design_matrix.values[train_index],is_expressed.values[train_index])
	lin_mod=LogisticRegression(C=1.0,penalty='l2',fit_intercept=False,max_iter=500).fit(design_matrix.values[train_index],is_expressed.values[train_index])

	lin_mod_pred=lin_mod.predict_proba(design_matrix.values[test_index])
	rf_mod_pred=rf_mod.predict_proba(design_matrix.values[test_index])


	prec_logreg,recall_logreg,thresh_logreg = precision_recall_curve(is_expressed.values[test_index],lin_mod_pred[:,1])
	prec_rf,recall_rf,thresh_rf = precision_recall_curve(is_expressed.values[test_index],rf_mod_pred[:,1])

	logreg_f1_scores = 2*recall_logreg*prec_logreg/(recall_logreg+prec_logreg)
	rf_f1_scores = 2*recall_rf*prec_rf/(recall_rf+prec_rf)


	identify_expression_avg_precision_logreg = average_precision_score(is_expressed.values[test_index],lin_mod_pred[:,1])
	identify_expression_avg_precision_rf= average_precision_score(is_expressed.values[test_index],rf_mod_pred[:,1])
	avg_precision_score_random = (is_expressed.values[test_index]).mean()


	max_f1_score_logreg=np.nanmax(logreg_f1_scores)
	max_f1_score_rf=np.nanmax(rf_f1_scores)

	prec_at_max_f1_logreg = prec_logreg[np.nanargmax(logreg_f1_scores)]
	prec_at_max_f1_rf = prec_rf[np.nanargmax(rf_f1_scores)]

	recall_at_max_f1_logreg = recall_logreg[np.nanargmax(logreg_f1_scores)]
	recall_at_max_f1_rf = recall_rf[np.nanargmax(rf_f1_scores)]

	score_max_f1_threshold_logreg = thresh_logreg[np.nanargmax(logreg_f1_scores)]
	score_max_f1_threshold_rf = thresh_rf[np.nanargmax(rf_f1_scores)]


	fpr_logreg,tpr_logreg,thresholds_logreg = roc_curve(1-is_expressed.values[test_index],lin_mod_pred[:,0])
	identify_non_expression_auc_logreg=auc(fpr_logreg,tpr_logreg)

	fpr_rf,tpr_rf,thresholds_rf = roc_curve(1-is_expressed.values[test_index],rf_mod_pred[:,0])
	identify_non_expression_auc_rf=auc(fpr_rf,tpr_rf)

	five_perc_fpr_thresh_logreg = thresholds_logreg[np.where(fpr_logreg<=0.05)[0][-1]]
	true_pos_rate_at_five_perc_logreg = tpr_logreg[np.where(fpr_logreg<=0.05)[0][-1]]

	five_perc_fpr_thresh_rf = thresholds_rf[np.where(fpr_rf<=0.05)[0][-1]]
	true_pos_rate_at_five_perc_rf= tpr_rf[np.where(fpr_rf<=0.05)[0][-1]]
	
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'RandomClassifierAvgPrecScore']=avg_precision_score_random
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionAvgPrecScore']=identify_expression_avg_precision_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1']=max_f1_score_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1Thresh']=score_max_f1_threshold_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1Precision']=prec_at_max_f1_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1Recall']=recall_at_max_f1_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'AsymptomaticAUC']=identify_non_expression_auc_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'Asymptomatic5PercFPRThresh']=five_perc_fpr_thresh_logreg
	logreg_output_table.loc['CV_{0:d}'.format(i+1),'Asymptomatic5PercFPR_TPR']=true_pos_rate_at_five_perc_logreg


	rf_output_table.loc['CV_{0:d}'.format(i+1),'RandomClassifierAvgPrecScore']=avg_precision_score_random
	rf_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionAvgPrecScore']=identify_expression_avg_precision_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1']=max_f1_score_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1Thresh']=score_max_f1_threshold_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1Precision']=prec_at_max_f1_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'ExpressionMaxF1Recall']=recall_at_max_f1_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'AsymptomaticAUC']=identify_non_expression_auc_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'Asymptomatic5PercFPRThresh']=five_perc_fpr_thresh_rf
	rf_output_table.loc['CV_{0:d}'.format(i+1),'Asymptomatic5PercFPR_TPR']=true_pos_rate_at_five_perc_rf

	penetrance_score_table.loc[design_matrix.index[test_index],'RF_Score']=rf_mod_pred[:,1]
	penetrance_score_table.loc[design_matrix.index[test_index],'LogReg_Score']=lin_mod_pred[:,1]

	penetrance_score_table.loc[design_matrix.index[test_index],'RF_Expressed_Flag']=rf_mod_pred[:,1]>=score_max_f1_threshold_rf
	penetrance_score_table.loc[design_matrix.index[test_index],'LogReg_Expressed_Flag']=lin_mod_pred[:,1]>=score_max_f1_threshold_logreg

	penetrance_score_table.loc[design_matrix.index[test_index],'RF_Aymptomatic_Flag']=rf_mod_pred[:,0]>=five_perc_fpr_thresh_rf
	penetrance_score_table.loc[design_matrix.index[test_index],'LogReg_Aymptomatic_Flag']=lin_mod_pred[:,0]>=five_perc_fpr_thresh_logreg

logreg_output_table.to_csv('/Path/to/ML/Model/Results/StopGainLogRegPerformanceTable.txt',sep='\t')
rf_output_table.to_csv('/Path/to/ML/Model/Results/StopGainRandomForestPerformanceTable.txt',sep='\t')

# global models 
rf_mod = RandomForestClassifier(min_samples_leaf=5,min_samples_split=10,criterion='log_loss',n_estimators=500)
final_rf_model=rf_mod.fit(design_matrix.values,is_expressed.values)
dump(final_rf_model, '/Path/to/ML/Model/Models/StopGainRandomForestModel.pth') 

final_lin_mod=LogisticRegression(C=1.0,penalty='l2',fit_intercept=False,max_iter=500).fit(design_matrix.values,is_expressed.values)
dump(final_lin_mod, '/Path/to/ML/Model/Models/StopGainLogisticRegressionModel.pth') 

penetrance_score_table.to_csv('/Path/to/ML/Model/Results/StopGainPenetranceScores.txt',sep='\t')
penetrance_score_table.to_pickle('/Path/to/ML/Model/Results/StopGainPenetranceScores.pth')

