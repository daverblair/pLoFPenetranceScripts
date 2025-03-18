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

num_splits=5


splice_table=pd.read_pickle('../Results/SpliceChangeExpressionData.pth')
splice_table=splice_table.loc[splice_table.IN_HQ_SET==True]

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

is_expressed = splice_table.BINARIZED_EXPRESSION

logreg_output_table = pd.DataFrame([],columns=['RandomClassifierAvgPrecScore','ExpressionAvgPrecScore','ExpressionMaxF1','ExpressionMaxF1Thresh','ExpressionMaxF1Precision','ExpressionMaxF1Recall','AsymptomaticAUC','Asymptomatic5PercFPRThresh','Asymptomatic5PercFPR_TPR'],index=['CV_{0:d}'.format(i+1) for i in range(num_splits)])

rf_output_table = pd.DataFrame([],columns=['RandomClassifierAvgPrecScore','ExpressionAvgPrecScore','ExpressionMaxF1','ExpressionMaxF1Thresh','ExpressionMaxF1Precision','ExpressionMaxF1Recall','AsymptomaticAUC','Asymptomatic5PercFPRThresh','Asymptomatic5PercFPR_TPR'],index=['CV_{0:d}'.format(i+1) for i in range(num_splits)])

biologic_artifact_score_table = pd.DataFrame([],columns=['RF_Score','LogReg_Score','RF_Expressed_Flag','LogReg_Expressed_Flag','RF_Aymptomatic_Flag','LogReg_Aymptomatic_Flag'],index=splice_table.index)

np.random.seed(80720151)
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

	biologic_artifact_score_table.loc[design_matrix.index[test_index],'RF_Score']=rf_mod_pred[:,1]
	biologic_artifact_score_table.loc[design_matrix.index[test_index],'LogReg_Score']=lin_mod_pred[:,1]

	biologic_artifact_score_table.loc[design_matrix.index[test_index],'RF_Expressed_Flag']=rf_mod_pred[:,1]>=score_max_f1_threshold_rf
	biologic_artifact_score_table.loc[design_matrix.index[test_index],'LogReg_Expressed_Flag']=lin_mod_pred[:,1]>=score_max_f1_threshold_logreg

	biologic_artifact_score_table.loc[design_matrix.index[test_index],'RF_Aymptomatic_Flag']=rf_mod_pred[:,0]>=five_perc_fpr_thresh_rf
	biologic_artifact_score_table.loc[design_matrix.index[test_index],'LogReg_Aymptomatic_Flag']=lin_mod_pred[:,0]>=five_perc_fpr_thresh_logreg

logreg_output_table.to_csv('../Results/SpliceChangeLogRegPerformanceTable.txt',sep='\t')
rf_output_table.to_csv('../Results/SpliceChangeRandomForestPerformanceTable.txt',sep='\t')

# global models 
rf_mod = RandomForestClassifier(min_samples_leaf=5,min_samples_split=10,criterion='log_loss',n_estimators=500)
final_rf_model=rf_mod.fit(design_matrix.values,is_expressed.values)
dump(final_rf_model, '../Models/SpliceChangeRandomForestModel.pth') 

final_lin_mod=LogisticRegression(C=1.0,penalty='l2',fit_intercept=False).fit(design_matrix.values,is_expressed.values)
dump(final_lin_mod, '../Models/SpliceChangeLogisticRegressionModel.pth') 

biologic_artifact_score_table.to_csv('../Results/SpliceChangeArtifactScores.txt',sep='\t')
biologic_artifact_score_table.to_pickle('../Results/SpliceChangeArtifactScores.pth')

