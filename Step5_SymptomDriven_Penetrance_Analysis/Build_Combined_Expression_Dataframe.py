import pandas as pd
import pickle
import tqdm
from sklearn.utils import resample
from sklearn.metrics import precision_recall_curve,average_precision_score,confusion_matrix
from statsmodels.stats.contingency_tables import mcnemar
import numpy as np

def _stat_perf_better_than_random(observations,predictions,n_resamples=10000,alternative='greater'):
	output_table={'BaselineAvgPrec':pd.NA,'ScoreAvgPrec':pd.NA,'PValue':pd.NA}
	baseline = observations.mean()
	score_avg_prec= average_precision_score(observations,predictions)
	output_table['BaselineAvgPrec']=baseline
	output_table['ScoreAvgPrec']=score_avg_prec

	obs_diff = score_avg_prec-baseline
	num_at_least_as_extreme=0
	for i in range(n_resamples):
		new_obs = resample(observations,replace=False)
		new_baseline = new_obs.mean()
		new_score_avg_prec =  average_precision_score(new_obs,predictions)
		if alternative=='greater':
			if (new_score_avg_prec-new_baseline)>=obs_diff:
				num_at_least_as_extreme+=1
		elif alternative=='less':
			if (new_score_avg_prec-new_baseline)<=obs_diff:
				num_at_least_as_extreme+=1
		else:
			if np.abs(new_score_avg_prec-new_baseline)>=np.abs(obs_diff):
				num_at_least_as_extreme+=1
	pval=num_at_least_as_extreme/n_resamples
	output_table['PValue']=pval
	return output_table

HPO_FPR=20
ANNOT_RATE=100
chip_confounded_diseases = ['BOPS']
dis_annots=pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/HaploinsuffientDiseasesCollapsed_Annotations.pth')
coverage_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_CoverageTable.pth')
covariate_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_Covariates.pth')
carrier_class_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CarrierInfoFiles/'


global_penetrance_table={'S_ID':[],'DIS':[],'PASSES_ANNOT_FILTER':[],'PENETRANCE_PROB':[],'HAS_DX':[]}
with open('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/10_PenetranceEstimation/UKBB/LoFPenetranceTables/LoFPenetrance_HPO_FPR-{0:d}_ANNOT_RATE-{1:d}.pth'.format(HPO_FPR,ANNOT_RATE),'rb') as f:
	penetrance_data = pickle.load(f)
for dis in chip_confounded_diseases:
	 del penetrance_data[dis]

for key in tqdm.tqdm(penetrance_data.keys()):
	prob_pen_table = penetrance_data[key]
	dx_data = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/DisDxDatasets/{0:s}_DxDataset.pth'.format(key))
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(key),'rb') as f:
		carrier_class_data=pickle.load(f)
	simple_filter = carrier_class_data['lof_table'].index[(carrier_class_data['lof_table'].LOFTEE_CLASS=='HC')*(carrier_class_data['lof_table'].TX_ANNOT=='MANE_SELECT')]
	global_penetrance_table['S_ID']+=list(prob_pen_table.index)
	global_penetrance_table['DIS']+=[key]*len(prob_pen_table)

	filtered = pd.Series(np.zeros(prob_pen_table.shape[0],dtype=int),prob_pen_table.index)
	filtered[filtered.index.intersection(simple_filter)]=1
	global_penetrance_table['PASSES_ANNOT_FILTER']+=list(filtered.values)
	global_penetrance_table['PENETRANCE_PROB']+=list(prob_pen_table.values)
	global_penetrance_table['HAS_DX']+=list(dx_data.loc[prob_pen_table.index,'HasDisDx'])

global_penetrance_table=pd.DataFrame(global_penetrance_table)
global_penetrance_table['ONSET']=global_penetrance_table.DIS.apply(lambda x:dis_annots.loc[x]['EarliestOnset'])
global_penetrance_table['DIS_TYPE']=global_penetrance_table.DIS.apply(lambda x:dis_annots.loc[x]['DiseaseType'])

global_penetrance_table['FIRST_OBS_AGE'] = pd.Series((coverage_table.loc[global_penetrance_table.S_ID]['StartDate']-coverage_table.loc[global_penetrance_table.S_ID]['BirthDate']).apply(lambda x:x.days/365).values,index=global_penetrance_table.index)
global_penetrance_table['RECRUITMENT_AGE']=pd.Series(covariate_table.loc[global_penetrance_table.S_ID]['Recruitment_Age'].values,index=global_penetrance_table.index)
global_penetrance_table['LAST_OBS_AGE']=pd.Series((coverage_table.loc[global_penetrance_table.S_ID]['EndDate']-coverage_table.loc[global_penetrance_table.S_ID]['BirthDate']).apply(lambda x:x.days/365).values,index=global_penetrance_table.index)
global_penetrance_table['NUM_TOTAL_VISITS']=pd.Series(coverage_table.loc[global_penetrance_table.S_ID]['NumTotalVisits'].values,index=global_penetrance_table.index)


global_penetrance_table_w_dx = global_penetrance_table.loc[pd.isna(global_penetrance_table.HAS_DX)==False]


precision_scores,recall_scores,thresholds = precision_recall_curve(global_penetrance_table_w_dx.HAS_DX,global_penetrance_table_w_dx.PENETRANCE_PROB)
f1_scores = 2*recall_scores*precision_scores/(recall_scores+precision_scores)
max_f1_score = np.nanmax(f1_scores)
prob_at_max_f1 = thresholds[np.nanargmax(f1_scores)]
symptom_model_pvalue = _stat_perf_better_than_random(global_penetrance_table_w_dx.HAS_DX,global_penetrance_table_w_dx.PENETRANCE_PROB)

f1_thresh_dict = {'F1_Value':max_f1_score,'ExpressionProbF1Thresh':prob_at_max_f1,'PrecisionF1Thres':precision_scores}
with open('../Results/SymptomSetPenetranceThreshold.pth','wb') as f:
	pickle.dump(f1_thresh_dict,f)
avg_precision_symp_probs = average_precision_score(global_penetrance_table_w_dx.HAS_DX,global_penetrance_table_w_dx.PENETRANCE_PROB)

global_penetrance_table['BINARIZED_PENETRANCE'] = pd.Series((global_penetrance_table.PENETRANCE_PROB>=prob_at_max_f1)+((pd.isna(global_penetrance_table.HAS_DX)==False)*(global_penetrance_table.HAS_DX==1)),dtype=np.float64)
global_penetrance_table['MAX_PENETRANCE'] = pd.Series(np.maximum(global_penetrance_table.PENETRANCE_PROB.values,np.array((pd.isna(global_penetrance_table.HAS_DX)==False)*(global_penetrance_table.HAS_DX==1).values,dtype=float)),dtype=np.float64,index=global_penetrance_table.index)
global_penetrance_table.to_pickle('../Results/GlobalPenetranceTable.pth')
global_penetrance_table.to_csv('../Results/GlobalPenetranceTable.txt',sep='\t')

symp_assigned_dx=np.array(global_penetrance_table_w_dx.PENETRANCE_PROB>=prob_at_max_f1,dtype=np.float64)
symp_assignd_confusion_mat = confusion_matrix(['HasDis' if x==1 else 'NoDis' for x in global_penetrance_table_w_dx.HAS_DX],['HasDis' if x==1 else 'NoDis' for x in symp_assigned_dx],labels = ['NoDis','HasDis'],normalize='all')
confusion_mat_stats = mcnemar(confusion_matrix(['HasDis' if x==1 else 'NoDis' for x in global_penetrance_table_w_dx.HAS_DX],['HasDis' if x==1 else 'NoDis' for x in symp_assigned_dx],labels = ['NoDis','HasDis']))


with open('../Results/Dx_Symptom_SummaryStats.txt','w') as f:
	f.write('Total # of Disease-LoF Pairs: {0:d}\n'.format(global_penetrance_table.shape[0]))
	f.write('Total # of Disease-LoF Pairs w/Dx Codes: {0:d}\n'.format(global_penetrance_table_w_dx.shape[0]))
	f.write('Total # of Diseases w/Dx Codes: {0:d}\n'.format(global_penetrance_table_w_dx.DIS.unique().shape[0]))
	f.write('Total # of Diseases w/o Dx Codes: {0:d}\n'.format(len(set(global_penetrance_table.DIS).difference(global_penetrance_table_w_dx.DIS.unique()))))
	f.write('\n')
	f.write('Fraction Penetrant per Dx: {0:.4f}\n'.format(global_penetrance_table_w_dx.HAS_DX.sum()/global_penetrance_table_w_dx.shape[0]))
	f.write('Fraction Penetrant per Symptoms (Dx Available Only): {0:.4f}\n'.format(global_penetrance_table_w_dx.PENETRANCE_PROB.sum()/global_penetrance_table_w_dx.shape[0]))
	f.write('Fraction Penetrant per Symptoms (All Diseases): {0:.4f}\n'.format(global_penetrance_table.PENETRANCE_PROB.sum()/global_penetrance_table.shape[0]))
	f.write('\n')
	f.write('Symptom Probs Average Precision Score for Identifying Dx: {0:.4f}\n'.format(avg_precision_symp_probs))
	f.write('Symptom Probs Average Precision Score P-value (Randomization Test; N=10000): {0:.2g}\n'.format(symptom_model_pvalue['PValue']))
	f.write('Symptom Probs Max F1 for Identifying Dx: {0:.4f}\n'.format(max_f1_score))
	f.write('Symptom Probs Threshold at Max F1: {0:.4f}\n'.format(prob_at_max_f1))
	f.write('Confusion Matrix at Max F1:\n')
	f.write('\tNoDisDx\tHasDisDx\n')
	f.write('No Symptom Dx\t{0:.4f}\t{1:.4f}\n'.format(symp_assignd_confusion_mat[0,0],symp_assignd_confusion_mat[0,1]))
	f.write('Has Symptom Dx\t{0:.4f}\t{1:.4f}\n'.format(symp_assignd_confusion_mat[1,0],symp_assignd_confusion_mat[1,1]))
	f.write('\n')
	f.write('Symptom-based vs Diagnosis-based Performance at Max F1 (McNemars Test): ({0:.4f}, {1:.4g})\n'.format(confusion_mat_stats.statistic, confusion_mat_stats.pvalue))



