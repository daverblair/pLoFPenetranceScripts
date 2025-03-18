import pandas as pd 
import numpy as np
import tqdm
from scipy.stats import fisher_exact
from pyfirth.PyFirth import PyFirth
import tqdm
import pickle



disease_table = pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/CollapsedDxTable.pth')
covariate_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_Covariates.pth')
carrier_class_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CarrierInfoFiles/'
dx_data_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/DisDxDatasets/'
chip_confounded_diseases = ['BOPS']


target_columns = list(covariate_table.columns[1:])
covariates_fully_obs=covariate_table.index[pd.isna(covariate_table[target_columns]).sum(axis=1)==0]
covariate_table=covariate_table.loc[covariates_fully_obs][target_columns]
covariates_fully_obs=covariate_table.index[pd.isna(covariate_table[target_columns]).sum(axis=1)==0]
covariate_table=covariate_table.loc[covariates_fully_obs]

def CorrectedLogOdds(contingency_table):
	contingency_table=np.array(contingency_table,dtype=np.float64)+0.5
	odds_ratio=contingency_table[0,0]*contingency_table[1,1]/(contingency_table[0,1]*contingency_table[1,0])
	return np.log(odds_ratio)

output_table={'DisAbbrev':[],'ContingencyTable':[],'FisherLogOddsRatio':[],'FisherPVal':[],'LogisticRegressionResults':[]}
for dis_abbrev in tqdm.tqdm(disease_table.index.drop(chip_confounded_diseases)):
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		carrier_class_data=pickle.load(f)
	target_lof_carriers=carrier_class_data['lof_table'].index.intersection(covariate_table.index)
	lof_carrier_status = pd.Series(np.zeros(covariate_table.shape[0]),index=covariate_table.index)
	lof_carrier_status.loc[carrier_class_data['lof_table'].index.intersection(covariate_table.index)]=1
	dx_data=pd.read_pickle(dx_data_direc+'{0:s}_DxDataset.pth'.format(dis_abbrev))
	output_table['DisAbbrev']+=[dis_abbrev]

	if (len(target_lof_carriers)>0) and (dx_data['HasDisDx'].sum()>0):
		target_controls=carrier_class_data['non_carriers']

		num_dx_lof_carriers=dx_data.loc[target_lof_carriers]['HasDisDx'].sum()
		num_undx_lof_carriers=(dx_data.loc[target_lof_carriers]['HasDisDx']==0).sum()

		num_dx_lof_controls=dx_data.loc[target_controls]['HasDisDx'].sum()
		num_undx_lof_controls=(dx_data.loc[target_controls]['HasDisDx']==0).sum()


		contingency_table=np.array([[num_undx_lof_controls,num_dx_lof_controls],[num_undx_lof_carriers,num_dx_lof_carriers]])
		fisher_results=fisher_exact(contingency_table,alternative='greater')
		output_table['ContingencyTable']+=[contingency_table]
		output_table['FisherLogOddsRatio']+=[CorrectedLogOdds(contingency_table)]
		output_table['FisherPVal']+=[fisher_results.pvalue]


		local_design_matrix=covariate_table.copy().astype(float)
		local_design_matrix['HasLOF']=lof_carrier_status.astype(float)
		local_design_matrix['HasDx']=dx_data.loc[local_design_matrix.index]['HasDisDx'].astype(float)
		target_subjects = target_controls+list(target_lof_carriers)
		local_design_matrix=local_design_matrix.loc[dx_data.index[pd.isna(dx_data['HasDisDx'])==False].intersection(target_subjects)]

		if local_design_matrix.HasDx.sum()>10:
			firth_model=PyFirth(local_design_matrix,['Recruitment_Age','Birth_Sex']+['PC{0:d}'.format(i) for i in range(1,17)]+['HasLOF'],'HasDx',hasconst=False)
			lof_output=firth_model.fit('HasLOF',convergence_limit=1e-6,step_limit=100)
			output_table['LogisticRegressionResults']+=[[lof_output['ParamTable'].loc['HasLOF']['BETA'],lof_output['ParamTable'].loc['HasLOF']['SE'],lof_output['PVal']]]
		else:
			output_table['LogisticRegressionResults']+=[[np.nan,np.nan,np.nan]]
	else:
		output_table['ContingencyTable']+=[np.array([[np.nan,np.nan],[np.nan,np.nan]])]
		output_table['FisherLogOddsRatio']+=[np.nan]
		output_table['FisherPVal']+=[np.nan]
		output_table['LogisticRegressionResults']+=[[np.nan,np.nan,np.nan]]	

output_table=pd.DataFrame(output_table)
output_table.set_index('DisAbbrev',inplace=True)
output_table.to_pickle('../Results/DxCodeAssociationAnalysis.pth')
output_table_txt=output_table.copy()
output_table_txt['ContingencyTable']=output_table_txt['ContingencyTable'].apply(lambda x: {'CTRL_NoDX':x[0,0],'CTRL_wDX':x[0,1],'CARR_NoDX':x[1,0],'CARR_wDX':x[1,1]})
output_table_txt['LogisticRegressionResults']=output_table_txt['LogisticRegressionResults'].apply(lambda x: {'Beta':x[0],'SE':x[1],'PVal':x[2]})
output_table_txt.to_csv('../Results/DxCodeAssociationAnalysis.txt',sep='\t')



