import pandas as pd 
import numpy as np
import tqdm
import statsmodels.api as sm
from statsmodels.robust import mad
from scipy.stats import rankdata,norm
import tqdm
import pickle

def BrunnerMunzel_wVariance(pd_x,pd_y,alternative='two-sided'):
	no_nans = pd_x.index[pd.isna(pd_x)==False]
	x = pd_x.loc[no_nans].values

	no_nans=pd_y.index[pd.isna(pd_y)==False]
	y = pd_y.loc[no_nans].values


	nx = len(x)
	ny = len(y)

	rankc = rankdata(np.concatenate((x, y)))
	rankcx = rankc[0:nx]
	rankcy = rankc[nx:nx+ny]
	rankcx_mean = np.mean(rankcx)
	rankcy_mean = np.mean(rankcy)
	rankx = rankdata(x)
	ranky = rankdata(y)
	rankx_mean = np.mean(rankx)
	ranky_mean = np.mean(ranky)

	Sx = np.sum(np.power(rankcx - rankx - rankcx_mean + rankx_mean, 2.0))
	Sx /= nx - 1
	Sy = np.sum(np.power(rankcy - ranky - rankcy_mean + ranky_mean, 2.0))
	Sy /= ny - 1
	bm_stat = -1*nx * ny * (rankcy_mean - rankcx_mean)
	bm_se = (nx + ny) * np.sqrt(nx * Sx + ny * Sy)
	bm_z = bm_stat/bm_se

	if alternative=='two-sided':
		p_val = min(norm(0,1).sf(bm_z),norm(0,1).cdf(bm_z))*2.0
	elif alternative == "greater":
		p_val = norm(0,1).sf(bm_z)
	else:
		p_val = norm(0,1).cdf(bm_z)
	return {'BrunMunzStat':bm_stat,'BrunMunzStatSE':bm_se,'P-value':p_val}


hpo_omop_fpr=20



disease_table = pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/CollapsedDxTable.pth')
covariate_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_Covariates.pth')
coverage_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_CoverageTable.pth')
carrier_class_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CarrierInfoFiles/'
chip_confounded_diseases = ['BOPS']

target_columns = list(covariate_table.columns[1:])
covariates_fully_obs=covariate_table.index[pd.isna(covariate_table[target_columns]).sum(axis=1)==0]
covariate_table=covariate_table.loc[covariates_fully_obs][target_columns]
covariates_fully_obs=covariate_table.index[pd.isna(covariate_table[target_columns]).sum(axis=1)==0]
covariate_table=covariate_table.loc[covariates_fully_obs]
covariate_table=covariate_table.loc[coverage_table.index[pd.isna(coverage_table.CoverageYears)==False].intersection(covariate_table.index)]

local_design_matrix=covariate_table.copy().astype(float)
local_design_matrix=sm.add_constant(local_design_matrix,prepend=True)
score_matrix=pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/8_OutlierScoreEstimation/UKBB/OutlierMatrices/OutlierScoreMatrix_HPO_FPR-{0:d}_ANNOT_RATE-100.pth'.format(hpo_omop_fpr))


output_table = pd.DataFrame([],columns=['PheRS_Median_pLoF_Carriers','PheRS_MAD_pLoF_Carriers','PheRS_Median_Controls','PheRS_MAD_Controls','BrunMunz_Stat','BrunMunz_SE','BrunMunz_Pvalue'],index= disease_table.index.drop(chip_confounded_diseases))
for i,dis_abbrev in enumerate(tqdm.tqdm(disease_table.index.drop(chip_confounded_diseases))):
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		carrier_class_data=pickle.load(f)
	target_lof_carriers=carrier_class_data['lof_table'].index.intersection(covariate_table.index)
	if (len(target_lof_carriers)>0) and (pd.isna(score_matrix[dis_abbrev]).sum()!=score_matrix.shape[0]):
		target_controls=carrier_class_data['non_carriers']
		target_subjects = target_controls+list(target_lof_carriers)
		covariate_model_results = sm.OLS(score_matrix.loc[local_design_matrix.index][dis_abbrev].values,local_design_matrix).fit()
		score_residuals = covariate_model_results.resid
		lof_median,lof_mad = score_residuals.loc[target_lof_carriers].median(),mad(score_residuals.loc[target_lof_carriers])
		control_median,control_mad = score_residuals.loc[target_controls].median(),mad(score_residuals.loc[target_controls])
		brunmunz_stats=BrunnerMunzel_wVariance(score_residuals.loc[target_lof_carriers],score_residuals.loc[target_controls],alternative='greater')

		output_table.loc[dis_abbrev,'PheRS_Median_pLoF_Carriers']=lof_median
		output_table.loc[dis_abbrev,'PheRS_MAD_pLoF_Carriers']=lof_mad
		output_table.loc[dis_abbrev,'PheRS_Median_Controls']=control_median
		output_table.loc[dis_abbrev,'PheRS_MAD_Controls']=control_mad
		output_table.loc[dis_abbrev,'BrunMunz_Stat']=brunmunz_stats['BrunMunzStat']
		output_table.loc[dis_abbrev,'BrunMunz_SE']=brunmunz_stats['BrunMunzStatSE']
		output_table.loc[dis_abbrev,'BrunMunz_Pvalue']=brunmunz_stats['P-value']
output_table.to_pickle('../Results/BrunnerMunzelResultsTable.pth')
output_table.to_csv('../Results/BrunnerMunzelResultsTable.txt',sep='\t')
