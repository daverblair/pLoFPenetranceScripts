import pandas as pd 
import numpy as np
import tqdm
from scipy.stats import beta
import tqdm
import pickle

PSEUDO_COUNT=0.5
disease_table = pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/CollapsedDxTable.pth')
covariate_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CovariateDatasets/UKBB_Covariates.pth')
carrier_class_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/CarrierInfoFiles/'
dx_data_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/DisDxDatasets/'
chip_confounded_diseases = ['BOPS']

dx_data_direc='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/5_Clinical_Datasets/UKBB/Datasets/DisDxDatasets/'

output_table={'DisAbbrev':[],'CarrierCount':[],'FilteredCarrierCount':[],'PPV':[],'FilteredVariant_PPV':[]}
for dis_abbrev in tqdm.tqdm(disease_table.index.drop(chip_confounded_diseases)):
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		carrier_class_data=pickle.load(f)
	target_lof_carriers=carrier_class_data['lof_table'].index.intersection(covariate_table.index)
	simple_filter = carrier_class_data['lof_table'].index[(carrier_class_data['lof_table'].LOFTEE_CLASS=='HC')*(carrier_class_data['lof_table'].TX_ANNOT=='MANE_SELECT')]
	target_lof_carriers_filtered=simple_filter.intersection(covariate_table.index)

	dx_data=pd.read_pickle(dx_data_direc+'{0:s}_DxDataset.pth'.format(dis_abbrev))
	output_table['DisAbbrev']+=[dis_abbrev]
	if (len(target_lof_carriers)>0) and (dx_data['HasDisDx'].sum()>0):

		num_dx_lof_carriers=dx_data.loc[target_lof_carriers]['HasDisDx'].sum()
		num_undx_lof_carriers=(dx_data.loc[target_lof_carriers]['HasDisDx']==0).sum()

		penetrance=beta(num_dx_lof_carriers+PSEUDO_COUNT,num_undx_lof_carriers+PSEUDO_COUNT).mean()
		penetrance_95=beta(num_dx_lof_carriers+PSEUDO_COUNT,num_undx_lof_carriers+PSEUDO_COUNT).interval(confidence=0.95)

		output_table['PPV']+=[{'Mean':penetrance,'95_CI':penetrance_95}]
		output_table['CarrierCount']+=[num_dx_lof_carriers+num_undx_lof_carriers]

		num_dx_lof_carriers_filtered=dx_data.loc[target_lof_carriers_filtered]['HasDisDx'].sum()
		num_undx_lof_carriers_filtered=(dx_data.loc[target_lof_carriers_filtered]['HasDisDx']==0).sum()

		penetrance_filtered=beta(num_dx_lof_carriers_filtered+PSEUDO_COUNT,num_undx_lof_carriers_filtered+PSEUDO_COUNT).mean()
		penetrance_95_filtered=beta(num_dx_lof_carriers_filtered+PSEUDO_COUNT,num_undx_lof_carriers_filtered+PSEUDO_COUNT).interval(confidence=0.95)

		output_table['FilteredVariant_PPV']+=[{'Mean':penetrance_filtered,'95_CI':penetrance_95_filtered}]
		output_table['FilteredCarrierCount']+=[num_dx_lof_carriers_filtered+num_undx_lof_carriers_filtered]

	else:
		output_table['PPV']+=[{'Mean':np.nan,'95_CI':(np.nan,np.nan)}]
		output_table['FilteredVariant_PPV']+=[{'Mean':np.nan,'95_CI':(np.nan,np.nan)}]
		output_table['CarrierCount']+=[np.nan]
		output_table['FilteredCarrierCount']+=[np.nan]

output_table=pd.DataFrame(output_table)
output_table.set_index('DisAbbrev',inplace=True)
output_table.to_pickle('../Results/DxCodePPV.pth')
output_table.to_csv('../Results/DxCodePPV.txt',sep='\t')

