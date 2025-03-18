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
chip_confounded_diseases = ['BOPS']

combined_expression_table = pd.read_pickle('../Results/GlobalPenetranceTable.pth')

output_table={'DisAbbrev':[],'CarrierCount':[],'FilteredCarrierCount':[],'HasDxData':[],'PPV':[],'FilteredVariant_PPV':[]}

for dis_abbrev in tqdm.tqdm(disease_table.index.drop(chip_confounded_diseases)):
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		carrier_class_data=pickle.load(f)
	target_lof_carriers=carrier_class_data['lof_table'].index.intersection(covariate_table.index)
	simple_filter = carrier_class_data['lof_table'].index[(carrier_class_data['lof_table'].LOFTEE_CLASS=='HC')*(carrier_class_data['lof_table'].TX_ANNOT=='MANE_SELECT')]
	target_lof_carriers_filtered=simple_filter.intersection(covariate_table.index)

	output_table['DisAbbrev']+=[dis_abbrev]

	sub_table = combined_expression_table.loc[combined_expression_table.DIS==dis_abbrev].set_index('S_ID')
	if (len(target_lof_carriers)>0):
		if (pd.isna(sub_table['HAS_DX']).sum()!= sub_table.shape[0]):
			output_table['HasDxData']+=[True]
		else:
			output_table['HasDxData']+=[False]

		dis_specific_combined_expression = sub_table.loc[sub_table.DIS==dis_abbrev]['MAX_PENETRANCE']

		num_exp_lof_carriers=dis_specific_combined_expression.loc[target_lof_carriers].sum()
		num_unexp_lof_carriers=target_lof_carriers.shape[0]-num_exp_lof_carriers

		expression=beta(num_exp_lof_carriers+PSEUDO_COUNT,num_unexp_lof_carriers+PSEUDO_COUNT).mean()
		expression_95=beta(num_exp_lof_carriers+PSEUDO_COUNT,num_unexp_lof_carriers+PSEUDO_COUNT).interval(confidence=0.95)

		output_table['PPV']+=[{'Mean':expression,'95_CI':expression_95}]
		output_table['CarrierCount']+=[num_exp_lof_carriers+num_unexp_lof_carriers]


		num_exp_lof_carriers_filtered=dis_specific_combined_expression.loc[target_lof_carriers_filtered].sum()
		num_unexp_lof_carriers_filtered=target_lof_carriers_filtered.shape[0]-num_exp_lof_carriers_filtered

		expression_filtered=beta(num_exp_lof_carriers_filtered+PSEUDO_COUNT,num_unexp_lof_carriers_filtered+PSEUDO_COUNT).mean()
		expression_95_filtered=beta(num_exp_lof_carriers_filtered+PSEUDO_COUNT,num_unexp_lof_carriers_filtered+PSEUDO_COUNT).interval(confidence=0.95)

		output_table['FilteredVariant_PPV']+=[{'Mean':expression_filtered,'95_CI':expression_95_filtered}]	
		output_table['FilteredCarrierCount']+=[num_exp_lof_carriers_filtered+num_unexp_lof_carriers_filtered]	
	else:
		output_table['PPV']+=[{'Mean':np.nan,'95_CI':(np.nan,np.nan)}]
		output_table['FilteredVariant_PPV']+=[{'Mean':np.nan,'95_CI':(np.nan,np.nan)}]
		output_table['CarrierCount']+=[np.nan]
		output_table['HasDxData']+=[np.nan]
		output_table['FilteredCarrierCount']+=[np.nan]
		


output_table=pd.DataFrame(output_table)
output_table.set_index('DisAbbrev',inplace=True)
output_table.to_pickle('../Results/CombinedDxSymptomsPPV.pth')
output_table.to_csv('../Results/CombinedDxSymptomsPPV.txt',sep='\t')

