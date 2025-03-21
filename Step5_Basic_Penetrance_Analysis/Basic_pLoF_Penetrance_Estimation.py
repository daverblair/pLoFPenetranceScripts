import pandas as pd 
import numpy as np
import tqdm
from scipy.stats import beta
import tqdm
import pickle

# #This script computes basic penetrance estimates for pLoFs using a simple beta-binomial model. It requires the following input information:
# 1) Table of basic disease information (disease_table)
# 2) Table containing the biobank covariates (age at enrollment, sex, ancestry PCs; must be indexed by subject ID)
# 3) Path to the pLoF-carrier tables created using Step3_Build_Carrier_and_NonCarrier_Datasets/BuildCarrier_NonCarrier_Datasets.py in Step4_CreateClinicalDatasets directory
# 4) Table of rare disease diagnoses, with one row per biobank subject and a simple binary indicator ('HasDisDx') that indicates whether the subject carries a disease-specific diagnosis code for this disorder. See CollapsedDxTable.txt in Auxillary_Data/CollapsedDxTable.txt for list of disease-specific diagnostic codes. These tables can be generated with BuildDiseaseDxTables.py in Step4_CreateClinicalDatasets directory

PSEUDO_COUNT=0.5


disease_table = pd.read_pickle('/Path/to/Auxillary_Data/CollapsedDxTable.pth')
covariate_table = pd.read_csv('/Path/to/Covariate/Data/Covariates.txt',sep='\t')
covariate_table.set_index('Subject_ID',inplace=True)
carrier_class_direc='Path/to/CarrierDatasets/CarrierInfoFiles/'
dx_data_direc='/Path/to/Haploinsufficient/Disease/Dx/Tables/'

# Bohring Opitz was dropped from the analysis, as it's associated gene carries many pLoFs due to CHiP, confounding the analysis
chip_confounded_diseases = ['BOPS']


output_table={'DisAbbrev':[],'CarrierCount':[],'FilteredCarrierCount':[],'Penetrance':[],'FilteredVariant_Penetrance':[]}
for dis_abbrev in tqdm.tqdm(disease_table.index.drop(chip_confounded_diseases)):
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		carrier_class_data=pickle.load(f)

	# drop carriers that are not in the covariate table due to missing data
	target_lof_carriers=carrier_class_data['lof_table'].index.intersection(covariate_table.index)

	# simple filter, must be in MANE Select and have HC LOFTEE flag
	simple_filter = carrier_class_data['lof_table'].index[(carrier_class_data['lof_table'].LOFTEE_CLASS=='HC')*(carrier_class_data['lof_table'].TX_ANNOT=='MANE_SELECT')]
	target_lof_carriers_filtered=simple_filter.intersection(covariate_table.index)

	dx_data=pd.read_pickle(dx_data_direc+'{0:s}_DxDataset.pth'.format(dis_abbrev))
	output_table['DisAbbrev']+=[dis_abbrev]
	if (len(target_lof_carriers)>0) and (dx_data['HasDisDx'].sum()>0):

		num_dx_lof_carriers=dx_data.loc[target_lof_carriers]['HasDisDx'].sum()
		num_undx_lof_carriers=(dx_data.loc[target_lof_carriers]['HasDisDx']==0).sum()

		penetrance=beta(num_dx_lof_carriers+PSEUDO_COUNT,num_undx_lof_carriers+PSEUDO_COUNT).mean()
		penetrance_95=beta(num_dx_lof_carriers+PSEUDO_COUNT,num_undx_lof_carriers+PSEUDO_COUNT).interval(confidence=0.95)

		output_table['Penetrance']+=[{'Mean':penetrance,'95_CI':penetrance_95}]
		output_table['CarrierCount']+=[num_dx_lof_carriers+num_undx_lof_carriers]

		num_dx_lof_carriers_filtered=dx_data.loc[target_lof_carriers_filtered]['HasDisDx'].sum()
		num_undx_lof_carriers_filtered=(dx_data.loc[target_lof_carriers_filtered]['HasDisDx']==0).sum()

		penetrance_filtered=beta(num_dx_lof_carriers_filtered+PSEUDO_COUNT,num_undx_lof_carriers_filtered+PSEUDO_COUNT).mean()
		penetrance_95_filtered=beta(num_dx_lof_carriers_filtered+PSEUDO_COUNT,num_undx_lof_carriers_filtered+PSEUDO_COUNT).interval(confidence=0.95)

		output_table['FilteredVariant_Penetrance']+=[{'Mean':penetrance_filtered,'95_CI':penetrance_95_filtered}]
		output_table['FilteredCarrierCount']+=[num_dx_lof_carriers_filtered+num_undx_lof_carriers_filtered]

	else:
		output_table['Penetrance']+=[{'Mean':np.nan,'95_CI':(np.nan,np.nan)}]
		output_table['FilteredVariant_Penetrance']+=[{'Mean':np.nan,'95_CI':(np.nan,np.nan)}]
		output_table['CarrierCount']+=[np.nan]
		output_table['FilteredCarrierCount']+=[np.nan]

output_table=pd.DataFrame(output_table)
output_table.set_index('DisAbbrev',inplace=True)
output_table.to_pickle('/Path/to/Basic/Penetrance/Results/BasicPenetrance.pth')
output_table.to_csv('/Path/to/Basic/Penetrance/Results/BasicPenetrance.txt',sep='\t')

