import pandas as pd 
import numpy as np 
import pickle
from sklearn.utils import resample
from scipy.stats import kruskal,spearmanr


#This script creates three data tables for machine learning penetrance model training or validation (depending on the dataset). It requires:
# 1) The expression measurement dataset (see Build_Global_Disease_Expression_Table.py)
# 2) The flags for low EHR data coverage (see IdentifyLowCoverageSubjects.py)
# 3) All of the carrier tables (see BuildCarrier_NonCarrier_Datasets.py)
# The script yields 3 tables, one for each pLoF class. 
# Note, the in HQ set was used to filter training data for UKBB model fitting. 

combined_penetrance_table = pd.read_pickle('/Path/to/Symptom/Penetrance/Results/GlobalExpressionTable.pth')
ehr_flags = pd.read_pickle('/Path/to/Coverage/Results/LowEHRDataFlags.pth')
combined_penetrance_table=combined_penetrance_table.loc[ehr_flags.LowEHRDataFlag==False]
carrier_class_direc='Path/to/CarrierDatasets/CarrierInfoFiles/'
allowed_diseases = combined_penetrance_table.DIS.unique()
chip_confounded_diseases = ['BOPS']

stop_gain_table={'INDEX':[],'SUBJECT_ID':[],'VARIANT_ID':[],'DIS':[],'ONSET':[],'LOF_TYPE':[],'SYMBOL':[],'CLINVAR_ANNOT':[],'CLINVAR_RATING':[],'TX_ANNOT':[],'LOFTEE_CLASS':[],'CADD':[],'PRED_NMD_ESCAPE':[],'FRAC_AA_IMPACTED':[],'POSSIBLE_MET_RESCUE':[],'PROB_EXPRESSION':[],'BINARIZED_EXPRESSION':[],'IN_HQ_SET':[]}
frameshift_table={'INDEX':[],'SUBJECT_ID':[],'VARIANT_ID':[],'DIS':[],'ONSET':[],'LOF_TYPE':[],'SYMBOL':[],'CLINVAR_ANNOT':[],'CLINVAR_RATING':[],'TX_ANNOT':[],'LOFTEE_CLASS':[],'CADD':[],'IN_LAST_CODING_EXON':[],'FRAC_AA_IMPACTED':[],'POSSIBLE_MET_RESCUE':[],'PROB_EXPRESSION':[],'BINARIZED_EXPRESSION':[],'IN_HQ_SET':[]}
splice_table={'INDEX':[],'SUBJECT_ID':[],'VARIANT_ID':[],'DIS':[],'ONSET':[],'LOF_TYPE':[],'SYMBOL':[],'CLINVAR_ANNOT':[],'CLINVAR_RATING':[],'TX_ANNOT':[],'LOFTEE_CLASS':[],'CADD':[],'SPLICE_TYPE':[],'SPLICE_AI_SCORE':[],'FRAC_AA_IMPACTED':[],'SPLICE_RESCUE_FLAGS':[],'PROB_EXPRESSION':[],'BINARIZED_EXPRESSION':[],'IN_HQ_SET':[]}



for dis_abbrev in allowed_diseases:
	with open(carrier_class_direc+'{0:s}_CarrierInfo.pth'.format(dis_abbrev),'rb') as f:
		carrier_class_data=pickle.load(f)
	lof_table = carrier_class_data['lof_table']
	for idx,subject in combined_penetrance_table.loc[combined_penetrance_table.DIS==dis_abbrev]['S_ID'].items():
		if lof_table.loc[subject,'LOF_TYPE']=='STOP_GAINED':
			stop_gain_table['INDEX']+=[idx]
			stop_gain_table['SUBJECT_ID']+=[subject]
			stop_gain_table['VARIANT_ID']+=[lof_table.loc[subject,'VARIANT']]
			stop_gain_table['DIS']+=[dis_abbrev]
			stop_gain_table['ONSET']+=[combined_penetrance_table.loc[idx,'ONSET']]
			stop_gain_table['LOF_TYPE']+=['STOP_GAINED']
			stop_gain_table['SYMBOL']+=[lof_table.loc[subject,'SYMBOL']]
			stop_gain_table['CLINVAR_ANNOT']+=[lof_table.loc[subject,'CLINVAR_ANNOT']]
			stop_gain_table['CLINVAR_RATING']+=[lof_table.loc[subject,'CLINVAR_RATING']]			
			stop_gain_table['TX_ANNOT']+=[lof_table.loc[subject,'TX_ANNOT']]
			stop_gain_table['LOFTEE_CLASS']+=[lof_table.loc[subject,'LOFTEE_CLASS']]
			stop_gain_table['CADD']+=[lof_table.loc[subject,'CADD']]
			stop_gain_table['PRED_NMD_ESCAPE']+=[lof_table.loc[subject,'PRED_NMD_ESCAPE'] if pd.isna(lof_table.loc[subject,'PRED_NMD_ESCAPE'])==False else 'FALSE']
			stop_gain_table['FRAC_AA_IMPACTED']+=[lof_table.loc[subject,'FRAC_AA_IMPACTED'] if pd.isna(lof_table.loc[subject,'FRAC_AA_IMPACTED'])==False else 0.0]
			stop_gain_table['POSSIBLE_MET_RESCUE']+=[lof_table.loc[subject,'MET_RESCUE']]
			stop_gain_table['PROB_EXPRESSION']+=[combined_penetrance_table.loc[idx][['PENETRANCE_PROB','HAS_DX']].max()]
			stop_gain_table['BINARIZED_EXPRESSION']+=[combined_penetrance_table.loc[idx]['BINARIZED_PENETRANCE']]
			stop_gain_table['IN_HQ_SET']+=[True if ((combined_penetrance_table.loc[idx]['BINARIZED_PENETRANCE']==1) or (combined_penetrance_table.loc[idx]['PENETRANCE_PROB']==0.0)) else False]

		if lof_table.loc[subject,'LOF_TYPE']=='FRAMESHIFT':
			frameshift_table['INDEX']+=[idx]
			frameshift_table['SUBJECT_ID']+=[subject]
			frameshift_table['VARIANT_ID']+=[lof_table.loc[subject,'VARIANT']]
			frameshift_table['DIS']+=[dis_abbrev]
			frameshift_table['ONSET']+=[combined_penetrance_table.loc[idx,'ONSET']]
			frameshift_table['LOF_TYPE']+=['FRAMESHIFT']
			frameshift_table['SYMBOL']+=[lof_table.loc[subject,'SYMBOL']]
			frameshift_table['CLINVAR_ANNOT']+=[lof_table.loc[subject,'CLINVAR_ANNOT']]
			frameshift_table['CLINVAR_RATING']+=[lof_table.loc[subject,'CLINVAR_RATING']]
			frameshift_table['TX_ANNOT']+=[lof_table.loc[subject,'TX_ANNOT']]
			frameshift_table['LOFTEE_CLASS']+=[lof_table.loc[subject,'LOFTEE_CLASS']]
			frameshift_table['CADD']+=[lof_table.loc[subject,'CADD']]
			frameshift_table['IN_LAST_CODING_EXON']+=[lof_table.loc[subject,'LAST_CODING_EXON']]
			frameshift_table['FRAC_AA_IMPACTED']+=[lof_table.loc[subject,'FRAC_AA_IMPACTED'] if pd.isna(lof_table.loc[subject,'FRAC_AA_IMPACTED'])==False else 0.0]
			frameshift_table['POSSIBLE_MET_RESCUE']+=[lof_table.loc[subject,'MET_RESCUE']]
			frameshift_table['PROB_EXPRESSION']+=[combined_penetrance_table.loc[idx][['PENETRANCE_PROB','HAS_DX']].max()]
			frameshift_table['BINARIZED_EXPRESSION']+=[combined_penetrance_table.loc[idx]['BINARIZED_PENETRANCE']]
			frameshift_table['IN_HQ_SET']+=[True if ((combined_penetrance_table.loc[idx]['BINARIZED_PENETRANCE']==1) or (combined_penetrance_table.loc[idx]['PENETRANCE_PROB']==0.0)) else False]
			
		if lof_table.loc[subject,'LOF_TYPE']=='SPLICE_CHANGE':
			splice_table['INDEX']+=[idx]
			splice_table['SUBJECT_ID']+=[subject]
			splice_table['VARIANT_ID']+=[lof_table.loc[subject,'VARIANT']]
			splice_table['DIS']+=[dis_abbrev]
			splice_table['ONSET']+=[combined_penetrance_table.loc[idx,'ONSET']]
			splice_table['LOF_TYPE']+=['SPLICE_CHANGE']
			splice_table['SYMBOL']+=[lof_table.loc[subject,'SYMBOL']]
			splice_table['CLINVAR_ANNOT']+=[lof_table.loc[subject,'CLINVAR_ANNOT']]
			splice_table['CLINVAR_RATING']+=[lof_table.loc[subject,'CLINVAR_RATING']]		
			splice_table['TX_ANNOT']+=[lof_table.loc[subject,'TX_ANNOT']]
			splice_table['LOFTEE_CLASS']+=[lof_table.loc[subject,'LOFTEE_CLASS']]
			splice_table['CADD']+=[lof_table.loc[subject,'CADD']]
			splice_table['SPLICE_TYPE']+=[lof_table.loc[subject,'SPLICE_TYPE']]
			splice_table['SPLICE_AI_SCORE']+=[lof_table.loc[subject,'PRIMARY_SPLICE_SCORE']]
			splice_table['FRAC_AA_IMPACTED']+=[lof_table.loc[subject,'FRAC_AA_IMPACTED'] if pd.isna(lof_table.loc[subject,'FRAC_AA_IMPACTED'])==False else 0.0]
			splice_table['SPLICE_RESCUE_FLAGS']+=[lof_table.loc[subject,'SPLICE_RESCUE_EVENTS']]
			splice_table['PROB_EXPRESSION']+=[combined_penetrance_table.loc[idx][['PENETRANCE_PROB','HAS_DX']].max()]
			splice_table['BINARIZED_EXPRESSION']+=[combined_penetrance_table.loc[idx]['BINARIZED_PENETRANCE']]
			splice_table['IN_HQ_SET']+=[True if ((combined_penetrance_table.loc[idx]['BINARIZED_PENETRANCE']==1) or (combined_penetrance_table.loc[idx]['PENETRANCE_PROB']==0.0)) else False]


stop_gain_table=pd.DataFrame(stop_gain_table)
stop_gain_table.set_index('INDEX',inplace=True)
stop_gain_table.to_pickle('/Path/to/ML/Model/Data/StopGainExpressionData.pth')

frameshift_table=pd.DataFrame(frameshift_table)
frameshift_table.set_index('INDEX',inplace=True)
frameshift_table.to_pickle('/Path/to/ML/Model/Data/FrameshiftExpressionData.pth')

splice_table=pd.DataFrame(splice_table)
splice_table.set_index('INDEX',inplace=True)
splice_table.to_pickle('/Path/to/ML/Model/Data/Results/SpliceChangeExpressionData.pth')


shared_columns = ['SUBJECT_ID','VARIANT_ID', 'DIS','ONSET','LOF_TYPE', 'SYMBOL', 'CLINVAR_ANNOT', 'CLINVAR_RATING','TX_ANNOT','LOFTEE_CLASS', 'CADD','FRAC_AA_IMPACTED','PROB_EXPRESSION','BINARIZED_EXPRESSION','IN_HQ_SET']

all_variant_types = pd.concat([stop_gain_table[shared_columns],frameshift_table[shared_columns],splice_table[shared_columns]],axis=0)
all_variant_types=all_variant_types.loc[all_variant_types.IN_HQ_SET==True]
all_variant_types_mane_select_only = all_variant_types.loc[all_variant_types.TX_ANNOT=='MANE_SELECT']
all_training_variants = all_variant_types.VARIANT_ID.unique()
with open('/Path/to/ML/Model/Data/HQTrainingVariants.txt','w') as f:
	for v_id in all_training_variants:
		f.write(v_id+'\n')
