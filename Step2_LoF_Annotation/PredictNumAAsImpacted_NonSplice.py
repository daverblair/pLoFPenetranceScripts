import pandas as pd
import numpy as np
import tqdm
import sys
sys.path.append('/Users/davidblair/Desktop/Research/VariantAnnotation/AuxFunctions')
from SequenceAnalysisFunctions import *


pc_transcript_table=pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/HaploinsuffientGenes_CompleteTranscriptInfo.pth')
full_sequences=pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/HaploinsuffientGenes_FullSequences.pth')
variant_table = pd.read_pickle('../VariantData/AllHaploLOFVariants/AllHaploLOFVariants_NoScores.pth')
normalized_varaint_csq=variant_table['CONSEQUENCE'].apply(lambda x:NormalizeVariantCsq(x.split('&')))
pc_transcript_table['AA_LENGTH']=full_sequences.loc[pc_transcript_table.index]['AA'].apply(lambda x: len(x))



output={'VARIANT_ID':[],'SYMBOL':[],'CSQ_CLASS':[],'NUM_AA_IMPACTED':[],'FRAC_AA_IMPACTED':[],'FLAG':[]}
for variant_id in tqdm.tqdm(variant_table.index):
	variant_class=normalized_varaint_csq.loc[variant_id]
	symbol=variant_table.loc[variant_id].SYMBOL

	if variant_class!='SPLICE_CHANGE':
		impacted_tx=variant_table.loc[variant_id,'TX_ID']
		num_aa_lost,frac_aa_lost,flag=PredictAAImpacted_NoSplice(variant_class,variant_table.loc[variant_id,'PROT_POS'],variant_table.loc[variant_id,'REF_AA'],variant_table.loc[variant_id,'ALT_AA'],full_sequences.loc[impacted_tx]['AA'],pc_transcript_table.loc[impacted_tx]['TX_EXONS'],pc_transcript_table.loc[impacted_tx]['5_UTR'],pc_transcript_table.loc[impacted_tx]['3_UTR'],pc_transcript_table.loc[impacted_tx]['STRAND'],met_rescue_threshold=0.9)

		output['VARIANT_ID']+=[variant_id]
		output['SYMBOL']+=[symbol]
		output['CSQ_CLASS']+=[variant_class]
		output['NUM_AA_IMPACTED']+=[num_aa_lost]
		output['FRAC_AA_IMPACTED']+=[frac_aa_lost]
		output['FLAG']+=[flag]

	else:
		output['VARIANT_ID']+=[variant_id]
		output['SYMBOL']+=[symbol]
		output['CSQ_CLASS']+=[variant_class]
		output['NUM_AA_IMPACTED']+=[pd.NA]
		output['FRAC_AA_IMPACTED']+=[pd.NA]
		output['FLAG']+=[pd.NA]		


output=pd.DataFrame(output)
output.set_index('VARIANT_ID',inplace=True)
output.to_pickle('../VariantData/AllHaploLOFVariants/Scores/NumAAsImpacted_NonSplice.pth')
output.to_csv('../VariantData/AllHaploLOFVariants/Scores/NumAAsImpacted_NonSplice.txt',sep='\t',index=True)



