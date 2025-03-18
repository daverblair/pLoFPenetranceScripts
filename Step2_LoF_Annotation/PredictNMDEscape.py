import pandas as pd
import numpy as np
import tqdm
import gget
import sys
sys.path.append('/Users/davidblair/Desktop/Research/VariantAnnotation/AuxFunctions')
from SequenceAnalysisFunctions import *


pc_transcript_table=pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/HaploinsuffientGenes_CompleteTranscriptInfo.pth')
variant_table = pd.read_pickle('../VariantData/AllHaploLOFVariants/AllHaploLOFVariants_NoScores.pth')
normalized_varaint_csq=variant_table['CONSEQUENCE'].apply(lambda x:NormalizeVariantCsq(x.split('&')))



output={'VARIANT_ID':[],'SYMBOL':[],'CSQ_CLASS':[],'POSSIBLE_NMD_ESCAPE':[],'NMD_ESCAPE_FLAG':[]}
for variant_id in tqdm.tqdm(variant_table.index):
	variant_class=normalized_varaint_csq.loc[variant_id]
	symbol=variant_table.loc[variant_id].SYMBOL
	tx_id=variant_table.loc[variant_id].TX_ID
	if variant_class=='STOP_GAINED':
		alt_id,alt_info=AlternateVariantIndex(variant_id)
		nmd_escape,flag=PredictEscapeNMD(alt_info['POS'],pc_transcript_table.loc[tx_id]['TX_EXONS'],pc_transcript_table.loc[tx_id]['5_UTR'],pc_transcript_table.loc[tx_id]['3_UTR'],pc_transcript_table.loc[tx_id]['STRAND'])
	else:
		nmd_escape=pd.NA
		flag=pd.NA

	output['VARIANT_ID']+=[variant_id]
	output['SYMBOL']+=[symbol]
	output['CSQ_CLASS']+=[variant_class]
	output['POSSIBLE_NMD_ESCAPE']+=[nmd_escape]
	output['NMD_ESCAPE_FLAG']+=[flag]

output=pd.DataFrame(output)
output.set_index('VARIANT_ID',inplace=True)
output.to_pickle('../VariantData/AllHaploLOFVariants/Scores/PTC_NMDEscape.pth')
output.to_csv('../VariantData/AllHaploLOFVariants/Scores/PTC_NMDEscape.txt',sep='\t',index=True)


