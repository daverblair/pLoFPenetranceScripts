import pandas as pd
import numpy as np
import tqdm
import sys
sys.path.append('/Path/to/Auxillary_Functions/')
from SequenceAnalysisFunctions import NormalizeVariantCsq,ParseSpliceAI,ReturnPrimarySpliceInfo,PredictDonorGainRescue,PredictDonorLossRescue,PredictAcceptorGainRescue,PredictAcceptorLossRescue

#This script parses the output from SpliceAI, generating the genomic features described in the Supplementary Methods, including the fraction of amino acids impacted.


variant_table = pd.read_pickle('/Path/to/Variants/AllHaploLOFVariants_NoScores.pth')
spliceai_lof=ParseSpliceAI('/Path/to/SpliceAI/Results/AllHaploLOFVariants_spliceAI_Scores_500bp.vcf.gz')
normalized_varaint_csq=variant_table['CONSEQUENCE'].apply(lambda x:NormalizeVariantCsq(x.split('&')))
pc_transcript_table=pd.read_pickle('/Path/to/Auxillary_Data/HaploinsuffientGenes_CompleteTranscriptInfo.pth')
full_sequence_table=pd.read_pickle('/Path/to/Auxillary_Data/HaploinsuffientGenes_FullSequences.pth')

output={'VARIANT_ID':[],'SYMBOL':[],'CSQ_CLASS':[],'PRIMARY_SPLICE_TYPE':[],'PRIMARY_SPLICEAI_SCORE':[],'FRAC_CDS_IMPACTED':[],'RESCUE_FLAGS':[]}


for variant_id in tqdm.tqdm(variant_table.index):
	variant_class=normalized_varaint_csq.loc[variant_id]
	symbol=variant_table.loc[variant_id]['SYMBOL']
	output['VARIANT_ID']+=[variant_id]
	output['CSQ_CLASS']+=[variant_class]
	output['SYMBOL']+=[symbol]
	if (variant_class=='SPLICE_CHANGE'):
		tx_id = variant_table.loc[variant_id,'TX_ID']
		all_tx_spliceai=spliceai_lof.loc[[variant_id]]
		target_splice_ai_info = all_tx_spliceai.loc[all_tx_spliceai.TX_ID==tx_id]['SPLICE_AI'].iloc[0]
		splice_type,max_splice_score,frac_coding_nuc_impacted=ReturnPrimarySpliceInfo(target_splice_ai_info,variant_table.loc[variant_id,'POS'],pc_transcript_table.loc[tx_id].STRAND,pc_transcript_table.loc[tx_id].TX_EXONS,pc_transcript_table.loc[tx_id]['5_UTR'],pc_transcript_table.loc[tx_id]['3_UTR'])
		output['PRIMARY_SPLICE_TYPE']+=[splice_type]
		output['PRIMARY_SPLICEAI_SCORE']+=[max_splice_score]
		output['FRAC_CDS_IMPACTED']+=[frac_coding_nuc_impacted]

		if splice_type=='DONOR_GAIN':
			flags=PredictDonorGainRescue(target_splice_ai_info,variant_table.loc[variant_id,'POS'],pc_transcript_table.loc[tx_id].STRAND,pc_transcript_table.loc[tx_id].TX_EXONS,pc_transcript_table.loc[tx_id]['5_UTR'],pc_transcript_table.loc[tx_id]['3_UTR'])

		if splice_type=='DONOR_LOSS':
			flags=PredictDonorLossRescue(target_splice_ai_info,variant_table.loc[variant_id,'POS'],pc_transcript_table.loc[tx_id].STRAND,pc_transcript_table.loc[tx_id].TX_EXONS,pc_transcript_table.loc[tx_id]['5_UTR'],pc_transcript_table.loc[tx_id]['3_UTR'],full_sequence_table.loc[tx_id]['AA'])


		if splice_type=='ACC_GAIN':
			flags=PredictAcceptorGainRescue(target_splice_ai_info,variant_table.loc[variant_id,'POS'],pc_transcript_table.loc[tx_id].STRAND,pc_transcript_table.loc[tx_id].TX_EXONS,pc_transcript_table.loc[tx_id]['5_UTR'],pc_transcript_table.loc[tx_id]['3_UTR'])

		if splice_type=='ACC_LOSS':
			flags=PredictAcceptorLossRescue(target_splice_ai_info,variant_table.loc[variant_id,'POS'],pc_transcript_table.loc[tx_id].STRAND,pc_transcript_table.loc[tx_id].TX_EXONS,pc_transcript_table.loc[tx_id]['5_UTR'],pc_transcript_table.loc[tx_id]['3_UTR'],full_sequence_table.loc[tx_id]['AA'])

		output['RESCUE_FLAGS']+=[flags]

	else:
		output['PRIMARY_SPLICE_TYPE']+=[pd.NA]
		output['PRIMARY_SPLICEAI_SCORE']+=[pd.NA]
		output['FRAC_CDS_IMPACTED']+=[pd.NA]
		output['RESCUE_FLAGS']+=[pd.NA]

output=pd.DataFrame(output)
output.set_index('VARIANT_ID',inplace=True)
output.to_pickle('/Path/to/ScoreOutput/SpliceAIScores.pth')
output.to_csv('/Path/to/ScoreOutput/SpliceAIScores.txt',sep='\t',index=True)