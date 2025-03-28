import pandas as pd
import numpy as np
import tqdm
import sys
sys.path.append('/Users/davidblair/Desktop/Research/VariantAnnotation/AuxFunctions')
from SequenceAnalysisFunctions import NormalizeVariantCsq,AlternateVariantIndex


# This script creates a table of CADD scores for the pLoFs found in each biobank. Note, the variant data from the VCF files was submitted to the CADD website (https://cadd.gs.washington.edu/score), and this script parses the resulting tsv file. 

variant_table = pd.read_pickle('Path/to/Variant/Data/AllHaploLOFVariants_NoScores.pth')
cadd_lof=pd.read_csv('/Path/to/VariantCADDScores/AllHaploLOFVariants_CADD_Scores.tsv',sep='\t',skiprows=[0])
cadd_lof['VARIANT_ID']=cadd_lof[['#Chrom','Pos','Ref', 'Alt']].apply(lambda x:'chr'+'_'.join([str(x[cn]) for cn in x.index]),axis=1)
cadd_lof.set_index('VARIANT_ID',inplace=True)

output={'VARIANT_ID':[],'SYMBOL':[],'CSQ_CLASS':[],'CADD':[]}
for variant_id in tqdm.tqdm(variant_table.index):
	try:
		cadd_score=cadd_lof.loc[variant_id]['PHRED']
	except KeyError:
		try:
			cadd_score=cadd_lof.loc[AlternateVariantIndex(variant_id)[0]]['PHRED']
		except KeyError:
			cadd_score=pd.NA
			print('Failed to find CADD score for Variant {0:s}'.format(variant_id))
	output['VARIANT_ID']+=[variant_id]
	output['SYMBOL']+=[variant_table.loc[variant_id].SYMBOL]
	output['CSQ_CLASS']+=[NormalizeVariantCsq(variant_table.loc[variant_id]['CONSEQUENCE'].split('&'))]
	output['CADD']+=[cadd_score]
output=pd.DataFrame(output)
output.set_index('VARIANT_ID',inplace=True)
output.to_pickle('/Path/to/ScoreOutput/CADDScores.pth')
output.to_csv('/Path/to/ScoreOutput/CADDScores.txt',sep='\t',index=True)