from cyvcf2 import VCF,Writer
import pandas as pd
import copy

def ParseLoFVariants(transcript,position,delta_aa):
	if len(delta_aa)>1:
		delta_aa=delta_aa.split('/')
		return (transcript,position,delta_aa[0],delta_aa[1])
	elif len(delta_aa)==1:
		return (transcript,position,delta_aa[0],delta_aa[0])
	else:
		return (transcript,position,'','')

#initialize new 
gene_info_table = pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/HaploinsuffientGenes_BasicGeneInfo.pth')
pc_transcript_table=pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/HaploinsuffientGenes_CompleteTranscriptInfo.pth')


lof_table={'CHROM':[],'POS':[],'SYMBOL':[],'GENE_ID':[],'ID':[],'TX_ID':[],'REF':[],'ALT':[],'CONSEQUENCE':[],'PROT_POS':[],'REF_AA':[],'ALT_AA':[],'LOFTEE_CLASS':[],'CLINVAR_CLNSIG':[],'CLINVAR_CLNREVSTAT':[]}
rev_gene_list=open('../VariantData/TargetVariants/GeneList_wLoF.txt','w')
for symbol in gene_info_table.index:
	
	print('Identifying LoF variants for gene {0:s} ...'.format(symbol))
	f=open('../VariantData/TargetVariants/{0:s}_HaploLOFVariants_TargetCoords.txt'.format(symbol),'w')
	current_vcf=VCF("../VariantData/Annotated_GeneVCFs/{0:s}_UKBB_Exome_Annotated_NoSamples.vcf.gz".format(symbol))
	header_id_dict={}
	for header_entry in current_vcf.header_iter():
		if header_entry.type=='INFO':
			header_id_dict[header_entry.info()['ID']]=header_entry.info()['Description']

	CSQ_Columns=header_id_dict['CSQ'].strip('"').split('Format:')[1].strip().split('|')
	total_LoFs=0
	for variant in current_vcf:
		info=variant.INFO.get('CSQ').split(',')
		if len(variant.ALT)>1:
			AF=list(variant.INFO.get('AF'))
		else:
			AF=[variant.INFO.get('AF')]
		num_features=int(len(info)/len(variant.ALT))
		num_variants=len(variant.ALT)
		id_list=variant.ID.split(';')

		is_lof=False
		for t_num in range(num_features):
			for v_num in range(num_variants):
				parsed_info=dict(zip(CSQ_Columns,info[t_num*num_variants+v_num].split('|')))
				if parsed_info['Feature_type']=='Transcript' and parsed_info['Feature'] in gene_info_table.loc[symbol]['ALL_PC_TXs']:
					if parsed_info['LoF'] in set(['HC','LC']):
						#LOF Variant!
						is_lof=True
						lof_info=ParseLoFVariants(parsed_info['Feature'],parsed_info['Protein_position'],parsed_info['Amino_acids'])
						lof_table['CHROM']+=[variant.CHROM[3:]]
						lof_table['POS']+=[int(variant.POS)]
						lof_table['ID']+=['{0:s}_{1:d}_{2:s}_{3:s}'.format(variant.CHROM,variant.POS,variant.REF,variant.ALT[v_num])]
						lof_table['REF']+=[variant.REF]
						lof_table['ALT']+=[variant.ALT[v_num]]
						lof_table['SYMBOL']+=[symbol]
						lof_table['CONSEQUENCE']+=[parsed_info['Consequence']]
						lof_table['PROT_POS']+=[lof_info[1]]
						lof_table['REF_AA']+=[lof_info[2]]
						lof_table['ALT_AA']+=[lof_info[3]]
						lof_table['GENE_ID']+=[gene_info_table.loc[symbol]['ENSEMBL_ID']]
						lof_table['TX_ID']+=[parsed_info['Feature']]
						lof_table['LOFTEE_CLASS']+=[parsed_info['LoF']]
						lof_table['CLINVAR_CLNSIG']+=[parsed_info['ClinVar_CLNSIG']]
						lof_table['CLINVAR_CLNREVSTAT']+=[parsed_info['ClinVar_CLNREVSTAT']]	
						total_LoFs+=1
		if is_lof==True:
			f.write('{0:s}\t{1:d}\n'.format(variant.CHROM,variant.POS))
	if total_LoFs>0:
		rev_gene_list.write('{0:s}\n'.format(symbol))
	f.close()
rev_gene_list.close()

variant_table=pd.DataFrame(lof_table)
duplicated_info = variant_table.loc[variant_table['ID'].duplicated(keep=False)].copy()
duplicated_ids = variant_table.loc[variant_table['ID'].duplicated(keep=False)]['ID'].unique()
duplicated_info.set_index('ID',inplace=True)


variant_table_collapsed=variant_table.loc[variant_table['ID'].duplicated(keep=False)==False].copy()
variant_table_collapsed.set_index('ID',inplace=True)

def _identify_tx_label(symbol,tx_id):
	if tx_id==gene_info_table.loc[symbol]['MANE_SELECT_TX']:
		return 'MANE_SELECT'
	elif tx_id in gene_info_table.loc[symbol]['MANE_SELECT_TX']:
		return 'MANE_PLUS_CLINICAL'
	else:
		return 'NONE'

variant_table_collapsed['TX_ANNOT']=variant_table_collapsed[['SYMBOL','TX_ID']].apply(lambda x: _identify_tx_label(x['SYMBOL'],x['TX_ID']),axis=1)

#most variants have duplicate effects. To isolate only one transcript-variant pair, kept the top pair according to the following hierarchy:
# 1) MANE select
# 2) MANE Plus clinical
# 3) Longest transcript (CDS)
for d_id in duplicated_ids:
	current_dups=duplicated_info.loc[d_id]
	tx_ids = set(current_dups.TX_ID.values)
	symbol = current_dups.SYMBOL.unique()[0]
	if gene_info_table.loc[symbol]['MANE_SELECT_TX'] in tx_ids:
		target_tx = gene_info_table.loc[symbol]['MANE_SELECT_TX']
		target_tx = current_dups.loc[current_dups.TX_ID==target_tx].copy()
		target_tx['TX_ANNOT'] = pd.Series(['MANE_SELECT'],index=[d_id])
		variant_table_collapsed=variant_table_collapsed._append(target_tx,ignore_index=False)
	elif len(gene_info_table.loc[symbol]['MANE_PLUS_CLINICAL'].intersection(tx_ids))>0:
		mane_plus_clinical = list(gene_info_table.loc[symbol]['MANE_PLUS_CLINICAL'].intersection(tx_ids))
		target_tx = mane_plus_clinical[pc_transcript_table.loc[mane_plus_clinical]['CDS_LENGTH'].argmax()]
		target_tx = current_dups.loc[current_dups.TX_ID==target_tx].copy()
		target_tx['TX_ANNOT'] = pd.Series(['MANE_PLUS_CLINICAL'],index=[d_id])
		variant_table_collapsed=variant_table_collapsed._append(target_tx,ignore_index=False)
	else:
		all_impacted_tx = list(pc_transcript_table.index.intersection(tx_ids))
		target_tx = all_impacted_tx[pc_transcript_table.loc[all_impacted_tx]['CDS_LENGTH'].argmax()]
		target_tx = current_dups.loc[current_dups.TX_ID==target_tx].copy()
		target_tx['TX_ANNOT'] = pd.Series(['NONE'],index=[d_id])
		variant_table_collapsed=variant_table_collapsed._append(target_tx,ignore_index=False)

variant_table_collapsed=variant_table_collapsed.loc[variant_table_collapsed.index.sort_values()]
variant_table_collapsed=variant_table_collapsed[['CHROM','POS','SYMBOL','GENE_ID','TX_ID','TX_ANNOT','REF','ALT','CONSEQUENCE','PROT_POS','REF_AA','ALT_AA','LOFTEE_CLASS','CLINVAR_CLNSIG','CLINVAR_CLNREVSTAT']]
variant_table_collapsed.to_pickle('../VariantData/AllHaploLOFVariants/AllHaploLOFVariants_NoScores.pth')
variant_table_collapsed.to_csv('../VariantData/AllHaploLOFVariants/AllHaploLOFVariants_NoScores.txt',sep='\t')

