import pandas as pd
import pickle
import numpy as np
import tqdm


def ConvertReviewToStars(rev_string):
	if rev_string=='practice_guideline':
		return 4
	elif rev_string=='reviewed_by_expert_panel':
		return 3
	elif rev_string=='criteria_provided&_multiple_submitters&_no_conflicts':
		return 2
	elif rev_string=='criteria_provided&_conflicting_interpretations':
		return 1
	elif rev_string=='criteria_provided&_single_submitter':
		return 1
	elif rev_string=='no_assertion_criteria_provided':
		return 1
	elif rev_string=='no_assertion_provided':
		return 1
	else:
		return 0

def NormalizeClinVarAnnots(annot):
	possible_clinvar_annots={}
	possible_clinvar_annots['']='Unannotated'
	possible_clinvar_annots['not_provided']='Unannotated'
	possible_clinvar_annots['risk_factor']='Uncertain'
	possible_clinvar_annots['Benign/Likely_benign']='Benign/Likely_Benign'
	possible_clinvar_annots['Benign']='Benign/Likely_Benign'
	possible_clinvar_annots['Likely_benign']='Benign/Likely_Benign'
	possible_clinvar_annots['Pathogenic/Likely_pathogenic']='Pathogenic/Likely_Pathogenic'
	possible_clinvar_annots['Conflicting_interpretations_of_pathogenicity']='Uncertain'
	possible_clinvar_annots['Likely_pathogenic']='Pathogenic/Likely_Pathogenic'
	possible_clinvar_annots['Pathogenic']='Pathogenic/Likely_Pathogenic'
	possible_clinvar_annots['Uncertain_significance']='Uncertain'
	possible_clinvar_annots['Conflicting_classifications_of_pathogenicity']='Uncertain'
	return possible_clinvar_annots[annot]

sample_count_directory='/Users/davidblair/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/SampleCounts/'

#Disease Dx Info
dx_code_file=pd.read_pickle('~/Desktop/Research/PDS_Project/Data/ClinGenHaploinsufficientDiseases/CollapsedDxTable.pth')

## Load carrier data table
carrier_file=pd.read_pickle('~/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/AllHaploLOFVariants/AllHaploLOFVariants_FilteredCarriers.pth')

## Score tables
basic_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/AllHaploLOFVariants/AllHaploLOFVariants_NoScores.pth')
cadd_table= pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/AllHaploLOFVariants/Scores/CADDScores.pth')
aa_lost_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/AllHaploLOFVariants/Scores/NumAAsImpacted_NonSplice.pth')
splice_ai_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/AllHaploLOFVariants/Scores/SpliceAIScores.pth')
nmd_escape_table = pd.read_pickle('/Users/davidblair/Desktop/Research/PDS_Project/Analysis/4_GenotypeAnnotation/UKBB/VariantData/AllHaploLOFVariants/Scores/PTC_NMDEscape.pth')


## Remove subjects who have NO diagnostic data whatsoever, which can bias results
coverage_data=pd.read_pickle('../Datasets/CovariateDatasets/UKBB_CoverageTable.pth')
allowed_subjects=coverage_data.index[pd.isna(coverage_data['CoverageYears'])==False]

## remove subjects filtered out of exome QC analysis
pass_exome_qc = list(map(str,pd.read_csv('/Users/davidblair/Desktop/Research/PDS_Project/Data/UKBB/ExomeSampleQC/All_UKKBB_QCSamples.txt',sep='\t',header=None)[0].values))
allowed_subjects=allowed_subjects.intersection(pass_exome_qc)

# ## all LOF Carriers
all_lof_carriers=set().union(*carrier_file['CARRIERS'].values).intersection(allowed_subjects)


for dis_abbrev in tqdm.tqdm(dx_code_file.index):
	genes = dx_code_file.loc[dis_abbrev]['SYMBOLS']

	sample_count_tables={}
	for gene in genes:
		sample_count_file=pd.read_csv(sample_count_directory+'{0:s}_PerSampleCounts.scount'.format(gene),sep='\t')
		sample_count_file.set_index('#IID',inplace=True)
		sample_count_tables[gene]=sample_count_file.copy()


	target_lof_variants = carrier_file.index[carrier_file.SYMBOL.apply(lambda x: x in genes)]

	target_lof_carriers=list(set().union(*carrier_file.loc[target_lof_variants]['CARRIERS'].values).intersection(allowed_subjects))
	target_no_calls=list(set().union(*carrier_file.loc[target_lof_variants]['NO_CALLS'].values).intersection(allowed_subjects))
	subject_to_variant_map={}
	for variant in target_lof_variants:
		for subject in carrier_file.loc[variant]['CARRIERS']:
			try:
				subject_to_variant_map[subject].add(variant)
			except KeyError:
				subject_to_variant_map[subject]=set([variant])

	lof_table = pd.DataFrame([],index=target_lof_carriers,columns=['SYMBOL','VARIANT','LOF_TYPE','TX_ANNOT','LOFTEE_CLASS','CADD','LAST_CODING_EXON','MET_RESCUE','FRAC_AA_IMPACTED','PRED_NMD_ESCAPE','SPLICE_TYPE','PRIMARY_SPLICE_SCORE','SPLICE_RESCUE_EVENTS','CLINVAR_ANNOT','CLINVAR_RATING'])
	
	for subject,variants in subject_to_variant_map.items():
		# exclude subjescts with multiple variants. This a rare but confusing situation. 
		if len(variants)==1:
			variant = variants.pop()
			lof_table.loc[subject,'SYMBOL']=basic_table.loc[variant]['SYMBOL']
			lof_table.loc[subject,'VARIANT']=variant
			lof_table.loc[subject,'LOF_TYPE']=aa_lost_table.loc[variant]['CSQ_CLASS']
			lof_table.loc[subject,'TX_ANNOT']=basic_table.loc[variant]['TX_ANNOT']
			lof_table.loc[subject,'LOFTEE_CLASS']=basic_table.loc[variant]['LOFTEE_CLASS']
			lof_table.loc[subject,'CADD']=cadd_table.loc[variant]['CADD']
			if lof_table.loc[subject,'LOF_TYPE']!='SPLICE_CHANGE':
				lof_table.loc[subject,'LAST_CODING_EXON']=True if (pd.isna(aa_lost_table.loc[variant,'FLAG'])==False and (aa_lost_table.loc[variant,'FLAG']=='LAST_CODING_EXON')) else False
				lof_table.loc[subject,'MET_RESCUE']=True if (pd.isna(aa_lost_table.loc[variant,'FLAG'])==False and (aa_lost_table.loc[variant,'FLAG']=='POSSIBLE_MET_RESCUE')) else False
				lof_table.loc[subject,'FRAC_AA_IMPACTED']=aa_lost_table.loc[variant,'FRAC_AA_IMPACTED']
			else:
				lof_table.loc[subject,'LAST_CODING_EXON']=True if 'LAST_CODING_EXON' in set([x[0] for x in splice_ai_table.loc[variant,'RESCUE_FLAGS']]) else False
				lof_table.loc[subject,'MET_RESCUE']=True if 'POSSIBLE_MET_RESCUE' in set([x[0] for x in splice_ai_table.loc[variant,'RESCUE_FLAGS']]) else False
				lof_table.loc[subject,'FRAC_AA_IMPACTED']=splice_ai_table.loc[variant,'FRAC_CDS_IMPACTED']

			if lof_table.loc[subject,'LOF_TYPE']=='STOP_GAINED':
				lof_table.loc[subject,'PRED_NMD_ESCAPE']=nmd_escape_table.loc[variant]['NMD_ESCAPE_FLAG']
			elif lof_table.loc[subject,'LOF_TYPE']=='SPLICE_CHANGE':
				lof_table.loc[subject,'SPLICE_TYPE']=splice_ai_table.loc[variant,'PRIMARY_SPLICE_TYPE']
				lof_table.loc[subject,'PRIMARY_SPLICE_SCORE']=splice_ai_table.loc[variant,'PRIMARY_SPLICEAI_SCORE']
				lof_table.loc[subject,'SPLICE_RESCUE_EVENTS']=splice_ai_table.loc[variant,'RESCUE_FLAGS']

			lof_table.loc[subject,'CLINVAR_ANNOT']=NormalizeClinVarAnnots(basic_table.loc[variant]['CLINVAR_CLNSIG'])
			lof_table.loc[subject,'CLINVAR_RATING']=ConvertReviewToStars(basic_table.loc[variant]['CLINVAR_CLNREVSTAT'])

	other_rare_variant_carriers=set()
	non_carriers=allowed_subjects.difference(all_lof_carriers).difference(target_no_calls).difference(lof_table.index)
	for gene in genes:
		all_carriers = sample_count_tables[gene].index[sample_count_tables[gene][['HOM_ALT_SNP_CT','HET_SNP_CT','DIPLOID_NONSNP_NONSYMBOLIC_CT']].sum(axis=1)>0]
		non_lof_carriers = all_carriers.difference(lof_table.index)
		other_rare_variant_carriers.update(non_lof_carriers.intersection(allowed_subjects).difference(all_lof_carriers).difference(target_no_calls))
	non_carriers=list(non_carriers.drop(other_rare_variant_carriers))
	other_rare_variant_carriers=list(other_rare_variant_carriers)

	dataset={}
	dataset['lof_table']=lof_table
	dataset['other_rare_variant_carriers']=other_rare_variant_carriers
	dataset['non_carriers']=non_carriers

	with open('../Datasets/CarrierInfoFiles/{0:s}_CarrierInfo.pth'.format(dis_abbrev),'wb') as f:
		pickle.dump(dataset,f)




