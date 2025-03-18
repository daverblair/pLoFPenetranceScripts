import pandas as pd
import numpy as np
from cyvcf2 import VCF,Writer
import tqdm

# This script identifies the carriers of every pLoF variant, assuming they are stored in a (sorted) VCF file. If using this script on All of Us data, make sure to delete the line flagged below 
master_table=pd.read_pickle('/Path/to/Output/From/IdentifyLoFs.py/AllHaploLOFVariants_NoScores.pth')


#output table
variant_table=pd.DataFrame([],index=master_table.index)
variant_table.loc[variant_table.index,'SYMBOL']=master_table.loc[variant_table.index]['SYMBOL']
variant_table['CARRIER_FREQ']=pd.Series([np.nan for x in range(len(variant_table))],index=variant_table.index)
variant_table['CALL_RATE']=pd.Series([np.nan for x in range(len(variant_table))],index=variant_table.index)
variant_table['CARRIERS']=pd.Series([set() for x in range(len(variant_table))],index=variant_table.index)
variant_table['NO_CALLS']=pd.Series([set() for x in range(len(variant_table))],index=variant_table.index)


#
vcf_file=VCF('/Path/to/VCF/Files/AllHaploLOFVariants_wSamples_Sorted.vcf.gz')
header_id_dict={}
for header_entry in vcf_file.header_iter():
	if header_entry.type=='INFO':
		header_id_dict[header_entry.info()['ID']]=header_entry.info()['Description']
samples=np.array(vcf_file.samples)

for parent_idx,row in tqdm.tqdm(variant_table.iterrows(),total=variant_table.shape[0]):
	chrom,pos,ref,alt=parent_idx.split('_')
	#find the variant(s) at the position of interest for the current pLoF
	for variant in vcf_file('{0:s}:{1:s}-{1:s}'.format(chrom,pos)):
		#make sure the ref and alt alleles match what we are looking for
		if (variant.REF==ref) and (alt in variant.ALT):
			#make sure the variant passes recommended best practice filter for exome variant calls in the UKBB 
			######## DELETE THE FOLLOWING LINE AND ADJUST TABS IF USING THIS ON ALL OF US) ########
			if (variant.call_rate>=0.99) and ((variant.gt_depths>=10).sum()/(variant.gt_depths.shape[0])>=0.9):

				is_snp=variant.is_snp
				desired_allele=variant.ALT.index(alt)+1
				carrier_index = np.array([True if desired_allele in call else False for call in variant.genotypes])
				no_call_index=np.array([True if -1 in call else False for call in  variant.genotypes])
				corrected_af=carrier_index.sum()/(variant.num_called*2)
				if corrected_af <= 0.01: # Under HWE, this corresponds to a carrier frequency of 2%, so will definitely capture all variants with carrier freq <=0.001
					carrier_depths = variant.gt_depths[carrier_index]
					carrier_balance = variant.gt_alt_freqs[carrier_index]
					carrier_qscores = variant.gt_quals[carrier_index]

					#variant level quality control
					if is_snp:
						allowed=np.logical_and(np.logical_and(carrier_depths>=7,carrier_balance>=0.15),carrier_qscores>=30)
					else:
						allowed=np.logical_and(np.logical_and(carrier_depths>=10,carrier_balance>=0.2),carrier_qscores>=30)

					carriers=samples[carrier_index][allowed]
					no_calls=np.union1d(samples[no_call_index],samples[carrier_index][allowed==False])

					carrier_freq=len(carriers)/(samples.shape[0]-len(no_calls))
					new_call_rate=(samples.shape[0]-len(no_calls))/samples.shape[0]
					#carrier freq/call rate filter
					if  new_call_rate>=0.99 and (carrier_freq<=0.001):
						variant_table.loc[parent_idx,'CARRIER_FREQ']=carrier_freq
						variant_table.loc[parent_idx,'CALL_RATE']=new_call_rate
						variant_table.loc[parent_idx,'CARRIERS'].update(set(carriers))
						variant_table.loc[parent_idx,'NO_CALLS'].update(set(no_calls))
					else:
						print('Variant {0:s} was dropped, as it failed inclusion criteria.'.format(parent_idx))
						variant_table.drop(index=[parent_idx],inplace=True)

variant_table=variant_table.loc[variant_table['CARRIERS'].apply(lambda x: True if len(x)>0 else False)]

# Write the results to pickle format and tab-delimitted text
variant_table.to_pickle('/Path/to/Output/Carrier/Data/AllHaploLOFVariants_FilteredCarriers.pth')


variant_table_txt=variant_table.copy()
variant_table_txt['CARRIERS']=variant_table_txt['CARRIERS'].apply(lambda x: ','.join(list(x)))
variant_table_txt['NO_CALLS']=variant_table_txt['NO_CALLS'].apply(lambda x: ','.join(list(x)))
variant_table_txt.to_csv('/Path/to/Output/Carrier/Data/AllHaploLOFVariants_FilteredCarriers.txt',sep='\t')
		
